import numpy as np
import os, sys, errno
import copy
import fitsio
import galsim
from collections import OrderedDict
from astropy.io import fits
import warnings

# Balrog files
import tile as Tile
import filters
import grid
import injector
import balinput
import tile as Tile

# import pudb

#-------------------------------------------------------------------------------
# Config class and related functions

# The following base class is useful for accessing allowed parameter values
# without constructing a full config
class BaseConfig(object):

    # Define currently allowed and types for various objects. For `allowed`, inputs
    # MUST be one of the types. For `supported`, no guarantees are made if input is
    # not one of the types.

    _allowed_bands = 'griz'

    # TODO: Allow Piff when available!
    _supported_psf_types = ['DES_PSFEx']#, 'Piff'}
    _psf_extensions = {'DES_PSFEx' : 'psfexcat.psf'}#, 'Piff' : 'something.piff'}

    _non_inj_input_types = ['power_spectrum', 'nfw_halo', 'des_psfex', '_get']

    _valid_pos_sampling = ['uniform', 'sahar', 'RectGrid', 'HexGrid', 'MixedGrid']

    # NOTE: only used for background subtracted runs
    _valid_background_types = ['BKG', 'BKG+CCD', 'BKG+RN', 'BKG+SKY']

class Config(BaseConfig):
    '''
    Balrog Simulation configuration object. Contains the GalSim config file as
    well as additional simulation parameters.
    '''

    def __init__(self, args):
        # Save command-line arguments; args is type Namespace
        self.args = args
        self.config_dir = args.config_dir
        self.geom_file = args.geom_file
        self.tile_list_file = args.tile_list
        self.tile_dir = args.tile_dir
        self.psf_dir = args.psf_dir
        self.output_dir = args.output_dir
        self.vb = args.verbose
        self.nproc = args.nproc

        # NOTE: Most type checking of previous command-line args
        # (tile_list, geom_file, etc.) is handled in 'read_bal_gs_config()'
        # to allow inputs in config file
        if self.config_dir is None: self.config_dir = ''
        self.config_dir = os.path.abspath(self.config_dir)

        # Keeps track of current tile number
        self.tile_list = Tile.load_tile_list(self.tile_list_file)
        self.set_tile_num(self.tile_list[0])

        self._read_gs_config()
        self._load_tile_geometry()

        # Keeps track of current tile number
        self.tile_list = Tile.load_tile_list(self.tile_list_file)
        self.set_curr_tilename(self.tile_list[0])

        self._load_input_catalogs()

        if self.extinct_objs is True:
            self._load_ext_factors()
        else:
            self.ext_factors = None

        return

    def _read_gs_config(self):
        # Process .yaml config file
        # TODO: Allow multiple config types (JSON)

        self.gs_config_file = os.path.join(self.config_dir, self.args.config_file)

        # NOTE: `self.config` will be a list of `OrderedDict`'s with length given by
        # the number of configurations specified in the file. To work with a specific
        # simulated image, say with `galsim.config.BuildImage()`, you must pass a
        # *single* config file (e.g. BuildImage(config[3]) ).
        self.gs_config = galsim.config.process.ReadYaml(self.gs_config_file)

        # For now, we will not accept multi-output yaml files. We will work with
        # the loaded config object and build a multi-output config OrderedDict ourselves.
        if len(self.gs_config) != 1:
            raise AttributeError('For now, multi-output yaml files are not accepted! The ',
                                 'multi-output is handled with a newly created Balrog file.')
        else:
            # Keep track of config length, just to be sure
            self.gs_config_len = 1

        # Make a copy of the config dict as it exists now.
        self.orig_gs_config = copy.deepcopy(self.gs_config)

        # Process Balrog-specific params
        self._read_balrog_gs_config()

        # Keep track of whether gs_config has been modified
        self.gs_config_modified = False

        return

    def _read_balrog_gs_config(self):
        '''
        This helper function reads in Balrog-specific input parameters and sets defaults as
        necessary (Note: The type-checking is done in the GS class `BalrogImageBuilder`.
        Additionally, it ensures that only 1 of 'n_objects' or 'object_density' is passed.)
        '''

        self.parse_command_args()

        if self.gs_config[0]['image']['type'] != 'Balrog':
            raise ValueError('GalSim image type must be `Balrog`!')

        # NOTE: All of the type checking is now done in `injector.py` as part of the custom
        # GalSim image class
        conf = self.gs_config[0]['image']
        self.gs_config[0]['image'].update(injector.parse_bal_image_inputs(conf,
                                                                          self.gs_config[0]))

        # The above function now guarantees the required options are checked, present,
        # & safe
        im = self.gs_config[0]['image']
        self.realizations = im['realizations']
        self.n_realizations = im['n_realizations']
        self.n_objects = im['n_objects']
        self.object_density = im['object_density']
        self.bands = im['bands']
        self.data_version = im['version']
        self.run_name = im['run_name']
        self.extinct_objs = im['extinct_objs']
        self.inj_objs_only = im['inj_objs_only']
        self.pos_sampling = im['pos_sampling']

        self.bindx = dict(zip(self.bands, range(len(self.bands))))

        return

    def parse_command_args(self):
        '''
        Parse inputs that may have been passed as command-line arguments.
        '''

        # NOTE: config_dir and verbose have already been parsed correctly
        args = {'tile_list':self.tile_list_file, 'geom_file':self.geom_file, 'tile_dir':self.tile_dir,
                'psf_dir':self.psf_dir, 'output_dir':self.output_dir, 'nproc':self.nproc}

        base = {'tile_list':'image', 'geom_file':'image', 'tile_dir':'image', 'psf_dir':'image',
                'output_dir':'output', 'nproc':'image'}

        req = ['tile_list', 'geom_file']
        opt = ['tile_dir', 'psf_dir', 'output_dir', 'nproc']

        for arg, arg_val in args.items():
            try:
                # Have to handle `output_dir` to be consistent with GalSim
                if arg == 'output_dir':
                    config_val = self.gs_config[0][base[arg]]['dir']
                else:
                    config_val = self.gs_config[0][base[arg]][arg]
                if arg_val and config_val != arg_val:
                    if isinstance(arg_val, str) and (os.path.isfile(arg_val) or os.path.isdir(arg_val)):
                        # The following works for both files and directories
                        if not os.path.samefile(arg_val, config_val):
                            raise ValueError('Command-line argument {}={} '.format(arg, arg_val) +
                                            'is inconsistent with config value of {}!'.format(config_val))
                    else:
                        raise ValueError('Command-line argument {}={} '.format(arg, arg_val) +
                                        'is inconsistent with config value of {}!'.format(config_val))
                val = config_val

            except (KeyError, TypeError):
                if arg_val is None and arg in req:
                    raise ValueError('Must pass {} in command line or config file!'.format(arg))
                val = arg_val

            # Do any needed special processing of arg or set default
            if arg == 'tile_list':
                self.tile_list_file = os.path.abspath(val)
            elif arg == 'geom_file':
                self.geom_file = os.path.abspath(val)
            elif arg == 'tile_dir':
                if val is None:
                    val = ''
                self.tile_dir = os.path.abspath(val)
            elif arg == 'psf_dir':
                if val is None:
                    val = 'psfs'
                # Do not want to append CWD to `psf_dir`
                self.psf_dir = val
            elif arg == 'output_dir':
                if val is None:
                    val = 'balrog_outputs/'
                self.output_dir = os.path.abspath(val)
            elif arg == 'nproc':
                if val is None:
                    val = 1
                self.nproc = val
            # Can add others if needed
            # elif ...

        return

    def _load_tile_geometry(self):
        '''
        TODO: Make more general to allow for non-DES tiles.
        '''

        self.geom = fitsio.read(self.geom_file)
        tile_names = self.geom['TILENAME']
        self.tile_names = np.array([tile_name.strip() for tile_name in tile_names])

        # Unique area bounds from DES coadd tile geometry file
        uramin, uramax = self.geom['URAMIN'], self.geom['URAMAX']
        udecmin, udecmax = self.geom['UDECMIN'], self.geom['UDECMAX']
        self.u_areas = np.array([uramin, uramax, udecmin, udecmax])

        return

    def _load_ext_factors(self):
        # NOTE: In the future, we may want to generalize this to a different file.
        # For now, easiest to place this in the geometry file.

        # For now, file assumed to have structure of ['g','r','i','z']
        file_bindx = dict(zip('griz', range(4)))

        # NOTE: Can also grab EXTMAG, GALLONG, GALLAT, and MEANEBV for more complex
        # extinction methods
        try:
            ext_factors = np.array([self.geom['EXTFACT'][:, file_bindx[b]]
                                                 for b in self.bands])
        except KeyError:
            raise AttributeError('Column `EXTFACT` required in geometry file for setting '
                                 'extinction factors!')

        self.ext_factors = dict(zip(self.geom['TILENAME'], ext_factors.T))

        return

    def _load_input_catalogs(self):
        '''
        Load any relevant info from the input catalog(s)
        '''

        input_types = self.return_input_names()
        self.input_types = {}

        self.input_cats = {}
        self.input_nobjects = {}
        self.input_indx = {}

        self.inj_types= {}

        # TODO: We should generalize to a zeropoint for each input type, but this has
        # a lot of complications on the chip level. See chip._set_flux_factor()
        self.input_zp = None

        # Don't want to modify self.gs_config during processing
        gs_config = copy.deepcopy(self.gs_config[0])

        for i, input_type in enumerate(input_types):
            tname = self.curr_tilename
            input_obj = balinput.build_bal_input(input_type, gs_config, indx=i, tilename=tname)
            self.input_types[input_type] = input_obj
            self.inj_types[input_type] = input_obj.inj_type
            self.input_indx[input_type] = i
            if isinstance(input_obj, balinput.InputCatalog):
                # Return the number of objects in input catalog that will be injected
                # NOTE: Calling `ProcessInputNObjects()` builds the input catalog
                # in a minimal way that determines the number of objects from the
                # original catalog that make it past all of the mask cuts, etc.
                # self.input_nobjects = galsim.config.ProcessInputNObjects(gs_config)

                # Now that we are saving truth tables, it is necessary to load in the entire
                # catalog
                self.input_cats[input_type] = input_obj.cat
                self.input_nobjects[input_type] = input_obj.nobjects

            if input_obj.input_zp is not None:
                if self.input_zp is None:
                    self.input_zp = input_obj.input_zp
                elif input_obj.input_zp != self.input_zp:
                    raise ValueError('Balrog inputs with different zeropoints '
                                        'are not yet supported!')

        return

    def return_input_names(self):
        '''
        Return names of input types used for Balrog run.

        NOTE: This has fundamentally changed from the original input type checking.
        We now keep track of input types that are explicitly *not* allowed for injection
        to allow for native GalSim types.
        '''

        inputs = OrderedDict(self.gs_config[0]['input'])
        input_types = []

        for it in inputs:
            if it in self._non_inj_input_types:
                # We only want to include actual injection inputs due to the list structure
                # used by Galsim's multiple 'gal' injections.
                print('Input type {} is not currently supported for injection. '.format(it) +
                      'This is ok for non-injection inputs, e.g. a shear power spectrum. '
                      'Skipping this type for injection.')
                continue
            if it not in galsim.config.input.valid_input_types:
                print('Input type {} is not a native GalSim input type. '.format(it) +
                      'Make sure you have written and registered a valid input type.')
            input_types.append(it)

        return input_types

    def set_realization(self, i):
        '''
        Sets parameters relevant to the current injection realization i.
        '''

        self.curr_real = i

        return

    def set_curr_tilename(self, tilename):
        self.curr_tilename = tilename

        return

    def reset_gs_config(self):
        '''
        This function resets the gs_config after a tile is run.
        NOTE: This should only be used if multiple tiles are processed
        on a single batch.
        '''

        if self.gs_config_modified is True:
            # New tile, so reset config to original input
            self.gs_config = copy.deepcopy(self.orig_gs_config)
            self.gs_config_modified = False

        # If gs_config_modified is False, then there is nothing to reset!

        return

#-------------------------
# Related config functions

def setup_config(args):
    '''
    For now, just create new config object. Can make more complex later if needed.
    '''

    config = Config(args)

    return config
