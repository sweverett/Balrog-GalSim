import numpy as np
import os, sys, errno
import copy
import fitsio
import galsim
from collections import OrderedDict
from astropy.io import fits
import warnings

# Balrog files
import filters
import grid

# import pudb

#-------------------------------------------------------------------------------
# Config class and related functions

class Config(object):
    '''
    Balrog Simulation configuration object. Contains the GalSim config file as
    well as additional simulation parameters.
    '''

    # Define currently allowed and types for various objects. For `allowed`, inputs
    # MUST be one of the types. For `supported`, no guarantees are made if input is
    # not one of the types.

    _allowed_bands = 'grizy'

    # TODO: Allow Piff when available!
    _supported_psf_types = ['DES_PSFEx']#, 'Piff'}
    _psf_extensions = {'DES_PSFEx' : 'psfexcat.psf'}#, 'Piff' : 'something.piff'}

    _supported_input_types = ['ngmix_catalog', 'des_star_catalog', 'cosmos_chromatic_catalog',
                            'meds_catalog']
    _supported_gal_types = ['ngmix_catalog', 'cosmos_chromatic_catalog', 'meds_catalog']
    _supported_star_types = ['des_star_catalog']

    # Process configuration file
    def __init__(self, args):
        # Save command-line arguments; args is type Namespace
        self.args = args
        self.config_dir = args.config_dir
        self.geom_file = args.geom_file
        self.tile_list = args.tile_list
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
        self.tile_num = 0

        # Process GalSim config file
        self._read_gs_config()
        # Process geometry file
        self._load_tile_geometry()
        # Process input catalog(s)
        self._load_input_catalogs()

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
        Additionally, it ensures that only 1 of 'n_galaxies' or 'gal_density' is passed.)
        '''

        # TODO: Some of this type checking should be moved into `injector.py`

        if self.gs_config[0]['image']['type'] != 'Balrog':
            raise ValueError('GalSim image type must be \'Balrog\'!')

        self.parse_command_args()

        # Process input 'realizations' & 'n_realizations':
        try:
            reals = self.gs_config[0]['image']['realizations']
            if type(reals) is int:
                self.realizations = np.array([reals])
            elif isinstance(reals, list):
                if all(isinstance(x, int) for x in reals):
                    self.realizations = reals
                else:
                    raise TypeError('Passed realizations must be integers!')
            elif isinstance(reals, dict):
                try:
                    max_real = reals['max']
                    min_real = reals['min']
                    if min_real >= max_real:
                        raise ValueError('Key `min_real` must be less than `max_real`!')
                    if all(isinstance(x, int) for x in [min_real, max_real]):
                        self.realizations = [x for x in range(min_real, max_real+1)]
                    else:
                        raise TypeError('The realization values `min` and `max` must be ints!')
                except KeyError as e:
                    print(e + '\nThe only valid keys for `realizations` are `min` and `max`!')
            # Now check if `n_realizations` was also passed
            # NOTE: If n_realizations>len(realizations), this indicates that a larger simulation
            # is being split up between multiple runs. The main impact is star generation for
            # Y3 DES star catalogs, as they cannot be shuffled in this case to ensure all stars
            # are injected w/o repeats.
            try:
                n_reals = self.gs_config[0]['image']['n_realizations']
                if isinstance(n_reals, int):
                    if n_reals < len(self.realizations):
                        raise ValueError('`n_realizations` cannot be smaller than len(realizations).')
                    self.n_realizations = n_reals
                else:
                    raise TypeError('The value `n_realizations` must be an int!')
            except KeyError:
                # In this case, assume that n_realizations=len(realizations)
                self.n_realizations = len(self.realizations)

        except KeyError:
            try:
                n_reals = self.gs_config[0]['image']['n_realizations']
                # If it has been passed, warn user but still use
                warnings.warn('DEPRECATED: `n_realizations` without `realizations` has been ' +
                              'deprecated. Please use argument `realizations` (or both) instead.')
                if isinstance(n_reals, int):
                    self.n_realizations = n_reals
                    self.realizations = [x for x in range(n_reals)]
                else:
                    raise TypeError('The value `n_realizations` must be an int!')
            except KeyError:
                # Default is to use realization 0
                warnings.warn('No realization passed; using default of 0.')
                self.n_realizations = 1
                self.realizations = np.array([0])

        # Process input 'n_galaxies':
        try:
            # Number of galaxies per tile
            self.n_galaxies = self.gs_config[0]['image']['n_galaxies']
        except KeyError:
            self.n_galaxies = None

        # Process input 'gal_density':
        try:
            # Assumes units in galaxies per arcmin^2
            self.gal_density = self.gs_config[0]['image']['gal_density']
            # TODO: Should allow for more unit input types! e.g.:
            # den_unit = self.gs_config[0]['image']['gal_density']['unit']
            # self.gal_density = convert_units(val=self.gal_density, u1 = den_unit, u2 = 'arcmin^2')
        except KeyError:
            self.gal_density = None

        # This is ensured by `BalrogImageBuilder`, but not bad to check it here
        assert self.n_galaxies is None or self.gal_density is None

        # If both 'gal_density' and 'gal_density' are None, use default case:
        if self.n_galaxies is None and self.gal_density is None:
            self.n_galaxies = 1000 # Keep it small for default!
            print('Warning: Neither n_galaxies nor gal_density was passed in config file. ' +
                  'Using default of {} galaxies per tile'.format(self.n_galaxies))

        # Process input 'bands'
        try:
            # Grab selected bands in config, if present
            self.bands = self.gs_config[0]['image']['bands']

            # Make sure there aren't any incorrect inputs
            for band in self.bands:
                if band not in self._allowed_bands:
                    raise ValueError('Passed band {} is not one of the allowed bands in {}'\
                                     .format(band, _allowed_bands))
        except KeyError:
            # By default, use griz
            print('Warning: No injection bands were passed in config. Using `griz` by default')
            self.bands = 'griz'

        # Process input 'version'
        try:
            self.data_version = self.gs_config[0]['image']['version']
        except KeyError:
            # Warn user, but assume y3v02 for now
            warnings.warn('Data version not passed in config! Assuming y3v02.')
            self.data_version = 'y3v02'

        # Process input 'run_name'
        try:
            rname = self.gs_config[0]['image']['run_name']
            if not isinstance(rname, basestring):
                raise ValueError('The input `run_name` must be a string!')
            self.run_name = rname
        except KeyError:
            # TODO: Maybe come up with sensible default run name?
            #       Current metadata should provide enough info for now.
            # self.run_name = 'None'
            self.run_name = None

        # Process input 'inj_objs_only'. This is used to test Balrog injections on blank images
        try:
            inj_objs_only = self.gs_config[0]['image']['inj_objs_only']
            if type(inj_objs_only) is bool:
                # Default is to include chip noise in injection
                self.inj_objs_only = {'value':inj_objs_only, 'noise':'CCD'}
            elif isinstance(inj_objs_only, dict):
                # Is likely an OrderedDict, so convert
                inj_objs_only = dict(inj_objs_only)
                self.inj_objs_only = {}
                keys = ['value', 'noise']
                valid_noise = ['CCD', 'BKG', 'BKG+CCD', 'BKG+RN', 'BKG+SKY', 'None', None]

                if 'noise' not in inj_objs_only:
                    # Default is no noise
                    inj_objs_only['noise'] = None
                for key, val in inj_objs_only.items():
                    if key not in keys:
                        raise ValueError('{} is not a valid key for `inj_objs_only`! '.format(key) +
                                         'You may only pass the keys {}'.format(keys))
                    if (key == 'noise') and (val not in valid_noise):
                        raise ValueError('{} is not a valid value for the noise field!'.format(val) +
                                         'You may only pass the values {}'.format(valid_noise))
                    self.inj_objs_only[key] = val
            else:
                raise ValueError('The field \'inj_objs_only\' must be set with a bool or dict!')
        except KeyError:
            # Most runs will add objects to existing images
            self.inj_objs_only = {'value':False, 'noise':None}

        # Process input 'pos_sampling'
        self.pos_sampling = {}
        valid_pos_sampling = grid._valid_pos_sampling
        valid_grid_types= grid._valid_grid_types
        default_gs = 20 # arcsec
        try:
            ps = self.gs_config[0]['image']['pos_sampling']
            if isinstance(ps, basestring):
                # Then the string is the input type
                if ps not in valid_pos_sampling:
                    raise ValueError('{} is not a valid position sampling method. '.format(ps) +
                                    'Currently allowed methods are {}'.format(valid_pos_sampling))

                if ps in valid_grid_types:
                    print('No grid spacing passed; using default of {} arcsecs'.format(default_gs))
                    self.pos_sampling = {'type' : ps, 'grid_spacing' : default_gs}
                else:
                    self.pos_sampling = {'type' : ps}
            elif isinstance(ps, dict):
                if 'type' not in ps.keys():
                    raise ValueError('If `pos_sampling` is passed as a dict, then must set a type!')
                if (ps in valid_grid_types) and ('grid_spacing' not in ps.keys()):
                    print('No grid spacing passed; using default of {} arcsecs'.format(default_gs))
                    ps['grid_spacing'] = default_gs
                keys = ['type', 'grid_spacing', 'rotate', 'offset', 'angle_unit']
                for key, val in ps.items():
                    if key not in keys:
                        raise ValueError('{} is not a valid key for `pos_sampling`! '.format(key) +
                                         'You may only pass the keys {}'.format(keys))
                    if key == 'grid_spacing':
                        if val < 0.0:
                            raise ValueError('grid_spacing of {} is invalid; '.format(val) +
                                             'must be positive!')
                    if key == 'rotate':
                        if isinstance(val, str):
                            if val.lower() != 'random':
                                raise ValueError('{} is not a valid grid rotation type!'.format(val) +
                                                 'Must be `Random` or a number.')
                        self.pos_sampling.update({'rotate':val})
                        try:
                            unit = ps['angle_unit'].lower()
                            if unit == 'deg':
                                if (val<0.0) or (val>360.0):
                                    raise ValueError('rotate value of {} deg is invalid!'.format(val))
                            if unit == 'rad':
                                if (val<0.0) or (val>2*np.pi):
                                    raise ValueError('rotate value of {} rad is invalid!'.format(val))
                                else:
                                    raise ValueError('angle_unit of {} is invalid! '.format(unit) +
                                                     'only can pass `deg` or `rad`.')
                            self.pos_sampling.update({'angle_unit':unit})
                        except KeyError:
                            # Default of rad
                            if (val<0.0) or (val>2*np.pi):
                                raise ValueError('rotate value of {} rad is invalid!'.format(val))
                            self.pos_sampling.update({'angle_unit':'rad'})
                        self.pos_sampling.update({'rotate':val})

                    if key == 'offset':
                        if isinstance(val, str):
                            if val.lower() != 'random':
                                raise ValueError('{} is not a valid grid offset!'.format(val) +
                                'must be `Random` or an array.')
                            self.pos_sampling.update({'offset':val})
                        elif isinstance(val, list):
                            assert len(val)==2
                            self.pos_sampling.update({'offset':val})
                        else:
                            raise TypeError('grid offset of type {} is invalid!'.format(type(val)))


                    self.pos_sampling[key] = val

        except KeyError:
            # Most runs use uniform sampling
            self.pos_sampling['type'] = 'uniform'
            self.pos_sampling['grid_spacing'] = None

        return

    def parse_command_args(self):
        '''
        Parse inputs that may have been passed as command-line arguments.
        '''

        # pudb.set_trace()

        # NOTE: config_dir and verbose have already been parsed correctly
        args = {'tile_list':self.tile_list, 'geom_file':self.geom_file, 'tile_dir':self.tile_dir,
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
                self.tile_list = os.path.abspath(val)
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

        # Load unique area bounds from DES coadd tile geometry file
        with fits.open(self.geom_file) as hdu_geom:
            self.geom = hdu_geom[1].data
            self.tile_names = self.geom['TILENAME']
            uramin, uramax = self.geom['URAMIN'], self.geom['URAMAX']
            udecmin, udecmax = self.geom['UDECMIN'], self.geom['UDECMAX']
            # Unique tile area
            self.u_areas = np.array([uramin, uramax, udecmin, udecmax])
            # pudb.set_trace()

        return

    def _load_input_catalogs(self):
        '''
        Load any relevant info from the input catalog(s)
        '''

        # Determine input type
        input_cat_types = self._determine_input_types()

        self.input_cats = {}
        self.input_nobjects = {}

        # Keep track of which index corresponds to gals vs stars.
        # NOTE: Only one input catalog of each type is currently allowed!
        self.input_indx = {}
        self.input_types = {} # Input catalogs
        self.inj_types = {} # Injection types from catalogs
        self.sim_gals = False
        self.sim_stars = False

        # Don't want to modify self.gs_config during processing
        gs_config = copy.deepcopy(self.gs_config[0])

        for i, input_type in enumerate(input_cat_types):
            # TODO: This section could be generalized by making new `StarCatalog` and
            # `GalaxyCatalog` classes that users can set the following methods for, and
            # simply call those methods from here.

            if input_type in self._supported_gal_types:
                # Check that a galaxy cat type hasn't already been set
                if self.sim_gals is True:
                    raise ValueError('Can\'t set multiple input galaxy catalogs!')
                else:
                    self.sim_gals = True
                    self.input_indx['gals'] = i
                    self.input_types['gals'] = input_type

                    # TODO: Can we grab the injection type from the registered GS catalog?
                    # pudb.set_trace()
                    if input_type in ['ngmix_catalog', 'meds_catalog']:
                        if input_type == 'ngmix_catalog':
                            from ngmix_catalog import ngmixCatalog as CATMOD
                            from ngmix_catalog import ngmixCatalogLoader as LOADMOD
                        if input_type == 'meds_catalog':
                            from meds_catalog import MEDSCatalog as CATMOD
                            from meds_catalog import MEDSCatalogLoader as LOADMOD

                        # As we are outside of the GalSim executable, we need to register
                        # the input type explicitly
                        galsim.config.RegisterInputType(input_type, LOADMOD(CATMOD, has_nobj=True))

                        base_cat = input_type.split()[0]
                        self.inj_types['gals'] = base_cat.split('_')[0]+ 'Galaxy'

                        # This avoids a printed warning, and sets up the input correctly
                        # as no bands are passed in bal_config
                        gs_config['input'][input_type]['bands'] = 'griz'

                        # Version-specific settings
                        if self.data_version == 'y3v02':
                            # TODO: See if we can load this, rather than set explicitly (but true for y3v02)
                            # Set input catalog zeropoint
                            self.input_zp = 30.0
                        else:
                            # In future, can add updated input parsing
                            raise ValueError('No input parsing defined for {} catalogs for ' +
                            'data version {}. (y3v02 for DES Y3A2)'.format(input_type, self.data_version))

                    elif input_type == 'cosmos_chromatic_catalog':
                        self.inj_types['gals'] = 'COSMOSChromaticGalaxy'

                        # As we are outside of the GalSim executable, we need to register
                        # the input type explicitly
                        import input_cosmos_chromatic as icc
                        import scene_chromatic as sc
                        galsim.config.RegisterInputType('cosmos_chromatic_catalog',
                                                        icc.COSMOSChromaticLoader(
                                                        sc.COSMOSChromaticCatalog, has_nobj=True))

                        # Process a few fields only present for chromatic COSMOS galaxies
                        # pudb.set_trace()
                        try:
                            filter_dir = self.gs_config[0]['input']['cosmos_chromatic_catalog']['filter_dir']
                            self.filter_dir = os.path.abspath(filter_dir)
                        except KeyError:
                            # Default is set to cwd in Filter()
                            self.filter_dir = None
                        try:
                            use_filter_tables = self.gs_config[0]['input']['cosmos_chromatic_catalog']['use_filter_tables']
                            if type(use_filter_tables) != bool:
                                raise TypeError('The type of `use_filter_tables` must be a bool! ' +
                                                'Was passed a {}'.format(type(use_filter_tables)))
                            self.use_filter_tables = use_filter_tables

                        except KeyError:
                            # Default depends on other parameters
                            if self.filter_dir is None:
                                warnings.warn('Neither `filter_dir` nor `use_filter_tables` ' +
                                                'passed for input cosmos_chromatic_catalog. ' +
                                                'Using a constant throughput for COSMOS galaxies.')
                                self.use_filter_tables = False
                            else:
                                self.use_filter_tables = True

                        self.filters = filters.Filters(self.bands,
                                                use_transmission_tables=self.use_filter_tables,
                                                filter_dir=self.filter_dir)

                        # TODO: Check COSMOS zeropoint!
                        # Set input catalog zeropoint
                        self.input_zp = 25.94
                        # self.input_zp = 30.0

            elif input_type == 'des_star_catalog':
                if self.data_version == 'y3v02':
                    # Check that a star cat type hasn't already been set
                    if self.sim_stars is True:
                        raise ValueError('Can\'t set multiple input star catalogs!')
                    else:
                        self.sim_stars = True
                        self.input_indx['stars'] = i
                        self.input_types['stars'] = input_type
                        # TODO: Can we grab the injection type from the registered GS catalog?
                        self.inj_types['stars'] = 'desStar'

                    # pudb.set_trace()

                    import des_star_catalog

                    # As we are outside of the GalSim executable, we need to register
                    # the input type explicitly
                    galsim.config.RegisterInputType('des_star_catalog',
                            des_star_catalog.desStarCatalogLoader(
                            des_star_catalog.desStarCatalog, has_nobj=True))

                    valid_model_types = des_star_catalog.return_valid_model_types(
                                            data_version=self.data_version)
                    try:
                        self.star_model = gs_config['input']['des_star_catalog']['model_type']
                    except KeyError as e:
                        raise KeyError('Must pass a model_type if using a des_star_catalog ' +
                                       'for Balrog injections! See `des_star_catalog.py` for details.')

                    if self.star_model not in valid_model_types:
                        raise ValueError('The selected model_type {} '.format(self.star_model) +
                                         'is not a valid des_star_catalog model type!\n ' +
                                         'Valid types for data version {} are: {}'.format(
                                             self.data_version, valid_model_types))

                    prefix = self.star_model.split('_')[0]
                    if prefix == 'Model':
                        self.base_percent = 100
                    elif prefix == 'Extra':
                        self.base_percent = int(self.star_model.split('_')[1])
                    else:
                        # Then something very strange happened! Should have been caught above
                        raise ValueError

                    # Add bands to suppress a warning.
                    gs_config['input'][input_type]['bands'] = 'griz'

                    # A tile name is also a required input, so grab the first one
                    tile_list = load_tile_list(self.tile_list, self)
                    first_tile = tile_list[0]
                    gs_config['input'][input_type]['tile'] = first_tile

                    # TODO: See if we can load this, rather than set explicitly (but true for y3v02)
                    # Set input catalog zeropoint
                    self.input_zp = 30.0

                else:
                    # In future, can add updated input parsing
                    raise ValueError('No input parsing defined for DES star catalogs for ' +
                    'data version {}'.format(self.data_version))

            else:
                # Add more types later!
                # warnings.warn('Input type {} not used as an injected object type. '.format(input_type) +
                #               'May still be used to modify injections, e.g. global shear parameters.')
                # pudb.set_trace()
                input_cat_types.remove(input_type)
                raise ValueError('For now, only {} can be used for injections!'.format(_supported_input_types))

        # Return the number of objects in input catalog that will be injected
        # NOTE: Calling `ProcessInputNObjects()` builds the input catalog
        # in a minimal way that determines the number of objects from the
        # original catalog that make it past all of the mask cuts, etc.
        # self.input_nobjects = galsim.config.ProcessInputNObjects(gs_config)

        # Now that we are saving truth tables, it is necessary to load in the entire
        # catalog
        # TODO: possible race condition here?
        galsim.config.ProcessInput(gs_config)

        # Grab needed info from the proxy catalog
        # for i, input_type in enumerate(input_cat_types):
        # for  input_type in enumerate(self.input_types)
        for input_type in self.input_types.values():
            # Only do this for injection types (gals, stars)
            # if i not in self.input_indx.values(): continue
            cat_proxy = gs_config['input_objs'][input_type][0] # Actually a proxy
            self.input_cats[input_type] = cat_proxy.getCatalog()
            self.input_nobjects[input_type] = cat_proxy.getNObjects()

            # Need to load in additional parametric catalog for MEDSCatalog truth table
            if input_type == 'meds_catalog':
                self.meds_param_catalog = cat_proxy.getParamCatalog()

        return

    def _determine_input_types(self):
        '''
        # TODO: Check that is is one of the supported input types!
        # TODO: While just a placeholder method anyway, would be good to check that the
        # input catalog actually *is* an ngmix catalog!
        '''

        inputs = OrderedDict(self.gs_config[0]['input'])
        self.input_types = []

        for it in inputs:
            if it not in self._supported_input_types:
                # We only want to include actual injection inputs due to the list structure
                # used by Galsim's multiple 'gal' injections.
                warnings.warn('Input type {} is not currently supported for injection. '.format(it) +
                              'This is ok for other inputs, e.g. a shear power spectrum.')
                continue
            self.input_types.append(it)

        return self.input_types

    def set_realization(self, i):
        '''
        Sets parameters relevant to the current injection realization i.
        '''

        self.curr_real = i

        return

    def set_tile_num(self, i):
        self.tile_num = i

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
