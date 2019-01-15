import fileio as io
import grid
import galsim
import os, sys, errno
import ntpath
import fitsio
import yaml
import numpy as np
from astropy import wcs
from copy import deepcopy
import warnings
import subprocess
import datetime

# Balrog files
from chip import Chip
import fileio as io
import mathutil as util
import balobject as balobj
import grid

import pudb

#-------------------------------------------------------------------------------
# Tile class and functions

class Tile(object):
    """
    A `Tile` is a square ~0.7x0.7 deg^2 subsection of the DES footprint.
    Tiles overlap by 2 arcmin, but Balrog galaxies are only injected in the unique
    footrpint area.
    """

    def __init__(self, tile_name, config, realization=0):

        self.tile_name = tile_name

        # Grab input types
        self.input_types = config.input_types
        self.input_indx = config.input_indx
        self.inj_types = config.inj_types

        # Useful to store a few directories in each tile so config isn't passed around
        self.output_dir = config.output_dir

        # Tile index in geometry file list
        tile_names = config.tile_names
        indx = np.where(tile_names == self.tile_name)[0]
        if len(indx) == 0:
            raise ValueError('Tile {} is not contained in the tile list!'.format(self.tile_name))
        elif len(indx) > 1:
            raise ValueError('Tile {} appears more than once in the tile list!'.format(self.tile_name))
        else:
            self.indx = int(indx)

        self._determine_unique_area(config)
        self._set_wcs(config)

        # NOTE: `n_objects` and `object_density` are now defined to be *per realization*
        # Set the number of objects injected per realization
        self.objs_per_real = {}
        for inpt in self.input_types:
            if config.n_objects[inpt] is not None:
                # Fixed number of objects regardless of tile size
                self.objs_per_real[inpt] = config.n_objects[inpt]
            elif config.object_density[inpt] is not None:
                self.objs_per_real[inpt] = round(self.u_area * config.object_density)
            else:
                # This should only happen for grid runs
                assert (config.pos_sampling[inpt]['type'] in grid.BaseGrid()._valid_grid_types) \
                    or (config.pos_sampling[inpt]['type'] in grid.BaseGrid()._valid_mixed_types)
                # This will be set during grid creation
                self.objs_per_real[inpt] = None

        # Set tile directory structure
        self.dir = os.path.abspath(os.path.join(config.tile_dir, self.tile_name))

        self._set_bands(config)

        self.set_realization(realization)

        # Setup new Balrog config file for chip to be called by GalSim executable
        self._setup_bal_config(config)

        self._load_zeropoints(config)

        # Load background images, if needed
        self._load_backgrounds(config)

        self._set_noise(config)

        self._create_chip_list(config)

        # Keep track if any injections into tile have been made
        self.has_injections = False

        return

    def _determine_unique_area(self, config):
        '''Set the tile unique area for injections.'''

        # In [ramin, ramax, decmin, decmax] format:
        # self.u_area = config.u_areas[:, self.indx]
        self.ramin, self.ramax, self.decmin, self.decmax = config.u_areas[:, self.indx]

        # Account for tiles that cross over 360/0 boundary
        if self.ramin > self.ramax:
            # NOTE: For consistency w/ catalog, leave them flipped from expectation
            # ramin/ramax will always mean 'left'/'right' boundary
            # self.ramin, self.ramax = self.ramax, self.ramin
            self.ra_boundary_cross = True
        else:
            self.ra_boundary_cross = False

        # convert dec to radians, but keep ra in deg
        d1, d2 = np.deg2rad([self.decmin, self.decmax])
        r1, r2 = self.ramin, self.ramax

        # Set flag if tile goes over ra=0/360 deg boundary & correct r2
        if self.ramin > self.ramax:
            self.ra_flag = True
            r2 = r2 + 360.0
        else: self.ra_flag = False

        # In deg^2
        a = (180.0 / np.pi) * (r2 - r1) * (np.sin(d2) - np.sin(d1))
        # Save in arcmin^2
        self.u_area = 3600.0 * a

        return

    def _set_wcs(self, config):
        '''
        Load WCS info for each tile from geometry file.
        '''

        crpix1, crpix2 = config.geom['CRPIX1'][self.indx], config.geom['CRPIX2'][self.indx]
        crval1, crval2 = config.geom['CRVAL1'][self.indx], config.geom['CRVAL2'][self.indx]
        ctype1, ctype2 = config.geom['CTYPE1'][self.indx], config.geom['CTYPE2'][self.indx]
        cd1_1, cd1_2 = config.geom['CD1_1'][self.indx], config.geom['CD1_2'][self.indx]
        cd2_1, cd2_2 = config.geom['CD2_1'][self.indx], config.geom['CD2_2'][self.indx]

        # Create WCS object
        self.wcs = wcs.WCS()
        self.wcs.wcs.crpix = [crpix1, crpix2]
        self.wcs.wcs.crval = [crval1, crval2]
        self.wcs.wcs.ctype = [ctype1, ctype2]
        self.wcs.wcs.cd = [[cd1_1, cd1_2], [cd2_1, cd2_2]]

        # Set pixel information
        self.pixel_scale = config.geom['PIXELSCALE'][self.indx]
        # These are (ra, dec) and both 10,000 for DES tiles
        self.Npix_x = config.geom['NAXIS1'][self.indx]
        self.Npix_y = config.geom['NAXIS2'][self.indx]

        return

    def _set_bands(self, config):
        '''
        For now, just set to 'griz'. May want to use a different subset in the future.
        '''

        # For convenience of later functions
        self.bands = config.bands

        # Will store directory locatations for each band
        self.band_dir = {}

        for band in self.bands:
            self.band_dir[band] = os.path.join(self.dir, 'nullwt-{}'.format(band))
            # Make band directory if it doesn't already exist
	    try:
	        os.makedirs(self.band_dir[band])
            except OSError as e:
                if e.errno == errno.EACCES:
                    # ok if directory already exists
                    # print('permission error')
                    pass
                elif e.errno == errno.EEXIST:
                    # ok if directory already exists
                    pass
                else:
                    raise e

        return

    def _setup_bal_config(self, config):
        '''
        Sets up the Balrog configuration file that will house chip-specific parameters for
        eventual GalSim executable call. The output bal_config will be a list of dictionaries
        that contain the simulation parameters; the first has tile-wide parameters will
        subsequent appended entries have chip-specific simulation parameters. The main
        difference between this implementation and the gs_config one is that bal_config will
        eventually be re-written to a yaml file to be run by the GalSim executable rather than
        manually using GalSim python functions.
        # TODO: Should allow other formats (e.g. JSON) in future.
        '''

        self.set_bal_config_name(config)

        # Create bal config directory if it does not already exist
        if not os.path.exists(self.bal_config_dir):
            try:
                os.makedirs(self.bal_config_dir)
            except OSError as e:
                # Check for race condition
                if e.errno != errno.EEXIST:
                    raise

        # Balrog config will simply append the gs_config entered as an input. Each chip
        # to be injected will get a new section of the multi-output yaml file.
        with open(config.gs_config_file) as f:
            list_doc = yaml.safe_load_all(f)
            self.bal_config = list(list_doc)

        # Keep track of length of the multi-output yaml config file; want to make sure that
        # any major changes are monitored
        # Should have been enforced in gs_config init; just to make sure!
        assert len(self.bal_config) == 1
        self.bal_config_len = 1

        # NOTE: The GalSim documentation states that it is more efficient to put `nproc` in global
        # output field when constructing lots of files, but we have not found that to be the case
        try:
            self.bal_config[0]['image'].update({'nproc':config.nproc})
        except KeyError:
            self.bal_config[0]['image'] = {'nproc':config.nproc}

        return

    def _load_zeropoints(self, config, s_begin=0, s_end=4):
        '''
        Construct {chip : zeropoint} dictionaries for each band using the following files:
        {tile}/lists/{tile}_{band}_nullwt-flist-{version}.dat
        `s_begin` and `s_end` are used to determine how to grab a chip name from the chip
        filename.
        '''

        zp_dir = os.path.join(config.tile_dir, self.tile_name, 'lists')
        self.zeropoint_files = {}
        self.zeropoints = {}

        for band in self.bands:
            zp_filename = '{}_{}_nullwt-flist-{}.dat'.format(self.tile_name, band, config.data_version)
            self.zeropoint_files[band] = os.path.join(zp_dir, zp_filename)

            # Will store zeropoints as {chip_name : zeropoint} for each band
            self.zeropoints[band] = {}

            with open(self.zeropoint_files[band]) as f:
                for line in f.readlines():
                    line_data = line.replace('\n', '').split(' ')
                    # Determine chip name from path
                    chip_file, zp = ntpath.basename(line_data[0]), float(line_data[1])
                    # NOTE: Easier to handle chip filename rather than actual name
                    # chip_name = '_'.join(chip_file.split('_')[s_begin:s_end])

                    # self.zeropoints[band][chip_name] = zp
                    self.zeropoints[band][chip_file] = zp

        return

    def _load_backgrounds(self, config, s_begin=0, s_end=4):
        '''
        Load any needed background images.
        # NOTE: For now, these are only used for grid test images.
        '''

        # TODO: Change this to always look for 'BKG' or 'BKG+' inside of noise model!
        if config.inj_objs_only['noise'] in config._valid_background_types:
            self.bkg_file_list = {}
            self.bkg_files = {}
            for band in config.bands:
                self.bkg_files[band] = {}
                self.bkg_file_list[band] = os.path.join(self.dir, 'lists',
                    '{}_{}_bkg-flist-{}.dat'.format(self.tile_name, band, config.data_version))
                with open(self.bkg_file_list[band]) as f:
                    for line in f.readlines():
                        line_data = line.replace('\n', '').split(' ')
                        chip_file, file_name = line_data[0], ntpath.basename(line_data[0])
                        # NOTE: Only here for blank testing!
                        chip_name = '_'.join(file_name.split('_')[s_begin:s_end])
                        self.bkg_files[band][chip_name] = chip_file
        else:
            self.bkg_file_list = None
            self.bkg_files = None

        return

    def _set_noise(self, config):
        '''
        Set any tile-wide noise properties from the config.
        '''

        # For now, can be 'CCD', 'BKG', 'BKG+CCD', 'BKG+RN', 'BKG+SKY', or None
        self.noise_model = config.inj_objs_only['noise']

        return

    def _create_chip_list(self, config):
        '''
        Given nullweight file locations, create chip lists for each band.
        '''

        # Store all chips in dictionary
        self.chip_list = {}
        self.chips = {}

        for band, b_dir in self.band_dir.items():
            self.chip_list[band] = []
            self.chips[band] = []
            # Get list of files in nullwt directory
            try:
                file_list = os.listdir(b_dir)
            except OSError:
                #
                file_list = None
                continue
            # Check that it is an appropriate fits file to be injected into
            for f in file_list:
                if self.is_chip_image(config, f): self.chip_list[band].append(f)

                # Grab chip zeropoint for this file
                # QUESTION: Should we allow chips that aren't in the .dat file and
                #           just assign a zp of 30?
                zp = self.zeropoints[band][f]

                # Add chip to list
                filename = os.path.join(b_dir, f)
                self.chips[band].append(Chip(filename, band, config, tile_name=self.tile_name,
                                        zeropoint=zp, tile=self))


        return


    def is_chip_image(self, config, chip_file):
        '''
        Checks if passed file is an appropriate chip injection image given data version.
        '''

        if config.data_version == 'y3v02':
            # TODO: Make more generic in future
            if chip_file.endswith('nullwt.fits'):
                return True
            else:
                return False

    def set_realization(self, i):
        '''
        Sets parameters relevant to the current injection realization i.
        '''

        self.curr_real = i

        # If necessary, add more relevant details

        return

    def generate_objects(self, config, realization):
        '''
        Generate a list of positions and indices for objects (stars, galaxies, etc.) for a
        given realization (starts counting at 0).
        '''

        # NOTE: While in principle the object generation could depend explicitly on
        # realization, for now all generation types can create realization positions
        # and indices all at once during construction
        if realization == 0:
            self.inj_cats = balobj.BalInjectionCatalogs(config)
            self.inj_cats.generate_catalogs(config, self, realization)

            # NOTE: Single-object injection is supported for only the first input type
            # (Only used for testing)
            for input_type in self.input_types.keys():
                if self.inj_cats[input_type].single_obj_injection is True:
                    # Remove `index` from global config (can cause issues w/
                    # other input types
                    self.bal_config[0]['gal'].pop('index', None)
                    break
        else:
            # See note above
            pass

        return

    def write_bal_config(self):
        '''
        Write appended balrog config to a yaml file.
        TODO: In future, allow more config types! JSON?
        '''

        # Writes final tile-wide GalSim config file for current realization
        with open(self.bal_config_file, 'w') as f:
            yaml.dump_all(self.bal_config, f)

        return

    def reset_bal_config(self, config):
        '''
        Return bal_config to default state. For details, look at setup_bal_config().
        '''

        # Make sure that bal_config already exists
        if (not hasattr(self, 'bal_config')) or (self.bal_config is None):
            raise AttributeError('Bal_config does not exist yet - only use `reset_bal_config()`',
                                 ' after initializing with `setup_bal_config()`!')

        # Checks have already been made, so simply load in file.
        with open(config.gs_config_file) as f:
            list_doc = yaml.safe_load_all(f)
            self.bal_config = list(list_doc)

        assert len(self.bal_config) == 1
        self.bal_config_len = 1

        # Update bal_config filename for current realization
        self.set_bal_config_name(config)

        # Reset injections
        self.has_injections = False

        # Place `nproc` in the config if it isn't there already
        try:
            self.bal_config[0]['image'].update({'nproc':config.nproc})
        except KeyError:
            self.bal_config[0]['image'] = {'nproc':config.nproc}

        return

    def set_bal_config_name(self, config):
        '''
        Sets the correct balrog config filename given the current realization.
        '''
        filename = 'bal_config_' + str(self.curr_real) + '_' + self.tile_name + '.yaml'
        self.bal_config_dir = config.output_dir + '/configs/'
        self.bal_config_file = self.bal_config_dir + filename

        return

    def add_gs_injection(self, config, chip, input_type, real):
        '''
        This function appends the global GalSim config with an additional simulation to
        be done using nullwt chip-specific infomation and Balrog injected object positions
        in image coordinates.
        '''

        assert input_type in config.input_types
        inj_type = config.inj_types[input_type]

        # Determine which Balrog objects are contained in chip, and
        # get their image coordinates
        inj_cat = self.inj_cats[input_type]
        in_chip, pos_im = chip.contained_in_chip(inj_cat.pos[real])
        inj_pos_im = pos_im[in_chip]
        inj_indx = inj_cat.indx[real][in_chip]
        assert len(inj_indx) == len(inj_pos_im)
        chip.set_nobjects(len(inj_pos_im), inj_type)

        # Skip chips with no injections, except for a few special cases
        if (len(inj_indx) == 0) and (config.inj_objs_only['value'] is False):
            if config.vb > 1:
                print('No objects of type {} were passed on injection '.format(inj_type) +
                'to chip {}. Skipping injection.'.format(chip.name))
                return

        # If this is the first injection for the chip, then set up new entry in bal_config
        if chip.types_injected == 0:
            self.add_bal_config_entry()

        # Injection index
        i = self.bal_config_len - 1

        Ninput = len(self.input_types)
        Ninject = len(inj_indx)

        #-----------------------------------------------------------------------------------------------
        # Now set up common config entries if chip's first injection

        if chip.types_injected == 0:

            # Setup 'image' field
            chip_file = chip.filename
            self.bal_config[i]['image'] = {
                'initial_image' : chip_file,
                'wcs' : { 'file_name' : chip_file },
            }

            # The field `nobjects` will be set to the sum of each injection type contained
            # in the chip
            icount = 0
            nobjs = ''
            for itype in config.inj_types.values():
                if icount == 0:
                    pfx = '$'
                else:
                    pfx = '+'
                nobjs += pfx + '@image.N_{}'.format(itype)
                # Must initialize to 0 as size-zero injections are skipped
                self.bal_config[i]['image']['N_{}'.format(itype)] = 0
                icount +=1

            self.bal_config[i]['image']['nobjects'] = nobjs

            # TODO: Eventually, we should move the noise models to be specific for
            # each BalObject
            # If noise is to be added, do it here
            if self.noise_model:
                if self.noise_model in ['CCD', 'BKG+CCD']:
                    self.bal_config[i]['image']['noise'] = {
                        'type' : 'CCD',
                        'sky_level_pixel' : chip.sky_sigma**2,
                        'gain' : float(np.mean(chip.gain)),
                        'read_noise' : float(np.mean(chip.read_noise))
                    }
                elif self.noise_model in ['BKG+RN', 'BKG+SKY']:
                    if self.noise_model == 'BKG+RN':
                        sigma = float(np.mean(chip.read_noise))
                    elif self.noise_model == 'BKG+SKY':
                        sigma = chip.sky_sigma
                        self.bal_config[i]['image']['noise'] = {
                            'type' : 'Gaussian',
                            'sigma' : sigma
                        }
                if 'BKG' in self.noise_model:
                    # Use chip background file as initial image instead
                    self.bal_config[i]['image'].update({'initial_image' : chip.bkg_file})
            # Can add more noise models here!
            # elif ...

            # These fields need nothing besides dict initialization
            self.bal_config[i]['input'] = {}
            self.bal_config[i]['gal'] = {}
            self.bal_config[i]['stamp'] = {}

            # Setup the PSF
            psf_file = chip.psf_filename
            if psf_file is not None:
                self.bal_config[i]['input'] = {
                    'des_psfex' : {
                        'file_name' : psf_file,
                        'image_file_name' : chip_file,
                    }
                }

            # Setup 'output' field
            out_file = io.return_output_fname(config.output_dir,
                                                'balrog_images',
                                                str(self.curr_real),
                                                config.data_version,
                                                self.tile_name,
                                                chip.band,
                                                chip.name)
            self.bal_config[i]['output'] = {'file_name' : out_file}

            # NOTE: Some input types require the field `bands` to be set, and so will fail
            # if we do not set it for this chip injection (even though it is never used
            # for zero injections)
            for inj, ninj in chip.nobjects.items():
                if (ninj == 0) and (inj_cat.needs_band is True):
                    try:
                        # Need the inverse map of {input_type : inj_type}
                        inpt = {v:k for k, v in self.inj_types.iteritems()}[inj]
                        self.bal_config[i]['input'].update({inpt : {'bands' : chip.band}})
                    except KeyError: pass

            # If multiple input types, add list setup
            if len(self.input_types) > 1:
                list_structure_gal = deepcopy(self.bal_config[0]['gal'])
                list_structure_im = deepcopy(self.bal_config[0]['gal'])
                self.bal_config[i]['gal'].update(list_structure_gal)
                self.bal_config[i]['image'].update({'image_pos' : list_structure_im})

                # Build index array to choose between inputs in the 'gal' & 'image' fields
                icount = 0
                lower = ''
                # The following enumeration will be in the correct order as it is
                # a list of dicts rather than a nested dict
                # NOTE: The final 'else' will result in an index 1 larger than the list,
                # so it will automatically error if something goes wrong
                for i, inpt in enumerate(x['type'] for x in list_structure_gal['items']):
                    if icount == 0:
                        upper = lower + '@image.N_{}'.format(inpt)
                        indx_eval = '{} if obj_num<{} else {}'.format(icount, upper, icount+1)
                    else:
                        upper = '{}'.format(lower) + '+' + '@image.N_{}'.format(inpt)
                        indx_eval += ' if obj_num>={} and obj_num<{} else {}'.format(lower,
                                                                                     upper,
                                                                                     icount+1)
                    icount += 1
                    lower = upper

                self.bal_config[i]['gal'].update({
                    'index' : {
                        'type':'Eval',
                        'str' : indx_eval }
                    })
                self.bal_config[i]['image']['image_pos'].update({
                    'index' : {
                        'type':'Eval',
                        'str' : indx_eval }
                    })

                # Clean up entries
                for j in range(Ninput): self.bal_config[i]['image']['image_pos']['items'][j] = {}

        #-----------------------------------------------------------------------------------------------
        # Set common entries independent of list structure

        self.bal_config[i]['image'].update({'N_{}'.format(inj_type) : chip.nobjects[inj_type]})
        # Any extra fields to be set for a given input are handled here
        inj_cat.setup_chip_config(config, self.bal_config, chip, i)

        #-----------------------------------------------------------------------------------------------
        # If only one input type, don't use list structure
        # TODO: Some of this should be streamlined between 1 vs. multi list structure
        if Ninput == 1:

            # Set object number and positions
            x, y = inj_pos_im[:,0].tolist(), inj_pos_im[:,1].tolist()
            nobjs = len(inj_pos_im)
            self.bal_config[i]['image'].update({
                'image_pos' : {
                    'type' : 'XY',
                    'x' : { 'type' : 'List', 'items' : x },
                    'y' : { 'type' : 'List', 'items' : y }
                }
            })

            # Set the injected objects' catalog indices and flux factor
            indices = inj_indx.tolist()
            ff = float(chip.flux_factor)
            self.bal_config[i]['gal'].update({
                'scale_flux' : ff,
                'index' : {
                    'type' : 'List',
                    'items' : indices
                }
            })

            # Any extra fields to be set for a given input are added here
            inj_cat.build_single_chip_config(config, self.bal_config, chip, i)

        #-----------------------------------------------------------------------------------------------
        # If multiple input types, use list structure

        else:
            # Get object index
            indx = self.input_indx[input_type]

            # Make sure injection type indices are consistent
            assert self.bal_config[0]['gal']['items'][indx]['type'] == inj_type

            # Set object indices and flux factor
            indices = inj_indx.tolist()
            ff = float(chip.flux_factor)
            self.bal_config[i]['gal']['items'][indx].update({
                'scale_flux' : ff,
                'index' : {
                    'type' : 'List',
                    'items' : indices
                }
            })

            # Set object positions
            x, y = inj_pos_im[:,0].tolist(), inj_pos_im[:,1].tolist()
            self.bal_config[i]['image']['image_pos']['items'][indx].update({
                'type' : 'XY',
                'x' : { 'type' : 'List', 'items' : x },
                'y' : { 'type' : 'List', 'items' : y }
            })

            # Any extra fields to be set for a given input can be added here.
            inj_cat.build_multi_chip_config(config, self.bal_config, chip, i, indx)

        chip.types_injected += 1

        # A few final checks...
        if Ninput == chip.types_injected:
            # NOTE: If all injections have been completed but 'nobjects' is still 0, then
            # add a dummy injection if 'inj_objs_only' is true. This is to ensure that the
            # background is correctly set in `injector.py` during the GalSim call, even
            # for chips with no injections.
            if chip.total_n_objects == 0:
                # This edge case should have only happened for 'inj_objs_only'
                assert config.inj_objs_only['value'] is True

                # Use current inj object for dummy
                chip.set_nobjects(1, inj_type)
                self.bal_config[i]['image'].update({'N_{}'.format(inj_type) : 1})

                if Ninput == 1:
                    self.bal_config[i]['image']['image_pos'].update({
                        'type' : 'XY',
                        'x' : { 'type' : 'List', 'items' : [0] },
                        'y' : { 'type' : 'List', 'items' : [0] }
                    })
                    self.bal_config[i]['gal'] = {
                        'type' : 'Gaussian',
                        'sigma' : 1,
                        'flux' : 0.0 # GalSim will run, but no effective image added
                    }

                else:
                    # Use list format
                    k = self.input_indx['gals']
                    self.bal_config[i]['image']['image_pos']['items'][k].update({
                        'type' : 'XY',
                        'x' : { 'type' : 'List', 'items' : [0] },
                        'y' : { 'type' : 'List', 'items' : [0] }
                    })
                    self.bal_config[i]['gal']['items'][k] = {
                        'type' : 'Gaussian',
                        'sigma' : 1,
                        'flux' : 0.0 # GalSim will run, but no effective image added
                    }

        if self.has_injections is False:
            self.has_injections = True

        return

    def add_bal_config_entry(self):
        self.bal_config.append({})
        self.bal_config_len += 1

    def run_galsim(self, vb=0):
        '''
        Run full GalSim executable of the modified gs_config. Will inject all Balrog
        galaxies contained in a given tile for a given realization in all chips in
        all bands.
        '''

        # If there are no galaxies to inject, simply save chip file to balrog directory
        if self.has_injections is False:
            # NOTE: Can't run galsim executable (besides being a waste of time) as
            # the base layer yaml config won't be valid without extra additions.
            return

        # A new GalSim config file for Balrog injections has been created and all simulations
        # can now be run simultaneously using all GalSim machinery
        # bashCommand = 'galsim {} -v 2 -l gs_logfile'.format(self.bal_config_file)
        bashCommand = 'galsim {} -v {}'.format(self.bal_config_file, vb)

        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

        if vb>0:
            for line in iter(process.stdout.readline, ''): print(line.replace('\n', ''))

        # Needed to get the return code from GalSim
        streamdata = process.communicate()[0]
        rc = process.returncode

        return rc

    def write_truth_catalog(self, config):
        '''
        Writes a fits file that contains the subset of the input catalog that has been
        injected into the current tile, as well as a few extra columns.
        '''

        outfiles = {}
        truth = {}

        real = self.curr_real
        base_outfile = os.path.join(config.output_dir, 'balrog_images', str(real),
                                    config.data_version, self.tile_name)

        for input_type, inj_type in config.inj_types.items():
            inpt = config.input_types[input_type]
            inj = self.inj_cats[input_type]
            outfiles[inj_type] = inj.get_truth_outfile(base_outfile, real)
            if inpt.parametric_cat is None:
                truth[inj_type] = inpt.cat[inj.indx[real]]
            else:
                # Parametric truth catalog previously loaded in `load_input_catalogs()`
                truth[inj_type] = inpt.parametric_cat[inj.indx[real]]

            self._write_new_positions(truth, inj)

        for inj_type, outfile in outfiles.items():
            try:
                with fitsio.FITS(outfile, 'rw', clobber=True) as truth_table:

                    truth_table.write(truth[inj_type])

                    # Fill primary HDU with simulation metadata
                    # hdr = fits.Header()
                    hdr = {}
                    hdr['run_name'] = config.run_name
                    hdr['config_file'] = config.args.config_file
                    hdr['geom_file'] = config.geom_file
                    hdr['tile_list'] = config.args.tile_list
                    hdr['config_dir'] = config.config_dir
                    hdr['tile_dir'] = config.tile_dir
                    hdr['output_dir'] = config.output_dir
                    hdr['psf_dir'] = config.psf_dir
                    hdr['inj_time'] = str(datetime.datetime.now())
                    hdr['inj_bands'] = config.bands
                    hdr['n_objects'] = str(config.n_objects)
                    hdr['object_density'] = str(config.object_density)
                    hdr['pos_sampling'] = str(config.pos_sampling)
                    hdr['realizations'] = config.realizations
                    hdr['curr_real'] = self.curr_real
                    hdr['data_version'] = config.data_version
                    hdr['inj_type'] = inj_type

                    truth_table[0].write_keys(hdr)

            except IOError:
                # Directory structure will not exist if galsim call failed
                print('Warning: Injection for tile {}, realization {} failed! '
                    'Skipping truth-table writing.'.format(self.tile_name, self.curr_real))
                return

        return

    def _write_new_positions(self, truth_cat, inj_cat):
        pos = inj_cat.pos[self.curr_real]

        # Position re-writing (including default behaviour) has been moved to class
        # methods in `balobject.py`
        inj_cat.write_new_positions(truth_cat, self.curr_real)

        return

    def copy_extensions(self, config):
        '''
        Copy the remaining initial chip image extensions (e.g. weight and mask maps)
        to the new balrog injected images.
        '''

        for band in self.bands:
            out_band_dir = os.path.join(self.output_dir, 'balrog_images', str(self.curr_real),
                                        config.data_version, self.tile_name, band)
            chips = self.chips[band]
            for chip in chips:
                orig_fits = chip.filename

                bal_fits = io.return_output_fname(config.output_dir,
                                                    'balrog_images',
                                                    str(self.curr_real),
                                                    config.data_version,
                                                    self.tile_name,
                                                    chip.band,
                                                    chip.name)

                # As we want to rewrite orig_fits with bal_fits, we pass the same name for combined_fits
                combined_fits = bal_fits

                # Combine balrog image with extra extentions
                io.combine_fits_extensions(combined_fits, bal_fits, orig_fits, config=config)

        return

#-------------------------
# Related Tile functions

def create_tiles(config):
    '''
    Create list of `Tile` objects given input args and configuration file.
    '''

    tile_list = load_tile_list(config.tile_list, vb=config.vb)

    # Will keep a list of desired tiles
    tiles = []

    for tile_name in tile_list:
        tiles.append(Tile(tile_name, config))

    return tiles

def load_tile_list(tile_list_file, vb=0):

    # TODO: Make this check more robust...
    if tile_list_file.lower().endswith('.csv'):
        tile_list = io.open_csv_list(tile_list_file)
    elif tile_list_file.lower().endswith('.txt'):
        tile_list = io.open_txt_list(tile_list_file)
    else:
        raise Exception('`tile_list` must be in a `.csv` or `.txt` file!')

    if vb > 0:
        print('Loaded {} tiles...'.format(len(tile_list)))

    return tile_list

