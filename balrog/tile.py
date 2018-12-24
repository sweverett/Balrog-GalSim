import fileio as io
import grid
import galsim
import os, sys, errno
import ntpath
import fitsio
import yaml
import numpy as np
from astropy import wcs
import warnings
import subprocess
import datetime

# Balrog files
from chip import Chip
import fileio as io
import mathutil as util
import grid

# import pudb

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
        self.gals_pos = None # injected galaxy position array
        self.gals_indx = None # injected galaxy index array (from input catalog)
        self.stars_pos = None # injected star position array
        self.stars_indx = None # injected star index array (from input catalog)

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

        # NOTE: `n_galaxies` and `gal_density` are now defined to be *per realization*
        # Set the number of galaxies injected per realization
        if config.n_galaxies:
            # Then fixed number of galaxies regardless of tile size
            self.gals_per_real = config.n_galaxies
        else:
            # Then gal_density was set; galaxy number depends on area
            self.gals_per_real = round(self.u_area * config.gal_density)

        # Set tile directory structure
        self.dir = os.path.abspath(os.path.join(config.tile_dir, self.tile_name))
        self._set_bands(config)

        # Set initial injection realization number
        self.set_realization(realization)

        # Setup new Balrog config file for chip to be called by GalSim executable
        self._setup_bal_config(config)

        # Load zeropoint list from file
        self._load_zeropoints(config)

        # Load background images, if needed
        self._load_backgrounds(config)

        # Set noise properties
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

        # pudb.set_trace()

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
        self.bindx = dict(zip(self.bands, range(len(self.bands))))

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

        # pudb.set_trace()

        return

    def _load_backgrounds(self, config, s_begin=0, s_end=4):
        '''
        Load any needed background images.
        # NOTE: For now, these are only used for grid test images.
        '''

        # TODO: Change this to always look for 'BKG' or 'BKG+' inside of noise model!
        if config.inj_objs_only['noise'] in ['BKG', 'BKG+CCD', 'BKG+RN', 'BKG+SKY']:
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
			# TODO: Only here for blank testing!
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

        # pudb.set_trace()
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
                # pudb.set_trace()
                self.chips[band].append(Chip(filename, band, config, tile_name=self.tile_name,
                                        zeropoint=zp, tile=self))

        # pudb.set_trace()

        return


    def is_chip_image(self, config, chip_file):
        '''
        Checks if passed file is an appropriate chip injection image given data version.
        '''

        if config.data_version == 'y3v02':
            # TODO: Come up with a more robust check!
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

        if config.sim_stars is True:
            self.generate_stars(config, realization)

        if config.sim_gals is True:
            self.generate_galaxies(config, realization)

        # TODO: Can add new simulation types later. Transients, maybe?

        return

    def generate_stars(self, config, realization):
        '''
        For now (Y3), the star catalogs (including positions) are pre-computed. So we just
        need to declare some variables for future use.
        '''

        # pudb.set_trace()

        if config.input_types['stars'] == 'des_star_catalog':
            input_type = 'des_star_catalog'
            # Set tile-wide star injection parameters
            self.bal_config[0]['input'][input_type].update({'tile' : self.tile_name})

            if config.data_version == 'y3v02':
                # The first time, divide up catalog randomly between realizations
                if realization == config.realizations[0]:
                    # Can't guarantee star count consistency, so use dicts
                    # (moved from lists as realizations are no longer guaranteed
                    # to be contiguous)
                    self.stars_indx = {}
                    self.stars_pos = {}
                    self.Nstars = {}

                    self.star_model = config.star_model

                    # Set up proxy catalog
                    gs_config = copy.deepcopy(config.orig_gs_config[0])
                    gs_config['input'][input_type]['bands'] = 'g' # Band doesn't matter
                    gs_config['input'][input_type]['tile'] = self.tile_name

                    # Don't need galaxy info, so remove for speed
                    try: del gs_config['input'][config.input_types['gals']]
                    except KeyError: pass

                    # Make proxy catalog
                    galsim.config.ProcessInput(gs_config)
                    cat_proxy = gs_config['input_objs'][input_type][0]

                    # pudb.set_trace()
                    # If all stars in Sahar's catalogs were guaranteed to be in the
                    # unique tile region, then we would simply do:
                    #
                    # Ns = cat_proxy.getNObjects()
                    #
                    # However, her catalogs are structured such that they contain
                    # stars in the unique region of other tiles. So we only grab the
                    # ones inside the unique region.
                    indices = cat_proxy.indices_in_region([self.ramin, self.ramax],
                                                          [self.decmin, self.decmax],
                                                          boundary_cross=self.ra_boundary_cross)
                    Ns = len(indices)

                    # Update star catalog
                    config.input_cats[input_type] = cat_proxy.getCatalog()
                    config.input_nobjects[input_type] = Ns

                    # Randomize star catalog order and split into approximately equal parts
                    # NOTE: If n_realizations > len(realizations), then DO NOT randomly
                    # shuffle stars, as they are being injected across multiple jobs.
                    Nr = config.n_realizations
                    if Nr == len(config.realizations):
                        rand.shuffle(indices)
                    # pudb.set_trace()
                    indices = [np.array(indices[i::Nr]) for i in range(Nr)]

                    # Grab star positions
                    ra = config.input_cats[input_type]['RA_new']
                    dec = config.input_cats[input_type]['DEC_new']
                    assert len(ra)==len(dec)

                    # Sahar's DES Y3 star catalogs are all pre-computed, so we can set
                    # needed values for all realizations now.
                    # pudb.set_trace()
                    for real in config.realizations:
                        j = int(np.where(real==np.array(config.realizations))[0])
                        inds = indices[j]
                        r, d = ra[inds], dec[inds]

                        self.stars_indx[real] = inds
                        self.stars_pos[real] = np.column_stack((r, d))
                        self.Nstars[real] = len(inds)

        return

    def generate_galaxies(self, config, realization):
        '''
        Generate a list of galaxy positions and indices for tile, for a given realization
        (starts counting at 0). While the generation could in principle depend on the
        realization number, most methods will do all of the generation in the in the initial
        one.
        '''

        input_type = config.input_types['gals']
        if input_type in ['ngmix_catalog', 'cosmos_chromatic_catalog', 'meds_catalog']:
            # input_type = 'ngmix_catalog'
            gal_type = config.input_types['gals']
            if config.data_version == 'y3v02':
                # Generate galaxy positions and indices if this is the first realization

                if realization == config.realizations[0]:
                    # Can't guarantee galaxy count consistency, so use dicts
                    self.gals_pos = {}
                    self.gals_indx = {}
                    self.Ngals = {}

                    Ng = config.input_nobjects[gal_type]
                    Nr = config.n_realizations

                    for real in config.realizations:
                        ngals = self.gals_per_real
                        self.Ngals[real] = ngals

                        # pudb.set_trace()
                        # Generate galaxy coordinates
                        ps = config.pos_sampling
                        if ps['type'] == 'uniform':
                            ra = util.sample_uniform_ra(self.ramin, self.ramax, self.gals_per_real,
                                                boundary_cross=self.ra_boundary_cross)
                            dec = util.sample_uniform_dec(self.decmin, self.decmax, self.gals_per_real,
                                                 unit='deg')
                            self.gals_pos[real] = np.column_stack((ra, dec))

                        elif (ps['type']=='RectGrid') or (ps['type']=='HexGrid'):

                            # pudb.set_trace()

                            gs = ps['grid_spacing']

                            # Rotate grid if asked
                            try:
                                r = ps['rotate']
                                if (isinstance(r, str)) and (r.lower() == 'random'):
                                    if ps['type'] == 'RectGrid':
                                        self.grid_rot_angle = np.rand.uniform(0., np.pi/2.)
                                    elif ps['type'] == 'HexGrid':
                                        self.grid_rot_angle = np.rand.uniform(0., np.pi/3.)
                                else:
                                    unit = ps['angle_unit']
                                    if unit == 'deg':
                                        if (r>=0.0) and (r<360.0):
                                            self.grid_rot_angle = float(r)
                                        else:
                                            raise ValueError('Grid rotation of {} '.format(r) +
                                                             'deg is not valid!')
                                    else:
                                        if (r>=0.0) and (r<2*np.pi):
                                            self.grid_rot_angle = float(r)
                                        else:
                                            raise ValueError('Grid rotation of {} '.format(r) +
                                                             'rad is not valid!')
                            except KeyError:
                                self.grid_rot_angle = 0.0

                            # Offset grid if asked
                            try:
                                o = ps['offset']
                                if (isinstance(o, str)) and (o.lower() == 'random'):
                                    self.grid_offset = [np.rand.uniform(-gs/2., gs/2.),
                                                        np.rand.uniform(-gs/2., gs/2.)]
                                else:
                                    if isinstance(o, list):
                                        self.grid_offset = list(o)
                                    else:
                                        raise ValueError('Grid offset of {} '.format(r) +
                                                         'is not an array!')
                            except KeyError:
                                self.grid_offset = [0.0, 0.0]

                            try:
                                self.angle_unit = ps['angle_unit']
                            except KeyError:
                                self.angle_unit = None

                            # Creates the rectangular grid given tile parameters and calculates the
                            # image / world positions for each object
                            if ps['type'] == 'RectGrid':
                                tile_grid = grid.RectGrid(gs, self.wcs, Npix_x=self.Npix_x,
                                                          Npix_y=self.Npix_y,
                                                          pixscale=self.pixel_scale,
                                                          rot_angle = self.grid_rot_angle,
                                                          angle_unit = self.angle_unit,
                                                          pos_offset = self.grid_offset)
                            elif ps['type'] == 'HexGrid':
                                tile_grid = grid.HexGrid(gs, self.wcs, Npix_x=self.Npix_x,
                                                         Npix_y=self.Npix_y,
                                                         pixscale=self.pixel_scale,
                                                         rot_angle = self.grid_rot_angle,
                                                         angle_unit = self.angle_unit,
                                                         pos_offset = self.grid_offset)

                            self.gals_pos[real] = tile_grid.pos

                            # NOTE: We ignore the inputted ngals and use the correct grid value
                            # instead (user was already warned)
                            ngals = np.shape(tile_grid.pos)[0]
                            if ngals != self.gals_per_real:
                                warnings.warn('The passed n_galaxies : {}'.format(self.gals_per_real) +
                                              ' does not match the {} needed'.format(ngals) +
                                              ' for {} with spacing {}.'.format(ps['type'], gs) +
                                              ' Ignoring input n_galaxies (Only for grids!).')
                            self.Ngals[real] = ngals

                        # Generate galaxy indices (in input catalog)
                        # NOTE: Nearly all runs will generate a random sample of indices. However,
                        # for some testing it would be nice to use an identical galaxy for all
                        # injections. In this case, the user can set a single index in the 'gal'
                        # section of the global config
                        # pudb.set_trace()
                        try:
                            orig_indx = config.gs_config[0]['gal']['index']
                            if type(orig_indx) is int:
                                # Need to find original index of catalog
                                gs_config = copy.deepcopy(config.gs_config[0])
                                # Add dummy band index (band doesn't matter)
                                gs_config['input'][input_type].update({'bands':'g'})
                                galsim.config.ProcessInput(gs_config)
                                cat_proxy = gs_config['input_objs'][input_type][0] # Actually a proxy
                                cat = cat_proxy.getCatalog()
                                if input_type == 'ngmix_catalog':
                                    indx = int(np.where(cat['id']==orig_indx)[0])
                                elif input_type == 'meds_catalog':
                                    # ID's consistent between bands
                                    b = cat_proxy.getBands()[0]
                                    indx = int(np.where(cat[b]['id']==orig_indx)[0])
                                indices = indx * np.ones(ngals, dtype='int16')
                                del cat_proxy
                                del cat
                                # Now remove `index` from global config (can cause issues w/
                                # other input types
                                self.bal_config[0]['gal'].pop('index', None)
                            else:
                                raise TypeError('Can only set a global galaxy index in the ' +
                                                'config if it is an integer!')
                        except KeyError:
                            indices = np.random.choice(xrange(Ng), size=ngals)
                            # indices = util.sample_uniform_indx(0, config.input_nobjects[gal_type], Ng)
                        self.gals_indx[real] = indices

        else:
            raise Exception('No `generate_galaxies()` implementation for input of type ' +
                            '{}!'.format(input_type))

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

    def add_gs_injection(self, config, chip, inj_indx, inj_pos_im, inj_type):
        '''
        This function appends the global GalSim config with an additional simulation to
        be done using nullwt chip-specific infomation and Balrog injected galaxy/star positions
        in image coordinates.
        '''

        # TODO: This function has gotten unwieldly as Balrog has gotten more complex. We should
        # restructure this in kind with new input object classes that have their own method for
        # adding a gs injection to the config

        assert len(inj_indx) == len(inj_pos_im)

        # Skip chips with no injections, except for a few special cases
        # pudb.set_trace()
        if (len(inj_indx) == 0) and (config.inj_objs_only['value'] is False):
            if config.vb > 1:
                print('No objects of type {} were passed on injection '.format(inj_type) +
                'to chip {}. Skipping injection.'.format(chip.name))
                return

        if (inj_type != 'gals') and (inj_type != 'stars'):
            raise ValueError('Currently, the only injection types allowed are \'gals\' and \'stars\'.')

        # pudb.set_trace()

        # If this is the first injection for the chip, then set up new entry in bal_config
        if chip.types_injected == 0:
            self.bal_config.append({})
            self.bal_config_len += 1

        # Injection index
        i = self.bal_config_len - 1

        Ninput = len(self.input_types)
        Ninject = len(inj_indx)
        input_type = self.input_types[inj_type]

        #-----------------------------------------------------------------------------------------------
        # Now set up common config entries if chip's first injection

        # TODO: Can move some of the des_star_catalog fields (like tile and model_type) here

        if chip.types_injected == 0:

            # Setup 'image' field
            chip_file = chip.filename
            nobjs = '$@image.Ngals+@image.Nstars'
            self.bal_config[i]['image'] = {
                'Ngals' : 0,
                'Nstars' : 0,
                'initial_image' : chip_file,
                'nobjects' : nobjs,
                'wcs' : { 'file_name' : chip_file },
            }

            # If noise is to be added, do it here
            # pudb.set_trace()
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

            # Setup 'input' field (nothing besides dict init, for now)
            self.bal_config[i]['input'] = {}

            # Setup 'gal' field (nothing besides dict init, for now)
            self.bal_config[i]['gal'] = {}

            # Setup 'stamp' field (nothing besides dict init, for now)
            self.bal_config[i]['stamp'] = {}

            # Setup the 'psf' field
            psf_file = chip.psf_filename
            if psf_file is not None:
                self.bal_config[i]['input'] = {
                    'des_psfex' : {
                        'file_name' : psf_file,
                        'image_file_name' : chip_file,
                    }
                }

            # Setup 'output' field
            out_file = os.path.join(self.output_dir, 'balrog_images', str(self.curr_real),
                                    config.data_version, self.tile_name, chip.band,
                                    '{}_balrog_inj.fits'.format(chip.name))
            self.bal_config[i]['output'] = {'file_name' : out_file}

            # NOTE: Some input types require the field `bands` to be set, and so will fail
            # if we do not set it for this chip injection (even though it is never technically
            # used)
            if chip.Ngals == 0:
                try:
                    gal_type = self.input_types['gals']
                    self.bal_config[i]['input'].update({gal_type : {'bands' : chip.band}})
                except KeyError: pass
            if chip.Nstars == 0:
                try:
                    star_type = self.input_types['stars']
                    self.bal_config[i]['input'].update({star_type : {'bands' : chip.band}})
                except KeyError: pass

            # If multiple input types, add list setup
            if len(self.input_types) > 1:
                # TODO: Find a cleaner way to do this!
                list_structure_gal = copy.deepcopy(self.bal_config[0]['gal'])
                list_structure_im = copy.deepcopy(self.bal_config[0]['gal'])
                self.bal_config[i]['gal'].update(list_structure_gal)
                # TODO: Determine if we need M Jarvis's version of below:
                self.bal_config[i]['gal'].update({
                    'index' : {
                        'type' : 'Eval',
                        'str' : '0 if obj_num % (@image.Ngals + @image.Nstars) < @image.Ngals else 1'}
                        # 'str' : '0 if obj_num < @image.Ngals else 1'}
                    })
                    # 'index' : '$0 if obj_num % (@image.Ngals+@image.Nstars) < @image.Ngals else 1'})
                self.bal_config[i]['image'].update({'image_pos' : list_structure_im})
                self.bal_config[i]['image']['image_pos'].update({
                    'index' : {
                        'type' : 'Eval',
                        'str' : '0 if obj_num % (@image.Ngals + @image.Nstars) < @image.Ngals else 1'}
                        # 'str' : '$0 if obj_num < @image.Ngals else 1'}
                    })
                    # 'index' : '$0 if obj_num % (@image.Ngals+@image.Nstars) < @image.Ngals else 1'})
                # Clean up entries
                for j in range(Ninput): self.bal_config[i]['image']['image_pos']['items'][j] = {}

        #-----------------------------------------------------------------------------------------------
        # Set common entries independent of list structure

        if inj_type == 'gals':
            self.bal_config[i]['image'].update({'Ngals' : chip.Ngals})

        elif inj_type == 'stars':
            self.bal_config[i]['image'].update({'Nstars' : chip.Nstars})

        else: raise ValueError('For now, only `gals` or `stars` are valid injection types!')

        # NOTE: Any extra fields to be set for a given input can be added here.
        if (inj_type=='gals') and (input_type=='cosmos_chromatic_catalog'):
            # Set the bandpass
            self.bal_config[i]['stamp'].update({
                'type' : 'COSMOSChromatic',
                'bandpass' : config.filters[chip.band].band_config
            })

        # pudb.set_trace()
        if input_type in ['ngmix_catalog', 'meds_catalog', 'des_star_catalog']:
            # Only load into memory the needed band catalog information
            self.bal_config[i]['input'].update({
                input_type : {'bands' : chip.band}
            })

        #-----------------------------------------------------------------------------------------------
        # If only one input type, don't use list structure
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

            # NOTE: Any extra fields to be set for a given input can be added here.

            # pudb.set_trace()
            if input_type == 'meds_catalog':
                # Only use meds/psf files for needed band
                b = self.bindx[chip.band]
                meds_file = [self.bal_config[0]['input'][input_type]['meds_files'][b]]
                psf_file = [self.bal_config[0]['input'][input_type]['psf_files'][b]]
                self.bal_config[i]['input'][input_type].update({
                    'meds_files' : meds_file,
                    'psf_files' : psf_file
                })

                # TODO: We should switch all other input gal types to do the same.
                # Set the injection band
                self.bal_config[i]['gal'].update({
                    'band' : chip.band
                })

        #-----------------------------------------------------------------------------------------------
        # If multiple input types, use list structure

        else:
            # Get object index
            indx = self.input_indx[inj_type]

            # Make sure injection type indices are consistent
            assert self.bal_config[0]['gal']['items'][indx]['type'] == self.inj_types[inj_type]

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

            # NOTE: Any extra fields to be set for a given input can be added here.

            if input_type == 'meds_catalog':
                # Only use meds/psf files for needed band
                b = self.bindx[chip.band]
                meds_file = [self.bal_config[0]['input']['items'][indx]['meds_files'][b]]
                psf_file = [self.bal_config[0]['input']['items'][indx]['psf_files'][b]]
                self.bal_config[i]['input']['items'][indx].update({
                    'meds_files' : meds_file,
                    'psf_files' : psf_file
                })

                # TODO: We should switch all other input gal types to do the same.
                # Set the injection band
                self.bal_config[i]['gal']['items'][indx].update({
                    'band' : chip.band
                })

        chip.types_injected += 1

        # A few final checks...
        if Ninput == chip.types_injected:
            # NOTE: If all injections have been completed but 'nobjects' is still 0, then
            # add a dummy injection if 'inj_objs_only' is true. This is to ensure that the
            # background is correctly set in `injector.py` during the GalSim call, even
            # for chips with no injections.
            if chip.Ngals + chip.Nstars == 0:
                # This edge case should have only happened for 'inj_objs_only'
                assert config.inj_objs_only['value'] is True

                # pudb.set_trace()

                chip.set_Ngals(1)
                self.bal_config[i]['image'].update({'Ngals' : 1})

                if Ninput == 1:
                    self.bal_config[i]['image']['image_pos'].update({
                        'type' : 'XY',
                        'x' : { 'type' : 'List', 'items' : [0] },
                        'y' : { 'type' : 'List', 'items' : [0] }
                    })
                    self.bal_config[i]['gal'] = {
                        'type' : 'Gaussian',
                        'sigma' : 2,
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
                        'sigma' : 2,
                        'flux' : 0.0 # GalSim will run, but no effective image added
                    }

        if self.has_injections is False:
            self.has_injections = True

        # pudb.set_trace()

        return

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
                                    config.data_version, self.tile_name,
                                    '{}_{}_balrog_truth_cat'.format(self.tile_name, real))

        if config.sim_gals is True:
            itype = config.input_types['gals']
            outfiles['gals'] = base_outfile + '_gals.fits'
            # pudb.set_trace()
            # Now update ra/dec positions for truth catalog
            if itype == 'ngmix_catalog':
                truth['gals'] = config.input_cats[itype][self.gals_indx[real]]
            elif itype == 'meds_catalog':
                # Parametric truth catalog previously loaded in `load_input_catalogs()`
                truth['gals'] = config.meds_param_catalog[self.gals_indx[real]]
            # elif ...

            # Tries multiple column keywords
            self._write_new_positions(truth, 'gals', itype)

        if config.sim_stars is True:
            itype = config.input_types['stars']
            outfiles['stars'] = base_outfile + '_stars.fits'
            truth['stars'] = config.input_cats[config.input_types['stars']][self.stars_indx[real]]
            # Now update ra/dec positions for truth catalog
            if config.input_types['stars'] == 'des_star_catalog':
                truth['stars']['RA_new'] = self.stars_pos[self.curr_real][:,0]
                truth['stars']['DEC_new'] = self.stars_pos[self.curr_real][:,1]
            else:
                # Tries multiple column keywords
                self._write_new_positions(truth, 'stars', itype)

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
                    if config.n_galaxies:
                        hdr['n_galaxies'] = config.n_galaxies
                    if config.gal_density:
                        hdr['gal_density'] = config.gal_density
                    hdr['realizations'] = config.realizations
                    hdr['curr_real'] = self.curr_real
                    hdr['data_version'] = config.data_version

                    key = 'input_type_{}'.format(inj_type)
                    hdr[key] = config.input_types[inj_type]

                    truth_table[0].write_keys(hdr)

            except IOError:
                # Directory structure will not exist if galsim call failed
                print('Warning: Injection for tile {}, realization {} failed! '
                    'Skipping truth-table writing.'.format(self.tile_name, self.curr_real))
                return

        return

    def _write_new_positions(self, truth_cat, inj_type, in_type):
        # pudb.set_trace()
        if inj_type == 'gals': pos = self.gals_pos
        elif inj_type == 'stars': pos = self.stars_pos

        try:
            truth_cat[inj_type]['ra'] = pos[self.curr_real][:,0]
            truth_cat[inj_type]['dec'] = pos[self.curr_real][:,1]
        except KeyError:
            try:
                truth_cat[in_type]['RA'] = pos[self.curr_real][:,0]
                truth_cat[in_type]['DEC'] = pos[self.curr_real][:,1]
            except KeyError as e:
                raise('Tried to write truth positions using column names of ra/dec; RA/DEC.',
                        'either rename position columns or add an entry for {}'.format(in_type),
                        'in `write_truth_catalog()`\n',e)

        return

    def copy_extensions(self, config):
        '''
        Copy the remaining initial chip image extensions (e.g. weight and mask maps)
        to the new balrog injected images.
        '''

        # pudb.set_trace()
        for band in self.bands:
            out_band_dir = os.path.join(self.output_dir, 'balrog_images', str(self.curr_real),
                                        config.data_version, self.tile_name, band)
            chips = self.chips[band]
            for chip in chips:
                orig_fits = chip.filename

                # Grabbed from `add_gs_injection()`
                bal_fits = os.path.join(self.output_dir, 'balrog_images', str(self.curr_real),
                                        config.data_version, self.tile_name, chip.band,
                                        '{}_balrog_inj.fits'.format(chip.name))

                # As we want to rewrite orig_fits with bal_fits, we pass the same name for combined_fits
                combined_fits = bal_fits

                # Combine balrog image with extra extentions
                combine_fits_extensions(combined_fits, bal_fits, orig_fits, config=config)

        return

#-------------------------
# Related Tile functions

def create_tiles(config):
    '''
    Create list of `Tile` objects given input args and configuration file.
    '''

    tile_list = load_tile_list(config.tile_list, config)

    # pudb.set_trace()

    # Will keep a list of desired tiles
    tiles = []

    for tile_name in tile_list:
        tiles.append(Tile(tile_name, config))

    return tiles

def load_tile_list(tile_list_file, config):

    # TODO: Make this check more robust...
    if tile_list_file.lower().endswith('.csv'):
        tile_list = io.open_csv_list(tile_list_file)
    elif tile_list_file.lower().endswith('.txt'):
        tile_list = io.open_txt_list(tile_list_file)
    # elif ... # TODO: add more types!
    else:
        raise Exception('`tile_list` must be in a `.csv` or `.txt` file!')

    if config.vb > 0:
        print('Loaded {} tiles...'.format(len(tile_list)))

    return tile_list

def combine_fits_extensions(combined_file, bal_file, orig_file, config=None):
    '''
    Combine other needed image extensions (e.g. weight and mask maps).
    Modified from Nikolay's original version in BalrogPipeline.py.
    '''

    # Read in simulated image and pre-injection extensions
    sciIm = fitsio.read(bal_file, ext=0)
    sciHdr = fitsio.read_header(orig_file, ext=0)
    wgtIm, wgtHdr = fitsio.read(orig_file, ext=1, header=True)
    wgt_me_Im, wgt_me_Hdr = fitsio.read(orig_file, ext=2, header=True)
    mskIm, mskHdr = fitsio.read(orig_file, ext=3, header=True)

    # If injecting on blank images, then set these to some sensible values
    # pudb.set_trace()
    if config:
        if config.inj_objs_only['value'] is True:
            # TODO: This needs to be generalized for different noise models!
	    if config.inj_objs_only['noise'] == 'BKG+SKY':
		noise = sciHdr['SKYSIGMA']
	    else:
		noise = np.mean([sciHdr['RDNOISEA'], sciHdr['RDNOISEB']])
            inv_sky_var = 1.0 / (noise**2)
            wgtIm.fill(inv_sky_var)
            wgt_me_Im.fill(inv_sky_var)
            mskIm.fill(0)

    # Write all 3 extensions in the output file
    # NOTE: Usually combined_file will be the same name as bal_file,
    # so we delete the original injection and replace it with the new extensions
    if os.path.exists(combined_file):
        os.remove(combined_file)

    fits = fitsio.FITS(combined_file,'rw')
    fits.write(sciIm, header=sciHdr)
    fits.write(wgtIm, header=wgtHdr)
    fits.write(wgt_me_Im, header=wgtHdr)
    fits.write(mskIm, header=mskHdr)

    # Put back EXTNAME into headers
    fits[0].write_key('EXTNAME', 'SCI', comment="Extension name")
    fits[1].write_key('EXTNAME', 'WGT', comment="Extension name")
    fits[2].write_key('EXTNAME', 'WGT_ME', comment="Extension name")
    fits[3].write_key('EXTNAME', 'MSK', comment="Extension name")

    return

