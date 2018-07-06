#####################################################################
#
# I should probably add a description for balrog!
#
#
# Spencer Everett
# UCSC
# 12/10/2017
#####################################################################

import numpy as np
import random as rand
import os, sys, errno
import cPickle as pickle
import warnings
import subprocess
import shutil
import ntpath
import copy
import time
import galsim
import csv
import yaml
import argparse
import datetime
import itertools
import fitsio
from collections import OrderedDict
from astropy.io import fits
from astropy import wcs
from astropy.table import Table

# Balrog files
import injector
import grid
import filters

# Use for debugging
# import pudb

#-------------------------------------------------------------------------------
# Important todo's:
# TODO: Add a cut on TdByTe!
# TODO: Implement error handling for galaxy injections / gsparams! (Working solution, but need
#       to look into more detail per Erin)
# TODO: Clean up evals in add_gs_injection()!
# TODO: Figure out injector.py parameter parsing issue!
# TODO: Check pixel origin for transformations!
# TODO: CHECK COSMOS ZEROPOINTS!

# Some extra todo's:
# TODO: Have code automatically check fitsio vs astropy
# TODO: More general geometry file inputs
# TODO: Add check for python path!
# TODO: Figure out permission issue!
# TODO: Make a log system!
# TODO: Add a single seed for each tile w/r/t noise realizations
# TODO: Multi-line print statements are a bit bugged

# Questions:
# None curently

#-------------------------------------------------------------------------------
# Define currently allowed and types for various objects. For `allowed`, inputs
# MUST be one of the types. For `supported`, no guarantees are made if input is
# not one of the types.

# QUESTION: Will we ever use `u`?
_allowed_bands = 'grizy'

# TODO: Allow Piff when available!
_supported_psf_types = ['DES_PSFEx']#, 'Piff'}
_psf_extensions = {'DES_PSFEx' : 'psfexcat.psf'}#, 'Piff' : 'something.piff'}

# TODO: Incorporate postage stamps!
_supported_input_types = ['ngmix_catalog', 'des_star_catalog', 'cosmos_chromatic_catalog']#, 'postage_stamps'}
_supported_gal_types = ['ngmix_catalog', 'cosmos_chromatic_catalog']
_supported_star_types = ['des_star_catalog']

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

        # TESTING: Can remove in future
        config.flux_factors[self.tile_name] = {}

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
            # self.u_area = deg2arcmin(((self.ramax+360.0) - self.ramin)) * deg2arcmin(self.decmax - self.decmin)
        else: self.ra_flag = False
            # self.u_area = deg2arcmin((self.ramax - self.ramin)) * deg2arcmin(self.decmax - self.decmin)

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

        # pudb.set_trace()

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
                        chip_name = '_'.join(file_name.split('_')[s_begin+1:s_end+1])
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
                    gs_config['input'][input_type]['bands'] = 'griz'
                    gs_config['input'][input_type]['tile'] = self.tile_name

                    # Don't need galaxy info, so remove for speed
                    del gs_config['input'][config.input_types['gals']]

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
        if input_type == 'ngmix_catalog' or input_type == 'cosmos_chromatic_catalog':
            input_type = 'ngmix_catalog'
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

                        # Generate galaxy coordinates
                        ps = config.pos_sampling
                        if ps['type'] == 'uniform':
                            ra = sample_uniform_ra(self.ramin, self.ramax, self.gals_per_real,
                                                boundary_cross=self.ra_boundary_cross)
                            dec = sample_uniform_dec(self.decmin, self.decmax, self.gals_per_real,
                                                 unit='deg')
                            self.gals_pos[real] = np.column_stack((ra, dec))

                            # Generate galaxy indices (in input catalog)
                            indices = np.array(rand.sample(xrange(Ng), ngals))
                            # indices = sample_uniform_indx(0, config.input_nobjects[gal_type], Ng)
                            self.gals_indx[real] = indices

                        elif (ps['type']=='RectGrid') or (ps['type']=='HexGrid'):

                            # pudb.set_trace()

                            gs = ps['grid_spacing']

                            # Creates the rectangular grid given tile parameters and calculates the
                            # image / world positions for each object
                            if ps['type'] == 'RectGrid':
                                tile_grid = grid.RectGrid(gs, self.wcs, Npix_x=self.Npix_x,
                                                    Npix_y=self.Npix_y, pixscale=self.pixel_scale)
                            elif ps['type'] == 'HexGrid':
                                # But for the future...
                                tile_grid = grid.HexGrid(gs, self.wcs, Npix_x=self.Npix_x,
                                                    Npix_y=self.Npix_y, pixscale=self.pixel_scale)

                            self.gals_pos[real] = tile_grid.pos

                            # NOTE: We ignore the inputted ngals and use the correct grid value
                            # instead (user was already warned)
                            ngals = np.shape(tile_grid.pos)[0]
                            if ngals != self.gals_per_real:
                                warnings.warn('The passed n_galaxies : {}'.format(self.gals_per_real) +
                                              ' does not match the {} needed'.format(ngals) +
                                              ' for {} with spacing {}.'.format(ps['type'], gs) +
                                              ' Ignoring input n_galaxies.')
                            self.Ngals[real] = ngals

                            # Generate galaxy indices (in input catalog)
                            indices = np.array(rand.sample(xrange(Ng), ngals))
                            self.gals_indx[real] = indices

                        # pudb.set_trace()
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
            if self.noise_model is not None:
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
                if self.noise_model in ['BKG', 'BKG+CCD', 'BKG+RN']:
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
            out_file = os.path.join(self.output_dir, 'balrog_images', str(self.curr_real), \
                                    self.tile_name, chip.band, '{}_balrog_inj.fits'.format(chip.name))
            self.bal_config[i]['output'] = {'file_name' : out_file}

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
            if (inj_type=='gals') and (input_type=='ngmix_catalog'):
                # Set the band for injection
                self.bal_config[i]['input'].update({
                    self.input_types.values()[0] : {'bands' : chip.band}
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

            # pudb.set_trace()

            # NOTE: Any extra fields to be set for a given input can be added here.
            if (inj_type=='gals') and (input_type=='ngmix_catalog'):
                # Set the band for injection
                self.bal_config[i]['input'].update({
                    input_type : {'bands' : chip.band}
                })

        chip.types_injected += 1

        # If all injections have been completed but 'nobjects' is still 0, then add
        # a dummy injection if 'inj_objs_only' is true. This is to ensure that the
        # background is correctly set in `injector.py` during the GalSim call, even
        # for chips with no injections.
        # pudb.set_trace()
        if Ninput == chip.types_injected:
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

        # pudb.set_trace()

        if vb>0:
            # while True:
            #     line = p.stdout.readline()
            #     sys.stdout.write(line)
            # #     if not line: break
            # p = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            for line in iter(process.stdout.readline, ''):
                # from signal import signal, SIGPIPE, SIG_DFL
                # signal(SIGPIPE, SIG_DFL)
                try:
                    # sys.stdout.write(line)
                    # print('CASE 1 !\n')
                    print(line.replace('\n', ''))
                # NOTE: Testing!
                except IOError as e:
                    if e.errno == errno.EPIPE:
                        print('WARNING! PIPE ERROR! (1)')
                        raise
                    else:
                        print('WARNING! PIPE ERROR! (2)')
                        raise
        else:
            # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()

        # TODO: Would be nice to do something with the output / errors in future.
        # maybe implement an additional log? ISSUE: 8

        return

    def write_truth_catalog(self, config):
        '''
        Writes a fits file that contains the subset of the input catalog that has been
        injected into the current tile, as well as a few extra columns.
        '''

        # pudb.set_trace()

        outfiles = {}
        truth = {}


        real = self.curr_real
        base_outfile = os.path.join(config.output_dir, 'balrog_images', str(real),
                       self.tile_name, '{}_{}_balrog_truth_cat'.format(self.tile_name, real))

        if config.sim_gals is True:
            outfiles['gals'] = base_outfile + '_gals.fits'
            truth['gals'] = config.input_cats[config.input_types['gals']][self.gals_indx[real]]
            # Now update ra/dec positions for truth catalog
            if config.input_types['gals'] == 'ngmix_catalog':
                truth['gals']['ra'] = self.gals_pos[self.curr_real][:,0]
                truth['gals']['dec'] = self.gals_pos[self.curr_real][:,1]
        if config.sim_stars is True:
            outfiles['stars'] = base_outfile + '_stars.fits'
            truth['stars'] = config.input_cats[config.input_types['stars']][self.stars_indx[real]]
            # Now update ra/dec positions for truth catalog
            if config.input_types['stars'] == 'des_star_catalog':
                truth['stars']['RA_new'] = self.stars_pos[self.curr_real][:,0]
                truth['stars']['DEC_new'] = self.stars_pos[self.curr_real][:,1]

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

    def copy_extensions(self, config):
        '''
        Copy the remaining initial chip image extensions (e.g. weight and mask maps)
        to the new balrog injected images.
        '''

        # pudb.set_trace()
        for band in self.bands:
            out_band_dir = os.path.join(self.output_dir, 'balrog_images', str(self.curr_real), \
                                            self.tile_name, band)
            chips = self.chips[band]
            for chip in chips:
                orig_fits = chip.filename

                # Grabbed from `add_gs_injection()`
                bal_fits = os.path.join(self.output_dir, 'balrog_images', str(self.curr_real), \
                                        self.tile_name, chip.band, '{}_balrog_inj.fits'.format(chip.name))

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
        tile_list = open_csv_list(tile_list_file)
    elif tile_list_file.lower().endswith('.txt'):
        tile_list = open_txt_list(tile_list_file)
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

#-------------------------------------------------------------------------------
# Chip

class Chip(object):
    '''
    DES chip object.
    # NOTE: While DES chip image data is saved as (RA,DEC), the orientation
    # of the chip image is roated ccw by 90 deg's.
    '''

    def __init__(self, filename, band, config, tile_name=None, zeropoint=30.0, tile=None):

        self.filename = filename
        self.fits_filename = ntpath.basename(filename)
        self.tile_name = tile_name # Name of parent tile, if given
        self.band = band
        self.zeropoint = zeropoint

        # Will be set later
        self.Ngals, self.Nstars = 0, 0

        # Keeps track of how many input types for this chip have been added
        # to bal_config by add_gs_injection()
        self.types_injected = 0

        self._set_name(config)
        self._set_psf(config)
        self._set_wcs()
        self._set_noise(config)
        self._set_flux_factor(config)
        self._set_bkg(config, tile)

        return

    def _set_name(self, config, s_begin=0, s_end=4):
        '''
        Determine and set chip name given which subset of fits filename
        (separated by `_`'s) to select.
        # NOTE: For `inmasked` images, should be 0:4.
        '''

        self.name = '_'.join(self.fits_filename.split('_')[s_begin:s_end])

        return

    def _set_psf(self, config):
        '''
        Set the psf type / configuration for the chip.
        '''
        # Can load in type from global GalSim config file
        try:
            # If present, grab config psf type
            self.psf_type = config.gs_config[0]['psf']['type']

            # Check if PSF type is supported
            if self.psf_type in _supported_psf_types:
                self.psf_extension = _psf_extensions[self.psf_type]
                # pudb.set_trace()
                # NOTE: Due to peculiarities of DES_PSFEx GalSim class, we cannot keep psf
                # dir and filename separate; must be combined for an absolute path. This is
                # due to the psf file and chip file being stored in different directories.
                self.psf_dir = os.path.join(config.tile_dir, self.tile_name, config.psf_dir)
                self.psf_filename = os.path.join(self.psf_dir, self.name + '_' + self.psf_extension)

            else:
                # Some basic GalSim psf types will still work, even if not technically supported
                print('Warning: PSF type input {} is not one of the currently supported Balrog' \
                      + 'types: {}'.format(_supported_psf_types))
                self.psf_extension = None
                self.psf_filename = None

        except TypeError:
            # NOTE: For now, we will allow no psf if desired. Will just warn instead
            # raise Exception('Must set a psf type in global GalSim config! Currently ' \
            #              + 'supported types are {}'.format(supported_psf_types))
            print('Warning: No psf type set in global GalSim config! Currently ' \
                         + 'supported types are {}'.format(_supported_psf_types))
            self.psf_extension = None
            self.psf_filename = None
        return

    def _set_wcs(self):
        '''
        Get corners (chip not perfectly aligned to RA/DEC coordinate system).
        Especially useful for galaxy position rejection.
        NOTE: The nullwt chips are not oriented in the standard way.
        In a typical (RA,DEC) projection space:

        DEC increasing up
        .
        .
        1-----------4
        -           -
        -           -
        2-----------3....RA increasing right

        In the DES nullwt chip orientation:

        RA increasing up
        .
        .
        4-------3
        -       -
        -       -
        -       -
        -       -
        1-------2....DEC decreasing right

        Because of the strange orientation, wcs_world2pix([ra, dec]) will
        return correct image coordinates but will be flipped from what we
        normally expect; i.e. IN_PIX: (x,y) ~= (DEC, RA).
        This will affect how we implement contained_in_chip().
        '''

        hdr = fits.getheader(self.filename)
        # Get chip WCS
        self.wcs = wcs.WCS(hdr)

        self.ramin, self.ramax = hdr['RACMIN'], hdr['RACMAX']
        self.decmin, self.decmax = hdr['DECCMIN'], hdr['DECCMAX']
        rc = [hdr['RAC1'], hdr['RAC2'], hdr['RAC3'], hdr['RAC4']]
        dc = [hdr['DECC1'], hdr['DECC2'], hdr['DECC3'], hdr['DECC4']]
        self.corners = zip(rc,dc)

        # Round to nearest pixel (very slight offset)
        #NOTE: Should always be (2048x4096), but in principle could allow
        # different sizes
        self.corners_im = np.round(self.wcs.wcs_world2pix(self.corners,1))

        # Set naxis_im ranges (not (RA,DEC) ranges due to NOTE above)
        #NOTE: should be able to retrieve from header NAXISi, but just
        # to be sure...
        self.naxis1_range = [np.min(self.corners_im[:,0]), np.max(self.corners_im[:,0])]
        self.naxis2_range = [np.min(self.corners_im[:,1]), np.max(self.corners_im[:,1])]

        # pudb.set_trace()

        return

    def _set_noise(self, config):
        '''
        If desired, grab noise values from the chip header.
        '''

        if config.data_version == 'y3v02':
            hdr = fitsio.read_header(self.filename)
            self.sky_var = [hdr['SKYVARA'], hdr['SKYVARB']]
            self.sky_sigma = hdr['SKYSIGMA']
            self.gain = [hdr['GAINA'], hdr['GAINB']]
            self.read_noise = [hdr['RDNOISEA'], hdr['RDNOISEB']]

        return

    def _set_flux_factor(self, config):
        '''
        Calculate and set the flux factor needed to consistently lay down fluxes from the
        input catalog given different image zeropoints.
        '''
        self.flux_factor = np.power(10.0, 0.4 * (self.zeropoint - config.input_zp))

        # TESTING: Can remove in future
        config.flux_factors[self.tile_name][self.name] = self.flux_factor

        return

    def _set_bkg(self, config, tile):
        '''
        Set chip background file, if needed for grid test.
        '''

        if config.inj_objs_only['noise'] in ['BKG', 'BKG+CCD', 'BKG+RN', 'BKG+SKY']:
            assert tile is not None
            self.bkg_file = tile.bkg_files[self.band][self.name]

        return

    def contained_in_chip(self, pos):
        '''
        For an input vector of (RA,DEC) positions, returns a boolean vector
        of whether each position is contained within the chip image.
        '''

        # pudb.set_trace()

        # Load image bounds
        n1min, n1max = self.naxis1_range[0], self.naxis1_range[1]
        n2min, n2max = self.naxis2_range[0], self.naxis2_range[1]

        # convert positions to chip image coordinates
        pos_im = self.wcs.wcs_world2pix(pos,1)
        n1_im, n2_im = pos_im.T

        # Check if position image coords are within chip corners
        n1_in = (n1_im > n1min) & (n1_im < n1max)
        n2_in = (n2_im > n2min) & (n2_im < n2max)
        # Both indices must be in chip
        in_chip = n1_in & n2_in

        return in_chip, pos_im

    # NOTE: Depreciated
    def set_Ngals(self, Ng):
        self.Ngals = Ng

    def set_Nstars(self, Ns):
        self.Nstars = Ns

    def save_without_injection(self, outfile):
        '''
        If there are no Balrog galaxies to inject in the chip area, then save
        copy of current chip image in the new Balrog image format.
        '''

        # Only want to save the first HDU of nullwt image
        # TODO: Switch to fitsio eventually!
        with fits.open(self.filename) as f:
            hdu0 = f[0]
            try:
                # TODO: This is for old version of astropy!
                #hdu0.writeto(outfile, overwrite=True)
                hdu0.writeto(outfile, clobber=True)
            except (IOError, OSError):
                path = os.path.dirname(outfile)
                # To deal with race condition...
                while True:
                    try:
                        os.makedirs(path)
                        break
                    except OSError as e:
                        if e.errno != os.errno.EEXIST:
                            raise e
                        # Wait a bit before trying again!
                        time.sleep(0.5)

                # Now directory is guaranteed to exist
                # TODO: This is for old version of astropy!
                #hdu0.writeto(outfile, overwrite=True)
                hdu0.writeto(outfile, clobber=True)

        return

#-------------------------------------------------------------------------------
# Galaxy

class Galaxy(object):
    '''
    #TODO: Do we need a galaxy class? (probably not)
    '''

    def __init__(self):
        pass

#-------------------------------------------------------------------------------
# Config

class Config(object):
    '''
    Balrog Simulation configuration object. Contains the GalSim config file as
    well as additional simulation parameters.
    '''

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

        # TESTING: Can remove in future
        self.flux_factors = {}

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
        # Add config directory path, if needed
        # if self.args.config_dir:
        #     self.gs_config_file = self.args.config_dir + self.gs_config_file

        # NOTE: `self.config` will be a list of `OrderedDict`'s with length given by
        # the number of configurations specified in the file. To work with a specific
        # simulated image, say with `galsim.config.BuildImage()`, you must pass a
        # *single* config file (e.g. BuildImage(config[3]) ).
        self.gs_config = galsim.config.process.ReadYaml(self.gs_config_file)

        # pudb.set_trace()

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
            # DEPRECATED: Parse `n_realizations` as only input for older configs
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
                if band not in _allowed_bands:
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
                valid_noise = ['CCD', 'BKG', 'BKG+CCD', 'BKG+RN', 'BKG+SKY', 'None']

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
        valid_pos_sampling = ['uniform', 'RectGrid', 'HexGrid']
        default_gs = 40. # arcsec
        try:
            ps = self.gs_config[0]['image']['pos_sampling']
            if isinstance(ps, basestring):
                # Then the string is the input type
                if ps not in valid_pos_sampling:
                    raise ValueError('{} is not a valid position sampling method. '.format(ps) +
                                    'Currently allowed methods are {}'.format(valid_pos_sampling))

                print('No grid spacing passed; using default of {} arcsecs'.format(default_gs))
                self.pos_sampling = {'type' : ps, 'grid_spacing' : default_gs}
            elif isinstance(ps, dict):
                if 'type' not in ps.keys():
                    raise ValueError('If `pos_sampling` is passed as a dict, then must set a type!')
                if 'grid_spacing' not in ps.keys():
                    print('No grid spacing passed; using default of {} arcsecs'.format(default_gs))
                    ps['grid_spacing'] = default_gs
                keys = ['type', 'grid_spacing']
                for key, val in ps.items():
                    if key not in keys:
                        raise ValueError('{} is not a valid key for `pos_sampling`! '.format(key) +
                                         'You may only pass the keys {}'.format(keys))
                    if key == 'grid_spacing':
                        assert val > 0.

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
                'psf_dir':self.psf_dir, 'output_dir':self.output_dir}

        base = {'tile_list':'image', 'geom_file':'image', 'tile_dir':'image', 'psf_dir':'image',
                'output_dir':'output'}

        req = ['tile_list', 'geom_file']
        opt = ['tile_dir', 'psf_dir', 'output_dir']

        for arg, arg_val in args.items():
            try:
                # Have to handle `output_dir` to be consistent with GalSim
                if arg == 'output_dir':
                    config_val = self.gs_config[0][base[arg]]['dir']
                else:
                    config_val = self.gs_config[0][base[arg]][arg]
                if arg_val and config_val != arg_val:
                    # The following works for both files and directories
                    if not os.path.samefile(arg_val, config_val):
                        raise ValueError('Command-line argument {}={} '.format(arg, arg_val) +
                                        'is inconsistent with config value of {}!'.format(config_val))
                val = config_val

            except KeyError:
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
        Load any relevant info from the input catalog(s) (for now just ngmix and desStars)
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

            if input_type in _supported_gal_types:
                # Check that a galaxy cat type hasn't already been set
                if self.sim_gals is True:
                    raise ValueError('Can\'t set multiple input galaxy catalogs!')
                else:
                    self.sim_gals = True
                    self.input_indx['gals'] = i
                    self.input_types['gals'] = input_type
                    # TODO: Can we grab the injection type from the registered GS catalog?
                    if input_type == 'ngmix_catalog':
                        import ngmix_catalog
                        self.inj_types['gals'] = 'ngmixGalaxy'

                        # As we are outside of the GalSim executable, we need to register
                        # the input type explicitly
                        galsim.config.RegisterInputType('ngmix_catalog', ngmix_catalog.ngmixCatalogLoader(
                            ngmix_catalog.ngmixCatalog, has_nobj=True))

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
                            raise ValueError('No input parsing defined for ngmix catalogs for ' +
                            'data version {}'.format(self.data_version))


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
        galsim.config.ProcessInput(gs_config)

        # Grab needed info from the proxy catalog
        # for i, input_type in enumerate(input_cat_types):
        # for  input_type in enumerate(self.input_types)
        # pudb.set_trace()
        for input_type in self.input_types.values():
            # Only do this for injection types (gals, stars)
            # if i not in self.input_indx.values(): continue
            cat_proxy = gs_config['input_objs'][input_type][0] # Actually a proxy
            self.input_cats[input_type] = cat_proxy.getCatalog()
            self.input_nobjects[input_type] = cat_proxy.getNObjects()

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
            if it not in _supported_input_types:
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

#-------------------------------------------------------------------------------
# Unit conversions and misc functions

def deg2arcmin(val):
    return 60.0 * val

def arcmin2deg(val):
    return val / 60.0

def sample_uniform_ra(r1, r2, N=None, boundary_cross=False):
    '''
    Sample N random RA values from r1 to r2, where r1<r2. If boundary_cross
    is set to True, will adjust inputs to allow sampling over 360/0 boundary.
    '''

    # Adjust r1, r2 if they cross the 360/0 ra boundary
    if boundary_cross is True:
        # NOTE: ramin/ramax are saved the same way as the geometry file catalog,
        # which is to say ramin/ramax are always 'left'/'right' boundary, even if
        # ramax < ramin across boundary.
        # This means that r1, while 'min', is larger and needs to be adjusted down.
        r1_shift = r1 - 360.0
        shifted_dist = sample_uniform(r1_shift, r2, N)
        shifted_dist[shifted_dist < 0.0] += 360.0
        return shifted_dist

    else:
        # NOTE: Equivalent to below DEC procedure for ra = P(r2-r1)+r1 for P in (0,1)
        return sample_uniform(r1, r2, N)

def sample_uniform_dec(d1, d2, N=None, unit='deg'):
    '''
    Sample N random DEC values from d1 to d2, accounting for curvature of sky.
    '''

    if N is None: N=1

    # Convert dec to radians if not already
    if unit is 'deg' or unit is 'degree':
        d1, d2 = np.deg2rad(d1), np.deg2rad(d2)

    # Uniform sampling from 0 to 1
    P = np.random.rand(N)

    # Can't use `sample_uniform()` as dec needs angular weighting
    delta = np.arcsin(P * (np.sin(d2) - np.sin(d1)) +np.sin(d1))

    # Return correct unit
    if unit is 'deg' or unit is 'degree':
        return np.rad2deg(delta)
    elif unit is 'rad' or unit is 'radian':
        return delta
    else:
        raise TypeError('Only deg/degree or rad/radian are allowed as inputs!')

def sample_uniform_indx(n1, n2, N=None):
    'Samples N random indices from n1 to n2'

    return np.random.randint(n1, high=n2, size=N)

def sample_uniform(v1, v2, N=None):
    'Samples N random values from v1 to v2'

    if N is None: N=1

    return np.random.uniform(low=v1, high=v2, size=N)

#-------------------------------------------------------------------------------
# Input parsing and file loading

# TODO: Some of these functions are not yet implemented. In future, it would be nice
# to have some more robust and standard IO method calls that also incorporate
# astropy.io vs fitsio depending on what the user has.

# def return_output_name(config, ):
    # pass

def open_file(filename):
    pass

def open_fits_file(filename):
    pass

def open_csv_list(filename):
    with open(filename) as file:
        # delimiter = ',' by default
        reader = csv.reader(file)
        lst = list(reader)
        # flatten any iterative lists (in case commas aren't used as delimiter)
        lst = list(itertools.chain.from_iterable(lst))
        return lst

def open_txt_list(filename):
    with open(filename) as file:
        return [line.strip() for line in file]

def setup_output_dir(config, tiles):
    out_dir = config.output_dir
    images = os.path.join(out_dir, 'balrog_images')
    configs = os.path.join(out_dir, 'configs')
    dirs = [out_dir, images, configs]

    # Setup parent-level directories
    for d in dirs:
        try:
	    os.makedirs(d)
        except OSError as e:
	    if e.errno != errno.EEXIST:
	        raise

    # Setup tiles, reals, and bands
    for tile in tiles:
	pass

    return

def parse_args():
    '''
    Parse command line arguments. Many of the following can be set in the `bal_config` file as well;
    these are kept here only for convenience with previous versions where the arguments were
    required. Passing a value through both methods is fine, as long as they are not inconsistent with
    one another.
    '''

    parser = argparse.ArgumentParser()
    # Global GalSim config file with options that will be applied to all injection images.
    parser.add_argument('config_file', help='.yaml or .json confg file that specifies the desired GalSim'
                        'simulation and injection.')
    # Optional argument for config file(s) directory location (if not .)
    parser.add_argument('-c', '--config_dir', help='Directory that houses the GalSim and related config files.')
    # Required argument for tile geometry file
    parser.add_argument('-g', '--geom_file', help='Tile geometry file (e.g. `Y3A2_COADDTILE_GEOM`)')
    # Optional argument for tile list
    parser.add_argument('-l', '--tile_list', help='.txt or .csv file containing list of tiles to be processed.')
    # Optional argument for DES tiles directory location
    parser.add_argument('-t', '--tile_dir', help='Directory of DES tiles containing single exposure images.')
    # Optional argument for a tile's psf directory name (if not contained in tile_dir/{TILE}/psfs)
    parser.add_argument('-p', '--psf_dir', help='Directory that houses individual psf files for DES chips.')
    # Optional argument for output directory (if not .)
    parser.add_argument('-o', '--output_dir', help='Directory that houses output Balrog images.')
    # Optional argument for verbose messages
    parser.add_argument('-v', '--verbose', action='store', nargs='?', default='0', const='1',
                        help='Turn on verbose mode for additional messages.')

    return parser.parse_args()

#-------------------------------------------------------------------------------

# Run top-level Balrog script
# TODO: In future, make this `RunDESBalrog`, or something similar to allow for other stacks
def RunBalrog():
    '''
    Main Balrog call.
    '''

    # Parse command line arguments
    args = parse_args()

    # Check verbosity
    vb = args.verbose
    if vb: print('Input arguments: {}'.format(args))

    # Set up Balrog and GalSim configuration
    if vb: print('Setting up configuration...')
    config = setup_config(args)

    # Determine which tiles will have injections
    if vb: print('Creating tiles...')
    tiles = create_tiles(config)

    # Set up output directory structure
    if vb: print('Setting up output directory...')
    setup_output_dir(config, tiles)

    # Now loop over all tiles slated for injection:
    # TODO: This could be parallelized in future
    for i, tile in enumerate(tiles):
        # pudb.set_trace()
        config.reset_gs_config()
        config.set_tile_num(i)

        for real in config.realizations:
            if vb: print('Injecting Tile {}; realization {}'.format(tile.tile_name, real))
            # Reset gs config for each new realization
            tile.set_realization(real)
            tile.reset_bal_config(config)
            tile.generate_objects(config, real)
            # Allow for different band injections in different tiles
            bands = tile.bands
            for band in bands:
                # pudb.set_trace()
                for chip in tile.chips[band]:
                    # Reset injection counter
                    chip.types_injected = 0
                    if config.sim_gals is True:
                        # Determine which Balrog galaxies are contained in chip, and
                        # get their image coordinates
                        in_chip, pos_im = chip.contained_in_chip(tile.gals_pos[real])
                        gals_pos_im = pos_im[in_chip]
                        gals_indx = tile.gals_indx[real][in_chip]
                        chip.set_Ngals(len(gals_pos_im))

                        # pudb.set_trace()
                        # If any galaxies are contained in chip, add gs injection to config list
                        # (or add dummy injection for a few cases)
                        if (chip.Ngals > 0) or (config.inj_objs_only['value'] is True):
                            tile.add_gs_injection(config, chip, gals_indx, gals_pos_im, inj_type='gals')

                    if config.sim_stars is True:
                        # Determine which Balrog stars are contained in chip, and
                        # get their image coordinates
                        in_chip, pos_im = chip.contained_in_chip(tile.stars_pos[real])
                        stars_pos_im = pos_im[in_chip]
                        stars_indx = tile.stars_indx[real][in_chip]
                        chip.set_Nstars(len(stars_pos_im))

                        # If any stars are contained in chip, add gs injection to config list
                        # (or add dummy injection for a few cases)
                        if (chip.Nstars > 0) or (config.inj_objs_only['value'] is True):
                        # if chip.Nstars > 0:
                            tile.add_gs_injection(config, chip, stars_indx, stars_pos_im, inj_type='stars')

                    if (config.inj_objs_only['value'] is False) and (chip.Ngals + chip.Nstars) == 0:
                        # Don't want to skip image for a blank run; need to blank out the image!
                        # NOTE: The first check isn't actually required, as a dummy injection is
                        # added for 'inj_objs_only'=True cases that have no injections. This is
                        # needed to correct the background image in `injector.py`
                        # TODO: Eventually use a return_output_name() function
                        outfile = os.path.join(config.output_dir, 'balrog_images', str(tile.curr_real),
                                tile.tile_name, chip.band, '{}_balrog_inj.fits'.format(chip.name))
                        chip.save_without_injection(outfile)

            # Once all chips in tile have had Balrog injections, run modified config file
            # with GalSim
            if vb: print('Writing Balrog config...')
            tile.write_bal_config()
            if vb: print('Running GalSim for tile...')
            tile.run_galsim(vb=vb)
            if vb: print('Copying extra image planes...')
            tile.copy_extensions(config)
            if vb: print('Truth Catalog...')
            tile.write_truth_catalog(config)

            # pudb.set_trace()

    # TESTING: Can remove in future
    outfile = os.path.join(config.output_dir, 'configs', 'tile_flux_factors.p')
    with open(outfile, 'wb') as f: pickle.dump(config.flux_factors, f)
    # pudb.set_trace()
    return

if __name__ == '__main__':
    ret = RunBalrog()
    sys.exit(ret)
