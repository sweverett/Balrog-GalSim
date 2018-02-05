#####################################################################
#
# I should probably add a description for balrog!
#
#
# Last updated: 1/11/2018
#
# Spencer Everett
# UCSC
# 12/10/2017
#####################################################################

import numpy as np
# import pudb
import os, sys, errno
import subprocess
import shutil
import ntpath
import copy
import time
import galsim
import csv
import yaml
import argparse
import itertools
# TODO: Have automatic check for astropy vs. fitsio!
# try: import fitsio
# except: import astropy
# from fitsio import FITS, FITSHDR, etc.
from astropy.io import fits
from astropy import wcs

# Balrog Galsim image type
import injector

#-------------------------------------------------------------------------------
# Urgent todo's:
<<<<<<< HEAD
# TODO: Correctly normalize galaxy injections! (May be correct as is, per Erin)
# TODO: Implement error handling for galaxy injections / gsparams!
# TODO: Use fitsio when available!
=======
# TODO: Correctly normalize galaxy injections! (maybe ok, per Erin?)
# TODO: Implement error handling for galaxy injections / gsparams!
# TODO: For files without galaxy injections, save only first image! 
>>>>>>> 388649f7ed89e5f3bf030c0c460c6de7f9909462

# Some extra todo's:
# TODO: Redistribute config.galaxies_remainder among first m<n reals, rather than all at end
# TODO: Make filename concatenation more robust! (Mostly done)
# TODO: Clean up old gs_config!! (Maybe allow different options?)
        # NOTE: gs_config implementation currently does not work as intended
# TODO: More general geometry file inputs
# TODO: Add check for python path!
# TODO: Make a log system!
# TODO: Should be able to handle passing geometry file and directories in config OR command line consistently!

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
_supported_input_types = ['ngmix_catalog']#, 'postage_stamps'}

#-------------------------------------------------------------------------------
# Tile class and functions

class Tile(object):
    """
    A `Tile` is a square ?x? deg^2 subsection of the DES footprint. (more details).
    Tiles overlap by 2 arcmin, but balrog galaxies are only injected in the unique
    footrpint area.
    """

    def __init__(self, tile_name, config, realization=0):

        self.tile_name = tile_name
        self.gals_pos = None # injected galaxy position array
        self.gals_indx = None # injected galaxy index array (from input catalog)

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

        # Set the number of galaxies injected per realization
        if config.n_galaxies:
            # Then fixed number of galaxies regardless of tile size
            self.gals_per_real = int(config.n_galaxies / config.n_realizations)
            self.gals_remainder = int(config.n_galaxies % config.n_realizations)
        else:
            # Then gal_density was set; galaxy number depends on area
            self.gals_per_real = round((self.u_area * config.gal_density) / (1.0 * config.n_realizations))

        # Set tile directory structure and chip list
        self.dir = os.path.join(config.tile_dir, self.tile_name)
        self._set_bands(config)

        self._create_chip_list(config)

        # Set initial injection realization number
        self.set_realization(realization)

        # Setup new Balrog config file for chip to be called by GalSim executable
        self._setup_bal_config(config)

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
        'Load WCS info for each tile from geometry file.'

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

        return

    def _set_bands(self, config):
        '''
        For now, just set to 'griz'. May want to use a different subset in the future.
        '''

        # Will store directory locatations for each band
        self.band_dir = {}

        try:
            # Grab selected bands in config, if present
            self.bands = config.gs_config[0]['input'][config.input_type]['bands']

            # Make sure there aren't any incorrect inputs
            for band in self.bands:
                if band not in _allowed_bands:
                    raise ValueError('Passed band {} is not one of the allowed bands in {}'\
                                     .format(band, _allowed_bands))

        except KeyError:
            # By default, use griz
            print('Warning: No injection bands were passed in config. Using `griz` by default')
            self.bands = 'griz'

        for band in self.bands:
            self.band_dir[band] = os.path.join(self.dir, 'nullwt-{}'.format(band))
            # Make band directory if it doesn't already exist
            os.makedirs(self.band_dir[band])


        return

    def _create_chip_list(self,config):
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
            #TODO: Check if the directory exists! (check if solution works)
            try:
                file_list = os.listdir(b_dir)
            except OSError:
                # This tile does not have chip images in this band
                file_list = None
                continue
            # Check that it is a nullwt fits file
            for f in file_list:
                #TODO: Check if fits file is actually a nullwt image
                # ischip = self.is_nullwt_chip(f)
                if f.endswith('nullwt.fits'): self.chip_list[band].append(f)

                # Add chip to list
                filename = os.path.join(b_dir, f)
                # pudb.set_trace()
                self.chips[band].append(Chip(filename, band, config, tile_name=self.tile_name))

        # pudb.set_trace()

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

        # Create new config file
        filename = 'bal_config_' + str(self.curr_real) + '_' + self.tile_name + '.yaml'
        self.bal_config_dir = config.output_dir + '/configs/'
        self.bal_config_file = self.bal_config_dir + filename

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

        # Keep track if balrog config has been modified from original
        # TODO/NOTE: Probably don't need for this implementation; delete soon.
        # self.bal_config_modified is False:

        return

    def set_realization(self, i):
        '''
        Sets parameters relevant to the current injection realization i.
        '''

        self.curr_real = i

        # TODO: More relevant functions?

        return

    def generate_galaxies(self, config, realization):
        '''
        Generate a list of galaxy positions for tile, for a given realization
        (starts counting at 0).
        '''

        # Initialize galaxy positions and indices if this is the first realization
        if self.gals_pos is None:
            # self.gals_pos = np.zeros(shape=(config.n_realizations, self.gals_per_real), dtype=tuple)

            if (config.n_galaxies is not None) and (realization == config.n_realizations):
                # TODO: Better to distribute remaining m<n galaxies to the first m realizations!
                # Add remaining galaxies on final realization tile
                ngals = self.gals_per_real + self.gals_remainder
            else:
                ngals = self.gals_per_real
            self.gals_pos = np.zeros(shape=(config.n_realizations, ngals, 2))
            self.gals_indx = np.zeros(shape=(config.n_realizations, ngals))

        # Generate galaxy coordinates
        #TODO: Could add more sampling methods than uniform
        ra = sample_uniform_ra(self.ramin, self.ramax, self.gals_per_real, boundary_cross=self.ra_boundary_cross)
        dec = sample_uniform_dec(self.decmin, self.decmax, self.gals_per_real, unit='deg')
        # self.gals_pos[realization] = zip(ra, dec)
        self.gals_pos[realization] = np.column_stack((ra, dec))

        # pudb.set_trace()

        # Generage galaxy indices (in input catalog)
        #TODO: could add more sampling methods than uniform
        indices = sample_uniform_indx(self.gals_per_real, config.input_nobjects)
        self.gals_indx[realization] = indices

        #TODO: Generate any other needed galaxy properties!

        return

    def write_bal_config(self):
        '''
        Write appended balrog config to a yaml file.
        # TODO: In future, allow more config types!
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

        # Reset injections
        self.has_injections = False

        return

    def add_gs_injection(self, chip, gals_indx, gals_pos_im):
        '''
        This function appends the global GalSim config with an additional simulation to
        be done using nullwt chip-specific infomation and Balrog injected galaxy positions
        in image coordinates.
        NOTE/TODO: This is in contrast to doing the injection through the Config object.
        In future, will likely remove the alternate approach (but may keep it for multiple
        injection options).
        '''

        # TODO: Should probably do some size checks on gals_pos_im.
        # For now, just do the basics

        # Skip injection if no positions are passed
        if len(gals_pos_im) == 0:
            # TODO: Print warning!
            return

        # pudb.set_trace()

        # Initialize new output section in balrog config
        self.bal_config.append({})
        self.bal_config_len += 1

        # Injection index
        i = self.bal_config_len - 1

        # Set initial image, chip wcs, and galaxy positions
        chip_file = chip.filename
        x, y = gals_pos_im[:,0].tolist(), gals_pos_im[:,1].tolist()
	nobjs = len(gals_pos_im)
        # pudb.set_trace()
        self.bal_config[i]['image'] = {
            'initial_image' : chip_file,
	    'nobjects' : nobjs,
            'wcs' : { 'file_name' : chip_file },
            'image_pos' : {
                'type' : 'XY',
                'x' : { 'type' : 'List', 'items' : x },
                'y' : { 'type' : 'List', 'items' : y }
            }
        }

        # Set the injected galaxies' catalog indices
        indices = gals_indx.tolist()
        self.bal_config[i]['gal'] = {
            'index' : {
                'type' : 'List',
                'items' : indices
            }
        }

        # Set the PSF input file, if present
        psf_file = chip.psf_filename
        if psf_file is not None:
            self.bal_config[i]['input'] = {
                'des_psfex' : {
                    # TODO: Remove when possible
                    # 'dir' : chip.psf_dir,
                    'file_name' : psf_file,
                    'image_file_name' : chip_file,
                    }
            }

        # Set the initial image and output filenames
        out_file = os.path.join(self.output_dir, 'balrog_images', str(self.curr_real), self.tile_name, \
                                chip.band, '{}_balrog_inj.fits'.format(chip.name))
        self.bal_config[i]['output'] = {'file_name' : out_file}

        if self.has_injections is False:
            self.has_injections = True

        # pudb.set_trace()

        return

    def run_galsim(self, vb=False):
        '''
        Run full GalSim executable of the modified gs_config. Will inject all Balrog
        galaxies contained in a given tile for a given realization in all chips in
        all bands.
        '''

        # If there are no galaxies to inject, simply save chip file to balrog directory
        if self.has_injections is False:
            # TODO: Actually implement! For now, we just won't write.
            # NOTE: Can't run galsim executable (besides being a waste of time) as
            # the base layer yaml config won't be valid without extra additions.
            return

        # A new GalSim config file for Balrog injections has been created and all simulations
        # can now be run simultaneously using all GalSim machinery
        #bashCommand = 'galsim {} -v 2 -l gs_logfile'.format(self.bal_config_file)
        bashCommand = 'galsim {}'.format(self.bal_config_file)

        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

        if vb:
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
        # maybe implement an additional log?

        return

#-------------------------
# Related Tile functions

def create_tiles(args, config):
    '''
    Create list of `Tile` objects given input args and configuration file.
    '''

    tile_list = load_tile_list(args.tile_list, config)

    # pudb.set_trace()

    # Will keep a list of desired tiles
    tiles = []

    for tile_name in tile_list:
        tiles.append(Tile(tile_name, config))

    return tiles

def load_tile_list(tile_list_file, config):
    #TODO: Allow many file types

    if tile_list_file.lower().endswith('.csv'):
        tile_list = open_csv_list(tile_list_file)
    elif tile_list_file.lower().endswith('.txt'):
        tile_list = open_txt_list(tile_list_file)
    # elif ... # TODO: add more types!
    else:
        raise Exception('`tile_list` must be in a `.csv` or `.txt` file!')

    if config.vb:
        print('Loaded {} tiles...'.format(len(tile_list)))

    return tile_list

#-------------------------------------------------------------------------------
# Chip

class Chip(object):
    '''
    DES chip object.
    # NOTE: While DES chip image data is saved as (RA,DEC), the orientation
    # of the chip image is roated ccw by 90 deg's.
    '''

    def __init__(self, filename, band, config, tile_name=None):

        self.filename = filename
        self.fits_filename = ntpath.basename(filename)
        self.tile_name = tile_name # Name of parent tile, if given
        self.band = band

        self._set_name(config)
        self._set_psf(config)
        self._set_wcs()

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
        # TODO: For now just PSFEx, but should allow for Piff (or others!) in future.
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

    def _set_wcs(self, config=None):
        hdr = fits.getheader(self.filename)
        # Get chip WCS
        self.wcs = wcs.WCS(hdr)

        # Get corners (chip not perfectly aligned to RA/DEC coordinate system)
        # Especially useful for galaxy position rejection.
        #TODO: Determine if we actually need this!
        #NOTE: The nullwt chips are not oriented in the standard way.
        # In a typical (RA,DEC) projection space:
        #
        # DEC increasing up
        # .
        # .
        # 1-----------4
        # -           -
        # -           -
        # 2-----------3....RA increasing right
        #
        # In the DES nullwt chip orientation:
        #
        # RA increasing up
        # .
        # .
        # 4-------3
        # -       -
        # -       -
        # -       -
        # -       -
        # 1-------2....DEC decreasing right
        #
        # Because of the strange orientation, wcs_world2pix([ra, dec]) will
        # return correct image coordinates but will be flipped from what we
        # normally expect; i.e. IN_PIX: (x,y) ~= (DEC, RA).
        # This will affect how we implement contained_in_chip().

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

    def save_without_injection(self, outfile):
        '''
        If there are no Balrog galaxies to inject in the chip area, then save
        copy of current chip image in the new Balrog image format.
        '''

        # Only want to save the first HDU of nullwt image
        # TODO: Switch to fitsio eventually!
        with fits.open(outfile) as f:
            hdu0 = f[0]
            try:
                hdu0.writeto(outfile, overwrite=True)
            except OSError:
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
                hdu0.writeto(outfile, overwrite=True)

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

    def __init__(self, args):
        # Process configuration file

        #TODO: Initialize all class attributes defined alter to None

        # Save command-line arguments. 'args' is a Namespace, can access fields
        # as 'arg.config_dir'
        self.args = args
        self.vb = args.verbose
        # self.input_cat_file = args.input_catalog # Can now set in config file
        self.tile_dir = args.tile_dir
        self.config_dir = args.config_dir
        self.psf_dir = args.psf_dir
        self.output_dir = args.output_dir

        # Set directories to empty string if none passed
        if self.tile_dir is None: self.tile_dir = ''
        if self.config_dir is None: self.config_dir = ''
        if self.psf_dir is None: self.psf_dir = 'psfs'
        if self.output_dir is None: self.output_dir = 'balrog_images/'
        # TODO: Can we make a version of this work?
        # for d in [self.tile_dir, self.config_dir, self.output_dir]:
        #     if d is None: d = ''
        # TODO: Allow multiple config file types (.yaml, .json, etc.)

        # Process GalSim config file
        self._read_gs_config()
        # Process geometry file
        self._load_tile_geometry()
        # Process input catalog
        self._load_input_catalog()

        return

    def _read_gs_config(self):
        # Process .yaml config file
        # TODO: Allow multiple config types

        self.gs_config_file = self.config_dir + self.args.config_file
        # Add config directory path, if needed
        # if self.args.config_dir:
        #     self.gs_config_file = self.args.config_dir + self.gs_config_file

        # NOTE: `self.config` will be a list of `OrderedDict`'s with length given by
        # the number of configurations specified in the file. To work with a specific
        # simulated image, say with `galsim.config.BuildImage()`, you must pass a
        # *single* config file (e.g. BuildImage(config[3]) ).
        self.gs_config = galsim.config.process.ReadYaml(self.gs_config_file)

        # pudb.set_trace()

        # Double NOTE! Making the switch to multi-output yaml files.
        # NOTE: For now, we will not accept multi-output yaml files. We will work with
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

        # Process input 'n_realizations':
        try:
            self.n_realizations = self.gs_config[0]['image']['n_realizations']
        except KeyError:
            # Default is 1 realization
            self.n_realizations = 1

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
            #TODO: Should allow for more unit input types! e.g.:
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

        # if not self.n_galaxies:
        #     # Default is 1000 galaxies per tile
        #     self.n_galaxies = 100
        # if not self.gal_density:
        #     # Default is 10 galaxies per arcmin^2
        #     self.gal_density = 10.0 # arcmin^2

        # Finally, calculate number of galaxy injects per *realization*
        #TODO: Revisit to see if possible. It appeaers that tile areas may
        # note be uniform. (Spoiler; they're not)
        # if self.n_galaxies:
        #     self.gals_per_real = self.n_galaxies / self.n_realizations
        # if self.gal_density:
        #     # By now, density must be in arcmin^2
        #     tile_area =
        # self.gals_per_real =

        return

    def _load_tile_geometry(self):
        '''
        TODO: make more general
        '''

        self.geom_file = self.args.geom_file

        # now a required parameter!
        # if self.geom_file is None:
        #     # TODO: Make reasonable assumption about file in local directory
        #     # (e.g. serach for 'Y{}A{}GEOM', or similar)
        #     raise Exception('For now, must pass a geometry file! (e.g. `Y3A2_COADDTILE_GEOM.fits`)')

        # Add geom dir path, if needed (must be same as config directory)
        if self.config_dir:
            self.geom_file = self.args.config_dir + self.geom_file

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

    def _load_input_catalog(self):
        '''
        Load any relevant info from the input catalog (for now just ngmix)
        '''

        # Determine input type
        # TODO: For now just ngmix, but will want postage stamps in future
        input_cat_type = self._determine_input_type()

        if input_cat_type is 'ngmix_catalog':
            import ngmix_catalog
            # As we are outside of the GalSim executable, we need to register
            # the input type explicitly
            galsim.config.RegisterInputType('ngmix_catalog', ngmix_catalog.ngmixCatalogLoader(
                ngmix_catalog.ngmixCatalog, has_nobj=True))

            # Return the number of objects in input catalog that will be injected
            # NOTE: Calling `ProcessInputNObjects()` builds the input catalog
            # in a minimal way that determines the number of objects from the
            # original catalog that make it past all of the mask cuts, etc.
            self.input_nobjects = galsim.config.ProcessInputNObjects(self.gs_config[0])

        else:
            # Add more types later!
            raise ValueError('For now, only ngmix catalogs can be used for injections!')

        return

    def _determine_input_type(self):
        '''#TODO: For now we just select ngmix_catalog. In future, we may want to expand
        this to handle postage stamp inputs, for example.
        '''
        # TODO: Check that is is one of the supported input types!
        # TODO: While just a placeholder method anyway, would be good to check that the
        # input catalog actually *is* an ngmix catalog!

        self.input_type = 'ngmix_catalog'

        return self.input_type

    def set_realization(self, i):
        '''
        Sets parameters relevant to the current injection realization i.
        NOTE/TODO: May remove in future to store locally in tiles instead.
        '''

        self.curr_real = i

        # TODO: More relevant functions?

        return

    def add_gs_injection(self, tile_name, chip, gals_pos_im):
        '''
        WARNING: Deprecated!!
        This function appends the global GalSim config with an additional simulation to
        be done using nullwt chip-specific infomation and Balrog injected galaxy positions
        in image coordinates.
        NOTE: The original inputted config structure is saved in
        self.original_gs_config.
        NOTE/TODO: This is only used if you modify the config dictionary directly. We may
        be switching to appended yaml files soon.
        '''

        chip_file = chip.filename
        # TODO: Eventually, give it a more sensible name
        out_file = self.output_dir + tile_name + '/' + chip.band + '/BALROG_' + chip.fits_filename

        if self.gs_config_modified is False:
            # If this is the first GalSim injection added,
            # then just modify the original config
            i = 0
            self.gs_config_modified = True
        else:
            # Add a new entry to the config w/ same values
            # as original gs config file
            self.gs_config.append(self.orig_gs_config[0])
            i = len(self.gs_config) - 1

        # Set initial image, chip wcs, and galaxy positions
        x, y = gals_pos_im[:,0].tolist(), gals_pos_im[:,1].tolist()
	nobjs = len(gals_pos_im)
        self.gs_config[i]['image'] = {
            'initial_image' : chip_file,
	    'nobjects' : nobjs,
            'wcs' : { 'file_name' : chip_file },
            'image_pos' : {
                'type' : 'XY',
                'x' : { 'type' : 'List', 'items' : x },
                'y' : { 'type' : 'List', 'items' : y }
            }
        }

        # Set the initial image and output filenames
        self.gs_config[i]['output']['file_name'] = out_file

        return

    def append_bal_config(self, tile_name, chip, gals_pos_im):
        '''
        An alternative to add_gs_injection() (that may eventually replace it).
        Sidesteps issues with running GalSim manually and allows the use of the
        GalSim executable with multi-output yaml files.
        # TODO: Would be good to allow more config input types (like .json) in
        future.
        '''

        # have base config as only document; make sure it's not a multi-output file beforehand
        bal_config = list(yaml.safe_load_all(self.bal_config_file))

        # Make sure that config length hasn't been modified incorrectly
        c_len = len(bal_config)
        if c_len != self.bal_config_len:
            raise Exception('Warning: appended config has been incorrectly modified! Should have',
                            '{} entries, but currently has {}'.format(c_len, self.bal_config_len))

        # for each chip, append new dictionary & fill it w/ relevant info
        # bal_config
        # after all chips, use yaml.dump_all(config, appended_file) to write multi-output yaml file
        # use galsim executible on appended file and run all simultaneously!

        return

    def reset_gs_config(self):
        '''
        This function resets the gs_config after a tile is run.
        NOTE: This should only be used if multiple tiles are processed
        on a single batch.
        '''

        # TODO: Test, & make more robust!
        if self.gs_config_modified is True:
            # New tile, so reset config to original input
            self.gs_config = self.orig_gs_config
            self.gs_config_modified = False

        # If gs_config_modified is False, then there is nothing to reset!

        return

    def run_galsim(self):
        '''
        Run full GalSim executable of the modified gs_config. Will inject all Balrog
        galaxies contained in a given tile for a given realization in all chips in
        all bands.
        '''

        # Import any modules if requested
        # ImportModules(config)

        # pudb.set_trace()

        # Need to register Balrog as a valid image type
        # BalrogImageBuilder = balrog.Balr
        # galsim.config.RegisterImageType('Balrog', injector.BalrogImageBuilder())

        for i in range(len(self.gs_config)):
            if self.vb: print('Injecting chip {} of {}'.format(i,len(self.gs_config)))
            galsim.config.Process(self.gs_config[i])

        return

#-------------------------
# Related config functions

def setup_config(args):
    '''
    # TODO: For now, just create new config object. Can make more complex later.
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
        # ramax < ramin accross boundary.
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

    # Can just use `sample_unfiorm()`; only here for future naming
    # conventions.

    return sample_uniform(n1, n2, N)

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

def parse_args():
    '''
    Parse command line arguments.
    TODO: There are likely more optional inputs that need to be added.
    '''

    parser = argparse.ArgumentParser()
    # Required argument for GalSim config file
    parser.add_argument('config_file', help='.yaml or .json confg file that specifies the desired GalSim'
                        'simulation and injection.')
    # Required argument for tile list
    parser.add_argument('tile_list', help='.txt or .csv file containing list of tiles to be processed.')
    # Required argument for tile geometry file
    parser.add_argument('geom_file', help='Tile geometry file (e.g. `Y3A2_COADDTILE_GEOM`)')
    # Optional argument for config file(s) directory location (if not .)
    parser.add_argument('-c', '--config_dir', help='Directory that houses the GalSim and related config files.')
    # Optional argument for DES tiles directory location
    parser.add_argument('-t', '--tile_dir', help='Directory of DES tiles containing single exposure images.')
    # Optional argument for a tile's psf directory name (if not contained in tile_dir/{TILE}/psfs)
    parser.add_argument('-p', '--psf_dir', help='Directory that houses individual psf files for DES chips.')
    # Optional argument for output directory (if not .)
    parser.add_argument('-o', '--output_dir', help='Directory that houses output Balrog images.')
    # Optional argument for verbose messages
    parser.add_argument('-v', '--verbose', action='store_true', help='Turn on verbose mode for additional messages.')


    #TODO: Might add this functionality later, but for now these parameters are specified in GalSim config file
    # Required argument for input catalog
    # parser.add_argument('input_catalog', help='(For now) input ngmix catalog whose galaxies will be injected.')
    # Optional argument for injection density (arcmin^2)
    # parser.add_argument('-r', '--realization_density', help='Galaxy injection density (in arcmin^2).')
    # Optional argument for total injection number N
    # parser.add_argument('-n,', '--number', help='Total number of galaxies to be injected into tile.')

    # TODO: Should be able to handle passing geometry file and directories in config OR command line consistently!

    return parser.parse_args()

#-------------------------------------------------------------------------------

# Run top-level Balrog script
def RunBalrog():
    '''
    Main Balrog call.
    #TODO: Write description!
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
    tiles = create_tiles(args, config)

    # Now loop over all tiles slated for injection:
    # TODO: This should be parallelized with `multiprocessing` in future
    for tile in tiles:
        # pudb.set_trace()
        config.reset_gs_config()

        if vb: print('Injecting Tile {}'.format(tile.tile_name))

        # TODO: This (maybe?) should be parallelized with `multiprocessing` in future
        for real in range(config.n_realizations):
            # Reset gs config for each new realization
            # TODO/config: decide on implementation!
            # config.reset_gs_config()
            # config.set_realization(i)
            tile.reset_bal_config(config)
            tile.set_realization(real)
            tile.generate_galaxies(config, real)
            # Allow for different band injections in different tiles
            bands = tile.bands
            for band in bands:
                # pudb.set_trace()
                for chip in tile.chips[band]:
                    # Determine which Balrog galaxies are contained in chip, and
                    # get their image coordinates
                    in_chip, pos_im = chip.contained_in_chip(tile.gals_pos[real])
                    gals_pos = tile.gals_pos[real][in_chip] # NOTE: is this line needed?
                    gals_pos_im = pos_im[in_chip]
                    gals_indx = tile.gals_indx[real][in_chip]

                    # pudb.set_trace()

                    # If any galaxies are contained in chip, add gs injection to
                    # config list
                    if len(gals_pos) > 0:
                        tile.add_gs_injection(chip, gals_indx, gals_pos_im)
                    else:
                        # TODO: Eventually use a return_output_name() function
                        outfile = os.path.join(config.output_dir, 'balrog_images', str(tile.curr_real),
                                   tile.tile_name, chip.band, '{}_balrog_inj.fits'.format(chip.name))
                        chip.save_without_injection(outfile)

            # Once all chips in tile have had Balrog injections, run modified config file
            # with GalSim
            if vb is True: print('Writing Balrog config...')
            tile.write_bal_config()
            # pudb.set_trace()
            if vb is True: print('Running GalSim for tile...')
            tile.run_galsim(vb=vb)

            # pudb.set_trace()

    # pudb.set_trace()

    return

if __name__ == '__main__':
    ret = RunBalrog()
    sys.exit(ret)

