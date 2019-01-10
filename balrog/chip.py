import numpy as np
import os, sys, errno
import ntpath
from astropy import wcs
from astropy.io import fits
import fitsio
import warnings

#-------------------------------------------------------------------------------
# Chip class and functions

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
        self.nobjects = {}
        self.total_n_objects = 0

        # TODO: Can probably get rid of this after update
        # Keeps track of how many input types for this chip have been added
        # to bal_config by add_gs_injection()
        self.types_injected = 0
        self.N_inj_types = len(config.inj_types)

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
            if self.psf_type in config._supported_psf_types:
                self.psf_extension = config._psf_extensions[self.psf_type]
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
        if config.input_zp is not None:
            self.flux_factor = np.power(10.0, 0.4 * (self.zeropoint - config.input_zp))
        else:
            self.flux_factor = 1.0

        return

    def _set_bkg(self, config, tile):
        '''
        Set chip background file, if needed for grid test.
        '''

        if config.inj_objs_only['noise'] in config._valid_background_types:
            assert tile is not None
            self.bkg_file = tile.bkg_files[self.band][self.name]

        return

    def contained_in_chip(self, pos):
        '''
        For an input vector of (RA,DEC) positions, returns a boolean vector
        of whether each position is contained within the chip image.
        '''

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

    def set_nobjects(self, Nobjs, inj_type):
        self.nobjects[inj_type] = Nobjs
        self.total_n_objects += Nobjs

        return

    def save_without_injection(self, outfile):
        '''
        If there are no Balrog objects to inject in the chip area, then save
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
