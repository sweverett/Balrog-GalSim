#----------------------------------------------------------------------------
# Classes for new input type `MEDSCatalog`; intended for use with Balrog.
# Similar to GalSim's native RealGalaxyCatalog class, but with structural
# differences due to how MEDS stores it's postage stamp images.
# Some methods and naming choices were made to be as consistent as possible
# between the classes
#
# Contributors:
# Spencer Everett (UCSC)
#----------------------------------------------------------------------------

import galsim
from galsim import GSObject
#from galsim import RealGalaxyCatalog
from galsim.config.gsobject import RegisterObjectType
from galsim.config.gsobject import SkipThisObject
from galsim._pyfits import pyfits
import logging
import os
import meds
import numpy as np

import pudb

class MEDSGalaxy(GSObject):
    """A class describing real galaxies from some training dataset.  Its underlying implementation
    uses a Convolution instance of an InterpolatedImage (for the observed galaxy) with a
    Deconvolution of another InterpolatedImage (for the PSF).

    This class uses a catalog describing galaxies in some training data (for more details, see the
    RealGalaxyCatalog documentation) to read in data about realistic galaxies that can be used for
    simulations based on those galaxies.  Also included in the class is additional information that
    might be needed to make or interpret the simulations, e.g., the noise properties of the training
    data.  Users who wish to draw RealGalaxies that have well-defined flux scalings in various
    passbands, and/or parametric representations, should use the COSMOSGalaxy class.

    Because RealGalaxy involves a Deconvolution, `method = 'phot'` is unavailable for the
    drawImage() function.

    Initialization
    --------------

        >>> real_galaxy = galsim.RealGalaxy(real_galaxy_catalog, index=None, id=None, random=False,
        ...                                 rng=None, x_interpolant=None, k_interpolant=None,
        ...                                 flux=None, pad_factor=4, noise_pad_size=0,
        ...                                 gsparams=None)

    This initializes `real_galaxy` with three InterpolatedImage objects (one for the deconvolved
    galaxy, and saved versions of the original HST image and PSF). Note that there are multiple
    keywords for choosing a galaxy; exactly one must be set.

    Note that tests suggest that for optimal balance between accuracy and speed, `k_interpolant` and
    `pad_factor` should be kept at their default values.  The user should be aware that significant
    inaccuracy can result from using other combinations of these parameters; more details can be
    found in http://arxiv.org/abs/1401.2636, especially table 1, and in comment
    https://github.com/GalSim-developers/GalSim/issues/389#issuecomment-26166621 and the following
    comments.

    If you don't set a flux, the flux of the returned object will be the flux of the original
    HST data, scaled to correspond to a 1 second HST exposure (though see the `area_norm`
    parameter below, and also caveats related to using the `flux` parameter).  If you want a flux
    appropriate for a longer exposure, or for a telescope with a different collecting area than HST,
    you can either renormalize the object with the `flux_rescale` parameter, or by using the
    `exptime` and `area` parameters to `drawImage`.

    Note that RealGalaxy objects use arcsec for the units of their linear dimension.  If you
    are using a different unit for other things (the PSF, WCS, etc.), then you should dilate
    the resulting object with `gal.dilate(galsim.arcsec / scale_unit)`.

    @param real_galaxy_catalog  A RealGalaxyCatalog object with basic information about where to
                            find the data, etc.
    @param index            Index of the desired galaxy in the catalog. [One of `index`, `id`, or
                            `random` is required.]
    @param id               Object ID for the desired galaxy in the catalog. [One of `index`, `id`,
                            or `random` is required.]
    @param random           If True, then select a random galaxy from the catalog.  If the catalog
                            has a 'weight' associated with it to allow for correction of selection
                            effects in which galaxies were included, the 'weight' factor is used to
                            remove those selection effects rather than selecting a completely random
                            object.
                            [One of `index`, `id`, or `random` is required.]
    @param rng              A random number generator to use for selecting a random galaxy
                            (may be any kind of BaseDeviate or None) and to use in generating
                            any noise field when padding.  This user-input random number
                            generator takes precedence over any stored within a user-input
                            CorrelatedNoise instance (see `noise_pad` parameter below).
                            [default: None]
    @param x_interpolant    Either an Interpolant instance or a string indicating which real-space
                            interpolant should be used.  Options are 'nearest', 'sinc', 'linear',
                            'cubic', 'quintic', or 'lanczosN' where N should be the integer order
                            to use. [default: galsim.Quintic()]
    @param k_interpolant    Either an Interpolant instance or a string indicating which k-space
                            interpolant should be used.  Options are 'nearest', 'sinc', 'linear',
                            'cubic', 'quintic', or 'lanczosN' where N should be the integer order
                            to use.  We strongly recommend leaving this parameter at its default
                            value; see text above for details.  [default: galsim.Quintic()]
    @param flux             Total flux, if None then original flux in image is adopted without
                            change.  Note that, technically, this parameter sets the flux of the
                            postage stamp image and not the flux of the contained galaxy.  These two
                            values will be strongly correlated when the signal-to-noise ratio of the
                            galaxy is large, but may be considerably different if the flux of the
                            galaxy is small with respect to the noise variations in the postage
                            stamp.  To avoid complications with faint galaxies, consider using the
                            flux_rescale parameter.  [default: None]
    @param flux_rescale     Flux rescaling factor; if None, then no rescaling is done.  Either
                            `flux` or `flux_rescale` may be set, but not both. [default: None]
    @param pad_factor       Factor by which to pad the Image when creating the
                            InterpolatedImage.  We strongly recommend leaving this parameter
                            at its default value; see text above for details.  [default: 4]
    @param noise_pad_size   If provided, the image will be padded out to this size (in arcsec)
                            with the noise specified in the real galaxy catalog. This is
                            important if you are planning to whiten the resulting image.  You
                            should make sure that the padded image is larger than the postage
                            stamp onto which you are drawing this object.
                            [default: None]
    @param area_norm        Area in cm^2 by which to normalize the flux of the returned object.
                            When area_norm=1 (the default), drawing with `drawImage` keywords
                            exptime=1 and area=1 will simulate an image with the appropriate number
                            of counts for a 1 second exposure with the original telescope/camera
                            (e.g., with HST when using the COSMOS catalog).  If you would rather
                            explicitly specify the collecting area of the telescope when using
                            `drawImage` with a `RealGalaxy`, then you should set area_norm equal to
                            the collecting area of the source catalog telescope when creating the
                            `RealGalaxy` (e.g., area_norm=45238.93416 for HST). [default: 1]
    @param gsparams         An optional GSParams argument.  See the docstring for GSParams for
                            details. [default: None]
    @param logger           A logger object for output of progress statements if the user wants
                            them.  [default: None]

    Methods
    -------

    There are no additional methods for RealGalaxy beyond the usual GSObject methods.
    """
    _req_params = {'band' : str}
    _opt_params = { 'icutout' : int,
                    'x_interpolant' : str,
                    'k_interpolant' : str,
                    'flux' : float,
                    'flux_rescale' : float,
                    'pad_factor' : float,
                    'noise_pad_size' : float,
                    'area_norm' : float
                  }
    _single_params = [ { 'index' : int , 'id' : str , 'random' : bool } ]
    _takes_rng = True

    _valid_bands = 'grizy'

    def __init__(self, meds_catalog, band, index=None, id=None, random=False,
                 rng=None, icutout=0, x_interpolant=None, k_interpolant=None, flux=None,
                 flux_rescale=None, pad_factor=4, noise_pad_size=0, area_norm=1.0, gsparams=None,
                 logger=None):

        if rng is None:
            rng = galsim.BaseDeviate()
        elif not isinstance(rng, galsim.BaseDeviate):
            raise TypeError("The rng provided to MEDSGalaxy constructor is not a BaseDeviate")
        self.rng = rng

        if band not in self._valid_bands:
            raise ValueError('Passed band %s is not valid ' %band +
                            '- only %s are currently accepted.' %self._valid_bands)
        if band not in meds_catalog.bands:
            raise ValueError('Passed band %s is not contained in loaded MEDS catalog.' % band)

        if flux is not None and flux_rescale is not None:
            raise TypeError("Cannot supply a flux and a flux rescaling factor!")

        # Get the index to use in the catalog
        if index is not None:
            if id is not None or random:
                raise AttributeError('Too many methods for selecting a galaxy!')
            # TODO:
            iobj = index
        elif id is not None:
            if random:
                raise AttributeError('Too many methods for selecting a galaxy!')
            print '\n id = %i' % id
            print '\n iobj = %i\n' % iobj
            iobj = meds_catalog.getIobjForID(id)
        elif random:
            ud = galsim.UniformDeviate(self.rng)
            iobj = int(meds_catalog.nobjects * ud())
            if hasattr(meds_catalog, 'weight'):
                # If weight factors are available, make sure the random selection uses the
                # weights to remove the catalog-level selection effects (flux_radius-dependent
                # probability of making a postage stamp for a given object).
                while ud() > meds_catalog.weight[iobj]:
                    # Pick another one to try.
                    iobj = int(meds_catalog.nobjects * ud())
        else:
            raise AttributeError('No method specified for selecting a galaxy!')
        if logger:
            logger.debug('MEDSGalaxy %d: Start MEDSGalaxy constructor.',iobj)
        # print('\niobj = {}!\n'.format(iobj))

        # Used for logging purposes
        ident = meds_catalog.ident[iobj]

        # Read in the galaxy, PSF images; for now, rely on pyfits to make I/O errors.
        self.gal_image = meds_catalog.getGalImage(iobj, icutout, band)
        if np.sum(self.gal_image.array) == 0:
            if logger:
                logger.debug('MEDSGalaxy %d: gal_image failed; total flux is 0.',iobj)
            raise SkipThisObject
        if logger:
            logger.debug('MEDSGalaxy %d: Got gal_image',iobj)

        self.psf_image = meds_catalog.getPSFImage(iobj, icutout, band)
        if np.sum(self.gal_image.array) == 0:
            if logger:
                logger.debug('MEDSGalaxy %d: psf_image failed; total flux is 0.',iobj)
            raise SkipThisObject
        if logger:
            logger.debug('MEDSGalaxy %d: Got psf_image',iobj)

        #self.noise = meds_catalog.getNoise(iobj, self.rng, gsparams)
        # We need to duplication some of the MEDSCatalog.getNoise() function, since we
        # want it to be possible to have the MEDSCatalog in another process, and the
        # BaseCorrelatedNoise object is not picklable.  So we just build it here instead.
        # noise_image, pixel_scale, var = meds_catalog.getNoiseProperties(iobj)
        # if logger:
        #     logger.debug('MEDSGalaxy %d: Got noise_image',iobj)

        self.catalog_files = { 'MEDS':meds_catalog.getMEDSFiles(),
                               'PSF' :meds_catalog.getPSFFiles()  }
        self.catalog = meds_catalog

        # TODO: Do we need any of this?
        # if noise_image is None:
        #     self.noise = galsim.UncorrelatedNoise(var, rng=self.rng, scale=pixel_scale,
        #                                           gsparams=gsparams)
        # else:
        #     ii = galsim.InterpolatedImage(noise_image, normalization="sb",
        #                                   calculate_stepk=False, calculate_maxk=False,
        #                                   x_interpolant='linear', gsparams=gsparams)
        #     self.noise = galsim.correlatednoise._BaseCorrelatedNoise(self.rng, ii, noise_image.wcs)
        #     self.noise = self.noise.withVariance(var)
        # if logger:
        #     logger.debug('MEDSGalaxy %d: Finished building noise',iobj)

        # Save any other relevant information as instance attributes
        self.index = iobj
        self.pixel_scale = float(meds_catalog.pixel_scale)
        self._x_interpolant = x_interpolant
        self._k_interpolant = k_interpolant
        self._pad_factor = pad_factor
        self._noise_pad_size = noise_pad_size
        self._flux = flux
        self._flux_rescale = flux_rescale
        self._area_norm = area_norm
        self._gsparams = gsparams

        # Convert noise_pad to the right noise to pass to InterpolatedImage
        if noise_pad_size:
            noise_pad = self.noise
        else:
            noise_pad = 0.

        # Build the InterpolatedImage of the PSF.
        self.original_psf = galsim.InterpolatedImage(
            self.psf_image, x_interpolant=x_interpolant, k_interpolant=k_interpolant,
            flux=1.0, gsparams=gsparams)
        if logger:
            logger.debug('MEDSGalaxy %d: Made original_psf', iobj)

        # Build the InterpolatedImage of the galaxy.
        # Use the stepk value of the PSF as a maximum value for stepk of the galaxy.
        # (Otherwise, low surface brightness galaxies can get a spuriously high stepk, which
        # leads to problems.)
        try :
            self.original_gal = galsim.InterpolatedImage(
                    self.gal_image, x_interpolant=x_interpolant, k_interpolant=k_interpolant,
                    pad_factor=pad_factor, noise_pad_size=noise_pad_size,
                    calculate_stepk=self.original_psf.stepk,
                    calculate_maxk=self.original_psf.maxk,
                    noise_pad=noise_pad, rng=self.rng, gsparams=gsparams)
        except RuntimeError as e:
            # TODO: Is this needed anymore?
            pass
            # self.original_gal = galsim.Gaussian(flux=0, sigma=1)
            # self.original_gal = None

        if logger:
            logger.debug('MEDSGalaxy %d: Made original_gal', iobj)

        # Only alter normalization if a change is requested
        if flux is not None or flux_rescale is not None or area_norm != 1:
            if flux_rescale is None:
                flux_rescale = 1.0
            flux_rescale /= area_norm
            if flux is not None:
                flux_rescale *= flux/self.original_gal.flux
            self.original_gal *= flux_rescale
            self.noise *= flux_rescale**2

        # Calculate the PSF "deconvolution" kernel
        psf_inv = galsim.Deconvolve(self.original_psf, gsparams=gsparams)

        # Initialize the SBProfile attribute
        self._conv = galsim.Convolve([self.original_gal, psf_inv], gsparams=gsparams)
        self._sbp = self._conv._sbp
        # self._sbp = self.original_gal._sbp
        if logger:
            logger.debug('MEDSGalaxy %d: Made gsobject', iobj)

        # TODO: Do we need this?
        # Save the noise in the image as an accessible attribute
        # self.noise = self.noise.convolvedWith(psf_inv, gsparams)

        if logger:
            logger.debug('MEDSGalaxy %d: Finished building MEDSGalaxy', iobj)

    # TODO: Can probably delete this!
    # @classmethod
    # def makeFromImage(cls, image, PSF, xi, **kwargs):
    #     """Create a RealGalaxy directly from image, PSF, and noise description.

    #     @param image  `Image` of the galaxy you want to simulate.
    #     @param PSF    `GSObject` representing the PSF of the galaxy image.  Note that this PSF
    #                   should include the response of the pixel convolution.
    #     @param xi     `CorrelatedNoise` or `UncorrelatedNoise` object characterizing the noise in
    #                   the input image.
    #     """
    #     noise_image = xi.drawImage()
    #     pixel_scale = noise_image.scale
    #     var = xi.getVariance()
    #     psf_image = PSF.drawImage(method='no_pixel')
    #     return RealGalaxy((image, psf_image, noise_image, pixel_scale, var))

    def __eq__(self, other):
        return (isinstance(other, MEDSGalaxy) and
                self.catalog == other.catalog and
                self.index == other.index and
                self._x_interpolant == other._x_interpolant and
                self._k_interpolant == other._k_interpolant and
                self._pad_factor == other._pad_factor and
                self._noise_pad_size == other._noise_pad_size and
                self._flux == other._flux and
                self._flux_rescale == other._flux_rescale and
                self._area_norm == other._area_norm and
                self._gsparams == other._gsparams)

    def __hash__(self):
        return hash(("galsim.MEDSGalaxy", self.catalog, self.index, self._x_interpolant,
                     self._k_interpolant, self._pad_factor, self._noise_pad_size, self._flux,
                     self._flux_rescale, self._area_norm, self._gsparams))

    def __repr__(self):
        s = 'galsim.MEDSGalaxy(%r, index=%r, '%(self.catalog, self.index)
        if self._x_interpolant is not None:
            s += 'x_interpolant=%r, '%self._x_interpolant
        if self._k_interpolant is not None:
            s += 'k_interpolant=%r, '%self._k_interpolant
        if self._pad_factor != 4:
            s += 'pad_factor=%r, '%self._pad_factor
        if self._noise_pad_size != 0:
            s += 'noise_pad_size=%r, '%self._noise_pad_size
        if self._flux is not None:
            s += 'flux=%r, '%self._flux
        if self._flux_rescale is not None:
            s += 'flux_rescale=%r, '%self._flux_rescale
        if self._area_norm != 1:
            s += 'area_norm=%r, '%self._area_norm
        s += 'rng=%r, '%self.rng
        s += 'gsparams=%r)'%self._gsparams
        return s

    def __str__(self):
        return 'balrog.MEDSGalaxy(index=%s)' % self.index

    def __getstate__(self):
        # The _sbp is picklable, but it is pretty inefficient, due to the large images being
        # written as a string.  Better to pickle the image and remake the InterpolatedImage.
        d = self.__dict__.copy()
        del d['_conv']
        del d['_sbp']
        return d

    def __setstate__(self, d):
        self.__dict__ = d
        psf_inv = galsim.Deconvolve(self.original_psf, gsparams=self._gsparams)
        self._conv = galsim.Convolve([self.original_gal, psf_inv], gsparams=self._gsparams)
        self._sbp = self._conv._sbp

class MEDSCatalog(object):
    """
    The MEDSCatalog class reads in and stores information about a specific training sample of
    realistic galaxies. We assume that all files containing the images (galaxies and PSFs) live in
    one directory; they could be individual files, or multiple HDUs of the same file.  Currently
    there is no functionality that lets this be a FITS data cube, because we assume that the object
    postage stamps will in general need to be different sizes depending on the galaxy size.

    Note that when simulating galaxies based on HST but using either realistic or parametric galaxy
    models, the COSMOSCatalog class may be more useful.  It allows the imposition of selection
    criteria and other subtleties that are more difficult to impose with RealGalaxyCatalog.

    While you could create your own catalog to use with this class, the typical use cases would
    be to use one of the catalogs that we have created and distributed.  There are three such
    catalogs currently, which can be use with one of the following initializations:

    1. A small example catalog is distributed with the GalSim distribution.  This catalog only
       has 100 galaxies, so it is not terribly useful as a representative galaxy population.
       But for simplistic use cases, it might be sufficient.  We use it for our unit tests and
       in some of the demo scripts (demo6, demo10, and demo11).  To use this catalog, you would
       initialize with

           >>> rgc = galsim.RealGalaxyCatalog('real_galaxy_catalog_23.5_example.fits',
                                              dir='path/to/GalSim/examples/data')

    2. There are two larger catalogs based on HST observations of the COSMOS field with around
       26,000 and 56,000 galaxies each with a limiting magnitude of F814W=23.5.  (The former is
       a subset of the latter.) For information about how to download these catalogs, see the
       RealGalaxy Data Download Page on the GalSim Wiki:

           https://github.com/GalSim-developers/GalSim/wiki/RealGalaxy%20Data

       Be warned that the catalogs are quite large.  The larger one is around 11 GB after unpacking
       the tarball.  To use one of these catalogs, you would initialize with

           >>> rgc = galsim.RealGalaxyCatalog('real_galaxy_catalog_23.5.fits',
                                              dir='path/to/download/directory')

    3. There is a catalog containing a random subsample of the HST COSMOS images with a limiting
       magnitude of F814W=25.2.  More information about downloading these catalogs can be found on
       the RealGalaxy Data Download page linked above.

    4. Finally, we provide a program that will download the large COSMOS sample for you and
       put it in the $PREFIX/share/galsim directory of your installation path.  The program is

           galsim_download_cosmos

       which gets installed in the $PREFIX/bin directory when you install GalSim.  If you use
       this program to download the COSMOS catalog, then you can use it with

           >>> rgc = galsim.RealGalaxyCatalog()

       GalSim knows the location of the installation share directory, so it will automatically
       look for it there.

    @param file_name  The file containing the catalog. [default: None, which will look for the
                      F814W<25.2 COSMOS catalog in $PREFIX/share/galsim.  It will raise an
                      exception if the catalog is not there telling you to run
                      galsim_download_cosmos.]
    @param sample     A keyword argument that can be used to specify the sample to use, i.e.,
                      "23.5" or "25.2".  At most one of `file_name` and `sample` should be
                      specified.
                      [default: None, which results in the same default as `file_name=None`.]
    @param dir        The directory containing the catalog, image, and noise files, or symlinks to
                      them. [default: None]
    @param preload    Whether to preload the header information.  If `preload=True`, the bulk of
                      the I/O time is in the constructor.  If `preload=False`, there is
                      approximately the same total I/O time (assuming you eventually use most of
                      the image files referenced in the catalog), but it is spread over the
                      various calls to getGalImage() and getPSFImage().  [default: False]
    @param logger     An optional logger object to log progress. [default: None]
    """
    _req_params = {'meds_files' : list, 'psf_files' : list, 'bands' : str}
    _opt_params = {'meds_dir' : str, 'psf_dir' : str, 'preload' : str,
                   'pixel_scale' : float, 'psf_hdu' : int, 'jacob_flip' : bool,
                   'snr_min' : float, 'snr_max' : float, 'flux_min' : float,
                   'flux_max' : float}
    _single_params = []
    _takes_rng = False

    _valid_bands = 'grizy'

    # _nobject_only is an intentionally undocumented kwarg that should be used only by
    # the config structure.  It indicates that all we care about is the nobjects parameter.
    # So skip any other calculations that might normally be necessary on construction.
    def __init__(self, meds_files, psf_files, bands, meds_dir=None, psf_dir=None, preload=True,
                 pixel_scale=0.263, psf_hdu=1, logger=None, _nobjects_only=False,
                 jacob_flip=False, snr_min=None, snr_max=None, flux_min=None, flux_max=None):

        bl = len(bands)
        if (bl != len(meds_files)) or (bl != len(psf_files)):
            raise ValueError('The number of passed bands (%s) does not match the number' % bl,
                             'of passed files!')
        self.bands = bands
        self.nbands = bl
        self.bindex = dict(zip('griz', range(self.nbands)))

        self.full_files = self._parse_files_dirs(meds_files, psf_files, meds_dir, psf_dir, bands)
        self.meds_files = self.full_files['MEDS']
        self.psf_files = self.full_files['PSF']
        # self.wcs_files = self.full_files['WCS']
        self.meds_dir = meds_dir
        self.psf_dir = psf_dir
        # self.wcs_dir = wcs_dir

        if (psf_hdu < 0) or (not isinstance(psf_hdu, int)):
            raise ValueError('`psf_hdu` must be a non-negative Int!')
        self.psf_hdu = psf_hdu

        if snr_min:
            if snr_min < 0.0:
                raise ValueError("The signal-to-noise ratio `snr` must be positive!")
            self.snr_min = snr_min
        else:
            self.snr_min = 0.0

        if snr_max:
            if snr_max < 0.0:
                raise ValueError("The signal-to-noise ratio `snr` must be positive!")
            elif snr_min and (snr_min > snr_max):
                raise ValueError("`snr_max` must be greater than `snr_min`!")
            self.snr_max = snr_max
        else:
            self.snr_max = None

        if flux_min:
            if flux_min < 0.0:
                raise ValueError("`flux_min` must be positive!")
            self.flux_min = flux_min
        else:
            self.flux_min = 0.0

        if flux_max:
            if flux_max < 0.0:
                raise ValueError("`flux_max` must be positive!")
            elif flux_min and (flux_min > flux_max):
                raise ValueError("`flux_max` must be greater than `flux_min`!")
            self.flux_max = flux_max
        else:
            self.flux_max = None

        self.jacob_flip = jacob_flip
        self._preload = preload

        self._setup_meds()

        # Get flags and create mask
        self._get_flags()
        self._make_mask()

        # Do mask cut
        self._mask_cut()

        # Exit early if that's all we needed. (set in `_setup_meds()`)
        if _nobjects_only: return

        # Already checked that IDs are consistent between bands, so choose the first one
        # Using 'ident' to be consistent with RealGalaxyCatalog
        b = self.bands[0]
        ident = self.cats[b]['id']

        # We want to make sure that the ident array contains all strings.
        # Strangely, ident.astype(str) produces a string with each element == '1'.
        # Hence this way of doing the conversion:
        self.ident = [ "%s"%val for val in ident ]

        # For now, not generalizing to noise files
        # TODO: Should we implement this in the future?
        self.noise_file_name = None

        self.pixel_scale = pixel_scale # Default of 0.263 is for DECam

        # TODO: Do we need this? self.variance = self.cat.field('noise_variance') # noise variance for image
        # self.band = self.cat.field('band') # bandpass in which apparent mag is measured, e.g., F814W

        # TODO: Do we need this?
        # Use uniform weighting in our case
        # self.weight = 1. / self.nobjects
        self._setup_weights()

        self.saved_noise_im = {}
        self.loaded_files = {}
        self.logger = logger

        # The pyfits commands aren't thread safe.  So we need to make sure the methods that
        # use pyfits are not run concurrently from multiple threads.
        from multiprocessing import Lock
        self.gal_lock = Lock()  # Use this when accessing gal files
        self.psf_lock = Lock()  # Use this when accessing psf files
        self.loaded_lock = Lock()  # Use this when opening new files from disk
        self.noise_lock = Lock()  # Use this for building the noise image(s) (usually just one)

        # Preload all files if desired
        if preload: self.preload()

        return

    def _parse_files_dirs(self, meds_files, psf_files, meds_dir, psf_dir, bands):

        # Convert `meds_files` and `psf_files` to dict's, since GalSim won't allow dict's as inputs
        meds_files = dict(zip(bands, meds_files))
        psf_files= dict(zip(bands, psf_files))

        # Check that the bands and files were inputted in the correct order (most files will
        # have band explicitly in filename as `_band_`)
        for ftype, files in {'MEDS':meds_files, 'PSF':psf_files}.items():
            for b in bands:
                if '_{}_'.format(b) not in files[b]:
                    for ob in bands.replace(b, ''):
                        if '_{}_'.format(ob) in files[band]:
                            raise IOError('''Inconsistent matching between inputted bands and {} files!
                                        '\'_{}_\' found in filename for band {}'''.format(ftype, ob, b))
                    warnings.warn('No \'_{}_\' string was found in {} filename for band {},'.format(b,b),
                                'but no explicit identifier for another band was found.',
                                'Please ensure that bands and input files are listed consistently.')

        valid_bands = MEDSCatalog._valid_bands
        full_meds_files = {}
        full_psf_files = {}
        # full_wcs_files = {}
        dirs = {}
        bands = []
        files = {'MEDS': meds_files, 'PSF': psf_files} #, 'WCS': wcs_files}
        full_files = {'MEDS': full_meds_files, 'PSF': full_psf_files} #, 'WCS': full_wcs_files}

        if (len(meds_files) != len(psf_files)): # or (len(psf_files) != len(wcs_files)):
            raise Exception('Must pass equal number of MEDS and PSF files!')

        for ftype, d in {'MEDS':meds_dir, 'PSF':psf_dir}.items(): #'WCS':wcs_dir}.items():

            if d is None:
                dirs[ftype] = ''
            else:
                if not os.path.isdir(d):
                    raise IOError('Passed %s dir of %s does not exist!' % (ftype, d) )
                dirs[ftype] = d

            for band, f in files[ftype].items():
                if band not in valid_bands:
                    raise ValueError('Passed band %s is not valid ' %band +
                                    '- only %s are currently accepted.' %valid_bands)
                bands.append(band)
                full_file = os.path.abspath(os.path.join(dirs[ftype], f))

                if not os.path.isfile(full_file):
                    raise IOError('%s file %s does not exist!' % (ftype, full_file))

                full_files[ftype][band] = full_file

                # if ftype is 'MEDS':
                #     full_meds_files[band] = full_file
                # elif ftype is 'PSF':
                #     full_psf_files[band] = full_file
                # elif ftype is 'WCS':
                #     full_wcs_files[band] = full_file

        return full_files

    def __del__(self):
        # Make sure to clean up pyfits open files if people forget to call close()
        self.close()

    def close(self):
        # Need to close any open files.
        # Make sure to check if loaded_files exists, since the constructor could abort
        # before it gets to the place where loaded_files is built.
        if hasattr(self, 'loaded_files'):
            for f in self.loaded_files.values():
                f.close()
        self.loaded_files = {}

    def getNObjects(self) : return self.nobjects

    def __len__(self): return self.nobjects

    def getMEDSFiles(self): return self.meds_files

    def getPSFFiles(self): return self.psf_files

    def getCatalog(self): return self.cats

    def getIobjForID(self, id):
        """Internal function to find which object index corresponds to the value
        ID in the ident field.
        """
        # Just to be completely consistent, convert id to a string in the same way we
        # did above for the ident array:
        id = "%s"%id
        if id in self.ident:
            return self.ident.index(id)
        else:
            raise ValueError('ID not found in list of IDs',id, self.ident)

    def _setup_meds(self):
        self.meds = {}
        self.cats = {}
        size = None

        for band in self.bands:
            m = meds.MEDS(self.meds_files[band])
            self.cats[band] = m.get_cat()
            if size:
                if m.size != size:
                    raise ValueError('The passed MEDS files do not have equal sizes!')
            size = m.size

            if self.preload:
                self.meds[band] = m
            else:
                self.meds[band] = None

            # TODO: Are we able to close m here?
            # m.close()

        # The input logger needs to know the original catalog size
        self.orig_nobjects = size

        # Object indices in original meds catalog
        self.orig_index = np.arange(self.orig_nobjects)

        # Check that all MEDS objects are consistent with one another
        b = self.bands[0]
        for band in self.bands:
            if (self.cats[b]['id'] != self.cats[band]['id']).any():
                raise Exception('The MEDS files IDs are not be consistent with one another!')

        return

    def _get_flags(self):
        # self.flags = {}
        # self.obj_flags = {}

        # for band in self.bands:
        #     # General flags
        #     self.flags[band] = self.cats[band]['flags']
        #     # self.obj_flags = self.catalog['obj_flags']

        # self.flags = np.ones((self.nbands, self.orig_nobjects))

        self.flags = np.array([self.cats[b]['flags'] for b in self.bands])

        return

    def _make_mask(self):

        # ones = np.ones(self.orig_nobjects, dtype=bool)
        # base = self.orig_nobjects * [True]
        # mask = zip(self.bands, lb*base)
        # mask = dict.fromkeys(self.bands, np.ones(self.orig_nobjects, dtype=bool))
        mask = np.ones((self.nbands, self.orig_nobjects), dtype=bool)

        for b, i in self.bindex.items():
            mask[i][self.flags[i] != 0] = False

            # Remove objects with flux_min < flux < flux_max
            flux = self.cats[b]['flux']
            mask[i][flux < self.flux_min] = False
            if self.flux_max:
                mask[i][flux > self.flux_max] = False

            # Remove objects with snr_min < S/N < snr_max
            snr = self.cats[b]['flux'] / self.cats[b]['flux_err']
            mask[i][snr < self.snr_min] = False
            if self.snr_max:
                mask[i][snr > self.snr_max] = False

        # for band in self.bands:
        #     # For now, remove objects with any flags present
        #     mask[band][self.flags[band] != 0] = False
        #     # mask[self.obj_flags !=0] = False
        #     # Extra flags for 'mof' catalogs

        #     # # Remove objects with snr_min < S/N < snr_max
        #     # mask[self.catalog[cp+'_s2n_r'] < self.snr_min] = False
        #     # if self.snr_max:
        #     #     mask[self.catalog[cp+'_s2n_r'] > self.snr_max] = False

        # other masks....

        self.mask = np.bitwise_or.reduce(mask, 0)

        return

    def _mask_cut(self):
        """Do mask cut defined in `make_mask()`."""

        for band in self.bands: self.cats[band] = self.cats[band][self.mask]

        self.orig_index = self.orig_index[self.mask]
        self.nobjects = len(self.orig_index)

        print('Ntotal: {}\nNobjects: {}'.format(self.orig_nobjects, self.nobjects))

        return

    def _setup_weights(self):
        # TODO: We may want to assign weights in the MEDS catalog metadata. For now, do nothing.
        # (Everything still work if self.weights is not set)
        pass

    # def _setup_psfs(self):
    #     self.psfs = {}

    #     for band in self.bands:
    #         self.psfs[band] = galsim.des.DES_PSFEx(self.psf_files[band], self.wcs_files[band])

    #     return

    # def preload(self):
    #         """Preload the files into memory.
    #         There are memory implications to this, so we don't do this by default.  However, it can be
    #         a big speedup if memory isn't an issue.
    #         """
    #         from ._pyfits import pyfits
    #         self.logger.debug('MEDSCatalog: start preload')
    #         for file_name in np.concatenate((self.gal_file_name , self.psf_file_name)):
    #             # numpy sometimes add a space at the end of the string that is not present in
    #             # the original file.  Stupid.  But this next line removes it.
    #             file_name = file_name.strip()
    #             if file_name not in self.loaded_files:
    #                 self.logger.debug('MEDSCatalog: preloading %s',file_name)
    #                 # I use memmap=False, because I was getting problems with running out of
    #                 # file handles in the great3 real_gal run, which uses a lot of rgc files.
    #                 # I think there must be a bug in pyfits that leaves file handles open somewhere
    #                 # when memmap = True.  Anyway, I don't know what the performance implications
    #                 # are (since I couldn't finish the run with the default memmap=True), but I
    #                 # don't think there is much impact either way with memory mapping in our case.
    #                 f = pyfits.open(file_name,memmap=False)
    #                 self.loaded_files[file_name] = f
    #                 # Access all the data from all hdus to force PyFits to read the data
    #                 for hdu in f:hdu.data

    def preload(self):
        """Preload the files into memory.

        There are memory implications to this, so we don't do this by default.  However, it can be
        a big speedup if memory isn't an issue.
        """

        # TODO: Implement!
        print('Warning - `preload` hasn\'t been implemented yet for MEDSCatalog.')
        return

    def _getFile(self, file_name, meds=False):
        from galsim._pyfits import pyfits
        if file_name in self.loaded_files:
            if self.logger:
                self.logger.debug('MEDSCatalog: File %s is already open',file_name)
            f = self.loaded_files[file_name]
        else:
            self.loaded_lock.acquire()
            # Check again in case two processes both hit the else at the same time.
            if file_name in self.loaded_files: # pragma: no cover
                if self.logger:
                    self.logger.debug('MEDSCatalog: File %s is already open',file_name)
                f = self.loaded_files[file_name]
            else:
                if self.logger:
                    self.logger.debug('MEDSCatalog: open file %s',file_name)
                if meds:
                    f = meds.MEDS(file_name)
                else:
                    f = pyfits.open(file_name,memmap=False)
                self.loaded_files[file_name] = f
            self.loaded_lock.release()
        return f

    def getBands(self): return self.bands

    def getGalImage(self, iobj, icutout, band):
        """
        Returns the galaxy at the given index and band as an Image object.
        """
        if self.logger:
            self.logger.debug('MEDSCatalog %d: Start getGalImage', iobj)
        if iobj >= self.nobjects:
            raise IndexError(
                'index %d given to getGalImage is out of range (0..%d)'
                % (iobj,self.nobjects-1))


        # NOTE: iobj is the index of the masked-catalog. However, it is
        # not trivial to mask the MEDS file. So we instead grab the original
        # index given iobj:
        iobj_orig = self.orig_index[iobj]

        if self.preload:
            array = self.meds[band].get_cutout(iobj_orig, icutout)
        else:
            f = self._getFile(self.meds_files[band], meds=True)
            # For some reason the more elegant `with gal_lock:` syntax isn't working for me.
            # It gives an EOFError.  But doing an explicit acquire and release seems to work fine.
            self.gal_lock.acquire()
            array = f.get_cutout(iobj_orig, icutout)
            self.gal_lock.release()

        im = galsim.Image(np.ascontiguousarray(array.astype(np.float64)),
                        scale=self.pixel_scale)
        return im

    def getPSF(self, iobj, icutout, band, wcs=None):
        """Returns the PSF at index `i` for a given `band` as a GSObject.
        Gets the PSF solution of the original MEDS images.
        """
        if self.logger:
            self.logger.debug('MEDSCatalog %d: Start getPSF', iobj)
        if iobj >= self.nobjects:
            raise IndexError(
                'index %d given to getPSF is out of range (0..%d)'
                % (iobj,self.nobjects-1))

        if wcs is None: wcs = self.getWCS(iobj, icutout, band)

        return galsim.des.DES_PSFEx(self.psf_files[band], wcs=wcs)

    def getPSFImage(self, iobj, icutout, band):
        """Returns the PSF of a given band as an Image object.
        """

        # Grab WCS and centroid position in image coordinates
        # for object
        wcs = self.getWCS(iobj, icutout, band)
        im_pos = wcs.origin

        psf = self.getPSF(iobj, icutout, band, wcs=wcs)
        pobj = psf.getPSF(im_pos)

        return pobj.drawImage()

    def getWCS(self, iobj, icutout, band):
        """Returns the WCS of an object for a given band and cutout.
        """

        # NOTE: iobj is the index of the masked-catalog. However, it is
        # not trivial to mask the MEDS file. So we instead grab the original
        # index given iobj:
        iobj_orig = self.orig_index[iobj]

        if self.preload:
            j = self.meds[band].get_jacobian(iobj_orig, icutout)
            x = self.meds[band][iobj_orig]['orig_col'][0]
            y = self.meds[band][iobj_orig]['orig_row'][0]
        else:
            f = self._getFile(self.psf_files[band], meds=True)
            # For some reason the more elegant `with psf_lock:` syntax isn't working for me.
            # It gives an EOFError.  But doing an explicit acquire and release seems to work fine.
            self.psf_lock.acquire()
            j = f.get_jacobian(iobj_orig, icutout)
            x = f[iobj_orig]['orig_col'][0]
            y = f[iobj_orig]['orig_row'][0]
            self.psf_lock.release()

        # x, y = j['col0'], j['row0'] # This is the im position *in the frame of the PS*
        dudx, dudy = j['dudcol'], j['dudrow']
        dvdx, dvdy = j['dvdcol'], j['dvdrow']

        # There was a bug in original version of MEDS maker that incorrectly flipped the
        # (x,y) image coordinates of the WCS origin. Will account for this if `jacob_flip`
        # is set to True
        if self.jacob_flip: x, y = y, x

        return galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=galsim.PositionD(x, y),
                                        world_origin=galsim.PositionD(0, 0))

    # def getPSFImage(self, ident, band, im_pos):
    #     """Returns the PSF of a given band as an Image object.
    #     Position `pos` expected in (TODO:??) coordinates
    #     """
    #     psf = getPSF(ident, band, im_pos)
    #     return psf.drawImage()

    # def getPSF(self, ident, band, im_pos):
    #     """Returns the PSF at index `i` for a given `band` as a GSObject.
    #     Gets the PSF solution of the original MEDS images.
    #     """
    #     if self.logger:
    #         self.logger.debug('MEDSCatalog %d: Start getPSFImage', ident)
    #     if iobj >= self.nobjects:
    #         raise IndexError(
    #             'index %d given to getPSFImage is out of range (0..%d)'
    #             % (iobj,self.nobjects-1))

    #     return self.psfs[band].getPSF(im_pos)

    def getNoiseProperties(self, i):
        """Returns the components needed to make the noise correlation function at index `i`.
           Specifically, the noise image (or None), the pixel_scale, and the noise variance,
           as a tuple (im, scale, var).
        """

        # TODO: implement!
        pass

        # if self.logger:
        #     self.logger.debug('RealGalaxyCatalog %d: Start getNoise',i)
        # if self.noise_file_name is None:
        #     im = None
        # else:
        #     if i >= len(self.noise_file_name):
        #         raise IndexError(
        #             'index %d given to getNoise is out of range (0..%d)'%(
        #                 i,len(self.noise_file_name)-1))
        #     if self.noise_file_name[i] in self.saved_noise_im:
        #         im = self.saved_noise_im[self.noise_file_name[i]]
        #         if self.logger:
        #             self.logger.debug('RealGalaxyCatalog %d: Got saved noise im',i)
        #     else:
        #         self.noise_lock.acquire()
        #         # Again, a second check in case two processes get here at the same time.
        #         if self.noise_file_name[i] in self.saved_noise_im:
        #             im = self.saved_noise_im[self.noise_file_name[i]]
        #             if self.logger:
        #                 self.logger.debug('RealGalaxyCatalog %d: Got saved noise im',i)
        #         else:
        #             from galsim._pyfits import pyfits
        #             with pyfits.open(self.noise_file_name[i]) as fits:
        #                 array = fits[0].data
        #             im = galsim.Image(np.ascontiguousarray(array.astype(np.float64)),
        #                               scale=self.pixel_scale[i])
        #             self.saved_noise_im[self.noise_file_name[i]] = im
        #             if self.logger:
        #                 self.logger.debug('RealGalaxyCatalog %d: Built noise im',i)
        #         self.noise_lock.release()

        # return im, self.pixel_scale[i], self.variance[i]

    def getNoise(self, i, rng=None, gsparams=None):
        """Returns the noise correlation function at index `i` as a CorrelatedNoise object.
           Note: the return value from this function is not picklable, so this cannot be used
           across processes.
        """

        # TODO: Implement!
        pass

        # im, scale, var = self.getNoiseProperties(i)
        # if im is None:
        #     cf = galsim.UncorrelatedNoise(var, rng=rng, scale=scale, gsparams=gsparams)
        # else:
        #     ii = galsim.InterpolatedImage(im, normalization="sb",
        #                                   calculate_stepk=False, calculate_maxk=False,
        #                                   x_interpolant='linear', gsparams=gsparams)
        #     cf = galsim.correlatednoise._BaseCorrelatedNoise(rng, ii, im.wcs)
        #     cf = cf.withVariance(var)
        # return cf

    def __repr__(self):
        return 'balrog.MEDSCatalog(%r)' % self.meds_files

    def __eq__(self, other):
        return (isinstance(other, MEDSCatalog) and
                self.meds_files == other.file_name and
                self.psf_files == other.psf_files and
                self.wcs_files == other.wcs_files and
                self.meds_dir == other.meds_dir and
                self.psf_dir == other.psf_dir and
                self.wcs_dir == other.wcs_dir)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(repr(self))

    def __getstate__(self):
        d = self.__dict__.copy()
        d['loaded_files'] = {}
        d['saved_noise_im'] = {}
        del d['gal_lock']
        del d['psf_lock']
        del d['loaded_lock']
        del d['noise_lock']
        return d

    def __setstate__(self, d):
        from multiprocessing import Lock
        self.__dict__ = d
        self.gal_lock = Lock()
        self.psf_lock = Lock()
        self.loaded_lock = Lock()
        self.noise_lock = Lock()
        pass

#####------------------------------------------------------------------------------------------------

from galsim.config.input import RegisterInputType, InputLoader

class MEDSCatalogLoader(galsim.config.InputLoader):
    """ The MEDSCatalog loader doesn't need anything special other than registration as a valid input type.
        These additions are only used for logging purposes.
    """

    def setupImage(self, MEDS_catalog, config, base, logger):
        # This method is blank for a general InputLoader, and a convenient place to put the logger
        if logger: # pragma: no cover
            # Only report as a warning the first time.  After that, use info.
            first = not base.get('_MEDSCatalogLoader_reported_as_warning',False)
            base['_MEDSCatalogLoader_reported_as_warning'] = True
            if first:
                log_level = logging.WARNING
            else:
                log_level = logging.INFO
            if 'input' in base:
                if 'meds_catalog' in base['input']:
                    cc = base['input']['meds_catalog']
                    if isinstance(cc,list): cc = cc[0]
                    out_str = ''
                    if 'meds_dir' in cc:
                        out_str += '\n  meds_dir = %s'%cc['meds_dir']
                    if 'meds_files' in cc:
                        out_str += '\n  meds_files = %s'%cc['meds_files']
                    if 'psf_dir' in cc:
                        out_str += '\n  psf_dir = %s'%cc['psf_dir']
                    if 'psf_files' in cc:
                        out_str += '\n  psf_files = %s'%cc['psf_files']
                    if out_str != '':
                        logger.log(log_level, 'Using user-specified MEDSCatalog: %s', out_str)
            logger.info("file %d: MEDS catalog has %d total objects.", base['file_num'],
                        MEDS_catalog.getNObjects())

# Need to add the MEDSCatalog class as a valid input_type.
galsim.config.RegisterInputType('meds_catalog', MEDSCatalogLoader(MEDSCatalog, has_nobj=True))

# There are two gsobject types that are coupled to this: MEDSGalaxy and MEDSGalaxyOriginal.

def _BuildMEDSGalaxy(config, base, ignore, gsparams, logger, param_name='MEDSGalaxy'):
    """@brief Build a MEDSGalaxy from the meds_catalog input item.
    """
    meds_cat = galsim.config.GetInputObj('meds_catalog', config, base, param_name)

    # Special: if index is Sequence or Random, and max isn't set, set it to nobjects-1.
    # But not if they specify 'id' or have 'random=True', which overrides that.
    if 'id' not in config:
        if 'random' not in config:
            galsim.config.SetDefaultIndex(config, meds_cat.getNObjects())
        else:
            if not config['random']:
                galsim.config.SetDefaultIndex(config, meds_cat.getNObjects())
                # Need to do this to avoid being caught by the GetAllParams() call, which will flag
                # it if it has 'index' and 'random' set (but 'random' is False, so medsly it's OK).
                del config['random']

    kwargs, safe = galsim.config.GetAllParams(config, base,
        req = MEDSGalaxy._req_params,
        opt = MEDSGalaxy._opt_params,
        single = MEDSGalaxy._single_params,
        ignore = ignore + ['num'])
    if gsparams: kwargs['gsparams'] = galsim.GSParams(**gsparams)

    kwargs['rng'] = galsim.config.GetRNG(config, base, logger, param_name)

    if 'index' in kwargs:
        index = kwargs['index']
        if index >= meds_cat.getNObjects() or index < 0:
            raise ValueError(
                "index=%s has gone past the number of entries in the MEDSCatalog"%index)
            # raise galsim.errors.GalSimConfigError(
            #     "index=%s has gone past the number of entries in the MEDSCatalog"%index)

    kwargs['meds_catalog'] = meds_cat
    logger.debug('obj %d: %s kwargs = %s',base.get('obj_num',0),param_name,kwargs)

    gal = MEDSGalaxy(**kwargs)
    # pudb.set_trace()

    return gal, safe


def _BuildMEDSGalaxyOriginal(config, base, ignore, gsparams, logger):
    """@brief Return the original image from a MEDSGalaxy using the real_catalog input item.
    """
    gal, safe = _BuildMEDSGalaxy(config, base, ignore, gsparams, logger,
                                   param_name='MEDSGalaxyOriginal')
    return gal.original_gal, safe


# Register these as valid gsobject types
RegisterObjectType('MEDSGalaxy', _BuildMEDSGalaxy, input_type='meds_catalog')
RegisterObjectType('MEDSGalaxyOriginal', _BuildMEDSGalaxyOriginal, input_type='meds_catalog')

def simReal(real_galaxy, target_PSF, target_pixel_scale, g1=0.0, g2=0.0, rotation_angle=None,
            rand_rotate=True, rng=None, target_flux=1000.0, image=None): # pragma: no cover
    """Deprecated method to simulate images (no added noise) from real galaxy training data.

    This function takes a RealGalaxy from some training set, and manipulates it as needed to
    simulate a (no-noise-added) image from some lower-resolution telescope.  It thus requires a
    target PSF (which could be an image, or one of our base classes) that represents all PSF
    components including the pixel response, and a target pixel scale.

    The default rotation option is to impose a random rotation to make irrelevant any real shears
    in the galaxy training data (optionally, the RNG can be supplied).  This default can be turned
    off by setting `rand_rotate = False` or by requesting a specific rotation angle using the
    `rotation_angle` keyword, in which case `rand_rotate` is ignored.

    Optionally, the user can specify a shear (default 0).  Finally, they can specify a flux
    normalization for the final image, default 1000.

    @param real_galaxy      The RealGalaxy object to use, not modified in generating the
                            simulated image.
    @param target_PSF       The target PSF, either one of our base classes or an Image.
    @param target_pixel_scale  The pixel scale for the final image, in arcsec.
    @param g1               First component of shear to impose (components defined with respect
                            to pixel coordinates), [default: 0]
    @param g2               Second component of shear to impose, [default: 0]
    @param rotation_angle   Angle by which to rotate the galaxy (must be an Angle
                            instance). [default: None]
    @param rand_rotate      Should the galaxy be rotated by some random angle?  [default: True;
                            unless `rotation_angle` is set, then False]
    @param rng              A BaseDeviate instance to use for the random selection or rotation
                            angle. [default: None]
    @param target_flux      The target flux in the output galaxy image, [default: 1000.]
    @param image            As with the GSObject.drawImage() function, if an image is provided,
                            then it will be used and returned.  [default: None, which means an
                            appropriately-sized image will be created.]

    @return a simulated galaxy image.
    """
    from .deprecated import depr
    depr('simReal', 1.5, '',
         'This method has been deprecated due to lack of widespread use.  If you '+
         'have a need for it, please open an issue requesting that it be reinstated.')
    # do some checking of arguments
    if not isinstance(real_galaxy, galsim.RealGalaxy):
        raise RuntimeError("Error: simReal requires a RealGalaxy!")
    if isinstance(target_PSF, galsim.Image):
        target_PSF = galsim.InterpolatedImage(target_PSF, scale=target_pixel_scale)
    if not isinstance(target_PSF, galsim.GSObject):
        raise RuntimeError("Error: target PSF is not an Image or GSObject!")
    if rotation_angle is not None and not isinstance(rotation_angle, galsim.Angle):
        raise RuntimeError("Error: specified rotation angle is not an Angle instance!")
    if (target_pixel_scale < real_galaxy.pixel_scale):
        import warnings
        message = "Warning: requested pixel scale is higher resolution than original!"
        warnings.warn(message)
    import math # needed for pi, sqrt below
    g = math.sqrt(g1**2 + g2**2)
    if g > 1:
        raise RuntimeError("Error: requested shear is >1!")

    # make sure target PSF is normalized
    target_PSF = target_PSF.withFlux(1.0)

    # rotate
    if rotation_angle is not None:
        real_galaxy = real_galaxy.rotate(rotation_angle)
    elif rotation_angle is None and rand_rotate:
        ud = galsim.UniformDeviate(rng)
        rand_angle = galsim.Angle(math.pi*ud(), galsim.radians)
        real_galaxy = real_galaxy.rotate(rand_angle)

    # set fluxes
    real_galaxy = real_galaxy.withFlux(target_flux)

    # shear
    if (g1 != 0.0 or g2 != 0.0):
        real_galaxy = real_galaxy.shear(g1=g1, g2=g2)

    # convolve, resample
    out_gal = galsim.Convolve([real_galaxy, target_PSF])
    image = out_gal.drawImage(image=image, scale=target_pixel_scale, method='no_pixel')

    # return simulated image
    return image


