#----------------------------------------------------------------------------
# Classes for new input type `ngmixCatalog`; intended for use with Balrog.
#
# Contributors:
# Spencer Everett (UCSC)
#----------------------------------------------------------------------------

import galsim
import galsim.config
import ngmix # Used for GMix -> GSObject conversion
import numpy as np
import logging
import warnings
from past.builtins import basestring # Python 2&3 compatibility

# TODO: Include noise, pixscale

class ngmixCatalog(object):
    # TODO: Update documentation w/ new parameters!
    """ Class that handles galaxy catalogs from ngmix. These are fits files with names typically
    of the form 'DES{####}{+/-}{####}-y{#}v{#}-{type}-{###}.fits'.

    `ngmix` is software written by Erin Sheldon.  If you want more detail about it,
    check out the github repo:

        https://github.com/esheldon/ngmix

    Much of this class as well as its corresponding loader/builder are designed by inspection of
    `des_psfex.py` and `scene.py`. Credit to the GalSim team.

    @param file_name       The file name to be read in, or a pyfits HDU in which case it is used
                           directly instead of being opened.
    @param dir             Optionally a directory name can be provided if the file_name does not
                           already include it.  (The image file is assumed to be in the same
                           directory.) (Default `dir = None`).  Cannot pass an HDU with this option.
    @param catalog_type    The type of the input ngmix catalog. Only those in `valid_catalog_types`
                           are currently supported. If none is passed, the type is attempted to be
                           inferred from the filename.
    @param bands           A string of the desired bands to simulate from (only griz allowed). For
                           example, selecting only the 'g' and 'r' bands would be done by setting
                           bands='gr'. If none are passed, the g-band is selected by default.
    @param snr_min         The lower allowed bound for signal-to-noise ratio (snr). Can be any
                           positive value, as long as it is smaller than `snr_max`. All objects with
                           negative snr are removed by default.
    @param snr_max         The upper allowed bound for snr. Unlikely to be used very often, but
                           included for completeness.
    @param t_frac          The cutoff used for object size (T) / object size error (T_err). All
                           objects below this cutoff will be removed. (Default: 0).
    @param _nobjects_only  This is only passed if GalSim wants to know how many input objects will
                           be used without processing the whole input catalog.
    """

    _req_params = { 'file_name' : str, 'bands' : str}
    _opt_params = { 'dir' : str, 'catalog_type' : str, 'snr_min' : float, 'snr_max' : float,
                    't_frac' : float, 't_min' : float, 't_max' : float, 'version' : str,
                    'de_redden' : bool, 'TdByTe' : float}
    _single_params = []
    _takes_rng = False

    # Only these ngmix catalog types currently supported
    # `gauss`: Single Gaussian fit
    # `cm`: Combined bulge+disk model
    # `bdf`: CM model with buldge-disk ratio fixed
    # NOTE: In principle, should be able to handle any type supported by ngmix
    # See Issue #79
    _valid_catalog_types = ['gauss','cm','bdf']

    # Only these color bands are currently supported for an ngmix catalog
    _valid_band_types = 'griz'

    # Dictionary of color band flux to array index in ngmix catalogs
    # TODO: This should be grabbed from the fits header rather than be hardcoded
    # See Issue #80
    _band_index = {'g' : 0, 'r' : 1, 'i' : 2, 'z' : 3}

    # NOTE: In previous versions, the catalog column name prefix didn't always
    # match the catalog type (e.g. 'mof' had a prefix of 'cm' for most columns).
    # This shouldn't be needed in the future but leaving for now
    _cat_col_prefix = {'gauss' : 'gauss', 'cm' : 'cm', 'bdf' : 'bdf'}

    def __init__(self, file_name, bands, dir=None, catalog_type=None, snr_min=None, snr_max=None,
                 t_frac=None, t_min=None, t_max=None, version=None, de_redden=False, TdByTe=None,
                 _nobjects_only=False):

        if dir:
            if not isinstance(file_name, basestring):
                raise ValueError("Cannot provide dir and an HDU instance!")
            import os
            file_name = os.path.join(dir,file_name)
        if not isinstance(file_name, basestring):
            raise ValueError("The input filename must be a string!")
        self.file_name = file_name

        if catalog_type is not None:
            if catalog_type in self._valid_catalog_types:
                if catalog_type not in self.file_name:
                    print("Inputted ngmix catalog type of `{}` ".format(catalog_type) +
                          "does not match filename, which is standard "
                          "for DES ngmix catalogs. Ensure this is correct.")
                self.cat_type = catalog_type
            else:
                raise ValueError("{} is not a currently supported ".format(catalog_type) +
                                 "ngmix catalog type!")
        else:
            # Attempt to determine catalog type from filename (this is generally true for DES ngmix
            # catalogs)
            match = 0
            for t in self._valid_catalog_types:
                if t in self.file_name:
                    match +=1
                    self.cat_type = t
            # Reject if multiple matches
            if match == 0:
                raise ValueError("No inputted ngmix catalog type, and no matches in filename! "
                                 "Please set a valid catalog type.")
            if match > 1:
                raise ValueError("No inputted ngmix catalog type, and multiple matches in filename!"
                                 " Please set a valid catalog type.")

        if TdByTe is not None:
            if self.cat_type != 'bdf':
                raise ValueError('Can only set a constant `TdByTe` for ngmix type `bdf`!')
            if TdByTe < 0:
                raise ValueError('TdByTe must be non-negative!')
        else:
            # This should almost always be 1
            TdByTe = 1.
        self._TdByTe = TdByTe

        # Catalog column name prefixes don't always match catalog type
        self.col_prefix = self._cat_col_prefix[self.cat_type]

        if isinstance(bands, basestring):
            # Strip all whitespace
            bands = bands.replace(" ", "")
            # More useful as a list of individual bands
            bands_list = list(bands)
            if set(bands_list).issubset(self._valid_band_types):
                self.bands = bands_list
            else:
                raise ValueError("The only valid color bands for a ngmix catalog are "
                                 "{}!".format(self._valid_band_types))
        else:
            # TODO: Wouldn't be a bad idea to allow a list of individual bands as well
            raise ValueError("Must enter desired color bands as a string! "
                             "(For example, `bands : \'gr\'`)")

        if snr_min is not None:
            if snr_min < 0.0:
                raise ValueError("The signal-to-noise ratio `snr` must be positive!")
            self.snr_min = snr_min
        else:
            self.snr_min = None

        if snr_max is not None:
            if snr_max < 0.0:
                raise ValueError("The signal-to-noise ratio `snr` must be positive!")
            elif snr_min and (snr_min > snr_max):
                raise ValueError("`snr_max` must be greater than `snr_min`!")
            self.snr_max = snr_max
        else:
            self.snr_max = None

        if t_min is not None:
            self.t_min = t_min
        else:
            self.t_min = None

        if t_max is not None:
            if (t_min is not None) and (t_min > t_max):
                raise ValueError("`t_max` must be greater than `t_min`!")
            self.t_max = t_max
        else:
            self.t_max = None

        if t_frac is not None:
            if t_frac < 0.0:
                raise ValueError("The allowed size/size_err fraction `t_frac` must be positive!")
            self.t_frac = t_frac
        else:
            self.t_frac = None

        if version is not None:
            if not isinstance(version, basestring):
                raise TypeError('`version` must be a string!`')
        self.version = version

        self.de_redden = de_redden

        self.read()

        return

    #----------------------------------------------------------------------------------------------

    def read(self):
        """Read in relevant catalog information"""

        from galsim._pyfits import pyfits

        if isinstance(self.file_name, basestring):
            # If a filename is passed:
            hdu_list = pyfits.open(self.file_name)
            model_fits = hdu_list[1]
        else:
            # If a fits HDU is directly passed:
            hdu_list = None
            model_fits = self.file_name

        self.catalog = model_fits.data

        # NB: As discussed in `scene.py`, there is a bug in the pyfits FITS_Rec class
        # that leads to memory leaks. The simplest workaround seems to be to convert
        # it to a regular numpy recarray.
        self.catalog = np.array(self.catalog, copy=True)

        # The input logger needs to know the original catalog size
        self.ntotal = len(self.catalog)

        if hdu_list: hdu_list.close()

        # Galaxy indices in original ngmix catalog
        self.orig_index = np.arange(self.ntotal)

        self.getFlags()
        self.makeMask()
        self.maskCut()

        if self.de_redden is True:
            self._check_reddening_factors()

        return

    def _check_reddening_factors(self):
        # TODO: Generalize in future, but we want something simple for Y3 Balrog
        cp = self.col_prefix

        # Might need more in future...
        # TODO: Make this consistent with whatever Brian decides!
        req_colnames = [cp+'_flux_deredden']

        for colname in req_colnames:
            if colname not in self.catalog.dtype.names:
                raise AttributeError('The column `{}` must be present in '.format(colname) +
                                     'the ngmix catalog if for de-reddening!')

        return

    #----------------------------------------------------------------------------------------------

    def getFlags(self):
        """Retrieve object flags, where implementation depends on catalog type."""

        # General flags
        self.flags = self.catalog['flags']

        # We don't want to cut on these explicitly anymore:

        #self.obj_flags = self.catalog['obj_flags']

        # ngmix catalog-specific flags
        #self.ngmix_flags = self.catalog[self.col_prefix+'_flags']

        # if self.cat_type == 'mof':
        #     # mof has additional flags
        #     self.mof_flags = self.catalog[self.col_prefix+'_mof_flags']

        return

    #----------------------------------------------------------------------------------------------

    def makeMask(self):
        """Add a masking procedure, if desired."""
        # TODO: Allow multiple masking procedures
        # TODO works as currently written, but should rewrite to conform to conventional masking
        # definition

        cp = self.col_prefix

        mask = np.ones(len(self.orig_index), dtype=bool)

        # For now, remove objects with any flags present
        mask[self.flags != 0] = False

        # No longer do these explicitly:
        # mask[self.obj_flags !=0] = False
        # mask[self.ngmix_flags !=0] = False
        # Extra flags for 'mof' catalogs
        # if self.cat_type == 'mof':
        #     mask[self.mof_flags != 0] = False

        # Remove any object with `|T/T_err|` < t_frac
        # (have to do absolute value as T can be negative)
        if self.t_frac is not None:
            T_fraction = self.catalog[cp+'_T'] / self.catalog[cp+'_T_err']
            mask[abs(T_fraction) < self.t_frac ] = False

        # Remove objects with snr_min < S/N < snr_max
        if self.snr_min is not None:
            mask[self.catalog[cp+'_s2n_r'] < self.snr_min] = False
        if self.snr_max is not None:
            mask[self.catalog[cp+'_s2n_r'] > self.snr_max] = False

        # Remove objects with size T outside of desired bounds
        if self.t_min is not None:
            mask[self.catalog[cp+'_T'] < self.t_min] = False
        if self.t_max is not None:
            mask[self.catalog[cp+'_T'] > self.t_max] = False

        self.mask = mask

        return

    #----------------------------------------------------------------------------------------------

    def maskCut(self):
        """Do mask cut defined in `makeMask()`."""
        self.catalog = self.catalog[self.mask]
        self.orig_index = self.orig_index[self.mask]
        self.nobjects = len(self.orig_index)

        return

    #----------------------------------------------------------------------------------------------

    def makeGalaxies(self, band, index=None, n_random=None, rng=None, gsparams=None):
        """
        Construct GSObjects from a list of galaxies in the ngmix catalog specified by `index`
        (or a randomly generated one). This is done using Erin's code.

        @param index            Index of the desired galaxy in the catalog for which a GSObject
                                should be constructed.  You can also provide a list or array of
                                indices, in which case a list of objects is returned. If None,
                                then a random galaxy (or more: see n_random kwarg) is chosen,
                                correcting for catalog-level selection effects if weights are
                                available. [default: None]
        @param n_random         The number of random galaxies to build, if 'index' is None.
                                [default: 1 (set below)]
        @param rng              A random number generator to use for selecting a random galaxy
                                (may be any kind of BaseDeviate or None) and to use in generating
                                any noise field when padding.  [default: None]
        @param gsparams         An optional GSParams argument.  See the docstring for GSParams for
                                details. [default: None]
        """

        if band not in self._valid_band_types:
            raise ValueError('Band {} is not a valid band type for an ngmix catalog!'.format(band))

        # Make rng if needed
        if index is None:
            if rng is None:
                rng = galsim.BaseDeviate()
            elif not isinstance(rng, galsim.BaseDeviate):
                raise TypeError("The rng provided to makeGalaxies is not a BaseDeviate")

        # Select random indices if necessary (no index given).
        if index is None:
            if n_random is None: n_random = 1
            index = self.selectRandomIndex(n_random, rng=rng)
        else:
            # n_random is set to None by default instead of 1 for this check
            if n_random is not None:
                warnings.warn("Ignoring input n_random, since indices were specified!")

        # Format indices correctly
        if hasattr(index, '__iter__'):
            indices = index
        else:
            indices = [index]

        # Now convert ngmix galaxy to GSObject, with details dependent on type
        # QSTN: Should we add a more general scheme for ngmix types?

        galaxies = []

        for index in indices:
            gal = self.ngmix2gs(index, band, gsparams)
            galaxies.append(gal)

        # Store the orig_index as gal.index
        for gal, idx in zip(galaxies, indices):
            gal.index = self.orig_index[idx]

        # Only return a list if there are multiple GSObjects
        if hasattr(index, '__iter__'):
            return galaxies
        else:
            return galaxies[0]

    #----------------------------------------------------------------------------------------------

    def ngmix2gs(self, index, band, gsparams=None):
        """
        This function handles the conversion of a ngmix galaxy to a GS object. The required
        conversion is different for each ngmix catalog type.
        @ param index       The ngmix catalog index of the galaxy to be converted.
        """

        # Column name prefix is not always the same as catalog type
        cp = self.col_prefix
        ct = self.cat_type

        # Grab galaxy shape and size (identical in all bands)
        T = self.catalog[cp+'_T'][index]
        g1, g2 = self.catalog[cp+'_g'][index]
        # We don't want to put these galaxies at their original locations, but could
        # grab them like this:
        # c1, c2 = self.catalog[cp+'_c'][index]

        # Convert from dict to actual GsParams object
        if gsparams is not None:
            gsp = galsim.GSParams(**gsparams)
        else:
            gsp = None

        # Grab current band flux
        if self.de_redden is True:
            flux_colname = cp + '_flux_deredden'
        else:
            flux_colname = cp + '_flux'

        flux = self.catalog[flux_colname][index][self._band_index[band]]

        # NOTE: It used to be the case that all Gaussian-Mixture parameters
        # were in the format of:
        # gm_pars = [centroid1, centroid2, g1, g2, T, flux]
        # (this is identical to ngmix `pars`, except that flux is a vector
        # of fluxes in all bands)
        # However, this now depends on the gmix type, so we have to wait
        # gm_pars = [0.0, 0.0, g1, g2, T, flux]

        # TODO: Implement this once we get a response back from Erin why
        # CM isn't included in this function
        # https://github.com/esheldon/ngmix/blob/master/ngmix/gmix.py#L39
        # This allows us to construct the given gmix type without knowing
        # gm.make_gmix_model(pars, model, **kw):

        # Build the appropriate Gaussian mixture for a given model
        
        # MEGAN changed these calls according to API update: 
        # https://github.com/esheldon/ngmix/blob/89b687b4adacd44b76ffcf13f442da655fc011fb/ngmix/gmix/gmix.py#L46
        if ct == 'gauss':
            # Uses 'simple' pars scheme
            gm_pars = [0.0, 0.0, g1, g2, T, flux]
            #gm = ngmix.gmix.GMixModel(gm_pars, 'gaussian')
            gm = ngmix.GMixModel(gm_pars, "gauss")

        elif ct == 'cm':
            fracdev = self.catalog[cp+'_fracdev'][index]
            TdByTe  = self.catalog[cp+'_TdByTe'][index]
            # Uses 'simple' pars scheme
            gm_pars = [0.0, 0.0, g1, g2, T, flux]
            #gm = ngmix.gmix.GMixCM(fracdev, TdByTe, gm_pars)
            gm = ngmix.GMixModel(gm_pars, "cm")

        elif ct == 'bdf':
            fracdev = self.catalog[cp+'_fracdev'][index]
            TdByTe = self._TdByTe
            # Uses different 'bdf' pars scheme
            gm_pars = [0.0, 0.0, g1, g2, T, fracdev, flux]
            #gm = ngmix.gmix.GMixBDF(pars=gm_pars, TdByTe=TdByTe)
            gm = ngmix.GMixModel(gm_pars, "bdf")

        # The majority of the conversion will be handled by `ngmix.gmix.py`
        gs_gal = gm.make_galsim_object(gsparams=gsp)

        return gs_gal

    #----------------------------------------------------------------------------------------------
    @staticmethod
    def _makeSingleGalaxy(ngmix_catalog, index, rng=None, gsparams=None):
        """ A static function that mimics the functionality of makeGalaxes() for single index.
        The only point of this class is to circumvent some serialization issues. This means it
        can be used through a proxy ngmixCatalog object, which is needed for the config layer.
        """

        # TODO: Write the static version of makeGalaxies!
        # NB: I have noticed an occasional memory issue with N>~200 galaxies. This may be related
        # to the serialization issues Mike talks about in scene.py
        pass

    #----------------------------------------------------------------------------------------------

    def selectRandomIndex(self, n_random=1, rng=None, _n_rng_calls=False):
        """
        Routine to select random indices out of the catalog.  This routine does a weighted random
        selection with replacement (i.e., there is no guarantee of uniqueness of the selected
        indices).  Weighting uses the weight factors available in the catalog, if any; these weights
        are typically meant to remove any selection effects in the catalog creation process.
        @param n_random     Number of random indices to return. [default: 1]
        @param rng          A random number generator to use for selecting a random galaxy
                            (may be any kind of BaseDeviate or None). [default: None]
        @returns A single index if n_random==1 or a NumPy array containing the randomly-selected
        indices if n_random>1.
        """

        # Set up the random number generator.
        if rng is None:
            rng = galsim.BaseDeviate()

        # QSTN: What is the weighting scheme for ngmix catalogs? Will need to adjust below code to
        # match (or exclude entierly)
        if hasattr(self.catalog, 'weight'):
            use_weights = self.catalog.weight[self.orig_index]
        else:
            print('Selecting random object without correcting for catalog-level selection effects.')
            use_weights = None

        # By default, get the number of RNG calls. Then decide whether or not to return them
        # based on _n_rng_calls.
        index, n_rng_calls = galsim.utilities.rand_with_replacement(
                n_random, self.nobjects, rng, use_weights, _n_rng_calls=True)

        if n_random>1:
            if _n_rng_calls:
                return index, n_rng_calls
            else:
                return index
        else:
            if _n_rng_calls:
                return index[0], n_rng_calls
            else:
                return index[0]

    #----------------------------------------------------------------------------------------------

    # The catalog type is referred to as a 'subtype' w.r.t BalInput types
    # (i.e. the BalInput type is 'ngmix_catalog' with subtype self.catalog_type)
    def getSubtype(self):
        return self.cat_type

    def getNObjects(self):
        # Used by input/logger methods
        return self.nobjects

    def getNTot(self):
        # Used by input/logger methods
        return self.ntotal

    def getCatalog(self):
        return self.catalog

    def getBands(self):
        return self.bands

    # TODO: Write remaining `get` methods once saved columns are determined

    #----------------------------------------------------------------------------------------------

    # Since makeGalaxies is a function, not a class, it needs to use an unconventional location
    # for defining certain config parameters.
    makeGalaxies._req_params = {'band': str}
    makeGalaxies._opt_params = {'index' : int,
                                'n_random': int}
    makeGalaxies._single_params = []
    makeGalaxies._takes_rng = True

class ngmixCatalogLoader(galsim.config.InputLoader):
    """ The ngmixCatalog loader doesn't need anything special other than registration as a valid
    input type. These additions are only used for logging purposes.
    """

    def setupImage(self, ngmix_catalog, config, base, logger):
        # This method is blank for a general InputLoader, and a convenient place to put the logger
        if logger: # pragma: no cover
            # Only report as a warning the first time.  After that, use info.
            first = not base.get('_ngmixCatalogLoader_reported_as_warning',False)
            base['_ngmixCatalogLoader_reported_as_warning'] = True
            if first:
                log_level = logging.WARNING
            else:
                log_level = logging.INFO
            if 'input' in base:
                if 'ngmix_catalog' in base['input']:
                    cc = base['input']['ngmix_catalog']
                    if isinstance(cc,list): cc = cc[0]
                    out_str = ''
                    if 'dir' in cc:
                        out_str += '\n  dir = %s'%cc['dir']
                    if 'file_name' in cc:
                        out_str += '\n  file_name = %s'%cc['file_name']
                    # TODO: Add any desired additional ngmix_catalog inputs
                    if out_str != '':
                        logger.log(log_level, 'Using user-specified ngmixCatalog: %s',out_str)
            logger.info("file %d: Ngmix catalog has %d total objects; %d passed initial cuts.",
                        base['file_num'], ngmix_catalog.getNTot(), ngmix_catalog.getNObjects())

# Need to add the ngmixCatalog class as a valid input_type.
galsim.config.RegisterInputType('ngmix_catalog', ngmixCatalogLoader(ngmixCatalog, has_nobj=True))

#####----------------------------------------------------------------------------------------------

def BuildNgmixGalaxy(config, base, ignore, gsparams, logger):
    """ Build a NgmixGalaxy type GSObject from user input."""

    ngmix_cat = galsim.config.GetInputObj('ngmix_catalog', config, base, 'NgmixGalaxy')

    # If galaxies are selected based on index, and index is Sequence or Random, and max
    # isn't set, set it to nobjects-1.
    if 'index' in config:
        galsim.config.SetDefaultIndex(config, ngmix_cat.getNObjects())

    # Grab necessary parameters
    req = ngmixCatalog.makeGalaxies._req_params
    opt = ngmixCatalog.makeGalaxies._opt_params
    single = ngmixCatalog.makeGalaxies._single_params
    ignore = ignore + ['num']

    kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt, single=single,
                                              ignore=ignore)

    # Guaranteed to be present as it is in _req_params
    band = kwargs['band']

    # Convert gsparams from a dict to an actual GSParams object
    if gsparams:
        kwargs['gsparams'] = galsim.GSParams(**gsparams)
    else:
        gsparams = None

    # This handles the case of no index passed in config file
    # Details are in ngmixCatalog
    rng = None
    if 'index' not in kwargs:
        rng = galsim.config.GetRNG(config, base, logger, 'NgmixGalaxy')
        kwargs['index'], n_rng_calls = ngmix_cat.selectRandomIndex(1, rng=rng, _n_rng_calls=True)

        # Make sure this process gives consistent results regardless of the number of processes
        # being used.
        if not isinstance(ngmix_cat, ngmixCatalog) and rng is not None:
            # Then ngmix_cat is really a proxy, which means the rng was pickled, so we need to
            # discard the same number of random calls from the one in the config dict.
            rng.discard(int(n_rng_calls))

    # Check that inputted/set index is valid
    index = kwargs['index']
    if index >= ngmix_cat.getNObjects():
        raise IndexError("%s index has gone past the number of entries in the catalog"%index)

    logger.debug('obj %d: NgmixGalaxy kwargs = %s',base.get('obj_num',0),kwargs)

    kwargs['ngmix_catalog'] = ngmix_cat

    # NOTE: This uses a static method of ngmixCatalog to save memory. Not needed for the moment, but
    # worth looking into when trying to save resources for large Balrog runs
    # ngmix_gal = ngmix_cat._makeSingleGalaxy(**kwargs)

    # Create GSObject galaxies from the ngmix catalog
    ngmix_gal = ngmix_cat.makeGalaxies(band, index=index, gsparams=gsparams)

    # The second item is "safe", a boolean that declares whether the returned value is
    # safe to save and use again for later objects (which is not the case for ngmixGalaxies).
    return ngmix_gal, False

# Register this builder with the config framework:
galsim.config.RegisterObjectType('ngmixGalaxy', BuildNgmixGalaxy, input_type='ngmix_catalog')

