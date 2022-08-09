#----------------------------------------------------------------------------
# Classes for new input type `desStarCatalog`; intended for use with Balrog,
# but should be extendable.
#
# Contributors:
# Spencer Everett (UCSC)
#----------------------------------------------------------------------------

import galsim
import galsim.config
import numpy as np
import os
import warnings
import logging
from past.builtins import basestring # Python 2&3 compatibility
from astropy.io import ascii

# Can use for debugging
# import pudb

class desStarCatalog(object):
    """ Class that handles Sahar's star catalogs for DES. These are cvs files with names typically
    of the form #TODO.

    # @param base_dir        Location of all tile-specific catalog directories.
    # @param model_type      A string of which star catalog model type to use
                             (e.g. 'Extra_10_percent_16.5-26.5').
    # @param tile            Which catalog tile to load.
    # @param bands           A string of the desired bands to simulate from (only griz allowed). For
    #                        example, selecting only the 'g' and 'r' bands would be done by setting
    #                        bands='gr'.
    # @param snr_min         The lower allowed bound for signal-to-noise ratio (snr). Can be any
    #                        positive value, as long as it is smaller than `snr_max`. All objects with
    #                        negative snr are removed by default.
    # @param file_type       Which file extension is used for the star catalogs (not needed yet, is csv)
    #                        between versions is different (not needed yet)
    # @param data_version    Specify which version of the catalog is being used if the processing
    #                        between versions is different (not needed yet)
    # @param zeropoint       Reference zeropoint used for star catalog (Default: 30)
    #                        objects below this cutoff will be removed. (Default: `t_frac=0.5`).
    # @param base_model      Used for Balrog simulations where multiple model types will be
    #                        needed accross many realizations of images. Must be a valid model type
    #                        and must be a divisor of the selected model type.
    #                        For example, a simulation with 10 realizations will need to use a model
    #                        type of `Extra_%_percent` where % is a multiple of 10 up to 100.
    #                        So the base_model would be 'Extra_10_percent'.
    # @param _nobjects_only  This is only passed if GalSim wants to know how many input objects will
                             be used without processing the whole input catalog.
    """

    _req_params = { 'base_dir' : str, 'model_type' : str, 'tile' : str, 'bands': str}
    # TODO: any others?
    _opt_params = { 'file_type' : str, 'data_version' : str, 'zeropoint' : float,
                    'base_model' : str}
    _single_params = []
    _takes_rng = False

    # TODO: Shouldn't be hard to hadd fits, but currently not needed
    _valid_file_types = ['cvs']

    # Data versions that have currently have implementations
    _valid_data_versions = ['des-pizza-slices-y6-v13'] # MEGAN CHANGED FROM 'y3v02'

    # From Sahar's catalog directory, currently in `/data/des20.b/data/sallam/Yanny-balrog/5sigma/`
    # NOTE: Gets set below after determining data version
    _valid_model_types = []

    # Only these color bands are currently supported for star injection
    _valid_band_types = 'griz'

    # Dictionary of color band flux to array index in star catalog
    _band_index = {'g' : 0, 'r' : 1, 'i' : 2, 'z' : 3}

    def __init__(self, base_dir, model_type, tile, bands, file_type=None, data_version=None,
                 zeropoint=None, base_model=None):

        if not os.path.isdir(base_dir):
            raise ValueError('{} is not a valid directory or does not exist!'.format(base_dir))
        self.base_dir = base_dir
        
        if not data_version:
            warnings.warn('No data version passed - assuming `des-pizza-slices-y6-v13`.') # MEGAN CHANGED FROM 'y3v02'
            data_version = 'des-pizza-slices-y6-v13' # MEGAN CHANGED FROM 'y3v02'
        if data_version not in self._valid_data_versions:
            raise ValueError('`{}` does not have an implementation built yet!'.format(
                self.data_version))
        self.data_version = data_version

        self._set_valid_model_types()

        # TODO: Would be nice to have a check for valid DES tile names!
        # tile check...
        self.tile = tile

        if model_type not in self._valid_model_types:
            raise ValueError('{} is not a valid model type! '.format(model_type) +
                             'Currently allowed types are {}'.format(self._valid_model_types))
        self.model_type = model_type

        # NOTE: This parameter is only used for Balrog simulations where multiple model types will be
        # needed accross many 'realizations' of images. Must also be a valid model type, and must be a
        # divisor of the selected model type.
        #
        # For example, a simulation with 10 realizations will need to use a model type of
        # `Extra_%_percent` where % is a multiple of 10 up to 100. So the base_model would be
        # `Extra_10_percent`.
        if base_model:
            if base_model not in self._valid_model_types:
                raise ValueError('Base model {} is not a valid model type! '.format(model_type) +
                                'Currently allowed types are {}'.format(self._valid_model_types))
            prefix = base_model.split('_')[0]
            if prefix == 'Model':
                self.base_percent = 100
            elif prefix == 'Extra':
                self.base_percent = base_model.split('_')[1]
            self.base_model = base_model
        else:
            self.base_model, self.base_percent = None, None

        if isinstance(bands, basestring):
            # Strip all whitespace
            bands = bands.replace(' ', '')
            # More useful as a list of individual bands
            bands_list = list(bands)
            if set(bands_list).issubset(self._valid_band_types):
                self.bands = bands_list
            else:
                raise ValueError("The only valid color bands for a des star catalog are \'griz'\'!")
        else:
            # TODO: Wouldn't be a bad idea to allow a list of individual bands as well
            raise ValueError('Must enter desired color bands as a string!' +
                                ' (For example, `bands : \'gr\'`')

        if file_type:
            if file_type not in self._valid_file_types:
                raise ValueError('{} is not a valid file type! '.format(file_type) +
                                 'Currently allowed types are {}'.format(self._valid_file_types))
        else:
            # Default is csv
            self.file_type = 'csv'

        if not zeropoint:
            if self.data_version == 'des-pizza-slices-y6-v13': #MEGAN CHANGED
                # Catalogs made w/ zp of 30
                self.zeropoint = 30.0
            else:
                # In future, may be different defaults
                warnings.warn('No zeropoint passed; using default value of 30.0')
                self.zeropoint = 30.0

        # self._setup_files()
        self._read_catalog()

        return

    #------------------------------------------------------------------------------------------------

    def _set_valid_model_types(self):
        '''
        Given data version, construct the allowed model types from Sahar's star catalogs.
        '''

        # Uses public function to allow external use
        #print("self.data_version", self.data_version) # MEGAN ADDED
        self._valid_model_types = return_valid_model_types(data_version=self.data_version)

        return

    def _read_catalog(self):
        '''
        Setup file directory structure given base directory and model type.
        Load in the star catalog for the given model type and tile name.
        '''

        if self.data_version == 'des-pizza-slices-y6-v13': #MEAGN CHANGED
            self.model_dir = os.path.join(self.base_dir, self.model_type)
            filename = 'Model_{}.{}'.format(self.tile, self.file_type)
            self.cat_file = os.path.join(self.model_dir, filename)

            if self.file_type == 'csv':
                # TODO: I believe GalSim has a built in CSV reader; worth looking at
                catalog = ascii.read(self.cat_file, format='csv')
                # NB: Copied due to possible memory leak issue; discussed in `scene.py`
                self.catalog = np.array(catalog, copy=True)

                # Total objects before cuts
                self.ntotal = len(self.catalog)

                # Star indices in original star catalog
                self.orig_index = np.arange(self.ntotal)

                # Get flags and create mask
                self._set_flags()
                self._make_mask()

                # Do mask cut
                self._mask_cut()

        # pudb.set_trace()

        return

    #------------------------------------------------------------------------------------------------

    def _set_flags(self):

        if self.data_version == 'des-pizza-slices-y6-v13': #MEGAN CHANGED
            # No flags in current catalogs
            self.flags = None

        return

    #------------------------------------------------------------------------------------------------

    def _make_mask(self):
        """Add a masking procedure, if desired."""
        # TODO: Allow multiple masking procedures
        # TODO works as currently written, but should rewrite to conform to conventional masking
        # definition

        mask = np.ones(len(self.orig_index), dtype=bool)

        if self.data_version == 'des-pizza-slices-y6-v13': #MEGAN CHANGED
            # No mask cuts in this version
            pass

        self.mask = mask

        return

    #------------------------------------------------------------------------------------------------

    def _mask_cut(self):
        """Do mask cut defined in `makeMask()`."""
        self.catalog = self.catalog[self.mask]
        self.orig_index = self.orig_index[self.mask]
        self.nobjects = len(self.orig_index)
        # print('Ntotal (stars): {}\nNobjects (stars): {}'.format(self.ntotal,self.nobjects))

        return

    #------------------------------------------------------------------------------------------------

    def make_stars(self, band, index=None, n_random=None, rng=None, gsparams=None):
        """
        Construct GSObjects from a list of stars in the des star catalog specified by `index`
        (or a randomly generated one).
        NOTE: If `image` type is set to `Balrog`, then there will be different behaviour. Will
        inject *all* stars whose positions are in the current image, *at* that position!

        @param index            Index of the desired star in the catalog for which a GSObject
                                should be constructed.  You can also provide a list or array of
                                indices, in which case a list of objects is returned. If None,
                                then a random star (or more: see n_random kwarg) is chosen,
                                correcting for catalog-level selection effects if weights are
                                available. [default: None]
        @param n_random         The number of random stars to build, if 'index' is None.
                                [default: 1 (set below)]
        @param rng              A random number generator to use for selecting a random star
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
                raise TypeError("The rng provided to make_stars is not a BaseDeviate")

        # Select random indices if necessary (no index given).
        if index is None:
            if n_random is None: n_random = 1
            index = self.selectRandomIndex(n_random, rng=rng)
        else:
            # n_random is set to None by default instead of 1 for this check
            if n_random is not None:
                import warnings
                warnings.warn("Ignoring input n_random, since indices were specified!")

        # Format indices correctly
        if hasattr(index, '__iter__'):
            indices = index
        else:
            indices = [index]

        stars = []

        for index in indices:
            star = self.star2gs(index, band, gsparams)
            stars.append(star)

        # Store the orig_index as star.index
        for star, idx in zip(stars, indices):
            star.index = self.orig_index[idx]

        # Only return a list if there are multiple GSObjects
        if hasattr(index, '__iter__'):
            return stars
        else:
            return stars[0]

    #------------------------------------------------------------------------------------------------

    def star2gs(self, index, band, gsparams=None):

        if (gsparams is not None) and (not isinstance(gsparams, galsim.GSParams)):
                    if isinstance(gsparams, dict):
                        # Convert to actual gsparams object
                        gsparams = galsim.GSParams(**gsparams)
                    else:
                        raise TypeError('Only `dict` and `galsim.GSParam` types allowed for '
                                        'gsparams; input has type of {}.'.format(type(gsparams)))
        # List of individual band GSObjects
        gsobjects = []

        # NOTE: Used to iterate over all desired bands here, but this feature has
        # been removed as it is not useful and led to the possiblity of unintended bugs
        # for band in self.bands:
        if self.data_version == 'des-pizza-slices-y6-v13': #MEGAN CHANGED
            # Needed for all mag calculations
            gmag = self.catalog['g_Corr'][index]

            # Grab current band magnitude
            if band == 'g':
                mag = gmag
            elif band == 'r':
                mag = gmag - self.catalog['gr_Corr'][index]
            elif band == 'i':
                mag = gmag - self.catalog['gr_Corr'][index] - self.catalog['ri_Corr'][index]
            elif band == 'z':
                mag = gmag - self.catalog['gr_Corr'][index] - self.catalog['ri_Corr'][index] \
                            - self.catalog['iz_Corr'][index]
            else:
                raise ValueError('Band {} is not an allowed band input '.format(band) +
                                    'for data_version of {}!'.format(self.data_version))

            # Now convert to flux
            flux = np.power(10.0, 0.4 * (self.zeropoint - mag))

            # (star cats calibrated at zp=30)
            # NOTE: Will just use `scale_flux` parameter in galsim config for now
            # flux_factor = ...

            # Stars are treated as a delta function, to be convolved w/ set PSF
            gs_star = galsim.DeltaFunction(flux)

        else:
            # Should already be checked by this point, but just to be safe:
            raise ValueError('There is no implementation for `star2gs` for data_version ' +
                                'of {}'.format(self.data_version))

        return gs_star

    #------------------------------------------------------------------------------------------------
    @staticmethod
    def _make_single_star(des_star_catalog, index, rng=None, gsparams=None):
        """ A static function that mimics the functionality of make_stars() for single index.
        The only point of this class is to circumvent some serialization issues. This means it can be used
        through a proxy desStarCatalog object, which is needed for the config layer.
        """

        # TODO: Write the static version of make_stars! (We don't need it for prototyping Balrog, however)
        pass

    #------------------------------------------------------------------------------------------------

    def selectRandomIndex(self, n_random=1, rng=None, _n_rng_calls=False):
        """
        Routine to select random indices out of the catalog.  This routine does a weighted random
        selection with replacement (i.e., there is no guarantee of uniqueness of the selected
        indices).  Weighting uses the weight factors available in the catalog, if any; these weights
        are typically meant to remove any selection effects in the catalog creation process.
        @param n_random     Number of random indices to return. [default: 1]
        @param rng          A random number generator to use for selecting a random star.
                            (may be any kind of BaseDeviate or None). [default: None]
        @returns A single index if n_random==1 or a NumPy array containing the randomly-selected
        indices if n_random>1.
        """

        # Set up the random number generator.
        if rng is None:
            rng = galsim.BaseDeviate()

        # QSTN: What is the weighting scheme for des star catalogs? Will need to adjust below code
        # to match (or exclude entierly)
        if hasattr(self.catalog, 'weight'):
            use_weights = self.catalog.weight[self.orig_index]
        else:
            import warnings
            warnings.warn('Selecting random object without correcting for catalog-level selection effects.')
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

    #------------------------------------------------------------------------------------------------

    def indices_in_region(self, rlims, dlims, boundary_cross=False):
        '''
        Returns the indices of all stars contained within the ra/dec limits.
        `boundary_cross` is passed as r1 > r2 in `Y3A2_COADDTILE_GEOM.fits`
        if the tile passes over the 0/360 degree boundary.
        '''

        # pudb.set_trace()

        r1, r2 = rlims[0], rlims[1]
        d1, d2 = dlims[0], dlims[1]

        if not boundary_cross:
            assert r1 < r2
            indices = np.where( (self.catalog['RA_new']>r1) &
                                (self.catalog['RA_new']<r2) &
                                (self.catalog['DEC_new']>d1) &
                                (self.catalog['DEC_new']<d2) )[0]
        else:
            # Account for positions crossing accross RA=0/360
            assert r1 > r2
            indices = np.where(
                ( (self.catalog['RA_new']>=r1) & (self.catalog['RA_new']<360) ) | # r1<ra<360
                ( (self.catalog['RA_new']>= 0) & (self.catalog['RA_new']<=r2) )  & # 0=<ra<r2
                ( (self.catalog['DEC_new']>d1) & (self.catalog['DEC_new']<d2) ) )[0]

        return indices

    #------------------------------------------------------------------------------------------------

    def getNObjects(self):
        # Used by input/logger methods
        return self.nobjects

    def getNTot(self):
        # Used by input/logger methods
        return self.ntotal

    def getCatalog(self, indices=None):
        if indices is None:
            return self.catalog
        else:
            return self.catalog[indices]

    #------------------------------------------------------------------------------------------------

    # Since make_stars is a function, not a class, it needs to use an unconventional location for defining
    # certain config parameters.
    make_stars._req_params = {'band' : str}
    make_stars._opt_params = {'index' : int,
                              'n_random': int}
    make_stars._single_params = []
    make_stars._takes_rng = True

#####------------------------------------------------------------------------------------------------
# Helper Functions

def return_valid_model_types(data_version='des-pizza-slices-y6-v13'):  #MEGAN changed from 'y3v02'
    '''
    Useful to have this separate from _set_valid_model_types() for outside use.
    '''

    if data_version == 'des-pizza-slices-y6-v13': #MEGAN changed from 'y3v02'
        # There are the full-density models...
        valid_model_types = ['Model', 'Model_16.5-26.5', 'Model_16.5-27.5']
        # ...and the 'extra' density models, which are partitioned by percentage
        percents = np.arange(10, 110, 10, dtype=int)
        for per in percents:
            valid_model_types.append('Extra_{}_percent'.format(per))
            valid_model_types.append('Extra_{}_percent_16.5-26.5'.format(per))
            valid_model_types.append('Extra_{}_percent_16.5-27.5'.format(per))

    else:
        raise ValueError('{} does not have an implementation built yet!'.format(data_version))

    return valid_model_types

#####------------------------------------------------------------------------------------------------

class desStarCatalogLoader(galsim.config.InputLoader):
    """
    The desStarCatalog loader doesn't need anything special other than registration as a valid input type.
    These additions are only used for logging purposes.
    """

    def setupImage(self, des_star_catalog, config, base, logger):
        # This method is blank for a general InputLoader, and a convenient place to put the logger
        if logger: # pragma: no cover
            # Only report as a warning the first time.  After that, use info.
            first = not base.get('_desStarCatalogLoader_reported_as_warning',False)
            base['_desStarCatalogLoader_reported_as_warning'] = True
            if first:
                log_level = logging.WARNING
            else:
                log_level = logging.INFO
            if 'input' in base:
                if 'des_star_catalog' in base['input']:
                    cc = base['input']['des_star_catalog']
                    if isinstance(cc,list): cc = cc[0]
                    out_str = ''
                    if 'base_dir' in cc:
                        out_str += '\n  dir = %s'%cc['base_dir']
                    if 'model_type' in cc:
                        out_str += '\n  model_type = %s'%cc['model_type']
                    if out_str != '':
                        logger.log(log_level, 'Using user-specified desStarCatalog: %s',out_str)
            logger.info('file %d: DES star catalog has %d total objects; %d passed initial cuts.',
                        base['file_num'], des_star_catalog.getNTot(),des_star_catalog.getNObjects())

# Need to add the desStarCatalog class as a valid input_type.
galsim.config.RegisterInputType('des_star_catalog', desStarCatalogLoader(desStarCatalog,
                                                                         has_nobj=True))

#####------------------------------------------------------------------------------------------------

def build_desStar(config, base, ignore, gsparams, logger):
    '''
    Build a desStar type GSObject from user input.
    NOTE: If `image` type is set to `Balrog`, then there will be different behaviour. Will
    inject *all* stars whose positions are in the current image, *at* that position!
    '''

    des_star_cat = galsim.config.GetInputObj('des_star_catalog', config, base, 'des_star_catalog')

    # If stars are selected based on index, and index is Sequence or Random, and max
    # isn't set, set it to nobjects-1.
    if 'index' in config:
        galsim.config.SetDefaultIndex(config, des_star_cat.getNObjects())

    # Grab necessary parameters
    req = desStarCatalog.make_stars._req_params
    opt = desStarCatalog.make_stars._opt_params
    single = desStarCatalog.make_stars._single_params
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
    # Details are in desStarCatalog
    rng = None
    if 'index' not in kwargs:
        rng = galsim.config.GetRNG(config, base, logger, 'DES_Star')
        kwargs['index'], n_rng_calls = des_star_cat.selectRandomIndex(1, rng=rng, _n_rng_calls=True)

        # Make sure this process gives consistent results regardless of the number of processes
        # being used.
        if not isinstance(des_star_cat, desStarCatalog) and rng is not None:
            # Then des_star_cat is really a proxy, which means the rng was pickled, so we need to
            # discard the same number of random calls from the one in the config dict.
            rng.discard(int(n_rng_calls))

    # Check that inputted/set index is valid
    index = kwargs['index']
    if index >= des_star_cat.getNObjects():
        raise IndexError("%s index has gone past the number of entries in the catalog"%index)

    logger.debug('obj %d: DES_Star kwargs = %s',base.get('obj_num',0),kwargs)

    kwargs['des_star_catalog'] = des_star_cat

    # NOTE: This uses a static method of desStarCatalog to save memory. Not needed for the moment, but
    # worth looking into when trying to save resources for large Balrog runs
    # star = des_star_catalog._make_single_star(**kwargs)

    # Create GSObject stars from the star catalog
    des_stars = des_star_cat.make_stars(band, index=index, gsparams=gsparams)

    # The second item is "safe", a boolean that declares whether the returned value is
    # safe to save and use again for later objects (which is not the case for des_stars).
    return des_stars, False

# Register this builder with the config framework:
galsim.config.RegisterObjectType('desStar', build_desStar, input_type='des_star_catalog')
