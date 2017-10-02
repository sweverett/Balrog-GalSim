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
from past.builtins import basestring # Python 2&3 compatibility
import pdb

# We need ngmix.gmix.GMixCM, but this isn't included in the default `__init__.py`
# for ngmix. TODO: Ask if this can be added!
# import sys
# sys.path.append('/home/spencer/Documents/Software/ngmix/ngmix')
# from gmix import GMixCM

class ngmixCatalog(object):
    """ Class that handles galaxy catalogs from ngmix. These are usually stored as *?? files.

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
    """

    _req_params = { 'file_name' : str }
    _opt_params = { 'dir' : str, 'catalog_type' : str, 'bands': str}
    _single_params = []
    _takes_rng = False

    # Only these ngmix catalog types currently supported
    # `gauss`: Single Gaussian fit
    # `cm`: Combined bulge+disk model
    # `mof`: CM model with multi-object fitting
    # TODO: In principle, should be able to handle any type supported by ngmix
    _valid_catalog_types = ['gauss','cm','mof']

    # Only these color bands are currently supported for an ngmix catalog
    _valid_band_types = ['g','r','i','z']

    # Dictionary of color band flux to array index in ngmix catalogs
    _band_index = {'g' : 0, 'r' : 1, 'i' : 2, 'z' : 3}

    # The catalog column name prefix doens't always match the catalog type (e.g. 'mof' has a prefix
    # of 'cm' for most columns). Set this for each new supported catalog type.
    _cat_col_prefix = {'gauss' : 'gauss', 'cm' : 'cm', 'mof' : 'cm'}

    def __init__(self, file_name, dir=None, catalog_type=None, bands=None):

        if dir:
            if not isinstance(file_name, basestring):
                raise ValueError("Cannot provide dir and an HDU instance!")
            import os
            file_name = os.path.join(dir,file_name)
        self.file_name = file_name

        if catalog_type:
            if catalog_type in self._valid_catalog_types:
                if catalog_type not in filename:
                    import warnings
                    warnings.warning("Inputted ngmix catalog type of {} does not match filename, which is standard for DES ngmix catalogs. Ensure this is correct.".format(catalog_type))
                self.cat_type = catalog_type
            else:
                raise ValueError("{} is not a currently supported ngmix catalog type!".format(catalog_type))
        else:
            # Attempt to determine catalog type from filename (this is generally true for DES ngmix
            # catalogs)
            match = 0
            for type in self._valid_catalog_types:
                if type in self.file_name:
                    match +=1
                    self.cat_type = type
            # Reject if multiple matches
            if match == 0:
                raise ValueError("No inputted ngmix catalog type, and no matches in filename! Please set a valid catalog type.")
            if match > 1:
                raise ValueError("No inputted ngmix catalog type, and multiple matches in filename! Please set a valid catalog type.")

        # Catalog column name prefixes don't always match catalog type (e.g. 'cm' is still used for many 'mof' columns)
        self.col_prefix = self._cat_col_prefix[self.cat_type]

        if bands:
            if isinstance(bands, basestring):
                # Strip all whitespace
                bands = bands.replace(" ", "")
                # More useful as a list of individual bands
                bands_list = list(bands)
                if set(bands_list).issubset(self._valid_band_types):
                    self.bands = bands_list
                else:
                    raise ValueError("The only valid color bands for a ngmix catalog are \'griz'\'!")
            else:
                # TODO: Wouldn't be a bad idea to allow a list of individual bands as well
                raise ValueError("Must enter desired color bands as a string! (For example, `bands : \'gr\'`)")
        else:
            # Set the 'g' band to be the default behaviour, although warn user
            import warnings
            warnings.warn('No color band chosen; selecting \'g\'-band by default.')
            self.bands = ['g']

        self.read()

        return

    #------------------------------------------------------------------------------------------------

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

        # NB: As discussed in `scene.py`, there is a bug in the pyfits FITS_Rec class that leads to memory leaks.
        # The simplest workaround seems to be to convert it to a regular numpy recarray.
        self.catalog = np.array(self.catalog, copy=True)

        # The input logger needs to know the original catalog size
        self.ntotal = len(self.catalog)

        # TODO: Any others needed??

        # Close file!
        if hdu_list: hdu_list.close()

        # Galaxy indices in original ngmix catalog
        self.orig_index = np.arange(self.ntotal)

        # Create mask
        self.flags = self.catalog[self.cat_type+'_flags']
        mask = np.ones(len(self.orig_index), dtype=bool)
        mask = self.makeMask(mask)

        # Do mask cut
        self.catalog = self.catalog[mask]
        self.orig_index = self.orig_index[mask]
        self.nobjects = len(self.orig_index)

        return

    #------------------------------------------------------------------------------------------------

    def makeMask(self,mask):
        """Add a masking procedure, if desired."""
        # TODO: Implement any desired masking procedure
        return mask

    #------------------------------------------------------------------------------------------------

    def makeGalaxies(self, index=None, n_random=None, rng=None):
        """
        Construct GSObjects from a list of galaxies in the ngmix catalog specified by `index` (or a randomly generated one).
        This is done using Erin's code.

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
        """

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
                import warnings
                warnings.warn("Ignoring input n_random, since indices were specified!")

        # Format indices correctly
        if hasattr(index, '__iter__'):
            indices = index
        else:
            indices = [index]

        # print('\nIndices = {}'.format(indices))

        # Now convert ngmix galaxy to GSObject, with details dependent on type
        # TODO: Should we add a more general scheme for ngmix types?

        galaxies = []

        for index in indices:
            gal = self.ngmix2gs(index)
            galaxies.append(gal)

        # Store the orig_index as gal.index
        for gal, idx in zip(galaxies, indices):
            gal.index = self.orig_index[idx]

        # Only return a list if there are multiple GSObjects
        if hasattr(index, '__iter__'):
            return galaxies
        else:
            return galaxies[0]

    #------------------------------------------------------------------------------------------------

    def ngmix2gs(self,index):
        """
        This function handles the conversion of a ngmix galaxy to a GS object. The required conversion is
        different for each ngmix catalog type.
        @ param index       The ngmix catalog index of the galaxy to be converted.
        """

        # TODO: Make sure that the index usage is consistent with original/masked indices!!
        # TODO: In principle, I believe that we should be returning ChromaticObjects, rather
        #       than GSObjects! Look at details of implementation in future.

        # Column name prefix is not always the same as catalog type
        cp = self.col_prefix
        ct = self.cat_type

        # Grab galaxy shape and size (identical in all bands)
        T = self.catalog[cp+'_T'][index]
        g1, g2 = self.catalog[cp+'_g'][index]
        # We don't want to put these galaxies at their original locations!
        # c1, c2 = self.catalog[cp+'_c'][index]

        # List of individual band GSObjects
        gsobjects = []

        # Iterate over all desired bands
        for band in self.bands:
            # Grab current band flux
            flux = self.catalog[cp+'_flux'][index][self._band_index[band]]

            # Gaussian-Mixture parameters are in the format of:
            # gm_pars = [centroid1, centroid2, g1, g2, T, flux]
            # (this is identical to ngmix catalogs, except that flux is a vector
            # of fluxes in all color bands)
            #
            # TODO: Naively I would expect to put shape (g1,g2) in this array,
            # but I will follow what is done in ESheldon's repo `nsim/sigms.py`
            # https://github.com/esheldon/nsim/
            gm_pars = [0.0, 0.0, 0.0, 0.0, T, flux]

            if ct == 'gauss':
                import warning
                # This should be fixed very soon!
                warnings.warn("\nWARNING: GAUSS catalog conversion hasn't been implemented yet!\n")
                pass
            elif ct == 'cm':
                # The majority of the conversion will be handled by `ngmix/gmix.py`
                # Build the appropriate Gaussian mixture for a cm-model
                fracdev = self.catalog[cp+'_fracdev'][index]
                TdByTe = self.catalog[cp+'_TdByTe'][index]
                # print('\nfracdev, TdByTe = {}, {}\n'.format(fracdev, TdByTe))
                gm = ngmix.gmix.GMixCM(fracdev, TdByTe, gm_pars)

                # Convert to a GSObject
                gal = gm.make_galsim_object()
                gal = gal.shear(g1=g1, g2=g2)

                # TODO: Implement shearing of intrinsic shape
                # now shear it
                # if 'shear' in pars:
                #     shear=pars['shear']
                #     gal = gal.shear(g1=shear.g1, g2=shear.g2)

                gsobjects.append(gal)

            elif ct == 'mof':
                import warning
                # This should be fixed very soon!
                warnings.warn("\nWARNING: MOF catalog conversion hasn't been implemented yet!\n")
                pass
            else:
                # TODO: While this check has already been made, it is possible that a new valid catalog type
                # is added but whose GSObject conversion isn't implemented. Should come up with a more robust check!
                raise ValueError('The ngmix catalog type {} is not yet supported!'.format(ct))

        # TODO: We may want to use a ChromaticObject in the end, but for now we will just add the GSObjects together
        gs_gal = galsim.Add(gsobjects)

        return gs_gal

    #------------------------------------------------------------------------------------------------
    @staticmethod
    def _makeSingleGalaxy(ngmix_catalog, index, rng=None, gsparams=None):
        """ A static function that mimics the functionality of makeGalaxes() for single index.
        The only point of this class is to circumvent some serialization issues. This means it can be used
        through a proxy ngmixCatalog object, which is needed for the config layer.
        """

        # TODO: Write the static version of makeGalaxies! (We don't need it for prototyping Balrog, however)
        pass

    #------------------------------------------------------------------------------------------------

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

        # QSTN: What is the weighting scheme for ngmix catalogs? Will need to adjust below code to match (or exclude entierly)
        if hasattr(self.catalog, 'weight'):
            use_weights = self.catalog.weight[self.orig_index]
        else:
            import warnings
            warnings.warn('Selecting random object without correcting for catalog-level '
                          'selection effects.')
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

    def getNObjects(self):
        # Used by input/logger methods
        return self.nobjects

    def getNTot(self):
        # Used by input/logger methods
        return self.ntotal

    # TODO: Write remaining `get` methods once saved columns are determined

    #------------------------------------------------------------------------------------------------

    # Since makeGalaxies is a function, not a class, it needs to use an unconventional location for defining
    # certain config parameters.
    makeGalaxies._req_params = {}
    makeGalaxies._opt_params = { "index" : int,
                               "n_random": int
                               # TODO: Any other optional params?
                             }
    makeGalaxies._single_params = []
    makeGalaxies._takes_rng = True

#####------------------------------------------------------------------------------------------------

class ngmixCatalogLoader(galsim.config.InputLoader):
    """ The ngmixCatalog loader doesn't need anything special other than registration as a valid input type.
        These additions are only used for logging purposes.
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
galsim.config.RegisterInputType('ngmix_catalog', ngmixCatalogLoader(ngmixCatalog))

#####------------------------------------------------------------------------------------------------

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

    kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt, single=single, ignore=ignore)

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

    # NB: This uses a static method of ngmixCatalog to save memory. Not needed for the moment, but
    # worth looking into when trying to save resources for large Balrog runs
    # ngmix_gal = ngmix_cat._makeSingleGalaxy(**kwargs)

    # Create GSObject galaxies from the ngmix catalog
    ngmix_gal = ngmix_cat.makeGalaxies(index)

    # The second item is "safe", a boolean that declares whether the returned value is
    # safe to save and use again for later objects (which is not the case for ngmixGalaxies).
    return ngmix_gal, False

# Register this builder with the config framework:
galsim.config.RegisterObjectType('ngmixGalaxy', BuildNgmixGalaxy, input_type='ngmix_catalog')

