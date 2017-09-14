## Classes for new input type `ngmixCatalog`; intended for use with Balrog.
import galsim
import galsim.config
import numpy as np
import logging
from past.builtins import basestring # Python 2&3 compatibility
import pdb

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
    """

    _req_params = { 'file_name' : str }
    _opt_params = { 'dir' : str }
    _single_params = []
    _takes_rng = False

    def __init__(self, file_name=None, dir=None):

        if dir:
            if not isinstance(file_name, basestring):
                raise ValueError("Cannot provide dir and an HDU instance!")
            import os
            file_name = os.path.join(dir,file_name)
        self.file_name = file_name
        self.read()

    #------------------------------------------------------------------------------------------------

    def read(self):
        """ Read in relevant catalog information"""
        from galsim._pyfits import pyfits
        if isinstance(self.file_name, basestring):
            # If a filename is passed:
            hdu_list = pyfits.open(self.file_name)
            # QSTN: Is hdu 0 or 1 used for ngmix catalog? I'm guessing 1
            hdu = hdu_list[1]
        else:
            # If a fits HDU is directly passed:
            hdu_list = None
            hdu = self.file_name

        # TODO: Read in relevant ngmix catalog header info

        # TODO: Read in relevant ngmix catalog column info (handle a few obvious cases below)

        self.catalog = hdu.data

        # NB: As discussed in `scene.py`, there is a bug in the pyfits FITS_Rec class that leads to memory leaks.
        # The simplest workaround seems to be to convert it to a regular numpy recarray.
        self.catalog = np.array(self.catalog, copy=True)

        # The input logger needs to know the original catalog size
        self.ntotal = len(self.catalog)

        # Others ...

        # Close file!
        if hdu_list: hdu_list.close()

        # Galaxy indices in original ngmix catalog
        self.orig_index = np.arange(self.nobjects)
        mask = np.ones(len(self.orig_index), dtype=bool)

        # TODO: Do any desired masking procedure

        # Do mask cut
        self.orig_index = self.orig_index[mask]
        self.nobjects = len(self.orig_index)

        # TODO: Check for valid value types of inputs

        # TODO: Save relevant instance values

        return

    #------------------------------------------------------------------------------------------------

    def makeGalaxies(self, index=None, n_random=1, rng=None):
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
                                [default: 1]
        @param rng              A random number generator to use for selecting a random galaxy
                                (may be any kind of BaseDeviate or None) and to use in generating
                                any noise field when padding.  [default: None]
        """

        # Make rng if needed
        if index is None or gal_type == 'real':
            if rng is None:
                rng = galsim.BaseDeviate()
            elif not isinstance(rng, galsim.BaseDeviate):
                raise TypeError("The rng provided to makeGalaxy is not a BaseDeviate")

        # Select random indices if necessary (no index given).
        if index is None:
            if n_random is None: n_random = 1
            index = self.selectRandomIndex(n_random, rng=rng)
        else:
            if n_random is not None:
                import warnings
                warnings.warn("Ignoring input n_random, since indices were specified!")

        # Format indices correctly
        if hasattr(index, '__iter__'):
            indices = index
        else:
            indices = [index]

        # TODO: Implement Erin's ngmix->GSObject code

        return

    #------------------------------------------------------------------------------------------------

    @staticmethod
    def _makeSingleGalaxy(ngmix_catalog, index, rng=None, gsparams=None):
        """ A static function that mimics the functionality of makeGalaxy() for single index.
        The only point of this class is to circumvent some serialization issues. This means it can be used
        through a proxy ngmixCatalog object, which is needed for the config layer.
        """

        # TODO: Write the static version of makeGalaxy! (We don't need it for prototyping Balrog, however)
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

    def GetNTot(self):
        # Used by input/logger methods
        return self.ntotal

    # TODO: Write remaining `get` methods once saved columns are determined

    #------------------------------------------------------------------------------------------------

    # Since makeGalaxy is a function, not a class, it needs to use an unconventional location for defining
    # certain config parameters.
    makeGalaxy._req_params = {}
    makeGalaxy._opt_params = { "index" : int,
                               "n_random": int
                               # QSTN: Any other optional params?
                             }
    makeGalaxy._single_params = []
    makeGalaxy._takes_rng = True

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
    req = galsim.ngmixCatalog.makeGalaxy._req_params
    opt = galsim.ngmixCatalog.makeGalaxy._opt_params
    single = galsim.ngmixCatalog.makeGalaxy._single_params
    ignore = ignore + ['num']

    kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt,
                                              single=single ignore=ignore)

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
        if not isinstance(ngmix_cat, galsim.ngmixCatalog) and rng is not None:
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
    ngmix_gal = ngmix_cat.makeGalaxes(index)

    # The second item is "safe", a boolean that declares whether the returned value is
    # safe to save and use again for later objects (which is not the case for ngmixGalaxies).
    return ngmix_gal, False

# Register this builder with the config framework:
galsim.config.RegisterObjectType('ngmixGalaxy', BuildNgmixGalaxy, input_type='ngmix_catalog')

