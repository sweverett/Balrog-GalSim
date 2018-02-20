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
# import pudb

# TODO: Include noise, pixscale
# TODO: Understand T=0 cases! (shouldn't be)
# TODO: Handle case of gparams=None

class ngmixCatalog(object):
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
                           objects below this cutoff will be removed. (Default: `t_frac=0.5`).
    @param _nobjects_only  This is only passed if GalSim wants to know how many input objects will
                           be used without processing the whole input catalog.
    """

    _req_params = { 'file_name' : str }
    _opt_params = { 'dir' : str, 'catalog_type' : str, 'bands': str, 'snr_min' : float,
                    'snr_max' : float, 't_frac' : float, 't_min' : float, 't_max' : float}
    _single_params = []
    _takes_rng = False

    # Only these ngmix catalog types currently supported
    # `gauss`: Single Gaussian fit
    # `cm`: Combined bulge+disk model
    # `mof`: CM model with multi-object fitting
    # NOTE: In principle, should be able to handle any type supported by ngmix
    _valid_catalog_types = ['gauss','cm','mof']

    # Only these color bands are currently supported for an ngmix catalog
    _valid_band_types = ['g','r','i','z']

    # Dictionary of color band flux to array index in ngmix catalogs
    _band_index = {'g' : 0, 'r' : 1, 'i' : 2, 'z' : 3}

    # The catalog column name prefix doens't always match the catalog type (e.g. 'mof' has a prefix
    # of 'cm' for most columns). Set this for each new supported catalog type.
    _cat_col_prefix = {'gauss' : 'gauss', 'cm' : 'cm', 'mof' : 'cm'}

    def __init__(self, file_name, dir=None, catalog_type=None, bands=None, snr_min=None, snr_max=None,
                 t_frac=None, t_min=None, t_max=None, _nobjects_only=False):

        if dir:
            if not isinstance(file_name, basestring):
                raise ValueError("Cannot provide dir and an HDU instance!")
            import os
            file_name = os.path.join(dir,file_name)
        if not isinstance(file_name, basestring):
            raise ValueError("The inputted filename must be a string!")
        self.file_name = file_name

        if catalog_type:
            if catalog_type in self._valid_catalog_types:
                if catalog_type not in self.file_name:
                    import warnings
                    warnings.warn( ("Inputted ngmix catalog type of `{}` does not match filename, which is standard "
                                     "for DES ngmix catalogs. Ensure this is correct.".format(catalog_type) ) )
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
                raise ValueError("`t_max` must be greater than `t_min`!")
            self.snr_max = snr_max
        else:
            self.snr_max = None

        if t_min:
            if t_min < 0.0:
                raise ValueError("The object size cutoff `t_min` must be positive!".format())
            self.t_min = t_min
        else:
            self.t_min = 0.0

        if t_max:
            if t_max < 0.0:
                raise ValueError("The object size cutoff `t_max` must be positive!".format())
            elif t_min and (t_min > t_max):
                raise ValueError("`t_max` must be greater than `t_min`!")
            self.t_max = t_max
        else:
            self.t_max = None

        if t_frac:
            if t_frac < 0.0:
                raise ValueError("The allowed size/size_err fraction `t_frac` must be positive!")
            self.t_frac = t_frac
        else:
            # Default of t_frac = 0.5, but warn user that some objects will be removed.
            import warnings
            warnings.warn("No t_frac cutoff was chosen; removing objects with `T/T_err` < 0.5.")
            self.t_frac = 0.5

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

        # Close file!
        if hdu_list: hdu_list.close()

        # Galaxy indices in original ngmix catalog
        self.orig_index = np.arange(self.ntotal)

        # Get flags and create mask
        self.getFlags()
        self.makeMask()

        # Do mask cut
        self.maskCut()

        # pudb.set_trace()

        return

    #------------------------------------------------------------------------------------------------

    def getFlags(self):
        """Retrieve object flags, where implementation depends on catalog type."""
        self.flags = self.catalog[self.col_prefix+'_flags']

        # TODO: Check for additional flags
        if self.cat_type == 'mof':
            self.flags_mof = self.catalog[self.col_prefix+'_mof_flags']

        return

    #------------------------------------------------------------------------------------------------

    def makeMask(self):
        """Add a masking procedure, if desired."""
        # TODO: Allow multiple masking procedures
        # TODO works as currently written, but should rewrite to conform to conventional masking
        # definition

        cp = self.col_prefix

        mask = np.ones(len(self.orig_index), dtype=bool)

        # For now, remove objects with any flags present
        mask[self.flags != 0] = False
        # Extra flags for 'mof' catalogs
        if self.cat_type == 'mof':
            mask[self.flags_mof != 0] = False

        # Remove any object with `T/T_err` < t_frac
        T_fraction = self.catalog[cp+'_T'] / self.catalog[cp+'_T_err']
        mask[T_fraction < self.t_frac ] = False

        # Remove objects with snr_min < S/N < snr_max
        mask[self.catalog[cp+'_s2n_r'] < self.snr_min] = False
        if self.snr_max:
            mask[self.catalog[cp+'_s2n_r'] > self.snr_max] = False

        # Remove objects with size T outside of desired bounds
        mask[self.catalog[cp+'_T'] < self.t_min] = False
        if self.t_max:
            mask[self.catalog[cp+'_T'] > self.t_max] = False

        # TODO: Quick fix for the moment, should figure out rigorous cutoff
        # Looks like we likely want SNR > 10 and T/T_err>0.5; leaving for the moment
        # but should be deleted soon
        # cut = 1e-5
        # mask[self.catalog[self.col_prefix+'_T'] < cut] = False

        self.mask = mask

        return

    #------------------------------------------------------------------------------------------------

    def maskCut(self):
        """Do mask cut defined in `makeMask()`."""
        self.catalog = self.catalog[self.mask]
        self.orig_index = self.orig_index[self.mask]
        self.nobjects = len(self.orig_index)
        print('Ntotal: {}\nNobjects: {}'.format(self.ntotal,self.nobjects))

        return

    #------------------------------------------------------------------------------------------------

    def makeGalaxies(self, index=None, n_random=None, rng=None, gsparams=None):
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
        @param gsparams         An optional GSParams argument.  See the docstring for GSParams for
                                details. [default: None]
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
            gal = self.ngmix2gs(index,gsparams)
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

    def ngmix2gs(self, index, gsparams):
        """
        This function handles the conversion of a ngmix galaxy to a GS object. The required conversion is
        different for each ngmix catalog type.
        @ param index       The ngmix catalog index of the galaxy to be converted.
        """

        # TODO: Make sure that the index usage is consistent with original/ indices!!

        # Column name prefix is not always the same as catalog type
        cp = self.col_prefix
        ct = self.cat_type

        # Grab galaxy shape and size (identical in all bands)
        # T = self.catalog[cp+'_T'][index]
        # TODO: Check this conversion! (Not listed in catalogs, but appears to work.)
        # unit_conv = 60.0 # Convert between arcmin & arcsec
        unit_conv = 1.0 # If no conversion is needed
        T = unit_conv*self.catalog[cp+'_T'][index]
        g1, g2 = self.catalog[cp+'_g'][index]
        # We don't want to put these galaxies at their original locations!
        # c1, c2 = self.catalog[cp+'_c'][index]

        # List of individual band GSObjects
        gsobjects = []

        # Convert from dict to actual GsParams object
        # TODO: Currently fails if gsparams isn't passed!
        if gsparams: gsp = galsim.GSParams(**gsparams)
        else: gsp = None
        # galsim.GSParams(**gsparams)

        # Only here for testing purposes
        # print('\nOBJ ID: {}\n'.format(self.catalog['id'][index]))
        # print('\nOBJ (RA,DEC): ({},{})\n'.format(self.catalog['ra'][index],self.catalog['dec'][index]))

        # Iterate over all desired bands
        for band in self.bands:
            # Grab current band flux
            flux = self.catalog[cp+'_flux'][index][self._band_index[band]]

            if ct == 'gauss':
                # Mimicking implementation in Erin's ngmix/gmix.py function `make_galsim_object()`
                # This avoids using the whole ngmix machinery for the simple Gaussian case, in which
                # the conversion is more straightforward.
                T_round = ngmix.moments.get_Tround(T, g1, g2)
                sigma_round = np.sqrt(T_round/2.0)
                gal = galsim.Gaussian(flux=flux, sigma=sigma_round, gsparams=gsp)
                gal = gal.shear(g1=g1, g2=g2)

            elif ct == 'cm' or ct == 'mof':
                # The majority of the conversion will be handled by `ngmix/gmix.py`

                # Gaussian-Mixture parameters are in the format of:
                # gm_pars = [centroid1, centroid2, g1, g2, T, flux]
                # (this is identical to ngmix catalogs, except that flux is a vector
                # of fluxes in all color bands)
                #
                # NOTE: Naively I would expect to put shape (g1,g2) in this array,
                # but I will follow what is done in ESheldon's repo `nsim/sigms.py`
                # https://github.com/esheldon/nsim/
                gm_pars = [0.0, 0.0, 0.0, 0.0, T, flux]

                # Build the appropriate Gaussian mixture for a cm-model
                # NOTE: we use a slightly modified version of the GMixCM class.
                fracdev = self.catalog[cp+'_fracdev'][index]
                TdByTe = self.catalog[cp+'_TdByTe'][index]
                # print('\nfracdev, TdByTe = {}, {}\n'.format(fracdev, TdByTe))
                gm = BalGMixCM(fracdev, TdByTe, gm_pars)

                # pudb.set_trace()

                # Convert to a GSObject

                gal = gm.make_galsim_object(gsp) # Uses customized ngmix code
                # gal = gm.make_galsim_object() # Uses native ngmix

                # print('gal._gsparams = {}'.format(gal._gsparams))
                # print('_gsparams.__dict__ = {}'.format(gal._gsparams.__dict__))

                if ct == 'mof':
                    # TODO: Add any extra processing for mof catalogs, if needed
                    pass

            else:
                # TODO: While this check has already been made, it is possible that a new valid catalog type
                # is added but whose GSObject conversion isn't implemented. Should come up with a more robust check!
                raise ValueError('The ngmix catalog type {} is not yet supported!'.format(ct))

            # Testing!!
            # TODO: Fix this!
            # print('(1) gal gsparams = {}'.format(gal.getGSParams()))
            # print('(2) gal._gsparams = {}'.format(gal._gsparams))
            # gal._gsparams = gsp
            # print('(3) gal gsparams = {}'.format(gal.getGSParams()))
            # print('(4) gal._gsparams = {}'.format(gal._gsparams))

            # Give intrinsic shape
            gal = gal.shear(g1=g1, g2=g2)

            # TODO: Implement shearing of intrinsic shape
            # now shear it
            # if 'shear' in pars:
            #     shear=pars['shear']
            #     gal = gal.shear(g1=shear.g1, g2=shear.g2)

            # Add galaxy in given band to list
            gsobjects.append(gal)

        # NOTE: If multiple bands were passed, the fluxes are simply added together.
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
        # NB: I have noticed an occasional memory issue with N>~200 galaxies. This may be related to the serialization
        # issues Mike talks about in scene.py
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

    def getNObjects(self):
        # Used by input/logger methods
        return self.nobjects

    def getNTot(self):
        # Used by input/logger methods
        return self.ntotal

    def getCatalog(self):
        return self.catalog

    # TODO: Write remaining `get` methods once saved columns are determined

    #------------------------------------------------------------------------------------------------

    # Since makeGalaxies is a function, not a class, it needs to use an unconventional location for defining
    # certain config parameters.
    makeGalaxies._req_params = {}
    makeGalaxies._opt_params = { "index" : int,
                               "n_random": int
                             }
    makeGalaxies._single_params = []
    makeGalaxies._takes_rng = True

#####------------------------------------------------------------------------------------------------

class BalGMixCM(ngmix.gmix.GMixCM):
    '''
    This is a slightly modified version of Erin Sheldon's GmixCM class from ngmix.gmix. Rather than
    updating his repo for a few small changes, we overwrite the desired methods here.
    '''

    def make_galsim_object(self, gsparams=None):
        '''
        Make a galsim representation for the gaussian mixture.
        NOTE: The only difference from Erin's original code is the type checking and passing
        of gsparams. This is especially useful for increasing fft sizes.
        '''

        # Accept optional gsparams for GSObject construction
        if gsparams and ( type(gsparams) is not galsim._galsim.GSParams ):
            if type(gsparams) is dict:
                # Convert to actual gsparams object
                gsparams = galsim.GSParams(**gsparams)
            else:
                raise TypeError('Only `dict` and `galsim.GSParam` types allowed for'
                                'gsparams; input has {} type.'.format(type(gsparams)))

        data = self._get_gmix_data()

        row, col = self.get_cen()

        gsobjects = []
        for i in xrange(len(self)):
            flux = data['p'][i]
            T = data['irr'][i] + data['icc'][i]
            # assert(T!=0.0)
            if T==0:
                # TODO: Explore this issue more! T is being stored as nonzero -
                #       what is causing this?
                continue
            e1 = (data['icc'][i] - data['irr'][i])/T
            e2 = 2.0*data['irc'][i]/T

            rowshift = data['row'][i]-row
            colshift = data['col'][i]-col

            g1,g2 = ngmix.shape.e1e2_to_g1g2(e1,e2)

            Tround = ngmix.moments.get_Tround(T, g1, g2)
            sigma_round = np.sqrt(Tround/2.0)

            gsobj = galsim.Gaussian(flux=flux, sigma=sigma_round,gsparams=gsparams)

            gsobj = gsobj.shear(g1=g1, g2=g2)
            gsobj = gsobj.shift(colshift, rowshift)

            gsobjects.append( gsobj )

        gs_obj = galsim.Add(gsobjects)

        return gs_obj

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
galsim.config.RegisterInputType('ngmix_catalog', ngmixCatalogLoader(ngmixCatalog, has_nobj=True))

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
    # try:
    ngmix_gal = ngmix_cat.makeGalaxies(index=index,gsparams=gsparams)
    # except:
        # print('\nNOT WORKING!!!\n')

    # The second item is "safe", a boolean that declares whether the returned value is
    # safe to save and use again for later objects (which is not the case for ngmixGalaxies).
    return ngmix_gal, False

# Register this builder with the config framework:
galsim.config.RegisterObjectType('ngmixGalaxy', BuildNgmixGalaxy, input_type='ngmix_catalog')

