import galsim
import os

class COSMOSChromaticCatalog(galsim.scene.COSMOSCatalog):
    '''
    Modified COSMOSCatalog class to allow chromatic objects when using configs.
    '''
    _req_params = {}
    _opt_params = { 'file_name' : str, 'sample' : str, 'dir' : str,
                    'preload' : bool, 'use_real' : bool,
                    'exclusion_level' : str, 'min_hlr' : float, 'max_hlr' : float,
                    'min_flux' : float, 'max_flux' : float, 'use_filter_tables' : bool,
                    'filter_dir' : str
                  }
    _single_params = []
    _takes_rng = False

    def __init__(self, file_name=None, sample=None, image_dir=None, dir=None, preload=False,
                 noise_dir=None, use_real=True, exclusion_level='marginal', min_hlr=0, max_hlr=0.,
                 min_flux=0., max_flux=0., _nobjects_only=False, exclude_bad=None,
                 exclude_fail=None, use_filter_tables=False, filter_dir=None):

        # Most type checking of use_filter_tables and filter_dir currently done in balrog_injection.py
        if not os.path.isdir(filter_dir):
            raise ValueError('Input {} for filter_dir is not a directory!'.format(filter_dir))

        return super(COSMOSChromaticCatalog, self).__init__(file_name=file_name,
                                                            sample=sample,
                                                            image_dir=image_dir,
                                                            dir=dir,
                                                            preload=preload,
                                                            noise_dir=noise_dir,
                                                            use_real=use_real,
                                                            exclusion_level=exclusion_level,
                                                            min_hlr=min_hlr,
                                                            max_hlr=max_hlr,
                                                            exclude_fail=exclude_fail)

    def makeGalaxy(self, index=None, gal_type=None, chromatic=False, noise_pad_size=5,
                   deep=False, sersic_prec=0.05, rng=None, n_random=None, gsparams=None):

        self.gal_type = gal_type

        return super(COSMOSChromaticCatalog, self).makeGalaxy(index=index,
                                      gal_type=gal_type,
                                      chromatic=chromatic,
                                      noise_pad_size=noise_pad_size,
                                      deep=deep,
                                      sersic_prec=sersic_prec,
                                      rng=rng,
                                      n_random=n_random,
                                      gsparams=gsparams)

    def _makeSingleGalaxy(cosmos_catalog, index, gal_type, noise_pad_size=5, deep=False,
                          rng=None, sersic_prec=0.05, gsparams=None, chromatic=False,
                          use_real=True):
        # A static function that mimics the functionality of COSMOSCatalog.makeGalaxy()
        # for single index and chromatic=False.
        # The only point of this class is to circumvent some pickling issues when using
        # config objects with type : COSMOSGalaxy.  It's a staticmethod, which means it
        # cannot use any self attributes.  Just methods.  (Which also means we can use it
        # through a proxy COSMOSCatalog object, which we need for the config layer.)

        if not use_real:
            if gal_type is None:
                gal_type = 'parametric'
            elif gal_type != 'parametric':
                raise ValueError("Only 'parametric' galaxy type is allowed when use_real == False")
        else:
            if gal_type is None:
                gal_type = 'real'

        if gal_type not in ['real', 'parametric']:
            raise ValueError("Invalid galaxy type %r"%gal_type)
        # We'll set these up if and when we need them.
        self._bandpass = None
        self._sed = None

        # Make rng if we will need it.
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

        if hasattr(index, '__iter__'):
            indices = index
        else:
            indices = [index]

        # Check whether this is a COSMOSCatalog meant to represent real or parametric objects, then
        # call the appropriate helper routine for that case.
        if gal_type == 'real':
            if chromatic:
                raise RuntimeError("Cannot yet make real chromatic galaxies!")
            gal_list = self._makeReal(indices, noise_pad_size, rng, gsparams)
        else:
            # If no pre-selection was done based on radius or flux, then we won't have checked
            # whether we're using the old or new catalog (the latter of which has a lot of
            # precomputations done).  Just in case, let's check here, though it does seem like a bit
            # of overkill to emit this warning each time.
            if 'hlr' not in self.param_cat.dtype.names:  # pragma: no cover
                import warnings
                warnings.warn(
                    'You seem to have an old version of the COSMOS parameter file. '+
                    'Please run `galsim_download_cosmos -s %s` '%self.use_sample+
                    'to re-download the COSMOS catalog '+
                    'and take advantage of pre-computation of many quantities..')

            gal_list = self._makeParametric(indices, chromatic, sersic_prec, gsparams)

        # If trying to use the 23.5 sample and "fake" a deep sample, rescale the size and flux as
        # suggested in the GREAT3 handbook.
        if deep:
            if self.use_sample == '23.5':
                # Rescale the flux to get a limiting mag of 25 in F814W when starting with a
                # limiting mag of 23.5.  Make the galaxies a factor of 0.6 smaller and appropriately
                # fainter.
                flux_factor = 10.**(-0.4*1.5)
                size_factor = 0.6
                gal_list = [ gal.dilate(size_factor) * flux_factor for gal in gal_list ]
            elif self.use_sample == '25.2':
                import warnings
                warnings.warn(
                    'Ignoring `deep` argument, because the sample being used already '+
                    'corresponds to a flux limit of F814W<25.2')
            else:
                import warnings
                warnings.warn(
                    'Ignoring `deep` argument, because the sample being used does not '+
                    'corresponds to a flux limit of F814W<23.5')

        # Store the orig_index as gal.index regardless of whether we have a RealGalaxy or not.
        # It gets set by _makeReal, but not by _makeParametric.
        # And if we are doing the deep scaling, then it gets messed up by that.
        # So just put it in here at the end to be sure.
        for gal, idx in zip(gal_list, indices):
            gal.index = self.orig_index[idx]
            if hasattr(gal, 'original'): gal.original.index = self.orig_index[idx]

        if hasattr(index, '__iter__'):
            return gal_list
        else:
            return gal_list[0]
    # def _makeSingleGalaxy(cosmos_catalog, index, gal_type, noise_pad_size=5, deep=False,
    #                       rng=None, sersic_prec=0.05, gsparams=None, chromatic=False):
    #     # A static function that mimics the functionality of COSMOSCatalog.makeGalaxy()
    #     # for single index and chromatic=False.
    #     # The only point of this class is to circumvent some pickling issues when using
    #     # config objects with type : COSMOSGalaxy.  It's a staticmethod, which means it
    #     # cannot use any self attributes.  Just methods.  (Which also means we can use it
    #     # through a proxy COSMOSCatalog object, which we need for the config layer.)

    #     if not cosmos_catalog.canMakeReal():
    #         if gal_type is None:
    #             gal_type = 'parametric'
    #         elif gal_type != 'parametric':
    #             raise ValueError("Only 'parametric' galaxy type is allowed when use_real == False")

    #     if gal_type not in ['real', 'parametric']:
    #         raise ValueError("Invalid galaxy type %r"%gal_type)

    #     if gal_type == 'real' and rng is None:
    #         rng = galsim.BaseDeviate()

    #     if gal_type == 'real':
    #         real_params = cosmos_catalog.getRealParams(index)
    #         gal = galsim.RealGalaxy(real_params, noise_pad_size=noise_pad_size, rng=rng,
    #                                 gsparams=gsparams)
    #     else:
    #         record = cosmos_catalog.getParametricRecord(index)
    #         gal = COSMOSChromaticCatalog._buildParametric(record, sersic_prec, gsparams, chromatic=chromatic)

    #     # If trying to use the 23.5 sample and "fake" a deep sample, rescale the size and flux as
    #     # suggested in the GREAT3 handbook.
    #     if deep:
    #         if self.use_sample == '23.5':
    #             # Rescale the flux to get a limiting mag of 25 in F814W when starting with a
    #             # limiting mag of 23.5.  Make the galaxies a factor of 0.6 smaller and appropriately
    #             # fainter.
    #             flux_factor = 10.**(-0.4*1.5)
    #             size_factor = 0.6
    #             gal = gal.dilate(size_factor) * flux_factor
    #         else:
    #             import warnings
    #             warnings.warn(
    #                 'Ignoring `deep` argument, because the sample being used already '+
    #                 'corresponds to a flux limit of F814W<25.2')

    #     # Store the orig_index as gal.index, since the above RealGalaxy initialization just sets it
    #     # as 0.  Plus, it isn't set at all if we make a parametric galaxy.  And if we are doing the
    #     # deep scaling, then it gets messed up by that.  If we have done some transformations, and
    #     # are also doing later transformation, it will take the `original` attribute that is already
    #     # there.  So having `index` doesn't help, and we also need `original.index`.
    #     gal.index = cosm
    #     if hasattr(gal, 'original'):
    #         gal.original.index = cosmos_catalog.getOrigIndex(index)

    #     return gal

    def getCatalog(self, gal_type='parametric'):
        if gal_type == 'real':
            return self.real_cat
        elif gal_type == 'parametric':
            return self.param_cat
        else:
            raise ValueError('gal_type must be real or parametric!')

    def getUseReal(self):
        return self.use_real

    # We add 'chromatic' as an optional parameter:
    makeGalaxy._opt_params = { "index" : int,
                               "gal_type" : str,
                               "noise_pad_size" : float,
                               "deep" : bool,
                               "sersic_prec": float,
                               "n_random": int,
                               "chromatic" : bool
                             }

    # The rest are the same as COSMOSCatalog:
    makeGalaxy._req_params = {}
    makeGalaxy._single_params = []
    makeGalaxy._takes_rng = True
