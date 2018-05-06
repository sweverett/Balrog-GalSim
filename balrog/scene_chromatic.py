import galsim

class COSMOSChromaticCatalog(galsim.scene.COSMOSCatalog):
    '''
    Modified COSMOSCatalog class to allow chromatic objects when using configs.
    '''

    def makeGalaxy(self, index=None, gal_type=None, chromatic=False, noise_pad_size=5,
                   deep=False, sersic_prec=0.05, rng=None, n_random=None, gsparams=None):

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
                          rng=None, sersic_prec=0.05, gsparams=None, chromatic=False):
        # A static function that mimics the functionality of COSMOSCatalog.makeGalaxy()
        # for single index and chromatic=False.
        # The only point of this class is to circumvent some pickling issues when using
        # config objects with type : COSMOSGalaxy.  It's a staticmethod, which means it
        # cannot use any self attributes.  Just methods.  (Which also means we can use it
        # through a proxy COSMOSCatalog object, which we need for the config layer.)

        if not cosmos_catalog.canMakeReal():
            if gal_type is None:
                gal_type = 'parametric'
            elif gal_type != 'parametric':
                raise ValueError("Only 'parametric' galaxy type is allowed when use_real == False")

        if gal_type not in ['real', 'parametric']:
            raise ValueError("Invalid galaxy type %r"%gal_type)

        if gal_type == 'real' and rng is None:
            rng = galsim.BaseDeviate()

        if gal_type == 'real':
            real_params = cosmos_catalog.getRealParams(index)
            gal = galsim.RealGalaxy(real_params, noise_pad_size=noise_pad_size, rng=rng,
                                    gsparams=gsparams)
        else:
            record = cosmos_catalog.getParametricRecord(index)
            gal = COSMOSChromaticCatalog._buildParametric(record, sersic_prec, gsparams, chromatic=chromatic)

        # If trying to use the 23.5 sample and "fake" a deep sample, rescale the size and flux as
        # suggested in the GREAT3 handbook.
        if deep:
            if self.use_sample == '23.5':
                # Rescale the flux to get a limiting mag of 25 in F814W when starting with a
                # limiting mag of 23.5.  Make the galaxies a factor of 0.6 smaller and appropriately
                # fainter.
                flux_factor = 10.**(-0.4*1.5)
                size_factor = 0.6
                gal = gal.dilate(size_factor) * flux_factor
            else:
                import warnings
                warnings.warn(
                    'Ignoring `deep` argument, because the sample being used already '+
                    'corresponds to a flux limit of F814W<25.2')

        # Store the orig_index as gal.index, since the above RealGalaxy initialization just sets it
        # as 0.  Plus, it isn't set at all if we make a parametric galaxy.  And if we are doing the
        # deep scaling, then it gets messed up by that.  If we have done some transformations, and
        # are also doing later transformation, it will take the `original` attribute that is already
        # there.  So having `index` doesn't help, and we also need `original.index`.
        gal.index = cosm
        if hasattr(gal, 'original'):
            gal.original.index = cosmos_catalog.getOrigIndex(index)

        return gal

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
