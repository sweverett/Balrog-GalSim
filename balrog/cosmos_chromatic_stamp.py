import galsim
import logging
import numpy as np

import pudb

class COSMOSChromaticStampBuilder(galsim.config.StampBuilder):
    def setup(self, config, base, xsize, ysize, ignore, logger):
        """
        Do the initialization and setup for building a postage stamp.

        In the base class, we check for and parse the appropriate size and position values in
        config (aka base['stamp'] or base['image'].

        Values given in base['stamp'] take precedence if these are given in both places (which
        would be confusing, so probably shouldn't do that, but there might be a use case where it
        would make sense).

        @param config       The configuration dict for the stamp field.
        @param base         The base configuration dict.
        @param xsize        The xsize of the image to build (if known).
        @param ysize        The ysize of the image to build (if known).
        @param ignore       A list of parameters that are allowed to be in config that we can
                            ignore here. i.e. it won't be an error if these parameters are present.
        @param logger       If given, a logger object to log progress.

        @returns xsize, ysize, image_pos, world_pos
        """
        # .. Do any custom setup you need to do.
        # Probably want to call the base class setup function to do the normal determination
        # of the size and position values.

        # Extra processing of 'bandpass' argument
        # Most needed type-checking is done in galsim.bandpass
        self._req_bp_fields = ['throughput', 'wave_type']
        self._opt_bp_fields = ['red_limit', 'blue_limit', 'zeropoint']
        try:
            bp = config['bandpass']
            for req in self._req_bp_fields:
                if req not in bp.keys():
                    raise ValueError('Must pass field {} for a bandpass object!'.format(req))
            # for opt in self._opt_bp_fields:
            #     if opt not in bp.keys():
            #         config['bandpass'][opt] = None
            for key in bp.keys():
                if key not in (self._req_bp_fields+self._opt_bp_fields):
                    raise ValueError('Field {} is not a valid entry for a bandpass!'.format(key))
        except KeyError:
            raise KeyError('`bandpass` is a required field for a COSMOSChromatic stamp!')

        extra_ignore = ignore + ['bandpass']
        return super(self.__class__, self).setup(config, base, xsize, ysize, extra_ignore, logger)

    def updateSkip(self, prof, image, method, offset, config, base, logger):
        """Before drawing the profile, see whether this object can be trivially skipped.

        The base method checks if the object is completely off the main image, so the
        intersection bounds will be undefined.  In this case, don't bother drawing the
        postage stamp for this object.

        @param prof         The profile to draw.
        @param image        The image onto which to draw the profile (which may be None).
        @param method       The method to use in drawImage.
        @param offset       The offset to apply when drawing.
        @param config       The configuration dict for the stamp field.
        @param base         The base configuration dict.
        @param logger       If given, a logger object to log progress.

        @returns whether to skip drawing this object.
        """

        # NOTE: There are currently unresolved issues with the image size checking of chromatic
        # objects. For now, we ignore any possible speed increases and skip the check.
        if isinstance(prof, galsim.ChromaticObject):
            return False

        if prof is not None and base.get('current_image',None) is not None:
            if image is None:
                prof = base['wcs'].toImage(prof, image_pos=base['image_pos'])
                # NOTE: Old version:
                # N = prof.getGoodImageSize(1.)
                if isinstance(prof, galsim.GSObject):
                    N = prof.getGoodImageSize(1.)
                elif isinstance(prof, galsim.ChromaticObject):
                    # TODO: Finish implementation
                    pass
                    # pudb.set_trace()
                    # # Find the suggested image size for each object given the choice of scale, and use the
                    # # maximum just to be safe.
                    # # print '\nprof.original = {}'.format(prof.original)
                    # print '\nprof.original.obj_list = {}'.format(prof.original.obj_list)
                    # # print '\nprof.objlist = {}'.format(prof.original.obj_list)
                    # obj_list = prof.original.obj_list
                    # possible_im_sizes = []
                    # for obj in obj_list:
                    #     print '\n obj : {}'.format(obj)
                    #     possible_im_sizes.append([ ob.getGoodImageSize(1.) for ob in obj])
                    # print 'possible_im_sizes : {}'.format(possible_im_sizes)
                    # N = np.max(possible_im_sizes)
                N += 2 + int(np.abs(offset.x) + np.abs(offset.y))
                bounds = galsim._BoundsI(1,N,1,N)
            else:
                bounds = image.bounds

            # Set the origin appropriately
            stamp_center = base['stamp_center']
            if stamp_center:
                bounds = bounds.shift(stamp_center - bounds.center)
            else:
                bounds = bounds.shift(base.get('image_origin',galsim.PositionI(1,1)) -
                                      galsim.PositionI(bounds.xmin, bounds.ymin))

            overlap = bounds & base['current_image'].bounds
            if not overlap.isDefined():
                logger.info('obj %d: skip drawing object because its image will be entirely off '
                            'the main image.', base['obj_num'])
                return True

        return False

    def draw(self, prof, image, method, offset, config, base, logger, **kwargs):
        """Draw the profile on the postage stamp image.
        This is a slightly modified version of `stamp.DrawBasic()` which allows drawing
        of chromatic objects.

        @param prof         The profile to draw.
        @param image        The image onto which to draw the profile (which may be None).
        @param method       The method to use in drawImage.
        @param offset       The offset to apply when drawing.
        @param config       The configuration dict for the stamp field.
        @param base         The base configuration dict.
        @param logger       If given, a logger object to log progress.

        @returns the resulting image
        """
        # ... draw prof onto the given image (making a new Image if necessary)
        if prof is None:
            return image
        else:
            logger = galsim.config.LoggerWrapper(logger)
            # Setup the kwargs to pass to drawImage
            # (Start with any additional kwargs given as extra kwargs to DrawBasic and add to it.)
            kwargs['image'] = image
            kwargs['offset'] = offset
            kwargs['method'] = method
            if 'wmult' in config and 'wmult' not in kwargs: # pragma: no cover
                kwargs['wmult'] = galsim.config.ParseValue(config, 'wmult', base, float)[0]
            if 'wcs' not in kwargs and 'scale' not in kwargs:
                kwargs['wcs'] = base['wcs'].local(image_pos = base['image_pos'])
            if method == 'phot' and 'rng' not in kwargs:
                kwargs['rng'] = galsim.config.GetRNG(config, base, logger, "method='phot'")

            # Check validity of extra phot options:
            max_extra_noise = None
            if 'n_photons' in config and 'n_photons' not in kwargs:
                if method != 'phot':
                    raise AttributeError('n_photons is invalid with method != phot')
                if 'max_extra_noise' in config:
                    logger.warning(
                        "Both 'max_extra_noise' and 'n_photons' are set in config dict, "+
                        "ignoring 'max_extra_noise'.")
                kwargs['n_photons'] = galsim.config.ParseValue(config, 'n_photons', base, int)[0]
            elif 'max_extra_noise' in config:
                max_extra_noise = galsim.config.ParseValue(config, 'max_extra_noise', base, float)[0]
                if method != 'phot' and max_extra_noise is not None:
                    raise AttributeError('max_extra_noise is invalid with method != phot')

            if 'poisson_flux' in config and 'poisson_flux' not in kwargs:
                if method != 'phot':
                    raise AttributeError('poisson_flux is invalid with method != phot')
                kwargs['poisson_flux'] = galsim.config.ParseValue(config, 'poisson_flux', base, bool)[0]

            if max_extra_noise is not None and 'max_extra_noise' not in kwargs:
                if max_extra_noise < 0.:
                    raise ValueError("image.max_extra_noise cannot be negative")
                if 'image' in base and 'noise' in base['image']:
                    noise_var = galsim.config.CalculateNoiseVariance(base)
                else:
                    raise AttributeError("Need to specify noise level when using max_extra_noise")
                if noise_var < 0.:
                    raise ValueError("noise_var calculated to be < 0.")
                max_extra_noise *= noise_var
                kwargs['max_extra_noise'] = max_extra_noise

            if logger.isEnabledFor(logging.DEBUG):
                # Don't output the full image array.  Use str(image) for that kwarg.
                alt_kwargs = dict([(k,str(kwargs[k]) if isinstance(kwargs[k],galsim.Image) else kwargs[k])
                                for k in kwargs])
                logger.debug('obj %d: drawImage kwargs = %s',base.get('obj_num',0), alt_kwargs)
                logger.debug('obj %d: prof = %s',base.get('obj_num',0),prof)
            try:
                # NOTE: Old version:
                # image = prof.drawImage(**kwargs)
                if isinstance(prof, galsim.GSObject):
                    image = prof.drawImage(**kwargs)
                elif isinstance(prof, galsim.ChromaticObject):
                    bp = {}
                    for key in (self._req_bp_fields+self._opt_bp_fields):
                        try:
                            bp[key] = config['bandpass'][key]
                        except KeyError:
                            bp[key] = None

                    bandpass = galsim.Bandpass(blue_limit=bp['blue_limit'], red_limit=bp['red_limit'],
                                            wave_type=bp['wave_type'], throughput=bp['throughput'],
                                            zeropoint=bp['zeropoint'])

                    image = prof.drawImage(bandpass=bandpass, **kwargs)

            except Exception as e: # pragma: no cover
                logger.debug('obj %d: prof = %r', base.get('obj_num',0), prof)
                raise
            return image

galsim.config.RegisterStampType('COSMOSChromatic', COSMOSChromaticStampBuilder())
