import galsim
import logging
import pdb
## Class for injecting simulated galaxies into pre-existing images.


class AddOnImageBuilder(galsim.config.image_scattered.ScatteredImageBuilder):
    
    def setup(self, config, base, image_num, obj_num, ignore, logger):
        ignore = ignore + ['initial_image']
        return super(AddOnImageBuilder, self).setup(config, base, image_num, obj_num, ignore, logger)

    
    def buildImage(self, config, base, image_num, obj_num, logger):
        im, cv = super(AddOnImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)
        initial_image_name = galsim.config.ParseValue(config,'initial_image', base, str)[0]
        initial_image = galsim.fits.read(initial_image_name)
        im += initial_image
        return im, cv

galsim.config.RegisterImageType('AddOn', AddOnImageBuilder())
# class StampBuilder(galsim.config.stamp.StampBuilder):

#     def setup(self, config, base, xsize, ysize, ignore, logger):
#         # ignore = ignore + [anything needed for ngmix c_model]
#         return super(StampBuilder, self).setup(config, base, xsize, ysize, ignore, logger)

#     def buildProfile(self, config, base, psf, gsparams, logger):
#         pass

class ngmixLoader(galsim.config.InputLoader):
    # Allow the user to not provide the image file.  In this case, we'll grab the wcs from the
    # config dict.
    def getKwargs(self, config, base, logger):
        req = { 'file_name' : str }
        opt = { 'dir' : str, 'catalog_file_name' : str }
        kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)

        if 'catalog_file_name' not in kwargs:
            if 'wcs' in base:
                kwargs['wcs'] = base['wcs']
            else:
                # Then we aren't doing normal config processing, so just use pixel scale = 1.
                kwargs['wcs'] = galsim.PixelScale(1.)

        return kwargs, safe

# First we need to add the class itself as a valid input_type.
galsim.config.RegisterInputType('des_ngmix', ngmixLoader(DES_Ngmix))

def buildDES_ngmix(config, base, ignore, gsparams, logger):
    """@brief Build a RealGalaxy type GSObject from user input.
    """
    des_ngmix = galsim.config.GetInputObj('des_ngmix', config, base, 'DES_Ngmix')

    req = { 'flux' : float , 'num' : int, 'image_pos' : galsim.PositionD }
    params, safe = galsim.config.GetAllParams(config, base, req=req, ignore=ignore)

    if 'image_pos' in params:
        image_pos = params['image_pos']
    elif 'image_pos' in base:
        image_pos = base['image_pos']
    else:
        raise ValueError("DES_Ngmix requested, but no image_pos defined in base.")

    # Convert gsparams from a dict to an actual GSParams object
    if gsparams: gsparams = galsim.GSParams(**gsparams)
    else: gsparams = None

    #psf = des_ngmix.getPSF(image_pos, gsparams=gsparams)
    # Because of serialization issues, the above call doesn't work.  So we need to
    # repeat the internals of getPSF here.
    # Also, this is why we have getSampleScale and getLocalWCS.  The multiprocessing.managers
    # stuff only makes available methods of classes that are proxied, not all the attributes.
    # So this is the only way to access these attributes.
    im = galsim.Image(des_ngmix.getPSFArray(image_pos))
    psf = galsim.InterpolatedImage(im, scale=des_ngmix.getSampleScale(), flux=1,
                                   x_interpolant=galsim.Lanczos(3), gsparams=gsparams)
    psf = des_ngmix.getLocalWCS(image_pos).toWorld(psf)

    if 'flux' in params:
        psf = psf.withFlux(params['flux'])

    # The second item here is "safe", a boolean that declares whether the returned value is
    # safe to save and use again for later objects.  In this case, we wouldn't want to do
    # that, since they will be at different positions, so the interpolated PSF will be different.
    return psf, False

# Register this builder with the config framework:
galsim.config.RegisterObjectType('DES_Ngmix', BuildDES_Ngmix, input_type='des_ngmix')

#----------------------------------------------------------------------------------------------
# Don't need an extra copy of this method; just for quick reference

def GetInputObj(input_type, config, base, param_name):
    """Get the input object needed for generating a particular value

    @param input_type   The type of input object to get
    @param config       The config dict for this input item
    @param base         The base config dict
    @param param_name   The type of value that we are trying to construct (only used for
                        error messages).
    """
    if 'input_objs' not in base or input_type not in base['input_objs']:
        raise ValueError("No input %s available for type = %s"%(input_type,param_name))

    if 'num' in config:
        num = galsim.config.ParseValue(config, 'num', base, int)[0]
    else:
        num = 0

    if num < 0:
        raise ValueError("Invalid num < 0 supplied for %s: num = %d"%(param_name,num))
    if num >= len(base['input_objs'][input_type]):
        raise ValueError("Invalid num supplied for %s (too large): num = %d"%(param_name,num))

    return base['input_objs'][input_type][num]
