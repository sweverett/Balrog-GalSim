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

class BalrogImageBuilder(AddOnImageBuilder):

    def setup(self, config, base, image_num, obj_num, ignore, logger):
        ignore = ignore + ['bands', 'n_realizations', 'n_galaxies', 'gal_density', 'version']
        # full_xsize, full_ysize = super(BalrogImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)
        full_xsize, full_ysize = super(BalrogImageBuilder, self).setup(config, base, image_num, obj_num, ignore, logger)

        # Note that we only type-check, not set default values; GalSim has no knowledge of these parameters and
        # default behaviour is handled in balrog.py
        params = galsim.config.GetAllParams(config, base, ignore=ignore)[0]
        for arg in ['n_realizations', 'n_galaxies']:
            if arg in params:
                if params[arg] < 0 or (type(params[arg]) is not int):
                    raise ValueError(arg + ' must be a positive integer!')

        if 'n_galaxies' in params and 'gal_density' in params:
            raise AttributeError("Only one of `n_galaxies` or `gal_density` is allowed; not both.")

        return full_xsize, full_ysize

    def buildImage(self, config, base, image_num, obj_num, logger):
        return super(BalrogImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)

galsim.config.RegisterImageType('Balrog', BalrogImageBuilder())

