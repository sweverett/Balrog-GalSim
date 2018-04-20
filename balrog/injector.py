import galsim
import logging

# Can use for debugging
# import pudb

## Class for injecting simulated galaxies into pre-existing images.


class AddOnImageBuilder(galsim.config.image_scattered.ScatteredImageBuilder):

    def setup(self, config, base, image_num, obj_num, ignore, logger):
        ignore = ignore + ['initial_image']
        return super(AddOnImageBuilder, self).setup(config, base, image_num, obj_num, ignore, logger)

    def buildImage(self, config, base, image_num, obj_num, logger):
        im, cv = super(AddOnImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)
        initial_image_name = galsim.config.ParseValue(config, 'initial_image', base, str)[0]
        initial_image = galsim.fits.read(initial_image_name)
        im += initial_image
        return im, cv

galsim.config.RegisterImageType('AddOn', AddOnImageBuilder())

# --------------------------------------------------------------------------------------------------

class BalrogImageBuilder(AddOnImageBuilder):

    def setup(self, config, base, image_num, obj_num, ignore, logger):
        extra_ignore = ignore + ['tile_list', 'geom_file', 'tile_dir', 'config_dir', 'psf_dir',
                                 'version', 'run_name', 'Ngals', 'Nstars', 'bands', 'n_galaxies',
                                 'n_realizations', 'gal_density', 'inj_objs_only', 'pos_sampling',
                                 'realizations']
        # full_xsize, full_ysize = super(BalrogImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)
        full_xsize, full_ysize = super(BalrogImageBuilder, self).setup(config, base, image_num,
                                                                       obj_num, extra_ignore, logger)

        # Note that we only type-check, not set default values; GalSim has no knowledge of these parameters and
        # default behaviour is handled in balrog.py
        opt = {'Ngals':int, 'Nstars':int, 'bands':str, 'n_realizations':int, 'n_galaxies':int,
                'gal_density':float, 'version':str, 'run_name':str, 'inj_objs_only':bool, 'pos_sampling':str}
        params = galsim.config.GetAllParams(config, base, opt=opt)[0]
        # TODO: Clean up and fix! Most params are not being passed as expected
        # print('config[inj_objs_only] = {}'.format(config['inj_objs_only']))
        # print('params = {}'.format(params))
        # print('config = {}'.format(config))
        # print('config[_get] = {}'.format(config['_get']))
        # print('base = {}'.format(base))
        for arg in ['n_realizations', 'n_galaxies']:
            if arg in params:
                if params[arg] < 0 or (type(params[arg]) is not int):
                    raise ValueError(arg + ' must be a positive integer!')

        if 'n_galaxies' in params and 'gal_density' in params:
            raise AttributeError('Only one of `n_galaxies` or `gal_density` is allowed; not both.')

        # It can be useful to test Balrog by injecting objects into blank images rather
        # than real ones. For now, this test does not generate any images itself but
        # rather ignores existing information during injection.
        if 'inj_objs_only' in params:
            if type(params['inj_objs_onlyl']) is not bool:
                raise ValueError('The field `inj_objs_only` must be set with a bool!')

        return full_xsize, full_ysize

    def buildImage(self, config, base, image_num, obj_num, logger):
        try:
            ioo = config['inj_objs_only']
            if (type(ioo) is bool) and (ioo is True):
                return super(AddOnImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)
            elif (isinstance(ioo, dict)) and (ioo['value'] is True):
                # Still want to use existing image if using BKG as noise
                if ioo['noise'] in ['BKG', 'BKG+noise']:
                    return super(BalrogImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)
                else:
                    return super(AddOnImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)
            else:
                # Default is to add on top of initial images
                return super(BalrogImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)

        except KeyError:
            # Default is to add on top of initial images
            return super(BalrogImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)

galsim.config.RegisterImageType('Balrog', BalrogImageBuilder())
