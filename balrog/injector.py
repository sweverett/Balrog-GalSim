import galsim
import logging

# Balrog files
import grid

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
        full_xsize, full_ysize = super(BalrogImageBuilder, self).setup(config, base, image_num,
                                                                    obj_num, extra_ignore, logger)

        # TODO: GetAllParams is failing for an unknown reason - possibly due to the inherited
        # constructor. Can do the type checking below for now
        # Note that we only type-check, not set default values; GalSim has no knowledge of these
        # parameters and default behaviour is handled in balrog_injection.py
        # print 'config=',config
        # print 'base=',base
        # req = {'bands'}
        # opt = {'Ngals':int, 'Nstars':int, 'n_realizations':int,
        #        'version':str, 'run_name':str, 'inj_objs_only':bool,
    #        'pos_sampling':str}
        # single = [{'n_galaxies':int, 'gal_density':float},
        #           {'n_realizations':int, 'realizations':[int, list, dict]}]
        # params = galsim.config.GetAllParams(config, base, req=req, opt=opt, single=single)[0]
        # print 'params=',params

        config = parse_bal_image_inputs(config)

        return full_xsize, full_ysize

    def buildImage(self, config, base, image_num, obj_num, logger):
        try:
            ioo = config['inj_objs_only']
            if (type(ioo) is bool) and (ioo is True):
                return super(AddOnImageBuilder, self).buildImage(config, base, image_num, obj_num, logger)
            elif (isinstance(ioo, dict)) and (ioo['value'] is True):
                # Still want to use existing image if changed to be BKG
                if (ioo['noise']) and ('BKG' in ioo['noise']):
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

def parse_bal_image_inputs(config):
    # Process input 'realizations' & 'n_realizations':
    try:
        reals = config['realizations']
        if isinstance(reals, int):
            config['realizations'] = np.array([reals])
        elif isinstance(reals, list):
            if not all(isinstance(x, int) for x in reals):
                raise TypeError('Passed realizations in list must be integers!')
        elif isinstance(reals, dict):
            try:
                max_real = reals['max']
                min_real = reals['min']
                if min_real >= max_real:
                    raise ValueError('Key `min_real` must be less than `max_real`!')
                if all(isinstance(x, int) for x in [min_real, max_real]):
                    config['realizations'] = [x for x in range(min_real, max_real+1)]
                else:
                    raise TypeError('The realization values `min` and `max` must be ints!')
            except KeyError as e:
                print(e + '\nThe only valid keys for `realizations` are `min` and `max`!')
        else:
            raise TypeError('`realizations` can only be an int, list, or dict!')

        # Now check if `n_realizations` was also passed
        # NOTE: If n_realizations>len(realizations), this indicates that a larger simulation
        # is being split up between multiple runs. The main impact is star generation for
        # Y3 DES star catalogs, as they cannot be shuffled in this case to ensure all stars
        # are injected w/o repeats.
        try:
            n_reals = config['n_realizations']
            if isinstance(n_reals, int):
                if n_reals < len(self.realizations):
                    raise ValueError('`n_realizations` cannot be smaller than len(realizations).')
            else:
                raise TypeError('The value `n_realizations` must be an int!')
        except KeyError:
            # In this case, assume that n_realizations=len(realizations)
            config['n_realizations'] = len(self.realizations)
    except KeyError:
        try:
            n_reals = config['n_realizations']
            # If it has been passed, warn user but still use
            warnings.warn('DEPRECATED: `n_realizations` without `realizations` has been '
                            'deprecated. Please use argument `realizations` (or both) instead. '
                            'Will assume realizations=1,2,...,n_realizations.')
            if isinstance(n_reals, int):
                config['realizations'] = [x for x in range(n_reals)]
            else:
                raise TypeError('The value `n_realizations` must be an int!')
        except KeyError:
            # Default is to use realization 0
            print('DEFAULT WARNING: No realization passed; using default of 0.')
            config['realizations'] = [0]
            config['n_realizations'] = 1

    for arg in ['n_realizations', 'n_galaxies', 'gal_density']:
        if arg not in config:
            config[arg] = None

    for arg in ['n_realizations', 'n_galaxies']:
        if config[arg] is not None:
            if (config[arg] < 0) or (not isinstance(config[arg], int)):
                raise ValueError(str(arg) + ' must be a positive integer!')

    # Process input 'gal_density':
    if config['gal_density'] is not None:
        gd = config['gal_density'] # arcmin^2
        if (not isinstance(gd, int)) and (not isinstance(gd, float)):
            raise TypeError('`gal_density` must be an int or float!')

    # Process input 'bands'
    try:
        bands = config['bands']
        for band in bands:
            if not isinstance(band, str):
                raise TypeError('Each passed band must be a string!')
    except KeyError:
        print 'config[bands]=',config['bands']
        raise KeyError('Must specify which bands are to be used for injection!')

    # Process input 'version'
    if 'version' not in config:
        # Warn user, but assume y3v02 for now
        print('DEFAULT WARNING: Data version not passed in config! Assuming y3v02.')
        config['data_version'] = 'y3v02'

    # Process input 'run_name'
    try:
        rname = config['run_name']
        if not isinstance(rname, basestring):
            raise ValueError('The input `run_name` must be a string!')
    except KeyError:
        # TODO: Maybe come up with sensible default run name?
        #       Current metadata should provide enough info for now.
        config['run_name'] = None

    # Process input 'inj_objs_only'. This is used to test Balrog injections on blank images
    try:
        inj_objs_only = config['inj_objs_only']
        if type(inj_objs_only) is bool:
            # Default is to include chip noise in injection
            config['inj_objs_only'] = {'value':inj_objs_only, 'noise':'CCD'}
        elif isinstance(inj_objs_only, dict):
            # Is likely an OrderedDict, so convert
            inj_objs_only = dict(inj_objs_only)
            ioo = {}
            keys = ['value', 'noise']
            valid_noise = ['CCD', 'BKG', 'BKG+CCD', 'BKG+RN', 'BKG+SKY', 'None', None]

            if 'noise' not in inj_objs_only:
                # Default is no noise
                inj_objs_only['noise'] = None
            for key, val in inj_objs_only.items():
                if key not in keys:
                    raise ValueError('{} is not a valid key for `inj_objs_only`!'.format(key)+
                                        ' You may only pass the keys {}'.format(keys))
                if (key == 'noise') and (val not in valid_noise):
                    raise ValueError('{} is not a valid value for the noise field!'.format(val)+
                                        ' You may only pass the values {}'.format(valid_noise))
                ioo[key] = val
            config['inj_objs_only'] = ioo
        else:
            raise ValueError('The field `inj_objs_only` must be set with a bool or dict!')
    except KeyError:
        # Most runs will add objects to existing images
        config['inj_objs_only'] = {'value':False, 'noise':None}

    # Process input 'pos_sampling'
    valid_pos_sampling = grid._valid_pos_sampling
    valid_grid_types = grid._valid_grid_types
    default_gs = 30 # arcsec
    try:
        ps = config['pos_sampling']
        if isinstance(ps, basestring):
            # Then the string is the input type
            if ps not in valid_pos_sampling:
                raise ValueError('{} is not a valid position sampling method. '.format(ps) +
                                'Currently allowed methods are {}'.format(valid_pos_sampling))

            if ps in valid_grid_types:
                print('No grid spacing passed; using default of {} arcsecs'.format(default_gs))
                config['pos_sampling'] = {'type' : ps, 'grid_spacing' : default_gs}
            else:
                config['pos_sampling'] = {'type' : ps}
        elif isinstance(ps, dict):
            if 'type' not in ps.keys():
                raise ValueError('If `pos_sampling` is passed as a dict, then must set a type!')
            if (ps in valid_grid_types) and ('grid_spacing' not in ps.keys()):
                print('No grid spacing passed; using default of {} arcsecs'.format(default_gs))
                ps['grid_spacing'] = default_gs
            keys = ['type', 'grid_spacing', 'rotate', 'offset', 'angle_unit']
            for key, val in ps.items():
                if key not in keys:
                    raise ValueError('{} is not a valid key for `pos_sampling`! '.format(key) +
                                        'You may only pass the keys {}'.format(keys))
                if key == 'grid_spacing':
                    if val < 0.0:
                        raise ValueError('grid_spacing of {} is invalid; '.format(val) +
                                            'must be positive!')
                if key == 'rotate':
                    if isinstance(val, str):
                        if val.lower() != 'random':
                            raise ValueError('{} is not a valid grid rotation type!'.format(val) +
                                                'Must be `Random` or a number.')
                    # config['pos_sampling'].update({'rotate':val})
                    try:
                        unit = ps['angle_unit'].lower()
                        if unit == 'deg':
                            if (val<0.0) or (val>360.0):
                                raise ValueError('rotate value of {} deg is invalid!'.format(val))
                        if unit == 'rad':
                            if (val<0.0) or (val>2*np.pi):
                                raise ValueError('rotate value of {} rad is invalid!'.format(val))
                            else:
                                raise ValueError('angle_unit of {} is invalid! '.format(unit) +
                                                    'only can pass `deg` or `rad`.')
                        # config['pos_sampling'].update({'angle_unit':unit})
                    except KeyError:
                        # Default of rad
                        if (val<0.0) or (val>2*np.pi):
                            raise ValueError('rotate value of {} rad is invalid!'.format(val))
                        # config['pos_sampling'].update({'angle_unit':'rad'})
                    # config['pos_sampling'].update({'rotate':val})

                if key == 'offset':
                    if isinstance(val, str):
                        if val.lower() != 'random':
                            raise ValueError('{} is not a valid grid offset!'.format(val) +
                            'must be `Random` or an array.')
                        # config['pos_sampling'].update({'offset':val})
                    elif isinstance(val, list):
                        assert len(val)==2
                        # config['pos_sampling'].update({'offset':val})
                    else:
                        raise TypeError('grid offset of type {} is invalid!'.format(type(val)))

                # config['pos_sampling'][key] = val

    except KeyError:
        # Most runs use uniform sampling
        config['pos_sampling']['type'] = 'uniform'
        config['pos_sampling']['grid_spacing'] = None

    # Must have *exactly* one of `n_galaxies` or `gal_density`, unless using a grid.
    if (config['n_galaxies'] is not None) and (config['gal_density'] is not None):
        raise ValueError('Only one of `n_galaxies` or `gal_density` is allowed; not both!')
    elif (config['n_galaxies'] is None) and (config['gal_density'] is None):
        if config['pos_sampling']['type'] not in grid._valid_grid_types:
            raise ValueError('Must pass one of `n_galaxies` or `gal_desnity` if not injecting '
                            'on a grid!')

    return config
