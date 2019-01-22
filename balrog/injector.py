import galsim
import logging
import numpy as np

# Balrog files
import grid
import config as Config

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
                                 'version', 'run_name', 'bands', 'n_objects', 'n_realizations',
                                 'object_density', 'inj_objs_only', 'pos_sampling', 'realizations',
                                 'extinct_objs', 'rotate_objs']

        # There are additionally the keywords `N_{inj_type}` that we want to ignore
        for key in config:
            if 'N_' in key:
                extra_ignore.append(key)

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
        # single = [{'n_objects':int, 'object_density':float},
        #           {'n_realizations':int, 'realizations':[int, list, dict]}]
        # params = galsim.config.GetAllParams(config, base, req=req, opt=opt, single=single)[0]
        # print 'params=',params

        config = parse_bal_image_inputs(config, base)

        return full_xsize, full_ysize

    def buildImage(self, config, base, image_num, obj_num, logger):
        try:
            ioo = config['inj_objs_only']
            if (type(ioo) is bool) and (ioo is True):
                return super(AddOnImageBuilder, self).buildImage(config,
                                                                 base,
                                                                 image_num,
                                                                 obj_num,
                                                                 logger)
            elif (isinstance(ioo, dict)) and (ioo['value'] is True):
                # Still want to use existing image if changed to be BKG
                if (ioo['noise']) and ('BKG' in ioo['noise']):
                    return super(BalrogImageBuilder, self).buildImage(config,
                                                                      base,
                                                                      image_num,
                                                                      obj_num,
                                                                      logger)
                else:
                    return super(AddOnImageBuilder, self).buildImage(config,
                                                                     base,
                                                                     image_num,
                                                                     obj_num,
                                                                     logger)
            else:
                # Default is to add on top of initial images
                return super(BalrogImageBuilder, self).buildImage(config,
                                                                  base,
                                                                  image_num,
                                                                  obj_num,
                                                                  logger)

        except KeyError:
            # Default is to add on top of initial images
            return super(BalrogImageBuilder, self).buildImage(config,
                                                              base,
                                                              image_num,
                                                              obj_num,
                                                              logger)

galsim.config.RegisterImageType('Balrog', BalrogImageBuilder())

def parse_bal_image_inputs(config, base):
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
                if n_reals < len(config['realizations']):
                    raise ValueError('`n_realizations` cannot be smaller than len(realizations).')
                if n_reals < 1:
                    raise ValueError('`n_realizations` must be a positive integer!')
            else:
                raise TypeError('The value `n_realizations` must be an int!')
        except KeyError:
            # In this case, assume that n_realizations=len(realizations)
            config['n_realizations'] = len(config['realizations'])
    except KeyError:
        try:
            n_reals = config['n_realizations']
            # If it has been passed, warn user but still use
            print('DEPRECATED: `n_realizations` without `realizations` has been '
                  'deprecated. Please use argument `realizations` (or both) instead. '
                  'Will assume realizations=1,2,...,n_realizations.')
            if isinstance(n_reals, int):
                if n_reals < 1:
                    raise ValueError('`n_realizations` must be a positive integer!')
                config['realizations'] = [x for x in range(n_reals)]
            else:
                raise TypeError('The value `n_realizations` must be an int!')
        except KeyError:
            # Default is to use realization 0
            print('DEFAULT WARNING: No realization passed; using default of 0.')
            config['realizations'] = [0]
            config['n_realizations'] = 1

    # Process inputs 'n_objects' and 'object_density'
    Ninput = 0
    input_list = []
    for inpt in base['input'].keys():
        if inpt not in Config.BaseConfig()._non_inj_input_types:
            Ninput += 1
            input_list.append(inpt)
    try:
        nobjs = config['n_objects']
        if isinstance(nobjs, int):
            if nobjs < 1:
                raise ValueError('`n_objects` must be a positive integer!')
            # If a single int is passed, all objects get that #
            config['n_objects'] = nobjs * np.ones(Ninput)
        elif isinstance(nobjs, dict):
            for itype, n in nobjs.items():
                if itype not in input_list:
                    raise ValueError('`n_objects` field {}'.format(itype) +
                                     'is not a passed input type!')
                if n < 1:
                    raise ValueError('All `n_objects` values must be a positive integer!')
            for inpt in input_list:
                if inpt not in nobjs.keys():
                    # This is ok if pos_sampling is grid for this input type - this
                    # is checked later
                    nobjs[inpt] = None
            config['n_objects'] = nobjs
        else:
            raise TypeError('`n_objects` must be passed as an int or a dict!')
    except KeyError:
        config['n_objects'] = {inpt : None for inpt in input_list}

    # TODO: We should be able to combine these two together if we loop over arg/val
    # *and* int / (float/int)
    try:
        density = config['object_density']
        if isinstance(density, (int, float)):
            if density < 0:
                raise ValueError('`object_density` must be positive!')
            # If a single int is passed, all objects get that #
            config['object_density'] = density * np.ones(Ninput)
        elif isinstance(density, dict):
            for itype, d in density.items():
                if itype not in input_list:
                    raise ValueError('`object_density` field {}'.format(itype) +
                                     'is not a passed input type!')
                if d < 0:
                    raise ValueError('All `object_density` values must be positive!')
            for inpt in input_list:
                if inpt not in density.keys():
                    # This is ok if pos_sampling is grid for this input type - this
                    # is checked later
                    density[inpt] = None
            config['object_density'] = density
        else:
            raise TypeError('`object_density` must be passed as a float or a dict!')
    except KeyError:
        config['object_density'] = {inpt : None for inpt in input_list}

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

    # Process input 'rotate_objs'
    # NOTE: Just the one for now, but could expand in future
    # (e.g. rotate all by 45 degrees)
    valid_rot_types = ['Random']
    try:
        rotate = config['rotate_objs']
        if isinstance(rotate, basestring):
            if rotate not in valid_rot_types:
                raise ValueError('`{}` is not a valid rotation type! '.format(rotate) +
                                 'For now, only {} are valid.'.format(valid_rot_types))
            # TODO: Only temporary; can change once we add more rotation types
            config['rotate_objs'] = True
        elif not isinstance(rotate, bool):
            raise ValueError('The input `rotate_objs` must be a string or bool!')
    except KeyError:
        # TODO: Maybe come up with sensible default run name?
        #       Current metadata should provide enough info for now.
        config['rotate_objs'] = False

    # Process input 'extinct_objs'
    try:
        ext = config['extinct_objs']
        if not isinstance(ext, bool):
            raise TypeError('The input `extinct_objs` must be a bool!')
    except KeyError:
        config['extinct_objs'] = False

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
    # TODO: Clean up after testing!
    bg = grid.BaseGrid()
    bc = Config.BaseConfig()
    valid_pos_sampling = bc._valid_pos_sampling
    valid_grid_types = bg._valid_grid_types
    valid_mixed_types = bg._valid_mixed_types
    default_gs = 30 # arcsec

    try:
        ps = config['pos_sampling']
    except KeyError:
        # Most runs use uniform sampling
        config['pos_sampling'] = {inpt : {
            'type':'uniform', 'grid_spacing':None
            } for inpt in input_list}
        ps = config['pos_sampling']

    if isinstance(ps, basestring):
        # Then the string is the input type
        if ps not in valid_pos_sampling:
            raise ValueError('{} is not a valid position sampling method. '.format(ps) +
                            'Currently allowed methods are {}'.format(valid_pos_sampling))

        if ps in valid_grid_types:
            print('No grid spacing passed; using default of {} arcsecs'.format(default_gs))
            config['pos_sampling'] = {inpt : {'type' : ps, 'grid_spacing' : default_gs}
                                        for inpt in input_list}
        else:
            config['pos_sampling'] = {inpt : {'type' : ps} for inpt in input_list}
    elif isinstance(ps, dict):
        def _parse_pos_sampling_dict(ps):
            bg = grid.BaseGrid()
            bc = Config.BaseConfig()
            valid_pos_sampling = bc._valid_pos_sampling
            valid_grid_types = bg._valid_grid_types
            valid_mixed_types = bg._valid_mixed_types

            # NOTE: We actually don't check for 'grid_type' as only 1
            # of the mixed input types needs it
            req_mixed_keys = ['inj_frac'] #, 'grid_type']
            if ps['type'] in valid_mixed_types:
                for key in req_mixed_keys:
                    if key not in ps:
                        raise AttributeError('`pos_sampling` must contain the key '
                                                '{} for a Mixed sampling type!'.format(key))

            keys = ['type', 'grid_type', 'grid_spacing', 'rotate', 'offset',
                    'angle_unit', 'inj_frac']
            for key, val in ps.items():
                if key not in keys:
                    raise ValueError('{} is not a valid key for `pos_sampling`! '.format(key) +
                                        'You may only pass the keys {}'.format(keys))
                if key == 'grid_type':
                    if ps['type'] not in valid_mixed_types:
                        raise ValueError('`grid_type` is only a valid key for `pos_sampling` '
                                            'when using a `MixedGrid`!')
                    # TODO: We should add consistency checks here to make sure that MixedGrid
                    # params between input types are consistent

                if key == 'inj_frac':
                    if ps['type'] not in valid_mixed_types:
                        raise ValueError('`inj_frac` is only a valid key for `pos_sampling` '
                                            'when using a `MixedGrid`!')
                    if not isinstance(val, float):
                        raise TypeError('`inj_frac` must be a float!')
                    if val <= 0 or val >= 1:
                        raise ValueError('`inj_frac` must be between 0 and 1!')

                if key == 'grid_spacing':
                    if val < 0.0:
                        raise ValueError('grid_spacing of {} is invalid; '.format(val) +
                                            'must be positive!')
                if key == 'rotate':
                    if isinstance(val, str):
                        if val.lower() != 'random':
                            raise ValueError('{} is not a valid grid rotation type!'.format(val) +
                                                'Must be `Random` or a number.')
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
                    except KeyError:
                        # Default of rad
                        if isinstance(val, (int, float)):
                            if (val<0.0) or (val>2*np.pi):
                                raise ValueError('rotate value of {} rad is invalid!'.format(val))

                if key == 'offset':
                    if isinstance(val, str):
                        if val.lower() != 'random':
                            raise ValueError('{} is not a valid grid offset!'.format(val) +
                            'must be `Random` or an array.')
                    elif isinstance(val, list):
                        assert len(val)==2
                    else:
                        raise TypeError('grid offset of type {} is invalid!'.format(type(val)))

        if 'type' in ps.keys():
            if (ps['type'] in valid_grid_types) and ('grid_spacing' not in ps.keys()):
                print('No grid spacing passed; using default of {} arcsecs'.format(default_gs))
                config['pos_sampling']['grid_spacing'] = default_gs
            # Same type for all inputs
            _parse_pos_sampling_dict(ps)
            ps_temp = config['pos_sampling']
            config['pos_sampling'] = {inpt : ps_temp for inpt in input_list}

        else:
            for inpt in input_list:
                if (ps[inpt]['type'] in valid_grid_types) \
                and ('grid_spacing' not in ps[inpt].keys()):
                    print('No grid spacing passed; using default of '
                    '{} arcsecs'.format(default_gs))
                    config['pos_sampling'][inpt]['grid_spacing'] = default_gs
                _parse_pos_sampling_dict(ps[inpt])

    # Must have *exactly* one of `n_objects` or `object_density` for each input, unless using a grid
    for inpt in input_list:
        if (config['n_objects'][inpt] is not None) and (config['object_density'][inpt] is not None):
            raise ValueError('Only one of `n_objects` or `object_density` is allowed for '
                             'input {}, not both!'.format(inpt))
        elif (config['n_objects'][inpt] is None) and (config['object_density'][inpt] is None):
            if (config['pos_sampling'][inpt]['type'] not in bg._valid_grid_types) \
            and (config['pos_sampling'][inpt]['type'] not in bg._valid_mixed_types):
                raise ValueError('Must pass one of `n_objects` or `object_density` for '
                                 'input {} if not injecting on a grid!'.format(inpt))

    return config
