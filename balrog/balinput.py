import galsim
import galsim.config.input as gsinput
import os
from copy import deepcopy

# Balrog files
import config as Config
import tile as Tile
import balobject
import filters

# import pudb

class BalInput(object):
    def __init__(self, input_type, gsconfig, indx=None, tilename=None):
        if not isinstance(input_type, basestring):
            raise TypeError('`input_type` must be a string!')
        self.input_type = input_type

        if not isinstance(gsconfig, dict):
            raise TypeError('`gsconfig` must be a dict!')
        self.gsconfig = deepcopy(gsconfig)
        self.config = self.gsconfig['input'][input_type]

        if (indx is not None) and (not isinstance(indx, int)):
            raise TypeError('`indx` must be a int!')
        self.indx = indx

        if (tilename is not None) and (not isinstance(tilename, str)):
            raise TypeError('`tilename` must be a str!')
        self.tilename = tilename

        try:
            self.data_version = self.gsconfig['image']['version']
        except KeyError:
            self.data_version= None

        # Not all inputs will have certain quantities
        self.input_zp = None
        self.needs_band = False
        self.is_de_reddened = False

        # ...

        return

    def update_tile(self, new_tile):
        assert isinstance(new_tile, Tile.Tile)
        self.tilename = new_tile.tile_name

        return

class InputCatalog(BalInput):

    def __init__(self, input_type, gsconfig, indx, tilename=None):
        super(InputCatalog, self).__init__(input_type, gsconfig, indx=None, tilename=tilename)

        return

    def _register_input_type(self, input_type, CATMOD, LOADMOD, has_nobj=False):
        galsim.config.RegisterInputType(input_type, LOADMOD(CATMOD, has_nobj=has_nobj))

        return

    def _setup_proxy_catalog(self, save_proxy=False, get_subtype=False):
        # Some input catalogs have additional subtypes, like ngmix_catalog
        # Sometimes need to have this information for truth catalog updates, but
        # can only do it here if save_proxy is False

        conf = deepcopy(self.gsconfig)

        # Need to remove other input types as they may not have been registered
        # by GalSim yet
        for inpt in conf['input']:
            if inpt != self.input_type:
                conf['input'].pop(inpt)

        galsim.config.ProcessInput(conf)

        try:
            cat_proxy = conf['input_objs'][self.input_type][0]
        except KeyError:
            # Internal field name changed in newer versions of GalSim
            cat_proxy = conf['_input_objs'][self.input_type][0]

        if save_proxy is True:
            self.proxy = cat_proxy
        else:
            self.proxy = None
        self.cat = cat_proxy.getCatalog()
        self.nobjects = cat_proxy.getNObjects()

        if get_subtype is True:
            self.sub_type = cat_proxy.getSubtype()
        else:
            self.sub_type = None

        # Need to load in additional parametric catalog for some input types
        try:
            self.parametric_cat = cat_proxy.getParamCatalog()
        except AttributeError:
            self.parametric_cat = None

        return

    def generate_inj_catalog(self, config, tile, realization, mixed, mixed_grid=None):
        inp_type = self.input_type
        inj_type = self.inj_type
        sub_type = self.sub_type

        self.inj_catalog = balobject.build_bal_inject_cat(inp_type,
                                                          inj_type,
                                                          sub_type,
                                                          tile,
                                                          needs_band=self.needs_band,
                                                          mixed=mixed)

        mixed_grid = self.inj_catalog.generate_objects(config, realization, mixed_grid=mixed_grid)

        return self.inj_catalog, mixed_grid

class DESInputCatalog(InputCatalog):

    def __init__(self, input_type, gsconfig, indx=None, tilename=None):
        super(DESInputCatalog, self).__init__(input_type, gsconfig, indx, tilename=tilename)

        if not (self.data_version == 'y3v02'):
            print('WARNING: Unrecognized data version of {} '.format(self.data_version) +
                  'Assumming a magnitude zeropoint of 30!')
            self.input_zp = 30.0

        self.needs_band = True

        try:
            self.is_de_reddened = self.gsconfig['input'][input_type]['de_redden']
        except:
            # Default is already set to False
            pass

        # This avoids a printed warning, and sets up the input correctly
        # as no bands are passed in bal_config
        self.gsconfig['input'][input_type]['bands'] = 'griz' # dummy bands

        return

class NGMIXInputCatalog(DESInputCatalog):
    def __init__(self, gsconfig, indx=None, tilename=None):
        super(NGMIXInputCatalog, self).__init__('ngmix_catalog', gsconfig, indx, tilename=tilename)

        from ngmix_catalog import ngmixCatalog as CATMOD
        from ngmix_catalog import ngmixCatalogLoader as LOADMOD

        # TODO: See if this is actually necessary as the corresopnding
        # catalog input may automatically register the objects
        self._register_input_type(self.input_type, CATMOD, LOADMOD, has_nobj=True)

        # TODO: Can we grab the injection type from the registered GS catalog?
        self.inj_type = 'ngmixGalaxy'

        # ngmix catalogs have additional subtypes, e.g. cm
        self._setup_proxy_catalog(get_subtype=True)

        return

class MEDSInputCatalog(DESInputCatalog):
    def __init__(self, gsconfig, indx=None, tilename=None):
        super(MEDSInputCatalog, self).__init__('meds_catalog', gsconfig, indx, tilename=tilename)

        from meds_catalog import MEDSCatalog as CATMOD
        from meds_catalog import MEDSCatalogLoader as LOADMOD

        # TODO: See if this is actually necessary as the corresopnding
        # catalog input may automatically register the objects
        self._register_input_type(self.input_type, CATMOD, LOADMOD, has_nobj=True)

        # TODO: Can we grab the injection type from the registered GS catalog?
        self.inj_type = 'MEDSGalaxy'

        self._setup_proxy_catalog()

        return

class DESStarInputCatalog(DESInputCatalog):
    def __init__(self, gsconfig, indx=None, tilename=None):
        # DES Star catalogs must be initialized with a tilename
        assert tilename is not None
        super(DESStarInputCatalog, self).__init__('des_star_catalog', gsconfig, indx,
                                                  tilename=tilename)

        from des_star_catalog import desStarCatalog as CATMOD
        from des_star_catalog import desStarCatalogLoader as LOADMOD

        # TODO: See if this is actually necessary as the corresopnding
        # catalog input may automatically register the objects
        self._register_input_type(self.input_type, CATMOD, LOADMOD, has_nobj=True)

        # TODO: Can we grab the injection type from the registered GS catalog?
        self.inj_type = 'desStar'

        if self.data_version == 'y3v02':
            import des_star_catalog
            valid_model_types = des_star_catalog.return_valid_model_types(
                                    data_version=self.data_version)
            try:
                self.star_model = self.config['model_type']
            except KeyError:
                raise KeyError('Must pass a model_type if using a des_star_catalog '
                               'for Balrog injections! See `des_star_catalog.py` for details.')

            if self.star_model not in valid_model_types:
                raise ValueError('The selected model_type {} '.format(self.star_model) +
                                  'is not a valid des_star_catalog model type!\n'
                                  'Valid types for data version {} are: {}'.format(
                                      self.data_version, valid_model_types))

            prefix = self.star_model.split('_')[0]
            if prefix == 'Model':
                self.base_percent = 100
            elif prefix == 'Extra':
                self.base_percent = int(self.star_model.split('_')[1])
            else:
                # Then something very strange happened! Should have been caught above
                raise ValueError

            self.gsconfig['input'][self.input_type]['tile'] = self.tilename

            # TODO: See if we can load this, rather than set explicitly (but true for y3v02)
            # Set input catalog zeropoint
            self.input_zp = 30.0

        else:
            # In future, can add updated input parsing
            raise ValueError('No input parsing defined for DES star catalogs for '
                             'data version {}'.format(self.data_version))

        self._setup_proxy_catalog()

        return

    def _setup_proxy_catalog(self):
        if self.gsconfig['image']['pos_sampling'][self.input_type]['type'].lower() == 'sahar':
            sp = True
        else:
            sp = False

        super(DESStarInputCatalog, self)._setup_proxy_catalog(save_proxy=sp)

        if sp is True:
            # NOTE: If all stars in Sahar's catalogs were guaranteed to be in the
            # unique tile region, then we would simply use the existing cat and nobjects.
            # However, her catalogs are structured such that they contain
            # stars in the unique region of other tiles. So we only grab the
            # ones inside the unique region.
            indices = self.proxy.indices_in_region([self.ramin, self.ramax],
                                                    [self.decmin, self.decmax],
                                                    boundary_cross=self.ra_boundary_cross)

            # Update star catalog
            self.cat = self.cat[indices]
            self.nobjects = len(indices)

        return

    def update_tile(self, new_tilename):
        super(DESStarInputCatalog, self).update_tile(new_tilename)

        # DES star catalogs are tile-specific, so load new ones!
        self.gsconfig['input'][self.input_type]['tile'] = self.tilename
        self._setup_proxy_catalog()

        return

# TODO: This should be filled with relevant construction info from Alex DW's udg_catalog class
class UDGInputCatalog(DESInputCatalog):
    def __init__(self, gsconfig, indx=None, tilename=None):
        super(UDGInputCatalog, self).__init__('udg_catalog', gsconfig, indx, tilename=tilename)

        return

class COSMOSInputCatalog(InputCatalog):
    def __init__(self, gsconfig, indx=None, tilename=None):
        super(COSMOSInputCatalog, self).__init__('cosmos_chromatic_catalog', gsconfig, indx,
                                                 tilename=tilename)

        from scene_chromatic import COSMOSChromaticCatalog as CATMOD
        from input_cosmos_chromatic import COSMOSChromaticLoader as LOADMOD

        # TODO: See if this is actually necessary as the corresopnding
        # catalog input may automatically register the objects
        self._register_input_type(self.input_type, CATMOD, LOADMOD, has_nobj=True)

        # TODO: Can we grab the injection type from the registered GS catalog?
        self.inj_type = 'COSMOSChromaticGalaxy'

        try:
            filter_dir = self.config['filter_dir']
            self.filter_dir = os.path.abspath(filter_dir)
        except KeyError:
            # Default is set to cwd in Filter()
            self.filter_dir = None
        try:
            use_filter_tables = self.config['use_filter_tables']
            if type(use_filter_tables) != bool:
                raise TypeError('The type of `use_filter_tables` must be a bool! '
                                'Was passed a {}'.format(type(use_filter_tables)))
            self.use_filter_tables = use_filter_tables

        except KeyError:
            # Default depends on other parameters
            if self.filter_dir is None:
                warnings.warn('Neither `filter_dir` nor `use_filter_tables` '
                              'passed for input cosmos_chromatic_catalog. '
                              'Using a constant throughput for COSMOS galaxies.')
                self.use_filter_tables = False
            else:
                self.use_filter_tables = True

        bands = gsconfig['image']['bands']
        self.filters = filters.Filters(bands,
                                       use_transmission_tables=self.use_filter_tables,
                                       filter_dir=self.filter_dir)

        # TODO: Check COSMOS zeropoint!
        self.input_zp = 25.94

        self._setup_proxy_catalog()

        return

def build_bal_input(input_type, gsconfig, indx=None, tilename=None):

    if input_type in BALROG_INPUT_TYPES:
        # User-defined input construction
        bal_input = BALROG_INPUT_TYPES[input_type](gsconfig, indx=indx, tilename=tilename)
    else:
        # Generic input construction
        if input_type not in gsinput.valid_input_types:
            raise ValueError('{} is not a native GalSim input type '.format(input_type) +
                'or a recognized Balrog input type. Make sure you have written '
                'and registered a valid GalSim input type')
        bal_input = BalInput(input_type, gsconfig, indx=indx, tilename=tilename)

    return bal_input

BALROG_INPUT_TYPES = {
    'ngmix_catalog' : NGMIXInputCatalog,
    'meds_catalog' : MEDSInputCatalog,
    'udg_catalog' : UDGInputCatalog,
    'des_star_catalog' : DESStarInputCatalog,
    'cosmos_chromatic_catalog' : COSMOSInputCatalog
    }
