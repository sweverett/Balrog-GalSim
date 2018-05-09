import galsim
import logging
from galsim.config import RegisterInputType, InputLoader

# This file adds input type cosmos_chromatic_catalog and gsobject typs COSMOSGalaxyChromatic.

# This file along with scene_chromatic only exist to add chromatic functionality to config processing.
import scene_chromatic as sc

class COSMOSChromaticLoader(InputLoader):
    def setupImage(self, cosmos_cat, config, base, logger):
        if logger: # pragma: no cover
            # Only report as a warning the first time.  After that, use info.
            first = not base.get('_COSMOSChromaticLoader_reported_as_warning',False)
            base['_COSMOSChromaticLoader_reported_as_warning'] = True
            if first:
                log_level = logging.WARNING
            else:
                log_level = logging.INFO
            if 'input' in base:
                if 'cosmos_catalog_chromatic' in base['input']:
                    cc = base['input']['cosmos_catalog_chromatic']
                    if isinstance(cc,list): cc = cc[0]
                    out_str = ''
                    if 'sample' in cc:
                        out_str += '\n  sample = %s'%cc['sample']
                    if 'dir' in cc:
                        out_str += '\n  dir = %s'%cc['dir']
                    if 'file_name' in cc:
                        out_str += '\n  file_name = %s'%cc['file_name']
                    if out_str != '':
                        logger.log(log_level, 'Using user-specified COSMOSCatalogChromatic: %s',out_str)
            logger.info("file %d: COSMOS chromatic catalog has %d total objects; %d passed initial cuts.",
                        base['file_num'], cosmos_cat.getNTot(), cosmos_cat.getNObjects())
            if 'gal' in base and 'gal_type' in base['gal']:
                if base['gal']['gal_type']=='parametric':
                    logger.log(log_level,"Using parametric galaxies.")
                else:
                    logger.log(log_level,"Using real galaxies.")

RegisterInputType('cosmos_chromatic_catalog', COSMOSChromaticLoader(sc.COSMOSChromaticCatalog))

def _BuildCOSMOSChromaticGalaxy(config, base, ignore, gsparams, logger):
    cosmos_cat = galsim.config.GetInputObj('cosmos_chromatic_catalog', config, base, 'COSMOSChromaticGalaxy')

    ignore = ignore + ['num']

    # Special: if galaxies are selected based on index, and index is Sequence or Random, and max
    # isn't set, set it to nobjects-1.
    if 'index' in config:
        galsim.config.SetDefaultIndex(config, cosmos_cat.getNObjects())

    kwargs, safe = galsim.config.GetAllParams(config, base,
        req = sc.COSMOSChromaticCatalog.makeGalaxy._req_params,
        opt = sc.COSMOSChromaticCatalog.makeGalaxy._opt_params,
        single = sc.COSMOSChromaticCatalog.makeGalaxy._single_params,
        ignore = ignore)
    if gsparams: kwargs['gsparams'] = galsim.GSParams(**gsparams)

    if 'gal_type' not in kwargs:
        if cosmos_cat.getUseReal(): kwargs['gal_type'] = 'real'
        else: kwargs['gal_type'] = 'parametric'

    rng = None
    if 'index' not in kwargs:
        rng = galsim.config.GetRNG(config, base, logger, 'COSMOSChromaticGalaxy')
        kwargs['index'], n_rng_calls = cosmos_cat.selectRandomIndex(1, rng=rng, _n_rng_calls=True)

        # Make sure this process gives consistent results regardless of the number of processes
        # being used.
        if not isinstance(cosmos_cat, sc.COSMOSChromaticCatalog) and rng is not None:
            # Then cosmos_cat is really a proxy, which means the rng was pickled, so we need to
            # discard the same number of random calls from the one in the config dict.
            rng.discard(int(n_rng_calls))

    # Even though gal_type is optional, it will have been set in the code above.  So we can at this
    # point assume that kwargs['gal_type'] exists.
    if kwargs['gal_type'] == 'real':
        if rng is None:
            rng = galsim.config.GetRNG(config, base, logger, 'COSMOSChromaticGalaxy')
        kwargs['rng'] = rng

    # NB. Even though index is officially optional, it will always be present, either because it was
    #     set by a call to selectRandomIndex, explicitly by the user, or due to the call to
    #     SetDefaultIndex.
    index = kwargs['index']
    if index >= cosmos_cat.getNObjects():
        raise IndexError(
            "%s index has gone past the number of entries in the catalog"%index)

    logger.debug('obj %d: COSMOSChromaticGalaxy kwargs = %s',base.get('obj_num',0),kwargs)

    # kwargs['cosmos_chromatic_catalog'] = cosmos_cat

    # Use a staticmethod of COSMOSCatalog to avoid pickling the result of makeGalaxy()
    # The RealGalaxy in particular has a large serialization, so it is more efficient to
    # make it in this process, which is what happens here.
    # gal = sc.COSMOSChromaticCatalog._makeSingleGalaxy(cosmos_cat,**kwargs)
    gal = cosmos_cat.makeGalaxy(**kwargs)

    return gal, safe

# Register this as a valid gsobject type
from galsim.config.gsobject import RegisterObjectType
RegisterObjectType('COSMOSChromaticGalaxy', _BuildCOSMOSChromaticGalaxy, input_type='cosmos_chromatic_catalog')
