import numpy as np
import galsim.config.input as gsinput
import os

# Balrog files
import mathutil as util
import grid

# import pudb

class BalObject(object):
    '''
    '''

    def __init__(self):
        pass

class BalInjectionCatalog(object):

    def __init__(self, input_type, inj_type, tile, needs_band=False):
        self.input_type = input_type
        self.inj_type = inj_type
        self.needs_band = needs_band

        self.pos = None
        self.indx = None
        self.nobjects = None

        self.truth_outfile = {}

        # Balrog catalogs are constructed for a given Tile and require
        # info for construction
        self.ramin, self.ramax = tile.ramin, tile.ramax
        self.decmin, self.decmax = tile.decmin, tile.decmax
        self.ra_boundary_cross = tile.ra_boundary_cross
        self.Npix_x, self.Npix_y = tile.Npix_x, tile.Npix_y
        self.pixel_scale = tile.pixel_scale
        self.wcs = tile.wcs
        self.tile_name = tile.tile_name

        # This is also set for a given Tile, as it may depend on its area
        self.objs_per_real = tile.objs_per_real[input_type]

        # Some tests will impose single-object injection - this is turned off by default
        # NOTE: Only supported by a few input types for now, identified below
        self.single_obj_injection = False

        # NOTE: In the future, it may be useful to save the actual input catalog objects
        # for each realization. For now it is just wasted memory
        # self.cat = None

        return

    def generate_objects(self, config, realization):
        # Generate positions and indices if this is the first realization
        if realization == config.realizations[0]:
            # Can't guarantee count consistency per real, so use dicts
            self.pos = {}
            self.indx = {}
            self.nobjects= {}

            input_type = self.input_type
            input_nobjects = config.input_nobjects[input_type]
            Nr = config.n_realizations

            for real in config.realizations:
                inj_nobjs = self.objs_per_real
                self.nobjects[real] = inj_nobjs

                # Generate object coordinates
                ps = config.pos_sampling
                if ps['type'] == 'uniform':
                    ra = util.sample_uniform_ra(self.ramin, self.ramax, config.objs_per_real,
                                        boundary_cross=self.ra_boundary_cross)
                    dec = util.sample_uniform_dec(self.decmin, self.decmax, config.objs_per_real,
                                        unit='deg')
                    self.pos[real] = np.column_stack((ra, dec))

                elif (ps['type']=='RectGrid') or (ps['type']=='HexGrid'):

                    gs = ps['grid_spacing']

                    # Rotate grid if asked
                    try:
                        r = ps['rotate']
                        if (isinstance(r, str)) and (r.lower() == 'random'):
                            if ps['type'] == 'RectGrid':
                                self.grid_rot_angle = np.rand.uniform(0., np.pi/2.)
                            elif ps['type'] == 'HexGrid':
                                self.grid_rot_angle = np.rand.uniform(0., np.pi/3.)
                        else:
                            unit = ps['angle_unit']
                            if unit == 'deg':
                                if (r>=0.0) and (r<360.0):
                                    self.grid_rot_angle = float(r)
                                else:
                                    raise ValueError('Grid rotation of {} '.format(r) +
                                                    'deg is not valid!')
                            else:
                                if (r>=0.0) and (r<2*np.pi):
                                    self.grid_rot_angle = float(r)
                                else:
                                    raise ValueError('Grid rotation of {} '.format(r) +
                                                    'rad is not valid!')
                    except KeyError:
                        self.grid_rot_angle = 0.0

                    # Offset grid if asked
                    try:
                        o = ps['offset']
                        if (isinstance(o, str)) and (o.lower() == 'random'):
                            self.grid_offset = [np.rand.uniform(-gs/2., gs/2.),
                                                np.rand.uniform(-gs/2., gs/2.)]
                        else:
                            if isinstance(o, list):
                                self.grid_offset = list(o)
                            else:
                                raise ValueError('Grid offset of {} '.format(r) +
                                                'is not an array!')
                    except KeyError:
                        self.grid_offset = [0.0, 0.0]

                    try:
                        self.angle_unit = ps['angle_unit']
                    except KeyError:
                        self.angle_unit = None

                    # Creates the rectangular grid given tile parameters and calculates the
                    # image / world positions for each object
                    if ps['type'] == 'RectGrid':
                        tile_grid = grid.RectGrid(gs, self.wcs, Npix_x=self.Npix_x,
                                                Npix_y=self.Npix_y,
                                                pixscale=self.pixel_scale,
                                                rot_angle = self.grid_rot_angle,
                                                angle_unit = self.angle_unit,
                                                pos_offset = self.grid_offset)
                    elif ps['type'] == 'HexGrid':
                        tile_grid = grid.HexGrid(gs, self.wcs, Npix_x=self.Npix_x,
                                                Npix_y=self.Npix_y,
                                                pixscale=self.pixel_scale,
                                                rot_angle = self.grid_rot_angle,
                                                angle_unit = self.angle_unit,
                                                pos_offset = self.grid_offset)

                    self.pos[real] = tile_grid.pos

                    # NOTE: We ignore the inputted nobjects and use the correct grid value
                    # instead (user was already warned)
                    inj_nobjs = np.shape(tile_grid.pos)[0]
                    self.nobjects[real] = inj_nobjs

                # Generate object indices (in input catalog)
                indices = np.random.choice(xrange(input_nobjects), size=inj_nobjs)
                self.indx[real] = indices

                # NOTE: This is where we could initialize the injection cat if needed
                # self.cat[real] = config.input_cats[inj_type][indices]

                return

    # TODO: Should really make this a static method instead
    def _check_for_single_obj_indx(self, config):
        # NOTE: Nearly all runs will generate a random sample of indices. However,
        # for some testing it would be nice to use an identical object for all
        # injections. In this case, the user can set a single index in the 'gal'
        # section of the global config
        # NOTE: Only currently suppported for sims with a single input type
        try:
            orig_indx = config.gs_config[0]['gal']['index']
            if type(orig_indx) is int:
                # Need to find original index of catalog
                gs_config = copy.deepcopy(config.gs_config[0])
                # Add dummy band index (band doesn't matter)
                gs_config['input'][input_type].update({'bands':'g'})
                galsim.config.ProcessInput(gs_config)
                cat_proxy = gs_config['input_objs'][input_type][0]
                cat = cat_proxy.getCatalog()
                # Specific catalog structures can then generate indices from
                # the proxy catalog
                return cat

            else:
                raise TypeError('Can only set a global object index in the ' +
                                'config if it is an integer!')
        except KeyError:
            return None

    def get_truth_outfile(self, base_outfile, real):
        truth_fname = '{}_{}_balrog_truth_cat_{}.fits'.format(self.tile_name, real, self.inj_type)
        self.truth_outfile[real] = os.path.join(base_outfile, truth_fname)
        return self.truth_outfile[real]

    def write_new_positions(self, truth_cat, realization):
        pos = self.pos[realization]

        # If nothing is set for a given custom input, try the obvious
        try:
            truth_cat[self.inj_type]['ra'] = pos[:,0]
            truth_cat[self.inj_type]['dec'] = pos[:,1]
        except KeyError:
            try:
                truth_cat[in_type]['RA'] = pos[:,0]
                truth_cat[in_type]['DEC'] = pos[:,1]
            except KeyError:
                raise('Tried to write truth positions using column names of ra/dec; RA/DEC. '
                      'Eeither rename position columns or overload `write_new_positions()` '
                      'for {}'.format(self.input_type))

        return

    def setup_chip_config(self, config, bal_config, chip, chip_indx):
        # Many injection types will requite nothing special in setup
        pass

    def build_single_chip_config(self, config, bal_config, chip, chip_indx):
        pass

    def build_multi_chip_config(self, config, bal_config, chip, chip_indx, input_indx):
        pass

class DESInjectionCatalog(BalInjectionCatalog):
    def __init__(self, input_type, inj_type, tile, needs_band):
        # All catalogs require band input
        assert needs_band is True
        super(DESInjectionCatalog, self).__init__(input_type, inj_type, tile, needs_band)

        return

    def setup_chip_config(self, config, bal_config, chip, chip_indx):
        # Only load into memory the needed band catalog information
        bal_config[chip_indx]['input'].update({
            self.input_type : {'bands' : chip.band}
        })

        return

class NGMIXInjectionCatalog(DESInjectionCatalog):
    def generate_objects(self, config, realization):
        super(NGMIXInjectionCatalog, self).generate_objects(config, realization)

        # NOTE: See `_check_for_single_obj_indx()` for why we sometimes do this for testing
        single_obj_cat = self._check_for_single_obj_indx(config)
        if single_obj_cat is not None:
            # Specific to ngmix_catalog structure:
            indx = int(np.where(cat['id']==orig_indx)[0])
            self.indx[realization] = indx * np.ones(self.nobjects[realization], dtype='int16')
            self.single_obj_injection = True

        return

class MEDSInjectionCatalog(DESInjectionCatalog):
    def generate_objects(self, config, realization):
        super(MEDSInjectionCatalog, self).generate_objects(config, realization)

        # NOTE: See `_check_for_single_obj_indx()` for why we sometimes do this for testing
        single_obj_cat = self._check_for_single_obj_indx(config)
        if single_obj_cat is not None:
            # Specific to meds_catalog structure:
            b = cat_proxy.getBands()[0] # ID's consistent between bands
            indx = int(np.where(cat[b]['id']==orig_indx)[0])
            self.indx[realization] = indx * np.ones(self.nobjects[realization], dtype='int16')
            self.single_obj_injection = True

        return

    def build_single_chip_config(self, config, bal_config, chip, chip_indx):
        # Only use meds/psf files for needed band
        b = config.bindx[chip.band]
        meds_file = [bal_config[0]['input'][self.input_type]['meds_files'][b]]
        psf_file = [bal_config[0]['input'][self.input_type]['psf_files'][b]]
        bal_config[chip_indx]['input'][self.input_type].update({
            'meds_files' : meds_file,
            'psf_files' : psf_file
        })

        bal_config[chip_indx]['gal'].update({
            'band' : chip.band
        })

        return

    def build_multi_chip_config(self, config, bal_config, chip, chip_indx, input_indx):
        # Only use meds/psf files for needed band
        b = config.bindx[chip.band]
        meds_file = [bal_config[0]['input']['items'][input_indx]['meds_files'][b]]
        psf_file = [bal_config[0]['input']['items'][input_indx]['psf_files'][b]]
        bal_config[chip_indx]['input']['items'][input_indx].update({
            'meds_files' : meds_file,
            'psf_files' : psf_file
        })

        bal_config[chip_indx]['gal']['items'][input_indx].update({
            'band' : chip.band
        })

        return

class DESStarInjectionCatalog(DESInjectionCatalog):
    def __init__(self, input_type, inj_type, tile, needs_band=False):
        super(DESStarInjectionCatalog, self).__init__(input_type, inj_type, tile, needs_band)

        # Might add something like this in the future, if we end up passing configs
        # during init...
        # if config.pos_sampling['type'] == 'sahar':
        #     self.sahar_pos = True
        # else:
        #     self.sahar_pos = False

        return

    def generate_objects(self, config, realization):
        if config.pos_sampling['type'] == 'sahar':
            # Sahar has pre-computed positions for her catalogs
            self._generate_sahar_coords()
            self.sahar_pos = True
        else:
            super(DESStarInjectionCatalog, self).generate_objects(config, realization)
            self.sahar_pos = False

        return

    def setup_chip_config(self, config, bal_config, chip, chip_indx):
        super(DESStarInjectionCatalog, self).setup_chip_config(config, bal_config, chip, chip_indx)

        # Only load into memory the needed band catalog information
        bal_config[0]['input'][self.input_type].update({'tile' : self.tile_name})

        return

    def _generate_sahar_coords(self, config, realization):
        '''
        For now (Y3), the star catalogs (including positions) are pre-computed. So we just
        need to declare some variables for future use.
        '''
        inp_type = self.input_type

        if self.data_version == 'y3v02':
            # The first time, divide up catalog randomly between realizations
            if realization == config.realizations[0]:
                # Can't guarantee star count consistency, so use dicts
                self.indx = {}
                self.pos = {}
                self.nobjects = {}

                # Randomize star catalog order and split into approximately equal parts
                # NOTE: If n_realizations > len(realizations), then DO NOT randomly
                # shuffle stars, as they are being injected across multiple jobs.
                Nr = config.n_realizations
                if Nr == len(config.realizations):
                    rand.shuffle(indices)
                indices = [np.array(indices[i::Nr]) for i in range(Nr)]

                # Grab star positions
                ra = config.input_cats[input_type]['RA_new']
                dec = config.input_cats[input_type]['DEC_new']
                assert len(ra)==len(dec)

                # Sahar's DES Y3 star catalogs are all pre-computed, so we can set
                # needed values for all realizations now.
                for real in config.realizations:
                    j = int(np.where(real==np.array(config.realizations))[0])
                    inds = indices[j]
                    r, d = ra[inds], dec[inds]

                    self.indx[real] = inds
                    self.pos[real] = np.column_stack((r, d))
                    self.nobjects[real] = len(inds)

        return

    def write_new_positions(self, truth_cat, realization):
        # Currently, all used DES star catalogs have Sahar's naming scheme anyway,
        # so this is check is not needed
        # if self.sahar_pos is True:
        truth_cat[self.inj_type]['RA_new'] = self.pos[realization][:,0]
        truth_cat[self.inj_type]['DEC_new'] = self.pos[realization][:,1]
        # else:
        #     super(DESStarInjectionCatalog, self).write_new_positions(truth_cat, realization)

        return

class COSMOSInjectionCatalog(BalInjectionCatalog):
    def setup_chip_config(self, config, bal_config, chip, chip_indx):
        # Set the bandpass
        bal_config[chip_indx]['stamp'].update({
            'type' : 'COSMOSChromatic',
            'bandpass' : config.filters[chip.band].band_config
        })

        return

# TODO: This should be filled with relevant construction info from Alex DW's udg_catalog class
class UDGInjectionCatalog(BalInjectionCatalog):
    pass

class Galaxy(object):
    '''
    # TODO: Do we need a galaxy class? (probably not)
    '''

    def __init__(self):
        pass

class Star(object):
    '''
    # TODO: Do we need a star class? (probably not)
    '''

    def __init__(self):
        pass

def build_bal_inject_cat(input_type, inj_type, tile, needs_band):
    if input_type in BALROG_INJECTION_TYPES:
        # User-defined injection catalog construction
        inject_cat = BALROG_INJECTION_TYPES[input_type](input_type, inj_type, tile, needs_band)
    else:
        # Generic injection catalog construction
        if input_type not in gsinput.valid_input_types:
            raise ValueError('{} is not a native GalSim input type '.format(input_type) +
                'or a recognized Balrog input type. Make sure you have written '
                'and registered a valid GalSim input type')
        inject_cat = BalInjectionCatalog(input_type, inj_type, tile, needs_band)

    return inject_cat

BALROG_INJECTION_TYPES = {
    'ngmix_catalog' : NGMIXInjectionCatalog,
    'meds_catalog' : MEDSInjectionCatalog,
    'udg_catalog' : UDGInjectionCatalog,
    'des_star_catalog' : DESStarInjectionCatalog,
    'cosmos_chromatic_catalog' : COSMOSInjectionCatalog
    }
