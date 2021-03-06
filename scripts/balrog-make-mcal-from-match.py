import numpy as np
from numba import njit
import fitsio
from argparse import ArgumentParser
from numpy.lib.recfunctions import merge_arrays, append_fields
from glob import glob
import os

# import pudb

from mcal_to_h5 import mcal_to_h5

parser = ArgumentParser()

parser.add_argument(
    'basedir',
    type=str,
    help='Directory location of Balrog tile outputs'
    )
parser.add_argument(
    'det_catalog',
    type=str,
    help='Filename of a Balrog detection catalog (created with a given matching criterion)'
    )
parser.add_argument(
    '--version',
    default=None,
    type=str,
    help='Mcal Balrog stack version'
)
parser.add_argument(
    '--outdir',
    default=None,
    type=str,
    help='Output location for stacked catalogs'
)
parser.add_argument(
    '--save_det_only',
    action='store_true',
    default=False,
    help='Only stack objects with a valid bal_id'
)
parser.add_argument(
    '--use_numba',
    action='store_true',
    default=False,
    help='Use numba functions when available'
)
parser.add_argument(
    '--real',
    default=0,
    type=int,
    help='band to download',
)
parser.add_argument(
    '--meds_conf',
    default='y3v02',
    type=str,
    help='MEDS configuration used'
)
parser.add_argument(
    '--max_shape',
    default=450000000,
    type=int,
    help='Maximum size of array'
    )
parser.add_argument(
    '--clobber',
    action='store_true',
    default=False,
    help='Use to overwrite existing match files'
)
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
)

# The following can be modified for given needs; this is simply a way of reducing
# the total number of saved columns / final product size
single_cols = ['id',
               'flags',
               'mask_frac',
               'ra',
               'dec',
               'nimage_tot',
               'nimage_use',
               'psfrec_g',
               'psfrec_T',
               'mcal_gpsf',
               'mcal_Tpsf',
               'gauss_flux',
               'gauss_flux_cov'
              ]

shear_types = ['_1p',
               '_1m',
               '_2p',
               '_2m'
              ]

shear_cols = ['mcal_g',
              'mcal_g_cov',
              'mcal_pars',
              'mcal_pars_cov',
              'mcal_T_r',
              'mcal_T_err',
              'mcal_s2n_r'
             ]

mcal_cols = single_cols + shear_cols + \
            [col+shear for shear in shear_types for col in shear_cols]

def write_stack(stack, outfile=None, clobber=False, save_det_only=False):
    assert stack is not None

    if os.path.exists(outfile):
        if clobber is True:
            os.remove(outfile)
        else:
            raise OSError('{} already exists!'.format(outfile))

        pass

    if save_det_only is True:
        # Can't do this earlier, as we didn't yet have a mapping from
        # bal_id to meas_id
        fitsio.write(outfile, stack[stack['bal_id'] >= 0])
    else:
        fitsio.write(outfile, stack)

    return

@njit
def numba_id_fill(bal_ids, meas_ids, cat_ids, cat_bal_id):
    assert len(bal_ids) == len(meas_ids)
    for i in range(len(bal_ids)):
        bid, mid = bal_ids[i], meas_ids[i]
        indx = np.where(mid == cat_ids)
        assert len(indx) == 1
        cat_bal_id[indx] = bid

    return cat_bal_id

@njit
def numba_id_fill2(det_in_tile, cat, cat_bal_id, mid, bid):
    indx = np.where(mid == cat['id'])
    assert len(indx) == 1
    cat_bal_id[indx] = bid
    return

if __name__ == "__main__":
    args = parser.parse_args()
    vb = args.vb
    conf = args.meds_conf
    real = args.real
    max_shape = args.max_shape
    clobber = args.clobber
    basedir = args.basedir
    outdir = args.outdir
    save_det_only = args.save_det_only
    use_numba = args.use_numba

    if outdir is None:
        outdir = ''
    else:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

    if args.version is None:
        vers = ''
    else:
        vers = '_v' + args.version

    # Grab needed truth info for Balrog matches
    if vb:
        print('Loading detection catalog cols...')

    det_filename = args.det_catalog
    det_cols = ['bal_id', 'meas_id', 'meas_tilename']
    det_cat = fitsio.read(det_filename, columns=det_cols)

    # OLD: Only worked if meas_id was part of bal_id
        # if vb:
        #     print('Recovering meas_id...')

        # # bal_id is defined in the following way:
        # # bal_id = '1' + realization + tilename (w/o 'DES' and +/- mapped to 1/0)
        #            + truth_cat_index
        # # As long as we stick to 10 realizations max (hah!), the following will work

        # # Doing this correctly with numpy is actually fairly complex, but this is still
        # # quick even with ~11 million rows
        # dt = [('meas_id', '>i8')]
        # # meas_id = np.array([i[11:] for i in det_cat['bal_id'].astype(str)], dtype=dt)
        # meas_id = np.array([int(i[11:]) for i in det_cat['bal_id'].astype(str)])
        # det_cat = append_fields(det_cat, 'meas_id', meas_id, usemask=False)

    # Grab all tiles from base directory
    tilepaths = glob(basedir+'/*/')
    tiles = [os.path.basename(os.path.normpath(tilepath)) for tilepath in tilepaths
                if 'DES' in tilepath]

    # We need to create 4 catalog types: [griz, riz] x [nbr, no-nbr]
    mcal_types = ('griz', True), ('griz', False), ('riz', True), ('riz', False)

    Nt = len(tiles)
    tiledir = {}
    mcal_files = {}
    for bands, nbrs in mcal_types:
        if vb:
            print('Starting run for bands: {}; nbrs: {}'.format(bands, nbrs))

        if nbrs is True:
            h5_outfile = 'balrog_mcal_stack-{}-{}-{}-mcal{}.h5'.format(conf, real, bands, vers)
        else:
            h5_outfile = 'balrog_mcal_stack-{}-{}-{}-noNB-mcal{}.h5'.format(conf, real, bands, vers)

        h5_outfile = os.path.join(outdir, h5_outfile)

        real = args.real
        size = 0
        nobjects = 0
        tiledir[(bands, nbrs)] = {}
        mcal_files[(bands, nbrs)] = {}

        for tile in tiles:
            tdir = os.path.join(basedir, tile)
            tiledir[(bands, nbrs)][tile] = tdir
            if nbrs is True:
                mcal_name = 'real_{}_{}-{}-{}-mcal.fits'.format(real, tile, conf, bands)
            else:
                mcal_name = 'real_{}_{}-{}-{}-noNB-mcal.fits'.format(real, tile, conf, bands)

            mcal_file = os.path.join(tdir, mcal_name)
            mcal_files[(bands, nbrs)][tile] = mcal_file

            # Need to know total number of objects for efficient memory stacking
            try:
                h = fitsio.read_header(mcal_file, ext=1)
                size += h['NAXIS2']

            except IOError:
                print('Error: file {} does not exist! Skipping tile.'.format(mcal_file))

        # ----------------------------------------------
        k = 0
        stack = None
        iter_end = 0
        for tile, mcal_file in mcal_files[(bands, nbrs)].items():
            k += 1
            if vb:
                print('Loading tile {} ({} of {})'.format(tile, k, Nt))

            # Grab detected meas_id's for balrog objects in this tile
            det_in_tile = det_cat[det_cat['meas_tilename']==tile]
            det_in_tile = det_in_tile[det_in_tile['meas_id']>0]

            try:
                cat = fitsio.read(mcal_file, columns=mcal_cols)
            except IOError as e:
                print('Following IO error occured:\n{}\nSkipping tile.'.format(e))
                continue

            lencat = len(cat)
            cat_bal_id = -1 * np.ones(lencat, dtype='i8')

            bal_ids = det_in_tile['bal_id'].astype('i8')
            meas_ids = det_in_tile['meas_id'].astype('i8')
            cat_ids = cat['id'].astype('i8')

            if use_numba is True:
                cat_bal_id = numba_id_fill(bal_ids, meas_ids, cat_ids, cat_bal_id)
            else:
                for det_obj in det_in_tile:
                    bid, mid = det_obj['bal_id'], det_obj['meas_id']
                    indx = np.where(mid == cat['id'])
                    assert len(indx) == 1
                    cat_bal_id[indx] = bid

            cat = append_fields(cat, 'bal_id', cat_bal_id, usemask=False)

            if stack is None:
                dt = cat.dtype
                # self.stack = Table(data=np.zeros(size, dtype=dt))
                stack = np.zeros(size, dtype=dt)

            stack[iter_end:iter_end+lencat] = cat[:]

            nobjects += lencat
            iter_end += lencat

        assert nobjects == size

        if vb:
            print('Writing stacked catalog...')
        fits_outfile = h5_outfile.replace('.h5', '.fits')
        write_stack(stack, outfile=fits_outfile, clobber=clobber,
                    save_det_only=save_det_only)

        if vb:
            print('Converting stacked Mcal to hdf5...')
        mcal_to_h5(fits_outfile, h5_outfile, bands, balrog=True, max_shape=max_shape, vb=vb)

#------------------------------------------------------------------------------------
# Old code, possibly useful in future:

    # temp = tilename (w/o 'DES' and +/- mapped to 1/0) + meas_id
    # temp = np.array([i[2:] for i in bid])
    # print t

    # tile_id, meas_id = temp[]
    # print bid
    # print 'type 1: ', type(bid)
    # print 'dtype 1: ', bid.dtype
    # temp  = np.char.lstrip(bid, '1')
    # print 'type 2: ',type(temp)
    # temp2 = np.char.lstrip(temp, '0')
    # print 'type 3: ',type(temp2)
    # print temp2[0:10]

