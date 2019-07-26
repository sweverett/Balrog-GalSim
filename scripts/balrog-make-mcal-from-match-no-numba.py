import numpy as np
import fitsio
from astropy.table import Table, Column, join, vstack
from argparse import ArgumentParser
from numpy.lib.recfunctions import merge_arrays, append_fields
from glob import glob
import os

from balrog_mcal_to_h5 import convert_mcal_to_h5

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
    '--gold_base',
    default=None,
    type=str,
    help='Base directory of gold products '
    'Use to add a few extra value-adds from the gold-like cats'
)
parser.add_argument(
    '--gold_subdir',
    default=None,
    type=str,
    help='Subdirectory of desired gold merged cat, if not in TILENAME'
)
# parser.add_argument(
#     '--keep_cache',
#     action='store_true',
#     default=False,
#     help='Keep cached individual value-added mcal fits'
# )
parser.add_argument(
    '--write_fits',
    action='store_true',
    default=False,
    help='Set to write a fits version of the full stack'
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
mcal_flats = {'mcal_T_r':'T',
              'mcal_T_err':'T_err',
              'mcal_s2n_r':'snr'
             }

mcal_vec2 = {'psfrec_g':'psf_e',
             'mcal_gpsf':'mcal_psf_e'
            }
mcal_vec2_ext = ['1','2']

mcal_vec3 = {'nimage_tot':'nimage_tot_',
             'nimage_use':'nimage_use_'
            }
mcal_vec3_ext = ['r','i','z']

mcal_vec4 = {'nimage_tot':'nimage_tot_',
             'nimage_use':'nimage_use_'
            }
single_cols = ['id',
               'flags',
               'mask_frac',
               'ra',
               'dec',
               'psfrec_T',
               'mcal_Tpsf',
               'nimage_tot',
               'nimage_use',
               'psfrec_g',
               'mcal_gpsf',
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

# In principle we'd match gold catalogs directly. Until then, we may
# want to include certain value-added cols from Gold-like cats in the
# matched ngmix catalogs
_gold_basename = 'TILENAME_merged.fits'
_gold_cols = ['COADD_OBJECT_ID',
              'TILENAME',
              'FLAGS_GOLD',
              'EXTENDED_CLASS_MOF',
              'EXTENDED_CLASS_SOF']

# TODO: Move somewhere more central in the future
_blacklisted_tiles = [
    'DES0000-0333',
    'DES0000-3706',
    'DES0002-0207',
    'DES0238-3457'
]

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
        # fitsio.write(outfile, stack[stack['bal_id'] >= 0])
        stack[stack['bal_id'] >= 0].write(outfile, overwrite=clobber)
    else:
        # fitsio.write(outfile, stack)
        stack.write(outfile, overwrite=clobber)

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
    # keep_cache = args.keep_cache
    write_fits = args.write_fits
    gold_base = args.gold_base
    gold_subdir = args.gold_subdir

    if outdir is None:
        outdir = ''
    else:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

    if args.version is None:
        vers = ''
    else:
        vers = '_v' + args.version

    if gold_base is not None:
        if gold_base == 'base':
            gold_base = basedir
        else:
            if not os.path.exists(gold_base):
                raise OSError('{} does not exist!'.format(gold_base))

    if gold_subdir is not None:
        if gold_base is None:
            raise ValueError('Can\'t set gold_subdir if gold_base isn\'t set!')

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

    cache_dir = os.path.join(outdir, 'cache')
    if not os.path.isdir(cache_dir):
        os.mkdir(cache_dir)

    for mcal_type in mcal_types:
        b, n = mcal_type[0], mcal_type[1]
        if n is True:
            nb = 'NB'
        else:
            nb = 'noNB'
        type_dir = os.path.join(cache_dir, b+'_'+nb)
        if not os.path.isdir(type_dir):
            os.mkdir(type_dir)

    Nt = len(tiles)
    blacklisted = []
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
            if tile in _blacklisted_tiles:
                if tile not in blacklisted:
                    blacklisted.append(tile)
                continue

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

        del h

        # ----------------------------------------------
        k = 0
        stack = None
        iter_end = 0
        # tile_cats = []
        for tile, mcal_file in mcal_files[(bands, nbrs)].items():
            k += 1
            if vb:
                print('Loading tile {} ({} of {})'.format(tile, k, Nt-len(blacklisted)))

            if tile in _blacklisted_tiles:
                print('Tile {} is in blacklist, skipping tile'.format(tile))
                continue

            if nbrs is True:
                nb = 'NB'
            else:
                nb = 'noNB'
            type_dir = os.path.join(cache_dir, bands+'_'+nb)

            # Grab detected meas_id's for balrog objects in this tile
            det_in_tile = det_cat[det_cat['meas_tilename']==tile]
            det_in_tile = det_in_tile[det_in_tile['meas_id']>0]

            try:
                cat = Table(fitsio.read(mcal_file, columns=mcal_cols))
            except IOError as e:
                print('Following IO error occured:\n{}\nSkipping tile.'.format(e))
                continue

            lencat = len(cat)
            cat_bal_id = -1 * np.ones(lencat, dtype='i8')

            bal_ids = det_in_tile['bal_id'].astype('i8')
            meas_ids = det_in_tile['meas_id'].astype('i8')
            cat_ids = cat['id'].astype('i8')

            for det_obj in det_in_tile:
                bid, mid = det_obj['bal_id'], det_obj['meas_id']
                indx = np.where(mid == cat['id'])
                assert len(indx) == 1
                cat_bal_id[indx] = bid

            # cat = append_fields(cat, 'bal_id', cat_bal_id, usemask=False)
            cat.add_column(Column(cat_bal_id, name='bal_id'))

            # Add certain gold value-adds if desired
            if gold_base is not None:
                gold_tbase = os.path.join(gold_base, tile)
                gold_filename = _gold_basename.replace('TILENAME', tile)
                if gold_subdir is None:
                    gold_subdir = ''
                gold_catfile = os.path.join(gold_tbase,
                                            gold_subdir,
                                            gold_filename)

                try:
                    gold_cat = Table(fitsio.read(gold_catfile,
                                                columns=_gold_cols))
                except IOError as e:
                    print('Following IO error occured:\n{}\nSkipping tile.'.format(e))
                    continue

                assert len(cat) == len(gold_cat)
                # Need the ID colnames to match
                gold_cat.rename_column('COADD_OBJECT_ID', 'id')
                cat = join(cat, gold_cat, keys='id', join_type='left')

            # tile_cats.append(cat)

            cache_filename = 'balrog_' + os.path.basename(mcal_file)
            cache_file = os.path.join(type_dir, cache_filename)

            if vb is True:
                print('Writing cache file...')
            cat.write(cache_file, overwrite=clobber)

            if write_fits is True:
                if stack is None:
                    dt = cat.dtype
                    stack = Table(np.full(size, -1, dtype=dt))

                stack[iter_end:iter_end+lencat] = cat[:]
                iter_end += lencat

            nobjects += lencat

        # if vb:
        #     print('Stacking catalogs...')
        # stack = vstack(tile_cats, join_type='exact')

        assert nobjects == size

        if write_fits is True:
            if vb:
                print('Writing stacked catalog...')
            fits_outfile = h5_outfile.replace('.h5', '.fits')
            write_stack(stack, outfile=fits_outfile, clobber=clobber,
                        save_det_only=save_det_only)

        if vb:
            print('Converting Mcal\'s to hdf5 stack...')
        convert_mcal_to_h5(type_dir, h5_outfile, bands)

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

