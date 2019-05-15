import numpy as np
import fitsio
from argparse import ArgumentParser
from glob import glob
import os

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
               'gauss_flux'
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
              'mcal_s2n_r',
              'psfrec_g',
              'mcal_gpsf',
             ]

mcal_cols = single_cols + shear_cols + \
            [col+shear for shear in shear_types for col in shear_cols]

# for col in mcal_cols:
#     print col

def write_stack(stack, outfile=None, clobber=False):
    assert stack is not None

    if os.path.exists(outfile):
        if clobber is True:
            os.remove(outfile)
        else:
            raise OSError('{} already exists!'.format(outfile))

    fitsio.write(outfile, stack)

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

    if outdir is None:
        outdir = ''
    else:
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

    if args.version is None:
        vers = ''
    else:
        vers = '-' + args.version

    # Grab needed truth info for Balrog matches
    if vb:
        print('Loading detection catalog cols...')

    det_filename = args.det_catalog
    det_cols = ['bal_id', 'meas_tilename']
    det_cat = fitsio.read(det_filename, cols=det_cols)

    if vb:
        print('Recovering meas_id...')

    # bal_id is defined in the following way:
    # bal_id = '1' + realization + tilename (w/o 'DES' and +/- mapped to 1/0) + meas_id
    # As long as we stick to 10 realizations max (hah!), the following will work

    # Doing this correctly with numpy is actually fairly complex, but this is still
    # quick even with ~11 million rows
    meas_id = np.array([i[11:] for i in det_cat['bal_id'].astype(str)])
    print meas_id
    print len(meas_id)

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
        iter_end = 0
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
        for tile, mcal_file in mcal_files[(bands, nbrs)].items():
            k += 1
            if vb:
                if k == 1:
                    extra = ' (first tile will be slow)'
                else:
                    extra = ''
                print('Loading tile {} ({} of {}){}'.format(tile, k, Nt, extra))

            try:
                cat = fitsio.read(mcal_file, cols=mcal_cols)
                lencat = len(cat)
            except IOError as e:
                print('Following IO error occured:\n{}\nSkipping tile.'.format(e))
                continue

            if stack is None:
                dt = cat.dtype
                # self.stack = Table(data=np.zeros(size, dtype=dt))
                stack = np.zeros(size, dtype=dt)

            stack[iter_end:iter_end+lencat] = cat[:]

            nobjects += lencat

        assert nobjects == size

        if vb:
            print('Writing stacked catalog...')
        fits_outfile = h5_outfile.replace('.h5', '.fits')
        write_stack(stack, outfile=fits_outfile, clobber=clobber)

        if vb:
            print('Converting stacked Mcal to hdf5...')
        mcal_to_h5(fits_outfile, h5_outfile, bands, max_shape=max_shape, vb=vb)

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

