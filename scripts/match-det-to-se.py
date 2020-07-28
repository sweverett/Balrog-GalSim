import fitsio
from astropy.table import Table
import numpy as np
from glob import glob
import os
import esutil.htm as htm
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'det_file',
    type=str,
    help='Balrog detection file to add SE cols to'
)
parser.add_argument(
    'gold_basedir',
    type=str,
    help='Base directory of Gold tile outputs'
)
parser.add_argument(
    'gold_subdir',
    type=str,
    help='Name of Gold subdir desired'
)
parser.add_argument(
    '--clobber',
    type=bool,
    action='set_true',
    help='Set to overwrite the output file (not original det file)'
)
parser.add_argument(
    '--test',
    type=bool,
    action='set_true',
    help='Set to only run through a few tiles'
)
parser.add_argument(
    '--vb',
    type=bool,
    action='set_true',
    help='Set for verbose priting'
)

def match(det_cat, gold_cat, depth=14, det_ratag='true_ra', det_dectag='true_dec', gold_ratag='RA',
          gold_dectag='DEC', match_radius=0.01 / 3600.)
    h = htm.HTM(depth)
    matcher = htm.Matcher(depth=depth,
                           ra=det_cat[det_ratag],
                           dec=det_cat[det_dectag])
    id_gold, id_det, dist = matcher.match(ra=gold_cat[gold_ratag],
                                     dec=gold_cat[gold_dectag],
                                     radius=match_radius)

    return id_det, id_gold

def main():
    args = parser.parse_args()
    det_file = args.det_file
    basedir = args.gold_basedir
    subdir = args.gold_subdir
    test = args.test
    vb = args.vb

    det = Table(fitsio.read(det_file))
    Ndet = len(det)

    gold_colnames = ['RA', 'DEC']

    default_val = -10000.
    for b in 'griz':
        col = 'MAG_AUTO_{}'.format(b.capitalize())
        gold_colnames.append(col)
        det[col] = default_val*np.ones(len(det))

        col = 'FLUX_RADIUS_{}'.format(b.capitalize())
        gold_colnames.append(col)
        det[col] = default_val*np.ones(len(det))

    tiles = glob(os.path.join(basedir, '/*'))
    tiles = [t in tiles if 'DES' in t]
    Nt = len(tiles)

    if test is True:
        dt = Nt // 5
        tiles = tiles[::dt]

    k = 0
    for tile in tiles:
        k += 1
        if vb is True:
            print('Matching tile {} ({} of {})'.format(tile, k, Nt))
        gold_fname = tile + '_merged.fits'
        gold_file = os.path.join(basedir, tile, gold_subdir, gold_fname)
        gold = Table(fitsio.read(gold_file, columns=gold_cols))

        id_det, id_gold = match(det, gold)
        for b in 'griz':
            for col in ['MAG_AUGO_{}'.format(b.capitalize()), 'FLUX_RADIUS_{}'.format(b.capitalize())]:
                det[col][id_det] = gold[col][id_gold]

    assert len(det) == Ndet

    # In principle should be an assertion, but want to test first
    for b in 'griz':
        for col in ['MAG_AUGO_{}'.format(b.capitalize()), 'FLUX_RADIUS_{}'.format(b.capitalize())]:
            N = len(det[col][det[col] == default_val])
            if N > 0:
                print('WARNING: There are {} rows still set as default for col {}!'.format(
                N, col))

    if vb is True:
        print('Writing new detection catalog...')
    outfile = det_file.replace('.fits', '_se_match.fits')
    det.write(outfile, overwrite=clobber)

    return

if __name__ == '__main__':
    main()
