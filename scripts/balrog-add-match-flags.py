import numpy as np
import fitsio
from astropy.table import Table, join
import esutil.htm as htm
import os
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'det_catalog',
    type=str,
    help='Filename of a Balrog detection catalog (created with a given matching criterion)'
    )
parser.add_argument(
    '--gold_catalog',
    default='/data/des90.a/data/yanny/gapc/Y3GaussAps0m.fits',
    type=str,
    help='Filename of a DES gold catalog to match to as a reference catalog.'
    )
parser.add_argument(
    '--det_bright_col',
    default='true_gap_riz_flux_deredden',
    type=str,
    help='Name of brightness col to compare by in the detection catalog'
    )
parser.add_argument(
    '--gold_bright_col',
    default='AVG_RIZ_GAP_FLUX_DEREDDEN',
    type=str,
    help='Name of brightness col to compare by in the gold catalog'
    )
parser.add_argument(
    '--depth',
    default=14,
    type=int,
    help='depth of HTM matcher'
    )
parser.add_argument(
    '--min_radius',
    default=0.5,
    type=float,
    help='Minimum match radius in arcsec'
    )
parser.add_argument(
    '--max_radius',
    default=2.0,
    type=float,
    help='Maximum match radius in arcsec'
    )
parser.add_argument(
    '--radius_step',
    default=0.25,
    type=float,
    help='Radius step in arcsec'
    )
parser.add_argument(
    '--det_ratag',
    default='true_ra',
    type=str,
    help='Name of detection ra column'
    )
parser.add_argument(
    '--det_dectag',
    default='true_dec',
    type=str,
    help='Name of detection dec column'
    )
parser.add_argument(
    '--gold_ratag',
    default='RA',
    type=str,
    help='Name of gold ra column'
    )
parser.add_argument(
    '--gold_dectag',
    default='DEC',
    type=str,
    help='Name of gold dec column'
    )
# parser.add_argument(
#     '--clobber',
#     default=False,
#     type=bool,
#     help='Set to True if you want to write to same file'
#     )
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
)

# def compute_flux_factors(fluxes, dereddened_fluxes):
#     '''expects np arrays or astropy tables'''
#     return dereddened_fluxes / flux_factors

def main():
    args = parser.parse_args()
    det_file = args.det_catalog
    gold_file = args.gold_catalog
    min_radius = args.min_radius / 3600. # arcsec2deg
    max_radius = args.max_radius / 3600. # arcsec2deg
    drad = args.radius_step / 3600. # arcsec2deg
    det_ratag = args.det_ratag
    det_dectag = args.det_dectag
    gold_ratag = args.gold_ratag
    gold_dectag = args.gold_dectag
    depth = args.depth
    # clobber = args.clobber
    vb = args.vb

    det_bright_col  = args.det_bright_col
    gold_bright_col = args.gold_bright_col

    outfile = det_file.replace('.fits', '_match_flags.fits')

    gold_cols = [gold_ratag, gold_dectag, gold_bright_col]

    if vb is True:
        print('Loading det...')
    det  = Table(fitsio.read(det_file))
    if vb is True:
        print('Loading gold...')
    gold = Table(fitsio.read(gold_file, columns=gold_cols))

    radii = np.arange(min_radius, max_radius+drad, drad)
    if vb is True:
        print('Matching the following radii: {:.2f}'.format(3600.*radii))

    for match_radius in radii:
        if vb is True:
            print('Initializing match_flag for radius {:.2f}...'.format(3600.*match_radius))
        match_flag = np.zeros(len(det), dtype='i4')

        if vb is True:
            print('Matching...')
        h = htm.HTM(depth)
        matcher = htm.Matcher(depth=depth, ra=det[det_ratag], dec=det[det_dectag])
        id_gold, id_det, dist = matcher.match(ra=gold[gold_ratag], dec=gold[gold_dectag],
                                            radius=match_radius)
        assert len(id_det) == len(id_gold)

        if vb is True:
            print('Done Matching')
            print('len(det) = {}'.format(len(det)))
            print('len(id_det) = {}'.format(len(id_det)))
            print('len(match_flag) = {}'.format(len(match_flag)))
            print('len(id_gold) = {}'.format(len(id_gold)))
            print('max(id_det) = {}'.format(np.max(id_det)))
            print('max(id_gold) = {}'.format(np.max(id_gold)))

        match_flag[id_det] = 1
        match_indices = zip(id_det, id_gold)

        Nmatches = len(id_det)

        if vb is True:
            print('There are {} injections near a gold object'.format(Nmatches))
            print('Checking brightness ratios...')
        k = 0
        for idet, igold in match_indices:
            k += 1
            if k%100 == 0:
                if vb is True: print('{} / {}'.format(k, Nmatches))
            if det[det_bright_col][idet] < gold[gold_bright_col][igold]:
               match_flag[idet] = 2

        det['match_flag_{}_asec'] = match_flag

    if vb is True:
        print('Writing to {}...'.format(outfile))
    det.write(outfile)

    return

if __name__ == '__main__':
    main()
