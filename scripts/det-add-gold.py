#!/usr/bin/env python
"""
Add FLAGS_GOLD to detection catalog.
Spencer Everett
"""

import sys
import fitsio
from astropy.table import Table, join
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'detfile',
    type=str,
    help='Detection catalog filename'
)
parser.add_argument(
    'matchfile',
    type=str,
    help='Matched catalog filename'
)
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
)

def main():

    args = parser.parse_args()
    det_file = args.detfile
    match_file = args.matchfile
    vb = args.vb

    if vb is True:
        print('Loading catalogs...')

    det_cat = Table(det_file)

    # Can add more later
    gold_cols = ['meas_FLAGS_GOLD', 'meas_EXTENDED_CLASS_MOF', 'meas_EXTENDED_CLASS_SOF']
    match_cat = Table(fitsio.read(match_file, columns=gold_cols))

    if vb is True:
        print('Joining catalogs...')

    new_det_cat = join(det_cat, match_cat, join_type='left', keys='bal_id')


    if vb is True:
        print('Writing new detection catalog...')

    outfile = det_file.replace('.fits', '_gold_added.fits')
    new_det_cat.write(outfile)

    if vb is True:
        print('Done!')

    return 0

if __name__ is "__main__":
    sys.exit(main())
