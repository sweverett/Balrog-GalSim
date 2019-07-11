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
    '--verbose',
    action='store_true',
    default=False,
    help='Set to print out more information'
)

def main():

    args = parser.parse_args()
    det_file = args.detfile
    match_file = args.matchfile

    det_cat = Table(det_file)

    # Can add more later
    gold_cols = ['FLAGS_GOLD', 'EXTENDED_CLASS_MOF', 'EXTENDED_CLASS_SOF']
    match_cat = Table(fitsio.read(match_file, columns=gold_cols))

    new_det_cat = join(det_cat, match_cat, join_type='left', keys='bal_id')

    outfile = det_file.replace('.fits', '_gadded.fits')
    new_det_cat.write(outfile)

    return 0

if __name__ is "__main__":
    return sys.exit(main())
