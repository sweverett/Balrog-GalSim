#!/usr/bin/env python
"""
Match Balrog catalogs and save matched catalogs,
along with detection and extinction information
"""

import healpy as hp
import numpy as np
import os
import fitsio
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'catfile',
    type=str,
    help='Detection catalog filename'
)
parser.add_argument(
    'basedir',
    type=str,
    help='Base directory of the healpix mask files'
)
parser.add_argument(
    '--footprint',
    default='y3a2_footprint_griz_1exp_v2.0.fits.gz',
    type=str,
    help='Filename of the footprint mask'
)
parser.add_argument(
    '--foreground',
    default='y3a2_foreground_mask_v2.1.fits.gz',
    type=str,
    help='Filename of the foreground mask'
)
parser.add_argument(
    '--badregions',
    default='y3a2_badregions_mask_v2.0.fits.gz',
    type=str,
    help='Filename of the badregions mask'
)
parser.add_argument(
    '--verbose',
    action='store_true',
    default=False,
    help='Set to print out more information'
)

def get_masks(mask_files, det_catalog, vb):
    ra  = det_catalog['true_ra']
    dec = det_catalog['true_dec']
    assert len(ra) == len(dec)

    if vb:
        print('{} rows in catalog'.format(len(ra)))

    if vb:
        print('Reading footprint mask...')
    hfoot=hp.read_map(mask_files['footprint'], nest=True)

    if vb:
        print('Reading foreground mask...')
    hfore=hp.read_map(mask_files['foreground'], nest=True)

    if vb:
        print('Reading badregions mask...')
    hbadr=hp.read_map(mask_files['badregions'], nest=True)

    deg2rad = np.pi / 180.0
    p = hp.ang2pix(4096, (90.0-dec)*deg2rad, ra*deg2rad, nest=True)

    foot_mask = hfoot[p]
    fore_mask = hfore[p]
    badr_mask = hbadr[p]

    return foot_mask, fore_mask, badr_mask

def main():
    args = parser.parse_args()
    vb = args.verbose
    basedir = args.basedir
    catfile = args.catfile

    if not os.path.isdir(args.basedir):
        raise IOError('{} does not exist!'.format(basedir))
    basedir = os.path.abspath(basedir)

    if not os.path.exists(catfile):
        raise IOError('{} does not exist!'.format(catfile))
    catfile = os.path.abspath(catfile)

    mask_files = {}
    mask_files['footprint']  = os.path.abspath(os.path.join(basedir, args.footprint))
    mask_files['foreground'] = os.path.abspath(os.path.join(basedir, args.foreground))
    mask_files['badregions'] = os.path.abspath(os.path.join(basedir, args.badregions))

    cat = fitsio.read(catfile, columns=['true_ra', 'true_dec'], ext=1)

    if vb:
        print('Grabbing masks...')
    foot_mask, fore_mask, badr_mask = get_masks(mask_files, cat, vb)

    fits = fitsio.FITS(catfile, 'rw')
    if vb:
        print('Writing footprint mask...')
    fits[1].insert_column('flags_footprint', foot_mask)

    if vb:
        print('Writing foreground mask...')
    fits[1].insert_column('flags_foreground', fore_mask)

    if vb:
        print('Writing badregions mask...')
    fits[1].insert_column('flags_badregions', badr_mask)

if __name__=="__main__":
    main()
