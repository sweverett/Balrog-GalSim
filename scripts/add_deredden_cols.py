#!/usr/bin/env python
"""
Calculate and add the dereddened flux & mag cols
for measured objects (true already has them)
"""

from argparse import ArgumentParser
import numpy as np
import fitsio
import os

parser=ArgumentParser()

parser.add_argument(
    'catfile',
    type=str,
    help='Filename of Balrog matched catalog'
)
parser.add_argument(
    'extfile',
    type=str,
    help='Filename of catalog that contains extinction information'
)
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
)

def main():
    args = parser.parse_args()
    vb = args.vb

    if vb:
        print('Reading in catalog info...')
    cat = fitsio.read(args.catfile, columns=['meas_cm_flux',
                                             'meas_cm_mag',
                                             'meas_tilename'])
    if vb:
        print('Reading in extinction info...')
    ext = fitsio.read(args.extfile, columns=['EXTFACT',
                                             'EXTMAG',
                                             'REDTILE'])

    # Create de-reddened columns
    flux_deredden = np.zeros(np.shape(cat['meas_cm_flux']))
    mag_deredden  = np.zeros(np.shape(cat['meas_cm_mag']))
    ext_fact      = np.zeros(np.shape(cat['meas_cm_flux']))
    ext_mag       = np.zeros(np.shape(cat['meas_cm_mag']))

    tiles = np.unique(cat['meas_tilename'])
    Nt = len(tiles)
    if vb:
        print '{} tiles to process'.format(Nt)
        print 'flux len: ', len(flux_deredden[flux_deredden==0])
        print 'mag len: ', len(mag_deredden[mag_deredden==0])

    k = 0
    for tile in tiles:
        k += 1
        if vb:
            print('Calculating de-reddened cols for tile {} ({} of {})'
                  .format(tile, k, Nt))
        indices = np.where(cat['meas_tilename'] == tile)
        cat_tile = cat[indices]
        exfact = ext[ext['REDTILE']==tile]['EXTFACT']
        exmag  = ext[ext['REDTILE']==tile]['EXTMAG']

        flux = cat_tile['meas_cm_flux'] / exfact
        mag  = cat_tile['meas_cm_mag'] - exmag

        flux_deredden[indices] = flux
        mag_deredden[indices] = mag
        ext_fact[indices] = exfact
        ext_mag[indices] = exmag

    if vb:
        print 'flux len: ', len(flux_deredden[flux_deredden==0])
        print 'mag len: ', len(mag_deredden[mag_deredden==0])

    fits = fitsio.FITS(args.catfile, 'rw')
    if vb:
        print('Writing dereddened fluxes...')
    fits[1].insert_column('meas_cm_flux_deredden', flux_deredden)
    if vb:
        print('Writing dereddened mags...')
    fits[1].insert_column('meas_cm_mag_deredden', mag_deredden)
    if vb:
        print('Writing extinction factors...')
    fits[1].insert_column('ext_fact', ext_fact)
    if vb:
        print('Writing extinction mags...')
    fits[1].insert_column('ext_mag', ext_mag)

    return

if __name__=="__main__":
    main()
