"""
Add EXTENDED object classifier columns to Balrog outputs.
For Y3, see:
https://cdcvs.fnal.gov/redmine/projects/des-y3/wiki/Y3_Extended_Classifier_v2
"""

import numpy as np
import os
import fitsio
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'catfile',
    type=str,
    help='ngmix or gold catalog to add classifier to'
)
parser.add_argument(
    '--ngmix_only',
    action='store_true',
    default=False,
    help='Use to only compute MOF & SOF EXTENDED classifiers.'
)
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
)

def add_ext_mof(cat, fits, vb=False):
    selection_1 = (cat['MOF_CM_T'] + 5. * cat['MOF_CM_T_ERR']) > 0.1
    selection_2 = (cat['MOF_CM_T'] + 1. * cat['MOF_CM_T_ERR']) > 0.05
    selection_3 = (cat['MOF_CM_T'] - 1. * cat['MOF_CM_T_ERR']) > 0.02
    ext_mof = selection_1.astype(int) + selection_2.astype(int) + selection_3.astype(int)

    # Special case
    ext_mof[cat['MOF_CM_T'] < -9000] = -9

    if vb:
        print('Writing EXTENDED_CLASS_MOF...')
    try:
        # If it hasn't been computed yet
        fits[1].insert_column('EXTENDED_CLASS_MOF', ext_mof)
    except:
        fits[1].write_column('EXTENDED_CLASS_MOF', ext_mof)

    return

def add_ext_sof(cat, fits, vb=False):
    selection_1 = (cat['SOF_CM_T'] + 5. * cat['SOF_CM_T_ERR']) > 0.1
    selection_2 = (cat['SOF_CM_T'] + 1. * cat['SOF_CM_T_ERR']) > 0.05
    selection_3 = (cat['SOF_CM_T'] - 1. * cat['SOF_CM_T_ERR']) > 0.02
    ext_sof = selection_1.astype(int) + selection_2.astype(int) + selection_3.astype(int)

    # Special case
    ext_sof[cat['SOF_CM_T'] < -9000] = -9

    if vb:
        print('Writing EXTENDED_CLASS_SOF...')
    try:
        # If it hasn't been computed yet
        fits[1].insert_column('EXTENDED_CLASS_SOF', ext_sof)
    except:
        fits[1].write_column('EXTENDED_CLASS_SOF', ext_sof)

    return

def add_ext_wavg(cat, fits, vb=False):
    selection_1 = (cat['WAVG_SPREAD_MODEL_I'] + 3. * cat['WAVG_SPREADERR_MODEL_I']) > 0.005
    selection_2 = (cat['WAVG_SPREAD_MODEL_I'] + 1. * cat['WAVG_SPREADERR_MODEL_I']) > 0.003
    selection_3 = (cat['WAVG_SPREAD_MODEL_I'] - 1. * cat['WAVG_SPREADERR_MODEL_I']) > 0.001
    ext_wavg = selection_1.astype(int) + selection_2.astype(int) + selection_3.astype(int)

    if vb:
        print('Writing EXTENDED_CLASS_WAVG...')
    try:
        # If it hasn't been computed yet
        fits[1].insert_column('EXTENDED_CLASS_WAVG', ext_wavg)
    except:
        fits[1].write_column('EXTENDED_CLASS_WAVG', ext_wavg)

    return

def add_ext_coadd(cat, fits, vb=False):
    selection_1 = (cat['SPREAD_MODEL_I'] + 3. * cat['SPREADERR_MODEL_I']) > 0.005
    selection_2 = (cat['SPREAD_MODEL_I'] + 1. * cat['SPREADERR_MODEL_I']) > 0.003
    selection_3 = (cat['SPREAD_MODEL_I'] - 1. * cat['SPREADERR_MODEL_I']) > 0.001
    ext_wavg = selection_1.astype(int) + selection_2.astype(int) + selection_3.astype(int)

    if vb:
        print('Writing EXTENDED_CLASS_COADD...')
    try:
        # If it hasn't been computed yet
        fits[1].insert_column('EXTENDED_CLASS_COADD', ext_coadd)
    except:
        fits[1].write_column('EXTENDED_CLASS_COADD', ext_coadd)

    return

def add_ext_mash_mof(cat, fits, vb=False):

    if vb:
        print('Writing EXTENDED_CLASS_MASH_MOF...')
    try:
        # If it hasn't been computed yet
        fits[1].insert_column('EXTENDED_CLASS_MASH_MOF', ext_mash_mof)
    except:
        fits[1].write_column('EXTENDED_CLASS_MASH_MOF', ext_mash_mof)

    return

def add_ext_mash_sof(cat, fits, vb=False):

    if vb:
        print('Writing EXTENDED_CLASS_MASH_SOF...')
    try:
        # If it hasn't been computed yet
        fits[1].insert_column('EXTENDED_CLASS_MASH_SOF', ext_mash_sof)
    except:
        fits[1].write_column('EXTENDED_CLASS_MASH_SOF', ext_mash_sof)

    return

def add_ext_models(cat, fits, ngmix_only, vb=False):

    if vb:
        print('Calculating EXTEND_CLASS_MOF...')
        add_ext_mof(cat, fits, vb=vb)

    if vb:
        print('Calculating EXTEND_CLASS_SOF...')
        add_ext_sof(cat, fits, vb=vb)

    if ngmix_only is not True:

        if vb:
            print('Calculating EXTEND_CLASS_WAVG...')
            add_ext_wavg(cat, fits, vb=vb)

        if vb:
            print('Calculating EXTEND_CLASS_COADD...')
            add_ext_coadd(cat, fits, vb=vb)

        if vb:
            print('Calculating EXTEND_CLASS_MASH_MOF...')
            add_ext_mash_mof(cat, fits, vb=vb)

        if vb:
            print('Calculating EXTEND_CLASS_MASH_SOF...')
            add_ext_mash_sof(cat, fits, vb=vb)

    return

def main():
    args = parser.parse_args()
    vb = args.vb
    catfile = args.catfile
    ngmix_only = args.ngmix_only

    if not os.path.exists(catfile):
        raise IOError('{} does not exist!'.format(catfile))
    catfile = os.path.abspath(catfile)

    if vb:
        print('Grabbing col data...')

    cols = ['mof_cm_T', 'mof_cm_T_err',
            'sof_cm_T', 'sof_cm_T_err']
    if ngmix_only is not True:
        cols += ['wavg_spread_model_i', 'wavg_spreaderr_model_i',
                 'spread_model_i', 'spreaderr_model_i']

    cat = fitsio.read(catfile, columns=cols, ext=1)

    # Used to write new cols
    fits = fitsio.FITS(catfile, 'rw')

    add_ext_models(cat, fits, ngmix_only=ngmix_only, vb=vb)

    return

if __name__=="__main__":
    main()
