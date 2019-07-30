#!/usr/bin/env python
'''
tile-add-gold-flags.py
This script will append a FLAGS_GOLD column according to Y3 Gold definitions
Usage: python tile-add-gold-flags.py -d [merged_files_dir]
Author: Nacho Sevilla (nsevilla@gmail.com)
Further edited by Spencer Everett
'''
import sys
import numpy as np
import fitsio
import os
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'merged_file',
    type=str,
    help='Merged & flattened gold-like catalog'
    )
parser.add_argument(
    '--mode',
    type=str,
    default='all',
    choices=['all', 'mof', 'sof'],
    help='Can choose to include only one of MOF and SOF in merging & flattening'
    )
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
    )

# Gold cols needed to create FLAGS_GOLD
flag_cols = ['MOF_FLAGS', 'SOF_FLAGS',
             'FLAGS_G', 'FLAGS_R', 'FLAGS_I', 'FLAGS_Z',
             'IMAFLAGS_ISO_G', 'IMAFLAGS_ISO_R', 'IMAFLAGS_ISO_I', 'IMAFLAGS_ISO_Z',
             # This needs to be added once we know where to add these cols to the
             # flatten & merge script (ask Nacho)
             # Doesn't matter for now as we can't compute flag32 anyway
             #'NEPOCHS_G', 'MAG_DETMODEL_I',
             'MAG_AUTO_G', 'MAG_AUTO_R', 'MAG_AUTO_I', 'MAG_AUTO_Z',
             'MAGERR_AUTO_G', 'MAGERR_AUTO_R', 'MAGERR_AUTO_I', 'MAGERR_AUTO_Z',
             'SOF_CM_MAG_G', 'SOF_CM_MAG_R', 'SOF_CM_MAG_I', 'SOF_CM_MAG_Z']

def remove_cols(remove_mode):
    for col in flag_cols:
        if remove_mode in col:
            flag_cols.remove(col)
    return

def main():

    args = parser.parse_args()
    merged_file = args.merged_file
    mode = args.mode

    if not os.path.isfile(merged_file):
        raise OSError('{} not found!'.format(merged_file))

    # While this should usually be FLAGS_GOLD, we have some Balrog-specific processing
    # where we will want to modify the name to be explicit to signify that the usual
    # definition has been changed
    if mode == 'all':
        flag_gold_colname = 'FLAGS_GOLD'
    elif mode == 'mof':
        flag_gold_colname = 'FLAGS_GOLD_MOF_ONLY'
        remove_cols('SOF')
    elif mode == 'sof':
        flag_gold_colname = 'FLAGS_GOLD_SOF_ONLY'
        remove_cols('MOF')

    skipflag = 0b0000000
    skipflag |= 0b100000

    print('Reading {}'.format(merged_file))
    ext = 1
    hdulist = fitsio.FITS(merged_file, mode='rw')
    # tdata = hdulist[1].read()
    tdata = fitsio.read(merged_file, columns=flag_cols, ext=1) # Less memory this way
    # cols = hdulist[ext].get_colnames()

    flag_gold = np.zeros(tdata.size, dtype=np.int32)

    if mode != 'sof':
        if (skipflag & 0b1) == False:
            print 'Processing MOF_FLAGS: flag 1'
            flag_gold[tdata['MOF_FLAGS'] != 0] += 1
    else:
        print('Can\'t process MOF_FLAGS: flag1 without MOF cols - skipping')

    if mode != 'mof':
        if (skipflag & 0b10) == False:
            print 'Processing SOF_FLAGS: flag 2'
            flag_gold[tdata['SOF_FLAGS'] != 0] += 2
        if (skipflag & 0b100) == False:
            print 'Processing SOF GAL_FIT_FAILURE: flag 4'
            flag_gold[(tdata['SOF_FLAGS'] == 1) |
                    (tdata['SOF_FLAGS'] > 2)] += 4
    else:
        print('Can\'t process SOF_FLAGS: flag2, flag4 without SOF cols - skipping')

    if (skipflag & 0b1000) == False:
        print 'Processing SExtractor FLAGS: flag 8'
        flag_gold[(tdata['FLAGS_G'] > 3) |
                  (tdata['FLAGS_R'] > 3) |
                  (tdata['FLAGS_I'] > 3) |
                  (tdata['FLAGS_Z'] > 3)] += 8
    if (skipflag & 0b10000) == False:
        print 'Processing IMAFLAGS: flag 16'
        flag_gold[(tdata['IMAFLAGS_ISO_G'] != 0) |
                  (tdata['IMAFLAGS_ISO_R'] != 0) |
                  (tdata['IMAFLAGS_ISO_I'] != 0) |
                  (tdata['IMAFLAGS_ISO_Z'] != 0)] += 6
    if (skipflag & 0b100000) == False:
        print 'Processing BBJ: flag 32'
        flag_gold[(tdata['NEPOCHS_G']==0) &
                  (tdata['MAGERR_AUTO_G'] < 0.05) &
                  ((tdata['MAG_DETMODEL_I'] - tdata['MAG_AUTO_I']) < -0.4)] += 32

    if mode != 'mof':
        if (skipflag & 0b1000000) == False:
            print 'Processing bright streak rejection: flag 64'
            flag_gold[(((tdata['MAGERR_AUTO_G'] < 0.01) &
                        (tdata['MAG_AUTO_G'] - tdata['SOF_CM_MAG_G'] < -0.5)) |
                    ((tdata['MAGERR_AUTO_R'] < 0.01) & (tdata['MAG_AUTO_R'] - tdata['SOF_CM_MAG_R'] < -0.5)) |
                    ((tdata['MAGERR_AUTO_I'] < 0.01) & (tdata['MAG_AUTO_I'] - tdata['SOF_CM_MAG_I'] < -0.5)) |
                    ((tdata['MAGERR_AUTO_Z'] < 0.01) & (tdata['MAG_AUTO_Z'] - tdata['SOF_CM_MAG_Z'] < -0.5))) &
                    ((tdata['MAG_AUTO_G']-tdata['MAG_AUTO_R'] < -1) | (tdata['MAG_AUTO_G']-tdata['MAG_AUTO_R'] > 4) |
                        (tdata['MAG_AUTO_R']-tdata['MAG_AUTO_I'] < -1) | (tdata['MAG_AUTO_R']-tdata['MAG_AUTO_I'] > 4) |
                        (tdata['MAG_AUTO_I']-tdata['MAG_AUTO_Z']< -1) | (tdata['MAG_AUTO_I']-tdata['MAG_AUTO_Z'] > 4))] += 64
    else:
        print('Can\'t process bright streak rejection: flag64 without SOF cols - skipping')

    print('Writing {}'.format(flag_gold_colname))
    hdulist[ext].insert_column(flag_gold_colname, flag_gold)

    return

if __name__ == '__main__':
    sys.exit(main())

