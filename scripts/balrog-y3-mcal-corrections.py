# The Y3 Balrog catalogs have inherited a few bugs or cut choices
# that have since been changed. We apply the following here:
#
# (1) Scale the mcal snr by 1/sqrt(2) [metacalibration bug]
# (2) Re-define size_ratio
# (3) Compute a new weight column based off of mcal properties
#
# Adapted from code by Daniel Gruen

import numpy as np
import h5py
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'bal_mcal_file',
    type=str,
    help='Balrog mcal file to add weights to'
    )
parser.add_argument(
    '--weight_file',
    type=str,
    default='/global/cscratch1/sd/troxel/cats_des_y3/y3_shape_w_grid_03_16_20_highsnr.txt',
    help='txt file that contains the weight grid')
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
    )


# Some grid params from Daniel
snmin=10
snmax=300
sizemin=0.5
sizemax=5
steps=20

def assign_loggrid(x, y, xmin=snmin, xmax=snmax, xsteps=steps, ymin=sizemin, ymax=sizemax, ysteps=steps):
    # return x and y indices of data (x,y) on a log-spaced grid that runs from [xy]min to [xy]max in [xy]steps

    x = np.maximum(x, xmin)
    x = np.minimum(x, xmax)

    y = np.maximum(y, ymin)
    y = np.minimum(y, ymax)

    logstepx = np.log10(xmax/xmin)/xsteps
    logstepy = np.log10(ymax/ymin)/ysteps

    indexx = (np.log10(x/xmin)/logstepx).astype(int)
    indexy = (np.log10(y/ymin)/logstepy).astype(int)

    indexx = np.minimum(indexx, xsteps-1)
    indexy = np.minimum(indexy, ysteps-1)

    return indexx, indexy

def mesh_average(quantity, indexx, indexy, steps, count):
    m = np.zeros((steps,steps))
    np.add.at(m,(indexx,indexy),quantity)
    m /= count
    return m

def apply_loggrid(x, y, grid, xmin=snmin, xmax=snmax, xsteps=steps, ymin=sizemin, ymax=sizemax, ysteps=steps):
    indexx,indexy = assign_loggrid(x, y, xmin, xmax, xsteps, ymin, ymax, ysteps)
    res = np.zeros(len(x))
    res = grid[indexx,indexy]
    return res

def main():
    args = parser.parse_args()
    bal_mcal_file = args.bal_mcal_file
    weight_file = args.weight_file
    vb = args.vb

    # Weight grid
    w = np.genfromtxt(weight_file)

    # Balrog mcal file (will overwrite some columns)
    mcal = h5py.File(bal_mcal_file, 'r+')

    shear_types = ['unsheared',
                   'sheared_1m',
                   'sheared_1p',
                   'sheared_2m',
                   'sheared_2p'
                   ]

    for stype in shear_types:
        if vb is True:
            print('Starting {}'.format(stype))

        base = 'catalog/' + stype + '/'

        # (1) Apply mcal snr correction
        if vb is True:
            print('Applying snr correction...')
        snr_orig = np.array(mcal[base+'snr'])
        snr_corr = snr_orig / np.sqrt(2)
        mcal[base+'snr'][:] = snr_corr

        # (2) Re-define size_ratio
        # in mcal_to_h5.py, it is defined as: mcal_T_r / psfrec_T
        # we are redefining it to be: mcal_T_r / mcal_Tpsf
        # NOTE: the column names mcal_T_r and mcal_Tpsf get renamed to
        # T and mcal_psf_T respectively in that same file
        # so we use that here: T / mcal_psf_T
        if vb is True:
            print('Re-defining size_ratio...')
        T = np.array(mcal['catalog/unsheared/T'])
        mcal_psf_T = np.array(mcal['catalog/unsheared/mcal_psf_T'])
        size_ratio_corr = T / mcal_psf_T
        mcal[base+'size_ratio'][:] = size_ratio_corr

        # (3) Compute relative weight column
        if vb is True:
            print('Computing weight column...')
        weights = apply_loggrid(snr_corr, size_ratio_corr, w)

        max_shape = 450000000
        lencat = len(weights)
        dtype = weights.dtype
        chunks = 1000000

        mcal.create_dataset(base+'weight', max_shape=(max_shape,), shape=(lencat,), dtype=dtype, chunks=(chunks,))
        mcal[base+'weight'][:] = weight

    return

if __name__ == '__main__':
    main()
