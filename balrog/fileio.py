import numpy as np
import fitsio
import os
import sys
import errno
import itertools
import csv
import warnings

# import pudb

#-------------------------------------------------------------------------------
# Input parsing and file loading

# TODO: Some of these functions are not yet implemented. In future, it would be nice
# to have some more robust and standard IO method calls that also incorporate
# astropy.io vs fitsio depending on what the user has.

# def return_output_name(config, ):
    # pass

def open_file(filename):
    pass

def open_fits_file(filename):
    pass

def open_csv_list(filename):
    with open(filename) as file:
        # delimiter = ',' by default
        reader = csv.reader(file)
        lst = list(reader)
        # flatten any iterative lists (in case commas aren't used as delimiter)
        lst = list(itertools.chain.from_iterable(lst))
        return lst

def open_txt_list(filename):
    with open(filename) as file:
        return [line.strip() for line in file]

def setup_output_dir(config, tiles):
    out_dir = config.output_dir
    images = os.path.join(out_dir, 'balrog_images')
    configs = os.path.join(out_dir, 'configs')
    dirs = [out_dir, images, configs]

    # Setup parent-level directories
    for d in dirs:
        try:
            os.makedirs(d)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise e

    # Setup tiles, reals, and bands
    for tile in tiles:
        pass

    return

def combine_fits_extensions(combined_file, bal_file, orig_file, config=None):
    '''
    Combine other needed image extensions (e.g. weight and mask maps).
    Modified from Nikolay's original version in BalrogPipeline.py.
    '''

    # Read in simulated image and pre-injection extensions
    sciIm = fitsio.read(bal_file, ext=0)
    sciHdr = fitsio.read_header(orig_file, ext=0)
    wgtIm, wgtHdr = fitsio.read(orig_file, ext=1, header=True)
    wgt_me_Im, wgt_me_Hdr = fitsio.read(orig_file, ext=2, header=True)
    mskIm, mskHdr = fitsio.read(orig_file, ext=3, header=True)

    # If injecting on blank images, then set these to some sensible values
    if config is not None:
        if config.inj_objs_only['value'] is True:
            # TODO: This needs to be generalized for different noise models!
            if config.inj_objs_only['noise'] == 'BKG+SKY':
                noise = sciHdr['SKYSIGMA']
            else:
                noise = np.mean([sciHdr['RDNOISEA'], sciHdr['RDNOISEB']])
            inv_sky_var = 1.0 / (noise**2)
            wgtIm.fill(inv_sky_var)
            wgt_me_Im.fill(inv_sky_var)
            mskIm.fill(0)

    # Write all 3 extensions in the output file
    # NOTE: Usually combined_file will be the same name as bal_file,
    # so we delete the original injection and replace it with the new extensions
    if os.path.exists(combined_file):
        os.remove(combined_file)

    fits = fitsio.FITS(combined_file,'rw')
    fits.write(sciIm, header=sciHdr)
    fits.write(wgtIm, header=wgtHdr)
    fits.write(wgt_me_Im, header=wgtHdr)
    fits.write(mskIm, header=mskHdr)

    # Put back EXTNAME into headers
    fits[0].write_key('EXTNAME', 'SCI', comment="Extension name")
    fits[1].write_key('EXTNAME', 'WGT', comment="Extension name")
    fits[2].write_key('EXTNAME', 'WGT_ME', comment="Extension name")
    fits[3].write_key('EXTNAME', 'MSK', comment="Extension name")

    return

def return_output_fname(output_dir, dir_name, real, version, tile_name, band, chip_name):
    return os.path.join(output_dir, dir_name, str(real), version, tile_name, band,
                        '{}_balrog_inj.fits'.format(chip_name))
