import os
import sys
import errno
import itertools
import csv
import warnings

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
	        raise

    # Setup tiles, reals, and bands
    for tile in tiles:
	pass

    return
