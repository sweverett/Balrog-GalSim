#!/usr/bin/env python
'''
balrog-make-gold.py
This script will create a 'Gold-like' Balrog catalog
for each tile by calling a variety of scripts to
flatten, merge, compute value-adds, etc.

Current steps per tile:
NB: This is a work in progress, as we can't match DESDM's
method exactly for Y3.
(1) Flatten & merge
(2) Compute FLAGS_GOLD
(3) Compute EXTENDED_CLASS_{}
(4) ?? (Calibrations are complicated for Balrog...)

Usage: python balrog-make-gold.py ....
Author: Spencer Everett (sweveret@ucsc.edu)
'''

import numpy as np
import os
import sys
import subprocess
import shlex
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'basedir',
    type=str,
    help='Location of Balrog tile outputs'
    )
parser.add_argument(
    'script_dir',
    type=str,
    help='Location of Balrog scripts'
    )
parser.add_argument(
    '--merged_subdir',
    type=str,
    default=None,
    help='Set if you want an additional subdirectory for merged outputs'
    )
parser.add_argument(
    '--ngmix_extended_only',
    action='store_true',
    default=False,
    help='Use to only compute MOF & SOF EXTENDED classifiers.'
    )
parser.add_argument(
    '--save_all',
    action='store_true',
    default=False,
    help='Set to save all flattened files'
    )
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
    )

# Names of relevant python scripts
flatten_merge_filename  = 'tile-flatten-merge.py'
gold_flags_filename     = 'tile-add-gold-flags.py'
extended_class_filename = 'tile-add-extended-cols.py'

def run_cmd(cmd, vb):
    if vb is True:
        print(cmd)

    cmd_args = shlex.split(cmd)
    process = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    if vb is True:
        for line in iter(process.stdout.readline, ''): print(line.replace('\n', ''))

    streamdata = process.communicate()[0]
    rc = process.returncode
    if rc !=0:
        sys.exit(rc)

    return

def main():
    args = parser.parse_args()
    basedir = args.basedir
    # script_dir = os.path.abspath(args.script_dir)
    script_dir = args.script_dir
    merged_subdir = args.merged_subdir
    save_all = args.save_all
    ngmix_only = args.ngmix_extended_only
    vb = args.vb

    if os.path.isdir(basedir) is not True:
        raise OSError('{} does not exist!'.format(basedir))

    # Grab all tiles from base directory
    tilepaths = glob(basedir+'/*/')
    tilepaths = [tilepath for tilepath in tilepaths if 'DES' in tilepath]
    tiles = [os.path.basename(os.path.normpath(tilepath)) for tilepath in tilepaths]

    if vb is True:
        print('--------------------------------------------------------------------')

    k = 0
    for tile, tilepath in zip(tiles, tilepaths):
        k += 1
        print('Tile {} ({} out of {})'.format(tile, k, len(tiles)))

        #--------------------------------------------------------------------------------
        # Flatten & merge catalogs
        if vb is True:
            print('Flattening & merging...')

        script_file = os.path.join(script_dir, flatten_merge_filename)
        data_dir = tilepath
        cmd = 'python {} {}'.format(script_file, data_dir)

        if merged_subdir is not None:
            out_dir  = os.path.join(data_dir, merged_subdir)
            cmd += ' --out_dir={}'.format(out_dir)
        else:
            out_dir = data_dir

        if save_all is True:
            cmd += ' --save_all'
        if vb:
            cmd += ' --vb'

        run_cmd(cmd, vb)

        #--------------------------------------------------------------------------------
        # Compute FLAGS_GOLD
        if vb is True:
            print('Computing FLAGS_GOLD...')

        script_file = os.path.join(script_dir, gold_flags_filename)
        merged_cat = os.path.join(out_dir, '{}_merged.fits'.format(tile))

        cmd = 'python {} {} '.format(script_file, merged_cat)
        if vb: cmd += ' --vb'

        run_cmd(cmd, vb)

        #--------------------------------------------------------------------------------
        # Compute EXTENDED CLASSes
        if vb is True:
            print('Computing EXTENDED_CLASS(es)')

        script_file = os.path.join(script_dir, extended_class_filename)
        cmd = 'python {} {}'.format(script_file, merged_cat)
        if ngmix_only is True:
            cmd += ' --ngmix_only'
        if vb is True:
            cmd += ' --vb'

        run_cmd(cmd, vb)

        if vb is True:
            print('--------------------------------------------------------------------')

    if vb is True:
        print('Done!')

    return 0

if __name__ == "__main__":
    sys.exit(main())
