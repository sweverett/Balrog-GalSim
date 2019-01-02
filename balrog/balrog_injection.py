#####################################################################
#
# I should probably add a description for balrog!
#
#
# Spencer Everett
# UCSC
# 12/10/2017
#####################################################################

import numpy as np
import os
import sys
import argparse

# Balrog files
import config as Config
import tile as Tile
import fileio as io

# Use for debugging
# import pudb

#-------------------------------------------------------------------------------
# Important todo's:
# TODO: Clean up evals in add_gs_injection()!
# TODO: Check pixel origin for transformations!

# Some extra todo's:
# TODO: Have code automatically check fitsio vs astropy
# TODO: More general geometry file inputs
# TODO: Add check for python path!
# TODO: Make a log system!
# TODO: Add a single seed for each tile w/r/t noise realizations
# TODO: Multi-line print statements are a bit bugged

# Questions:
# None curently


def parse_args():
    '''
    Parse command line arguments. Many of the following can be set in the `bal_config` file as well;
    these are kept here only for convenience with previous versions where the arguments were
    required. Passing a value through both methods is fine, as long as they are not inconsistent with
    one another.
    '''

    parser = argparse.ArgumentParser()
    # Global GalSim config file with options that will be applied to all injection images.
    parser.add_argument('config_file', help='.yaml or .json confg file that specifies the desired GalSim'
                        'simulation and injection.')
    # Optional argument for config file(s) directory location (if not .)
    parser.add_argument('-c', '--config_dir', help='Directory that houses the GalSim and related config files.')
    # Required argument for tile geometry file
    parser.add_argument('-g', '--geom_file', help='Tile geometry file (e.g. `Y3A2_COADDTILE_GEOM`)')
    # Optional argument for tile list
    parser.add_argument('-l', '--tile_list', help='.txt or .csv file containing list of tiles to be processed.')
    # Optional argument for DES tiles directory location
    parser.add_argument('-t', '--tile_dir', help='Directory of DES tiles containing single exposure images.')
    # Optional argument for a tile's psf directory name (if not contained in tile_dir/{TILE}/psfs)
    parser.add_argument('-p', '--psf_dir', help='Directory that houses individual psf files for DES chips.')
    # Optional argument for output directory (if not .)
    parser.add_argument('-o', '--output_dir', help='Directory that houses output Balrog images.')
    # Optional argument for GalSim # of processors
    parser.add_argument('-n', '--nproc', action='store', type=int,
                        help='Number of processors that GalSim will use when building files.')
    # Optional argument for verbose messages
    parser.add_argument('-v', '--verbose', action='store', nargs='?', default='0', const='1', type=int,
                        help='Turn on verbose mode for additional messages.')

    return parser.parse_args()

#-------------------------------------------------------------------------------

# Run top-level Balrog script
# TODO: In future, make this `RunDESBalrog`, or something similar to allow for other stacks
def RunBalrog():
    '''
    Main Balrog call.
    '''

    # Parse command line arguments
    args = parse_args()

    # Check verbosity
    vb = args.verbose
    if vb: print('Input arguments: {}'.format(args))

    # Set up Balrog and GalSim configuration
    if vb: print('Setting up configuration...')
    config = Config.setup_config(args)

    # Determine which tiles will have injections
    if vb: print('Creating tiles...')
    tiles = Tile.create_tiles(config)

    # Set up output directory structure
    if vb: print('Setting up output directory...')
    io.setup_output_dir(config, tiles)

    # Now loop over all tiles slated for injection:
    # TODO: This could be parallelized in future
    for i, tile in enumerate(tiles):
        # pudb.set_trace()
        config.reset_gs_config()
        config.set_tile_num(i)

        for real in config.realizations:
            if vb: print('Injecting Tile {}; realization {}'.format(tile.tile_name, real))
            # Reset gs config for each new realization
            tile.set_realization(real)
            tile.reset_bal_config(config)
            tile.generate_objects(config, real)
            # Allow for different band injections in different tiles
            bands = tile.bands
            for band in bands:
                # pudb.set_trace()
                for chip in tile.chips[band]:
                    # Reset injection counter
                    chip.types_injected = 0
                    if config.sim_gals is True:
                        # Determine which Balrog galaxies are contained in chip, and
                        # get their image coordinates
                        in_chip, pos_im = chip.contained_in_chip(tile.gals_pos[real])
                        gals_pos_im = pos_im[in_chip]
                        gals_indx = tile.gals_indx[real][in_chip]
                        chip.set_Ngals(len(gals_pos_im))

                        # pudb.set_trace()
                        # If any galaxies are contained in chip, add gs injection to config list
                        # (or add dummy injection for a few cases)
                        if (chip.Ngals > 0) or (config.inj_objs_only['value'] is True):
                            tile.add_gs_injection(config, chip, gals_indx, gals_pos_im, inj_type='gals')

                    if config.sim_stars is True:
                        # Determine which Balrog stars are contained in chip, and
                        # get their image coordinates
                        in_chip, pos_im = chip.contained_in_chip(tile.stars_pos[real])
                        stars_pos_im = pos_im[in_chip]
                        stars_indx = tile.stars_indx[real][in_chip]
                        chip.set_Nstars(len(stars_pos_im))

                        # If any stars are contained in chip, add gs injection to config list
                        # (or add dummy injection for a few cases)
                        if (chip.Nstars > 0) or (config.inj_objs_only['value'] is True):
                        # if chip.Nstars > 0:
                            tile.add_gs_injection(config, chip, stars_indx, stars_pos_im, inj_type='stars')

                    if (config.inj_objs_only['value'] is False) and (chip.Ngals + chip.Nstars) == 0:
                        # Don't want to skip image for a blank run; need to blank out the image!
                        # NOTE: The first check isn't actually required, as a dummy injection is
                        # added for 'inj_objs_only'=True cases that have no injections. This is
                        # needed to correct the background image in `injector.py`
                        # TODO: Eventually use a return_output_name() function
                        outfile = os.path.join(config.output_dir, 'balrog_images', str(tile.curr_real),
                                               config.data_version, tile.tile_name, chip.band,
                                               '{}_balrog_inj.fits'.format(chip.name))
                        chip.save_without_injection(outfile)

            # Once all chips in tile have had Balrog injections, run modified config file
            # with GalSim
            if vb: print('Writing Balrog config...')
            tile.write_bal_config()
            if vb: print('Running GalSim for tile...')
            rc = tile.run_galsim(vb=vb)
            if rc != 0: raise Exception('\nYou shall not pass! GalSim failed to complete successfully.')
            if vb: print('Copying extra image planes...')
            tile.copy_extensions(config)
            if vb: print('Truth Catalog...')
            tile.write_truth_catalog(config)

    return 0

if __name__ == '__main__':
    ret = RunBalrog()
    sys.exit(ret)
