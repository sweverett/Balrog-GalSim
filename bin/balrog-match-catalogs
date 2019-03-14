#!/usr/bin/env python
"""
Match Balrog catalogs and save matched catalogs,
along with detection and extinction information
"""

from argparse import ArgumentParser
import match
import os

parser=ArgumentParser()

parser.add_argument(
    'base',
    type=str,
    help='Base directory where tile directories are located'
)
parser.add_argument(
    '--conf',
    default='y3v02',
    type=str,
    help='MEDS configuration used'
)
parser.add_argument(
    '--real',
    default=0,
    type=int,
    help='band to download',
)
parser.add_argument(
    '--tile_list',
    default=None,
    type=list,
    help='List of tiles to stack (if not all)'
)
parser.add_argument(
    '--inj_type',
    default='gals',
    type=str,
    help='Injection type to match between catalogs (gals, stars, or both)'
)
parser.add_argument(
    '--match_radius',
    default=0.5,
    type=float,
    help='Match radius in arcsec'
)
parser.add_argument(
    '--outdir',
    default=None,
    type=str,
    help='Output location for stacked catalogs'
)
# parser.add_argument(
#     '--ext_file',
#     default=None,
#     type=str,
#     help='Location of extinction info file (to calculate dereddened quantities)'
# )
parser.add_argument(
    '--save_mags',
    action='store_true',
    default=False,
    help='Use to save truth mags to detection catalog'
)
parser.add_argument(
    '--clobber',
    action='store_true',
    default=False,
    help='Use to overwrite existing match files'
)
parser.add_argument(
    '--cache',
    action='store_true',
    default=False,
    help='Cache individual tile matches before stacking'
)
parser.add_argument(
    '--det_only',
    action='store_true',
    default=False,
    help='Set to only write out the detection catalog'
)
# parser.add_argument(
#     '--no_deredden',
#     action='store_false',
#     default=True,
#     help='Set to *not* calculate dereddened quantities for det catalog'
# )
parser.add_argument(
    '--verbose',
    action='store_true',
    default=False,
    help='Set to print out more information'
)
# TODO: Implement!
parser.add_argument(
    '--clean',
    action='store_true',
    help=('Remove all existing matched catalogs in base'),
)

def main():
    args = parser.parse_args()
    vb = args.verbose

    if args.clean is True:
        # Clean out existing matched catalogs
        # ...
        pass

    if args.outdir is None:
        outdir = os.getcwd()
    else:
        if not os.path.isdir(os.path.abspath(args.outdir)):
            raise ValueError('{} is not an existing directory!'.format(args.outdir))
        outdir = args.outdir

    # if args.ext_file is not None:
    #     ef = os.path.abspath(os.path.expanduser(args.ext_file))
    #     if not os.path.exists(ef):
    #         raise IOError('Passed extinction file {} does not exist!'.format(ef))
    #     # Need to save `true_cm_mag` in det catalog to compute dereddened quantities
    # ext_file = ef
    save_mags = args.save_mags

    if vb:
        print('Matching catalogs...')

    # NOTE: Can pass lots of more parameters here. Use outside of script if desired.
    matched_cats = match.MatchedCatalogs(args.base,
                                         meds_conf=args.conf,
                                         real=args.real,
                                         tile_list=args.tile_list,
                                         inj_type=args.inj_type,
                                         match_radius=args.match_radius/3600.0,
                                         vb=vb)

    if args.cache:
        # Write out individual matched catalogs
        if vb:
            print('Writing combined catalogs...')

        outdirs = []
        cat_outdir = os.path.join(outdir, 'matched_catalogs')
        outdirs.append(cat_outdir)

        if not args.det_only:
            det_outdir = os.path.join(outdir, 'det_catalogs')
            outdirs.append(det_outdirs)

        for d in outdirs:
            if not os.path.exists(d):
                os.mkdir(d)

        matched_cats.write_det_cats(outdir=det_outdir, clobber=args.clobber)

        if not args.det_only:
            matched_cats.write_combined_cats(outdir=cat_outdir, clobber=args.clobber)

    # Write out detection truth stack
    if vb:
        print('Writing detection truth stack...')
    det_outfile = matched_cats.write_truth_det_stack(outdir=outdir, clobber=args.clobber, save_mags=save_mags)

    # if ext_file:
    #     if vb:
    #         print('Calculating dereddened magnitudes...')
    #     calc_dereddened_mags(det_outfile, ext_file)

    if not args.det_only:
        # Write out a concatonated true & meas matched catalog
        if vb:
            print('Writing combined stack...')
        matched_cats.write_combined_stack(outdir=outdir, cache=args.cache, clobber=args.clobber)

    # To grab full truth catalog (detected or not), use the following:
    # full_truth_stack = matched_cats.get_full_true_stack()

    # To grab the matched truth & measured catalogs individually (matched by row),
    # use the following:
    # true_stack, meas_stack = matched_cats.get_matched_stack()

    # A few more examples:

    # matched_cats.write_stacks() # All (with matched cats separate)

    # matched_cats.write_full_truth_stack() # Individually

if __name__=="__main__":
    main()
