#!/usr/bin/env python
"""
Match Balrog catalogs and save matched catalogs,
along with detection and extinction information
"""

from argparse import ArgumentParser
import os

from balrog import match

parser=ArgumentParser()

parser.add_argument(
    'base',
    type=str,
    help='Base directory where tile directories are located'
)
parser.add_argument(
    '--match_type',
    default='default',
    choices=['default', 'mof_only', 'sof_only'],
    type=str,
    help='Set the type of MatchedCatalog created (NB: not the same as ngmix_type!)'
)
parser.add_argument(
    '--version',
    default=None,
    type=str,
    help='Stack version'
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
    choices=['gals', 'stars', 'both'],
    type=str,
    help='Injection type to match between catalogs (gals, stars, or both)'
)
parser.add_argument(
    '--ngmix_type',
    default='mof',
    choices=['mof', 'sof'],
    type=str,
    help='ngmix type to match (mof or sof). Defaults to mof'
)
parser.add_argument(
    '--ngmix_profile',
    default='bdf',
    type=str,
    help='ngmix profile type used as col prefix'
)
parser.add_argument(
    '--gold_base',
    default=None,
    type=str,
    help='Base directory of gold products '
    'Use to add a few extra value-adds from the gold-like cats'
)
parser.add_argument(
    '--gold_subdir',
    default=None,
    type=str,
    help='Subdirectory of desired gold merged cat, if not in TILENAME'
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
parser.add_argument(
    '--save_mags',
    default=True,
    type=bool,
    help='Use to save truth mags to detection catalog'
)
parser.add_argument(
    '--save_gap_flux',
    default=False,
    type=bool,
    help='Use to save avg riz Gaussian aperture fluxes in detection catalog'
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
parser.add_argument(
    '--vb',
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

def add_version(fname, version):
    if version is not None:
        fname = fname.replace('.fits', '_v{}.fits'.format(version))
    return fname

def main():
    args = parser.parse_args()
    vb = args.vb

    if args.clean is True:
        # Clean out existing matched catalogs
        # ...
        pass

    if args.outdir is None:
        outdir = os.getcwd()
    else:
        if not os.path.isdir(os.path.abspath(args.outdir)):
            os.mkdir(args.outdir)
        outdir = args.outdir

    if args.inj_type == 'stars':
        if args.ngmix_type is not None:
            # ngmix_type not relevant for star case
            args.ngmix_type = None

    if args.gold_base is not None:
        if args.gold_base == 'base':
            args.gold_base = args.base
        else:
            if not os.path.exists(args.gold_base):
                raise OSError('{} does not exist!'.format(args.gold_base))

    if args.gold_subdir is not None:
        if args.gold_base is None:
            raise ValueError('Can\'t set gold_subdir if gold_base isn\'t set!')

    save_mags = args.save_mags
    save_gap_flux = args.save_gap_flux

    if vb:
        print('Matching catalogs...')

    # NOTE: Can pass lots of more parameters here. Use outside of script if desired.
    matched_cats = match.MatchedCatalogs(args.base,
                                         match_type=args.match_type,
                                         meds_conf=args.conf,
                                         real=args.real,
                                         tile_list=args.tile_list,
                                         inj_type=args.inj_type,
                                         ngmix_type=args.ngmix_type,
                                         match_radius=args.match_radius/3600.0,
                                         extra_base=args.gold_base,
                                         extra_subdir=args.gold_subdir,
                                         prefix=args.ngmix_profile+'_',
                                         de_reddened=True,
                                         save_gap_flux=save_gap_flux,
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
            outdirs.append(det_outdir)

        for d in outdirs:
            if not os.path.exists(d):
                os.mkdir(d)

        matched_cats.write_det_cats(outdir=det_outdir, clobber=args.clobber)

        if not args.det_only:
            matched_cats.write_combined_cats(outdir=cat_outdir,
                                             clobber=args.clobber)

    # Write out detection truth stack
    if vb:
        print('Writing detection truth stack...')
    det_outfile = add_version('balrog_truth_det_stack.fits', args.version)
    det_outfile = matched_cats.write_truth_det_stack(outdir=outdir,
                                                     clobber=args.clobber,
                                                     save_mags=save_mags,
                                                     save_gap_flux=save_gap_flux,
                                                     outfile=det_outfile)

    if not args.det_only:
        # Write out a concatonated true & meas matched catalog
        if vb:
            print('Writing combined stack...')
        match_outfile = add_version('balrog_matched_catalog.fits', args.version)
        matched_cats.write_combined_stack(outfile=match_outfile,
                                          outdir=outdir,
                                          cache=args.cache,
                                          clobber=args.clobber)

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
