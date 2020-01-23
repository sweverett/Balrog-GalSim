from balrog.metacal import convert_mcal_to_h5
import os

parser = ArgumentParser()

parser.add_argument(
    'mcal_file',
    default=None,
    type=str,
    help='Filename of mcal stack to rebuild'
)
parser.add_argument(
    'outdir',
    default=None,
    type=str,
    help='Output location of mcal stack to rebuild'
)
parser.add_argument(
    '--mcal_type',
    default='riz-noNB',
    choices=['griz', 'griz-noNB', 'riz', 'riz-noNB'],
    type=str,
    help='Set the type of MCalCatalog created'
)
parser.add_argument(
    '--match_type',
    default='default',
    choices=['default', 'mof_only', 'sof_only'],
    type=str,
    help='Set the type of MatchedCatalog used (NB: not the same as ngmix_type!)'
)
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
)


if __name__ == "__main__":
    args = parser.parse_args()
    mcal_file = args.mcal_file
    outdir = args.outdir
    mcal_type = args.vb
    match_type = args.match_type
    vb = args.vb

    cache_dir = os.path.join(outdir, 'cache')

    s = mcal_type.split('-')
    if len(s) == 1:
        mcal_type = [(s[0], True)]
    else:
        mcal_type = [(s[0], False)]

    b, n = mcal_type[0], mcal_type[1]
    if n is True:
        nb = 'NB'
    else:
        nb = 'noNB'
    type_dir = os.path.join(cache_dir, b+'_'+nb)

    if vb:
        print('Rebuilding mcal cache to hdf5 stack...')
    convert_mcal_to_h5(type_dir, mcal_file, b, match_type=match_type)
