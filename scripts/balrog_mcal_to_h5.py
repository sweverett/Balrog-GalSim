from argparse import ArgumentParser

from balrog.metacal import convert_mcal_to_h5

parser = ArgumentParser()

parser.add_argument(
    'catdir',
    type=str,
    help='Directory location of cached mcal value-adds from Balrog'
    )
parser.add_argument(
    'outfile',
    type=str,
    help='Output filename for h5 stack'
    )
parser.add_argument(
    '--version_tag',
    default=None,
    type=str,
    help='Mcal Balrog stack version'
)
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
)


if __name__ == "__main__":
    args = parser.parse_args()
    convert_mcal_to_h5(args.catdir, args.outfile, version_tag=args.version_tag, vb=args.vb)
    sys.exit(0)
