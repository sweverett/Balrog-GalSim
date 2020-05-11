import fitsio
from argparse import ArgumentParser
import os
from glob import glob

parser = ArgumentParser()

parser.add_argument(
    'cache_dir',
    type=str,
    help='Directory of cache files'
    )
parser.add_argument(
    'outfile',
    type=str,
    help='Name of output stacked catalog'
    )
parser.add_argument(
    '-outdir',
    type=str,
    help='Directory of output file'
    )
parser.add_argument(
    '-ext',
    type=int,
    default=1,
    help='FITS extension to stack'
    )
parser.add_argument(
    '--clobber',
    action='store_true',
    default=False,
    help='Set to overwrite existing file'
    )
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
    )
parser.add_argument(
    '--test',
    action='store_true',
    default=False,
    help='Set to only stack a few files'
    )

def main():
    args = parser.parse_args()
    cache_dir = args.cache_dir
    outfile = args.outfile
    outdir = args.outdir
    ext = args.ext
    clobber = args.clobber
    vb = args.vb
    test = args.test

    if outdir is None:
        outdir = os.getcwd()
    elif os.path.exists(outdir) is False:
        raise OSError('{} does not exist!'.format(outdir))

    outfile = os.path.join(outdir, outfile)
    if os.path.exists(outfile):
        if clobber is True:
            os.remove(outfile)
        else:
            raise OSError('File already exists!')

    stack = fitsio.FITS(outfile, 'rw')

    catfiles = glob(os.path.join(cache_dir, 'DES*balrog*.fits'))
    N = len(catfiles)

    if test is True:
        if N >=5:
	    df = N / 5
	    catfiles = catfiles[::df]
	    N = len(catfiles)

    # Need to know total number of objects for efficient memory stacking
    size = 0
    if vb is True:
        print('Prepping stack...')
    # TODO: This can be removed depending on how efficient fitsio appending is
    for catfile in catfiles:
        h = fitsio.read_header(catfile, ext=ext)
        size += h['NAXIS2']
    print('Total stack size will be {}'.format(size))

    if vb is True:
        print('Stacking...')
    k = 0
    for catfile in catfiles:
        k += 1
        if vb is True:
            cfile = os.path.basename(catfile)
            print('Stacking {} ({} of {})'.format(cfile, k, N))
	if k == 1:
	    stack.write(fitsio.read(catfile))
	else:
            stack[ext].append(fitsio.read(catfile))

    assert stack[1].get_nrows() == size

    if vb is True:
        print('Done!')

if __name__ == '__main__':
    main()
