import fitsio
import numpy as np
import h5py
from astropy.table import Table, Column, vstack, join
import os
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'version',
    type=str,
    help='version of merged stack'
    )
parser.add_argument(
    '-run_list',
    nargs='+',
    type=str,
    help='list of run_names to merge (e.g. run2, run2a)'
    )
parser.add_argument(
    '-version_list',
    nargs='+',
    type=float,
    help='list of run versions corresponding to run_list (e.g. 1.4 1.4)'
    )
parser.add_argument(
    '--basedir',
    type=str,
    default='/data/des41.b/data/severett/Balrog/',
    help='Base directory of run outputs'
    )
parser.add_argument(
    '--ngmixer_type',
    type=str,
    default='sof',
    help='ngmixer photometry type to merge'
    )
parser.add_argument(
    '--outdir',
    default=None,
    type=str,
    help='Output location for merged catalogs'
    )
# parser.add_argument(
#     '--clobber',
#     action='store_true',
#     default=False,
#     help='Use to overwrite existing match files'
#     )
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
    )

def remove_duplicates(cat, vb=True):
    '''
    Balrog stack versions 1.4 and below have a small bug that
    seems to duplicate exactly 1 object, so check for these
    '''
    unq, unq_idx, unq_cnt = np.unique(cat['bal_id'],
                                      return_inverse=True,
                                      return_counts=True)
    Nunq = len(unq)
    Nobjs = len(cat)
    if Nunq != Nobjs:
        Ndups = Nobjs - Nunq
        dup_ids = unq[np.where(unq_cnt > 1)]
        if vb:
            print('Catalog has {} duplicate(s)!'.format(Ndups))
            print('Removing the following duplicates from catalog:')
            print(dup_ids)

        Nbefore = Nobjs
        for did in dup_ids:
            indx = np.where(cat['bal_id']==did)[0]

            L = len(indx)
            for i in range(L-1): # keep last one
                cat.remove_row(indx[i])

        Nobjs = len(cat)
        assert Nobjs == (Nbefore - Ndups)

        print('{} duplicates removed, catalog size now {}'.format(Ndups, Nobjs))

    return dup_ids

def compute_injection_counts(det_catalog):
    # `true_id` is the DF id
    unique, ucounts = np.unique(det_catalog['true_id'], return_counts=True)

    freq = Table()
    freq['true_id'] = unique
    freq['injection_counts'] = ucounts

    joined = join(det_catalog, freq, keys='true_id', join_type='left')

    assert len(det_catalog) == len(joined)

    return joined

def get_datasets(key, h):

    if key[-1] != '/':
        key += '/'

    out = []

    for name in h[key]:
        path = key + name

        if isinstance(h[path], h5py.Dataset):
            out += [path]
        else:
            out += get_datasets(path, h)

    return out

def stack_mcal_files(mcal_list, run_list, outfile, vb=True):
    '''
    expects a list of mcal h5 file objects (already loaded into memory)
    # Partially solved here:
    # https://stackoverflow.com/questions/49851046/merge-all-h5-files-using-h5py?noredirect=1&lq=1
    '''

    if os.path.exists(outfile):
        raise OSError('File already exists!')

    mcal_stack = h5py.File(outfile, 'w')

    k = 0
    Nmcal = len(mcal_list)
    cat_length = []
    total_length = 0
    datasets = None
    for mcal in mcal_list:
        k += 1
        mcal_datasets = get_datasets('/', mcal)
        if k == 1:
            datasets = mcal_datasets
        else:
            assert mcal_datasets == datasets

        l = len(mcal['catalog/unsheared/bal_id'])
        cat_length.append(l)
        total_length += l

    # Gets path name up until dataset
    # NOTE: This funky looking indexing grabs all parts of the path except
    # the final entry (colname) by reversing the text order, removing the first
    # path entry, and then reversing the text order again
    groups = list(set([i[::-1].split('/', 1)[1][::-1] for i in datasets]))
    groups = [i for i in groups if len(i) > 0]

    # sort groups based on depth
    idx    = np.argsort(np.array([len(i.split('/')) for i in groups]))
    groups = [groups[i] for i in idx]

    # create all groups that contain dataset that will be copied
    if vb is True:
        print('Creating groups & datasets in merged file')

    for group in groups:
        print('Creating group ', group)
        mcal_stack.create_group(group)

    for dataset in datasets:
        print('Creating dataset ', dataset)
        mcal_stack.create_dataset(dataset, maxshape=(total_length,), shape=(total_length,), dtype=mcal_list[0][dataset].dtype)

    # Add a few more explicitly
    mcal_stack.create_dataset('/catalog/unsheared/run_name', maxshape=(total_length,), shape=(total_length,), dtype='S10')

    # copy mcal catalogs
    if vb is True:
        print('Copying mcal catalogs...')
    iter_end = 0
    for l, mcal, rname in zip(cat_length, mcal_list, run_list):
        if vb is True:
            print('mcal = ', mcal)
        for path in datasets:
            # get group name
            group = path[::-1].split('/',1)[1][::-1]
            # minimum group name
            if len(group) == 0: group = '/'
            # copy data
            if vb is True:
                print('copying path {}...'.format(path))
            mcal_stack[path][iter_end:iter_end+l] = mcal[path][:]
        # Add a few more explicitly
        mcal_stack['/catalog/unsheared/run_name'][iter_end:iter_end+l] = Column(l*[rname], dtype='S10')

        iter_end += l

    return

def main():
    args = parser.parse_args()
    run_list = list(args.run_list)
    version_list = list(args.version_list)
    version = args.version
    basedir = args.basedir
    ngtype = args.ngmixer_type
    outdir = args.outdir
    vb = args.vb

    if outdir is None:
        outdir = os.getcwd()

    det_outdir = os.path.join(outdir, version, ngtype)
    det_outfile = os.path.join(det_outdir, 'balrog_detection_catalog_{}_y3-merged_v{}.fits'.format(
        ngtype, version))
    match_outdir = os.path.join(outdir, version, ngtype)
    match_outfile = os.path.join(match_outdir, 'balrog_matched_catalog_{}_y3-merged_v{}.fits'.format(
        ngtype, version))
    mcal_outdir = os.path.join(outdir, version, 'mcal')
    mcal_outfile = os.path.join(mcal_outdir, 'balrog_mcal_stack-y3v02-0-riz-noNB_y3-merged_v{}.h5'.format(
        version))

    Nobjects = {}
    det_list   = []
    match_list = []
    mcal_list  = []

    if vb is True:
        print('Merging the following Balrog runs: {}'.format(run_list))

    k = 0
    Nruns = len(run_list)
    Ndetections = {}
    Nmatches = {}
    for run, ver in zip(run_list, version_list):
        ver = str(ver)
        k += 1
        if vb is True:
            print('Starting {} v{} ({} of {})'.format(run, ver, k, Nruns))

        # Get the needed catalog file paths
        det_fname   = 'balrog_detection_catalog_sof_{}_v{}.fits'.format(run, ver)
        match_fname = 'balrog_matched_catalog_sof_{}_v{}.fits'.format(run, ver)
        mcal_fname  = 'balrog_mcal_stack-y3v02-0-riz-noNB-mcal_{}_v{}.h5'.format(run, ver)
        det_file   = os.path.join(basedir, run, 'stacked_catalogs', ver, ngtype, det_fname)
        match_file = os.path.join(basedir, run, 'stacked_catalogs', ver, ngtype, match_fname)
        mcal_file  = os.path.join(basedir, run, 'stacked_catalogs', ver, 'mcal', mcal_fname)

        if vb is True:
            print('Loading sof detection catalog...')
        det = Table(fitsio.read(det_file))
        Ndetections[run] = len(det)
        if vb is True:
            print('Loading sof matched catalog...')
        match = Table(fitsio.read(match_file))
        Nmatches[run] = len(match)
        if vb is True:
            print('Grabbing mcal catalog...')
        mcal = h5py.File(mcal_file, 'r')

        if vb is True:
            print('Adding new columns...')
        det.add_column(Column(len(det)*[run], name='run_name', dtype='S10'))
        match.add_column(Column(len(match)*[run], name='run_name', dtype='S10'))
        ## NOTE: Column added to stacked h5 file later in stack_mcal_files()

        det_list.append(det)
        match_list.append(match)
        mcal_list.append(mcal)

        del det
        del match
        del mcal

    # Stack catalogs
    if vb is True:
        print('Stacking detection catalogs...')
    det_stack = vstack(det_list)
    if vb is True:
        print('Stacking matched catalogs...')
    match_stack = vstack(match_list)

    # Compute injection counts
    det_stack = compute_injection_counts(det_stack)

    # Look for and remove any duplicate entries
    dup_ids1 = remove_duplicates(det_stack)
    dup_ids2 = remove_duplicates(match_stack)
    assert (dup_ids1 == dup_ids2).all()
    # NOTE: The bug that causes duplicates in the fits tables doesn't appear
    # in the mcal stack, so not needed there (which is good, as h5 files are
    # much more complicated to remove rows from)

    # Save merged catalogs
    if vb is True:
        print('Writing merged detection catalog...')
    det_stack.write(det_outfile)
    if vb is True:
        print('Writing merged matched catalog...')
    match_stack.write(match_outfile)
    if vb is True:
        print('Stacking & writing mcal catalog...')
    mcal_stack = stack_mcal_files(mcal_list, run_list, mcal_outfile, vb=vb)

    if vb is True:
        print('Finished!')

    return

if __name__ == '__main__':
    main()
