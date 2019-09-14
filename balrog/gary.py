import numpy as np
import fitsio
from astropy.table import Table, vstack, join
import os

class SetupGary(object):
    _valid_match_methods = ['none', 'unambiguous', 'ambiguous']
    _needs_data = ['unambiguous', 'ambiguous']

def setup_gary(meas_file, true_file, match_type=None, **kw):
    '''
    '''
    sg = SetupGary()
    if match_type == 'simple':
        return SimpleGary(meas_file, true_file, **kw)
    elif match_type == 'unambiguous':
        return UnambiguousGary(meas_file, true_file, **kw)
    elif match_type == 'ambiguous':
        return AmbiguousGary(meas_file, true_file, **kw)
    elif match_type is None:
        return AmbiguousGary(meas_file, true_file, **kw)
    else:
        raise ValueError('match_type={} is not a valid matching type!'.format(match_type) +
                         ' Must choose from {}'.format(sg._valid_match_methods))

class SimpleGary(object):
    _default_cols_simple = ['RA', 'DEC']
    _default_cols_simple_GOLD = ['RA', 'DEC']
    _default_cols_simple_ngmixer = ['ra', 'dec']

class Gary(SetupGary):
    '''
    Matching algorithm that accounts for ambiguous blends in a
    systematic way.

    Will put more details here after completion of the algorithm.

    Blending inclusion types:
    'none'        :
    'unambiguous' :
    'ambiguous'   :

    Ideas for merging w/ match.py:
    - Changing `Gary` to htm.Matcher format
    - incorporate some of the `col` functions into MatchedCatalog
    - incorporate match_method functions into MatchedCatalog
    '''

    _default_cols_ngmixer = ['ra', 'dec', 'cm_mag']
    _default_cols_data = ['RA', 'DEC', 'CM_MOF_MAG_I']

    # # TODO: Turn this on when we know how!
    # MATCH_BY_METHOD = {
    #     'simple' : gary._simple_match(),
    #     'unambiguous' : _unambiguous_match,
    #     'ambiguous' : _ambiguous_match
    # }

    def __init__(self, meas_file, true_file, data_file=None, match_method='ambiguous',
                 tilename=None, meas_cols=None, true_cols=None, data_cols=None):

        # # TODO: MOVE THIS TO CLASS VARIABLE WHEN WE UNDERSTAND HOW!
        self.MATCH_BY_METHOD = {
            'simple' : self._simple_match,
            'unambiguous' : self._unambiguous_match,
            'ambiguous' : self._ambiguous_match
        }

        files = {'meas':meas_file, 'true':true_file}
        if data_file is not None:
            files['data'] = data_file
        import pudb
        pudb.set_trace()
        for f in files.values():
            if not os.path.isfile(os.path.expanduser(f)):
                raise OSError('{} is not a file!'.format(f))

        self.meas_file = meas_file
        self.true_file = true_file
        self.data_file = data_file

        if match_method not in self._valid_match_methods:
            raise ValueError('{} is not a valid match method: {}'.format(
                             match_method, self._valid_match_methods))

        if match_method in self._needs_data:
            if data_file is None:
                raise ValueError('Must pass a data catalog file for ' +
                                 'match_method={}'.format(match_method))
            self.has_data = True
        else:
            if data_file is not None:
                raise ValueError('Do not pass in a Data catalog for ' +
                                 'match_method={}'.format(match_method))
            self.has_data = False

        self.match_method = match_method

        if tilename is not None:
            if not isinstance(tilename, str):
                raise TypeError('Tilename must be a str!')
        self.tilename = tilename

        self._set_load_cols(meas_cols, true_cols, data_cols)

        self._load_cats()

        self._match()

        #...

        return

    def _set_load_cols(self, meas_cols, true_cols, data_cols):
        self.cols = {}
        for name, cols in {'meas':meas_cols, 'true':true_cols, 'data':data_cols}.items():
            if cols is None:
                if self.match_method in self._needs_data:
                    self.cols[name] = self._default_cols_data
                else:
                    self.cols[name] = self._default_cols_simple
            else:
                self.cols[name] = cols

        return

    def _load_cats(self):
        import pudb
        pudb.set_trace()
        self.meas = Table(fitsio.read(self.meas_file, columns=self.cols['meas']))
        self.true = Table(fitsio.read(self.true_file, columns=self.cols['true']))

        if self.has_data is True:
            if self.tilename is None:
                # If no tilename has been set at construction time, then we
                # assume the file has not been indexed by tile
                self.data = pd.read_hdf(self.data_file)
            else:
                # This catalog is huge, so we assume it has been saved
                # in a H5 file indexed by tilename if it is set
                self.data = pd.read_hdf(self.data_file, key=self.tilename)
        else:
            self.data = None

        #...

        return

    def _match(self):
        # General setup here...

        self.MATCH_BY_METHOD[self.match_method]()

        return


    def _simple_match(self):
        print('Simple running...')
        return

    def _unambiguous_match(self):
        print('Unambiguous running...')
        return

    def _ambiguous_match(self):
        print('Ambiguous running...')
        return

def build_matched_catalog(match_type, *args, **kwargs):
    if match_type in MATCHED_CATALOG_TYPES:
        return MATCHED_CATALOG_TYPES[match_type](*args, **kwargs)
    else:
        raise ValueError('{} is not a valid matched catalog type. Allowed tyes are: {}'.format(
            match_type, MATCHED_CATALOG_TYPES.keys()))
