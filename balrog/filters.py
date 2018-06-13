import numpy as np
import os
# import pudb

# NOTE: Only DES bands supported at the moment, but easy to add more
_valid_bands = 'ugrizy'

class Filter(object):

    def __init__(self, band, transmission_table=None):

        if band not in _valid_bands:
            raise ValueError('Band {} does not have bandpass information stored yet '.format(band) +
                             'in `_setup_bandpass_config`! Valid bands are {}'.format(_valid_bands))
        self.band = band

        if transmission_table:
            if os.path.isfile(transmission_table):
                self.transmission_table = transmission_table
                self.has_transmission_table = True
            else:
                raise ValueError('Passed transmission lookup-table ' +
                                 '{} is not a file!'.format(transmission_table))
        else:
            self.transmission_table = None
            self.has_transmission_table = False

        self._setup_bandpass_config()

        return

    def _setup_bandpass_config(self):
        '''
        This sets up the config used by GalSim's bandpass object.
        Filter information taken from:
        http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=CTIO&gname2=DECam
        '''

        self.band_config = {}

        if self.has_transmission_table is True:
            # Red / blue limits are determined by the lookup table
            self.band_config['throughput'] = self.transmission_table
        else:
            # Taken from  http://www.ctio.noao.edu/noao/node/2232
            self.band_config['throughput'] = '0.85'

            # DES Bandpass information
            if self.band is 'u':
                # lambda_cen = 350. # nm
                # lambda_width = 75. # nm
                # self.band_config['blue_limit'] = lambda_cen - lambda_width
                # self.band_config['red_limit']  = lambda_cen + lambda_width
                self.band_config['blue_limit'] = 303.4 # nm
                self.band_config['red_limit']  = 403.9 # nm
            elif self.band is 'g':
                # lambda_cen = 475. # nm
                # lambda_width = 75. # nm
                # self.band_config['blue_limit'] = lambda_cen - lambda_width
                # self.band_config['red_limit']  = lambda_cen + lambda_width
                self.band_config['blue_limit'] = 392.9 # nm
                self.band_config['red_limit']  = 555.4 # nm
            elif self.band is 'r':
                # lambda_cen = 635. # nm
                # lambda_width = 150. # nm
                # self.band_config['blue_limit'] = lambda_cen - lambda_width
                # self.band_config['red_limit']  = lambda_cen + lambda_width
                self.band_config['blue_limit'] = 561.8 # nm
                self.band_config['red_limit']  = 726.0 # nm
            elif self.band is 'i':
                # lambda_cen = 775. # nm
                # lambda_width = 150. # nm
                # self.band_config['blue_limit'] = lambda_cen - lambda_width
                # self.band_config['red_limit']  = lambda_cen + lambda_width
                self.band_config['blue_limit'] = 698.3 # nm
                self.band_config['red_limit']  = 870.6 # nm
            elif self.band is 'z':
                # lambda_cen = 925. # nm
                # lambda_width = 150. # nm
                # self.band_config['blue_limit'] = lambda_cen - lambda_width
                # self.band_config['red_limit']  = lambda_cen + lambda_width
                self.band_config['blue_limit'] = 836.0 # nm
                self.band_config['red_limit']  = 1016.6 # nm
            elif self.band is 'y':
                # lambda_cen = 1000. # nm
                # lambda_width = 110. # nm
                # self.band_config['blue_limit'] = lambda_cen - lambda_width
                # self.band_config['red_limit']  = lambda_cen + lambda_width
                self.band_config['blue_limit'] = 940.0 # nm
                self.band_config['red_limit']  = 1080.5 # nm
            else:
                raise ValueError('Band {} does not have bandpass information stored '.format(band) +
                                'yet in `get_bandpass_config`!')

        self.band_config['wave_type'] = 'nm'

        return

class Filters(dict):
    '''
    Structure to hold all simulation filters.
    '''

    def __init__(self, bands, use_transmission_tables=True, filter_dir=None):

        # Band validity checked in Filter() construction
        self.bands = bands

        self.use_transmission_tables = use_transmission_tables
        self.transmission_tables = {}

        # pudb.set_trace()

        if filter_dir:
            if not os.path.isdir(filter_dir):
                raise ValueError('The passed filter_dir={} is not a directory!'.format(filter_dir))
            self.filter_dir = filter_dir
        else:
            # Default is to expect filters in current directory
            self.filter_dir = os.getcwd()

        for band in bands:
            if use_transmission_tables:
                table_name = transmission_table_name[band]
                self.transmission_tables[band] = os.path.abspath(
                    os.path.join(self.filter_dir, table_name))
            else:
                self.transmission_tables[band] = None

            self[band] = Filter(band, transmission_table=self.transmission_tables[band])

        return

transmission_table_name = {
    'u' : 'CTIO_DECam.u_filter.dat',
    'g' : 'CTIO_DECam.g_filter.dat',
    'r' : 'CTIO_DECam.r_filter.dat',
    'i' : 'CTIO_DECam.i_filter.dat',
    'z' : 'CTIO_DECam.z_filter.dat',
    'y' : 'CTIO_DECam.Y_filter.dat'
}
