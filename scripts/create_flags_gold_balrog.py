#!/usr/bin/env python
'''
create_flags_gold.py
This script will append a FLAGS_GOLD column according to Y3 Gold definitions
Usage: python create_flags_gold.py -d [merged_files_dir]
Author: Nacho Sevilla (nsevilla@gmail.com) 
'''
import sys
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy import coordinates as coo
import os
import glob
from optparse import OptionParser

def main():

    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--datadir",dest="data_dir",help="Directory where data is read from",default='/home/nsevilla/y3_balrog/uploads/')
    # Parse command line
    (options, args) = parser.parse_args()

    if not os.path.isdir(options.data_dir):
        print 'Path with data not found at',options.data_dir
        sys.exit(1)

    os.chdir(options.data_dir)
    objfiles = glob.glob('*merged.fits')
    skipflag = 0b0000000
    skipflag |= 0b100000

    for f,infilename in enumerate(objfiles):
        print 'Reading',infilename
        hdulist = fits.open(infilename)
        tdata = hdulist[1].data
        cols = tdata.columns

        flag_gold = np.zeros(tdata.size, dtype=np.int32)

        if (skipflag & 0b1) == False:
            print 'Processing MOF_FLAGS: flag 1'
            flag_gold[tdata['flags_mof'] != 0] += 1
        if (skipflag & 0b10) == False:
            print 'Processing SOF_FLAGS: flag 2'
            flag_gold[tdata['flags_sof'] != 0] += 2
        if (skipflag & 0b100) == False:
            print 'Processing SOF GAL_FIT_FAILURE: flag 4'
            flag_gold[(tdata['flags_sof'] == 1) | (tdata['flags_sof'] > 2)] += 4
        if (skipflag & 0b1000) == False:
            print 'Processing SExtractor FLAGS: flag 8'
            flag_gold[(tdata['FLAGS_G'] > 3) | (tdata['FLAGS_R'] > 3) | (tdata['FLAGS_I'] > 3) | (tdata['FLAGS_Z'] > 3)] += 8
        if (skipflag & 0b10000) == False:
            print 'Processing IMAFLAGS: flag 16'
            flag_gold[(tdata['IMAFLAGS_ISO_G'] != 0) | (tdata['IMAFLAGS_ISO_R'] != 0) | (tdata['IMAFLAGS_ISO_I'] != 0) | (tdata['IMAFLAGS_ISO_Z'] != 0)] += 6
        if (skipflag & 0b100000) == False:
            print 'Processing BBJ: flag 32'
            flag_gold[(tdata['NEPOCHS_G']==0) & (tdata['MAGERR_AUTO_G'] < 0.05) & 
               ((tdata['MAG_DETMODEL_I'] - tdata['MAG_AUTO_I']) < -0.4)] += 32
        if (skipflag & 0b1000000) == False:
            print 'Processing bright streak rejection: flag 64'
            flag_gold[(((tdata['MAGERR_AUTO_G'] < 0.01) & (tdata['MAG_AUTO_G'] - tdata['cm_mag_g_sof'] < -0.5)) | ((tdata['MAGERR_AUTO_R'] < 0.01) & (tdata['MAG_AUTO_R'] - tdata['cm_mag_r_sof'] < -0.5)) | ((tdata['MAGERR_AUTO_I'] < 0.01) & (tdata['MAG_AUTO_I'] - tdata['cm_mag_i_sof'] < -0.5)) | ((tdata['MAGERR_AUTO_Z'] < 0.01) & (tdata['MAG_AUTO_Z'] - tdata['cm_mag_z_sof'] < -0.5))) & ((tdata['MAG_AUTO_G']-tdata['MAG_AUTO_R'] < -1) | (tdata['MAG_AUTO_G']-tdata['MAG_AUTO_R'] > 4) | (tdata['MAG_AUTO_R']-tdata['MAG_AUTO_I'] < -1) | (tdata['MAG_AUTO_R']-tdata['MAG_AUTO_I'] > 4) | (tdata['MAG_AUTO_I']-tdata['MAG_AUTO_Z']< -1) | (tdata['MAG_AUTO_I']-tdata['MAG_AUTO_Z'] > 4))] += 64    

        print 'Creating new file'
        newcol = fits.Column(name='FLAGS_GOLD',format='I',array=flag_gold)
        #coldefs = fits.ColDefs([c1, c2, c3 ,c4, c5])
        tbhdu_upload = fits.BinTableHDU.from_columns(cols + newcol)
        outfilename = 'merged_wflags.fits' 
        print 'Writing',outfilename
        tbhdu_upload.writeto(outfilename,clobber=True)

if __name__ == '__main__':
    sys.exit(main())

