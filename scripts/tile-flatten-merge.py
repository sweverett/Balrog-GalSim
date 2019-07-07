#!/usr/bin/env python
'''
tile-flatten-merge.py
This script will create a file for upload containing the flattened data in the files from
the mof_dir directory (the MOF files) but only for those objects included in the corresponding
id file in coaddid_dir.
Usage: python tile-flatten-merge.py -d [mof_dir] -i [coaddid_dir] -o [out_dir]

Author: Nacho Sevilla (nsevilla@gmail.com) with some code from Erin Sheldon
Further edited by Spencer Everett to add some new features and make a less IO intensive version,
as we do not need intermediate flattened files for Balrog
Saw ~25% clock time improvement using `save_all=False`, as well as 60% less disk space
'''

import numpy as np
from astropy.io import fits
from astropy.table import Table, join
import os
import sys
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'data_dir',
    type=str,
    help='Directory where data is read from'
    )
parser.add_argument(
    '--out_dir',
    type=str,
    default=None,
    help='Directory where merged data is stored'
    )
parser.add_argument(
    '--save_all',
    action='store_true',
    default=False,
    help='Set to save all flattened files'
    )
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
    )

printout = False
bands = ['g','r','i','z']

onedim_cols = ['id','number','NUMBER','fofid','fof_id','coadd_run','flags','time_last_fit',
               'box_size','bdf_flags','obj_flags','mask_frac','masked_frac','psf_T','psfrec_T',
               'cm_flags','cm_T','cm_T_err','cm_T_s2n','cm_weight','cm_max_flags','cm_max_T',
               'cm_max_T_err','cm_max_T_s2n','cm_s2n_w','cm_chi2per','cm_dof','cm_flags_r',
               'cm_s2n_r','cm_T_r','cm_psf_T_r','cm_fracdev','cm_fracdev_noclip','cm_fracdec_err',
               'cm_TdByTe','cm_TdByTe_noclip','cm_mof_flags','cm_mof_num_itr','fofind','ra','dec',
               'bdf_nfev','bdf_s2n','bdf_T','bdf_T_err','bdf_T_ratio','bdf_fracdev',
               'bdf_fracdev_err','flagstr']

onedim_cols_addband = ['FLAGS','MAG_AUTO','X_IMAGE','Y_IMAGE','XWIN_IMAGE','YWIN_IMAGE',
                       'ERRAWIN_IMAGE','ERRBWIN_IMAGE','ERRTHETAWIN_IMAGE','ALPHAWIN_J2000',
                       'DELTAWIN_J2000','A_IMAGE','B_IMAGE','A_WORLD','B_WORLD','XMIN_IMAGE',
                       'XMAX_IMAGE','YMIN_IMAGE','YMAX_IMAGE','THETA_J2000','FLUX_RADIUS',
                       'IMAFLAGS_ISO','FLUX_AUTO','FLUXERR_AUTO','MAGERR_AUTO','KRON_RADIUS',
                       'BACKGROUND','THRESHOLD','FWHM_IMAGE','FWHM_WORLD']

bidim_cols = ['cm_pars_cov','cm_max_pars_cov','cm_g_cov','bdf_g_cov','bdf_pars_cov']

band_cols = ['nimage_tot','psf_flags','psf_flux','psf_flux_err','psf_flux_s2n','psf_flux_flags',
             'psf_mag','nimage_use','cm_flux_cov','cm_flux','cm_flux_s2n','cm_mag','cm_logsb',
             'cm_max_flux_cov','cm_max_flux','cm_max_flux_s2n','cm_max_mag','cm_max_logsb',
             'bdf_flux','bdf_mag','bdf_flux_err']

bidim_band_cols = ['cm_flux_cov','cm_max_flux_cov','cm_pars_cov',
                   'cm_max_pars_cov','cm_g_cov','bdf_flux_cov']

multi_cols = ['psf_g','psfrec_g','cm_pars','cm_g','cm_max_pars', 'cm_mof_abs_diff',
              'cm_mof_frac_diff','cm_mof_err_diff','bdf_pars','bdf_pars_err','bdf_g']

multi_cols_add = ['FLUX_APER','FLUXERR_APER','MAG_APER','MAGERR_APER']

def get_coldefs(descr, defs={}, bands=None, band_cols=None, bidim_cols=None):
    """
    Convert a numpy descriptor to a set of oracle
    column definitions

    array columns are converted to name_{dim1}_{dim2}...{dimn}

    parameters
    ----------
    descr: numpy type descriptor
        E.g. arr.dtype.descr
    defs: dict,optional
        A dict returning a list of field defs. It is keyed by field names from
        the array.  This can be used to over-ride the defaults, e.g. to use a
        different name or to over-ride conversions for arrays.
    """

    if defs is None:
        defs={}

    alldefs=[]
    def_template='%s not null'

    for d in descr:
        name=d[0]
        ot=get_oracle_type(d[1])

        if name in defs:
            alldefs += defs[name]
        elif len(d) == 2:
            # this is a scalar column... easy!
            defi=def_template % ot
            alldefs.append( (name,defi,d[1]) )
        else:
            dims=d[2]
            if not isinstance(dims,tuple):
                dims=(dims,)

            if (bands is not None
                    and band_cols is not None
                    and name in band_cols):
                names=get_band_arr_colnames(name,dims,bands,bidim_cols)
            else:
                names=get_arr_colnames(name,dims)

            for n in names:
                defi=def_template % (ot)
                alldefs.append( (n,defi,d[1]) )

    return alldefs

def get_arr_colnames(name, dims):
    """
    Get db names for an array, naming
        name_{num1}_{num2}...
    """
    ndim=len(dims)
    if ndim==1:
        names=get_arr1_colnames(name,dims)
    elif ndim==2:
        names=get_arr2_colnames(name,dims)
    else:
        raise ValueError("only support 1 and 2 d arrays")

    return names

def get_arr1_colnames(name, dims):
    """
    Get db names for an array, naming
        name_{num}
    """
    names=[]
    for n in xrange(1,dims[0]+1):
        names.append( '%s_%d' % (name,n) )

    return names

def get_arr2_colnames(name, dims):
    """
    Get db names for an array, naming
        name_{num1}_{num2}
    """
    names=[]
    for n1 in xrange(1,dims[0]+1):
        for n2 in xrange(1,dims[1]+1):
            names.append( '%s_%d_%d' % (name,n1,n2) )

    return names

def get_band_arr_colnames(name, dims, bands, bidim_cols):
    """
    Get db names for an array, naming
        name_{num1}_{num2}...
    """
    ndim=len(dims)
    if ndim==1 and (name not in bidim_cols):
        names=get_band_arr1_colnames(name,dims,bands)
    elif ndim==1 and (name in bidim_cols):
        names=get_band_arr2_colnames(name,[np.sqrt(dims),np.sqrt(dims)],bands)
    elif ndim==2:
        names=get_band_arr2_colnames(name,dims,bands)
    else:
        raise ValueError("only support 1 and 2 d arrays")

    return names

def get_band_arr1_colnames(name, dims, bands):
    """
    Get db names for an array, naming
        name_{num}
    """
    names=[]
    for i in xrange(dims[0]):
        n=bands[i]
        names.append( '%s_%s' % (name,n) )

    return names

def get_band_arr2_colnames(name, dims, bands):
    """
    Get db names for an array, naming
        name_{num1}_{num2}
    """
    names=[]
    for i1 in xrange(dims[0]):
        n1=bands[i1]
        for i2 in xrange(dims[1]):
            n2=bands[i2]
            names.append( '%s_%s_%s' % (name,n1,n2) )

    return names

def get_oracle_type(nt):
    if 'f4' in nt:
        ot='binary_float'
    elif 'f8' in nt:
        ot='binary_double'
    elif 'i1' in nt or 'u1' in nt:
        ot='number(3)'
    elif 'i2' in nt or 'u2' in nt:
        ot='number(5)'
    elif 'i4' in nt:
        ot='number(10)'
    elif 'i8' in nt:
        ot='number(19)'
    elif 'u8' in nt:
        ot='number(20)'
    elif 'S' in nt:
        slen=nt[2:]
        ot='varchar2(%s)' % slen
    else:
        raise ValueError("unsupported numpy type: '%s'" % nt)

    return ot

def get_fits_type(name):
    if name == "id":
        format = 'K'
    elif name.lower() == "number":
        format = 'J'
    elif 'nimage_tot' in name:
        format = 'J'
    elif name == "fofid":
        format = 'K'
    elif name == "fof_id":
        format = 'K'
    elif name == "ra":
        format = 'D'
    elif name == "dec":
        format = 'D'
    elif name == "XMIN_IMAGE":
        format = 'J'
    elif name == "YMIN_IMAGE":
        format = 'J'
    elif name == "XMAX_IMAGE":
        format = 'J'
    elif name == "YMAX_IMAGE":
        format = 'J'
    elif name == "X_IMAGE":
        format = 'E'
    elif name == "Y_IMAGE":
        format = 'E'
    elif name == "XWIN_IMAGE":
        format = 'D'
    elif name  == "YWIN_IMAGE":
        format = 'D'
    elif name == "ERRAWIN_IMAGE":
        format = 'E'
    elif name == "ERRBWIN_IMAGE":
        format = 'E'
    elif name == "ERRTHETAWIN_IMAGE":
        format = 'E'
    elif name == "ALPHAWIN_J2000":
        format = 'D'
    elif name == "DELTAWIN_J2000":
        format = 'D'
    elif name == "A_IMAGE":
        format = 'E'
    elif name == "B_IMAGE":
        format = 'E'
    elif name == "THETA_J2000":
        format = 'E'
    elif "WORLD" in name:
        format = 'E'
    elif "RADIUS" in name:
        format = 'E'
    elif "APER" in name:
        format = 'E'
    elif "AUTO" in name:
        format = 'E'
    elif "FWHM" in name:
        format = 'E'
    elif name == "THRESHOLD":
        format = 'E'
    elif name == "BACKGROUND":
        format = 'E'
    elif name == 'flagstr':
        format = 'A32'
    elif "flags" in name.lower():
        format = 'J'
    elif name == "time_last_fit":
        format = 'D'
    elif name == "box_size":
        format = 'I'
    elif "s2n" in name:
        format = 'D'
    elif "_T" in name:
        format = 'D'
    elif "_mag" in name: 
        format = 'D'
    elif "nimage_use" in name: 
        format = 'J'
    elif "frac" in name: 
        format = 'D'
    elif "_g" in name: 
        format = 'D'
    elif "pars" in name: 
        format = 'D'
    elif "_flux" in name: 
        format = 'D'
    elif "_logsb" in name: 
        format = 'D'
    elif "fracdev" in name:
        format = 'D'
    elif name == "cm_weight": 
        format = 'D'
    elif name == "cm_chi2per": 
        format = 'D'
    elif name == "cm_dof": 
        format = 'D'
    elif name == "cm_TdByTe": 
        format = 'E'
    elif name == "cm_mof_num_itr": 
        format = 'J'
    elif "cm_mof_abs_diff" in name: 
        format = 'D'
    elif "cm_mof_err_diff" in name: 
        format = 'D'
    elif name == "bdf_nfev":
        format = 'J'
    elif name == 'fofind': 
        format = 'K'
    else:
        print name,'not found'
        sys.exit(1)

    return format

def check_name(name,cols,prev,strip,printout=False):
    check = False
    n = 0
    colname = prev
    for enum,col in enumerate(cols):
        if printout:
            print col,name[0:len(name)-strip]
        if col == name[0:len(name)-strip]:
            check = True
            colname = col
            break

    return (check,colname)

def merge(args, tilename, filelist):
    data_dir = args.data_dir
    out_dir = args.out_dir
    save_all = args.save_all
    Nfiles = len(filelist)

    #--------------------------------------------------------------------------------
    # Nacho's original code: The below works if all files are saved out individually,
    # but this is pretty IO intensive and Balrog doesn't need the intermediate files
    if save_all is True:
        flatten(data_dir, filelist, out_dir, tilename, save_all=True)

        for f, fname in enumerate(filelist):
            flatname = os.path.splitext(fname)[0]+'_flat.fits'
            print 'Merging',flatname
            if f == 0:
                merged = Table.read(os.path.join(out_dir,flatname),format='fits')
            else:
                t = Table.read(os.path.join(out_dir,flatname),format='fits')
                # merged = join(t,merged,keys='NUMBER') # Original
                # There are multiple cols that are identical
                # (NUMBER, RA, DEC, etc.), so I think it is best to use all
                merged = join(t, merged)

    #--------------------------------------------------------------------------------

    # Can make merged table directly instead
    else:
        merged = flatten(data_dir, filelist, out_dir, tilename, save_all=False)

    if tilename == '':
        newfile = os.path.join(out_dir, 'merged.fits')
        merged.write(newfile, format='fits', overwrite=True)
        print 'Wrote', os.path.join(out_dir, 'merged.fits')
    else:
        newfile = os.path.join(out_dir, tilename+'_merged.fits')
        merged.write(newfile, format='fits', overwrite=True)
        print 'Wrote', os.path.join(out_dir, tilename+'_merged.fits')

    return

def flatten(data_dir, filelist, out_dir, tilename, save_all=False):

    for fid, data_file in enumerate(filelist):

        print 'Flattening',data_file
        data_hdu = fits.open(os.path.join(data_dir,data_file))
        data_tab = data_hdu[1].data

        defs = {}
        descr = data_tab.view(np.ndarray).dtype.descr
        alldefs = get_coldefs(descr,defs,bands,band_cols,bidim_band_cols)
        names = [d[0] for d in alldefs]
        formats = [d[2] for d in alldefs]
        cols = []
        prev_colname = ''
        prev_colname_add = '' 
        prev_colname_band = ''
        prev_colname_multi = ''
        prev_colname_multi_add = ''
        prev_colname_bidim = ''
        prev_colname_bidim_band = ''
        k = 0
        m = 0

        mofsofstr = ''
        if '-mof' in data_file:
            mofsofstr = 'MOF_'
        elif '-sof' in data_file:
            mofsofstr = 'SOF_'
        else:
            mofsofstr = ''

        for i in xrange(len(names)):
            nm = names[i].split('_')[len(names[i].split('_'))-1]
            strip = len(nm)
            if printout:
                print 'Checking',names[i],strip
            check_cols,colname = check_name(names[i],onedim_cols,prev_colname,0)
            check_band_cols,colname_band = check_name(names[i],band_cols,prev_colname_band,strip+1)
            check_bidim_cols,colname_bidim = check_name(names[i],bidim_cols,prev_colname_bidim,strip+1)
            check_bidim_band_cols,colname_bidim_band = check_name(names[i],bidim_band_cols,prev_colname_bidim_band,4)
            check_multi_cols,colname_multi = check_name(names[i],multi_cols,prev_colname_multi,strip+1)
            check_cols_add,colname_add = check_name(names[i],onedim_cols_addband,prev_colname_add,0)
            check_multi_cols_add,colname_multi_add = check_name(names[i],multi_cols_add,prev_colname_multi_add,strip+1)
            if printout:
                print check_cols,check_band_cols,check_bidim_cols,check_bidim_band_cols,check_multi_cols,check_cols_add,check_multi_cols_add
            if i == 0:
                n = 0
                m = 0
            if i > 0 and (prev_colname != colname or
                          colname_band != prev_colname_band or
                          colname_bidim != prev_colname_bidim or
                          colname_bidim_band != prev_colname_bidim_band or
                          colname_multi != prev_colname_multi):
                n = 0
                m = 0
                k = k+1
            if check_band_cols == True:
                cols.append(fits.Column(name=mofsofstr+names[i].upper(),format=get_fits_type(names[i]),array=data_tab[colname_band][:,n]))
                n = n + 1
                prev_colname_band = colname_band
            elif check_bidim_cols == True:
                cols.append(fits.Column(name=mofsofstr+names[i].upper(),format=get_fits_type(names[i]),array=data_tab[colname_bidim][:,n,m]))
                if n == len(data_tab[colname_bidim][0])-1:
                    n = 0
                    m = m + 1
                else:
                    n = n + 1
                prev_colname_bidim = colname_bidim
            elif check_bidim_band_cols == True:
                cols.append(fits.Column(name=mofsofstr+names[i].upper(),format=get_fits_type(names[i]),array=data_tab[colname_bidim_band][:,n,m]))
                if n == len(data_tab[colname_bidim_band][0])-1:
                    n = 0
                    m = m + 1
                else:
                    n = n + 1
                prev_colname_bidim_band = colname_bidim_band
            elif check_multi_cols == True:
                cols.append(fits.Column(name=mofsofstr+names[i].upper(),format=get_fits_type(names[i]),array=data_tab[colname_multi][:,n]))
                if n == len(data_tab[colname_multi][0])-1:
                    n = 0
                else:
                    n = n + 1
                prev_colname_multi = colname_multi
            elif check_multi_cols_add == True:
                idx = data_file.find('_cat')
                fileband = data_file[idx-1:idx]
                #print data_tab[colname_multi_add][:,n],n
                #print names[i]+'_'+fileband.upper(),colname_multi_add
                cols.append(fits.Column(name=mofsofstr+names[i].upper()+'_'+fileband.upper(),format=get_fits_type(names[i]),array=data_tab[colname_multi_add][:,n]))
                if n == len(data_tab[colname_multi_add][0])-1:
                    n = 0
                else:
                    n = n + 1
                prev_colname_multi_add = colname_multi_add
            elif check_cols_add == True:
                idx = data_file.find('_cat')
                fileband = data_file[idx-1:idx]
                newname = names[i].upper()+'_'+fileband.upper()
                cols.append(fits.Column(name=mofsofstr+newname,format=get_fits_type(names[i]),array=data_tab[names[i]]))
                prev_colname_add = colname_add
            else:
                if names[i] == "id":
                    newname = "COADD_OBJECT_ID"
                elif names[i] == "number":
                    newname = "NUMBER"
                elif names[i] == "ra":
                    newname = "RA"
                elif names[i] == "dec":
                    newname = "DEC"
                else:
                    newname = mofsofstr+names[i].upper()

                cols.append(fits.Column(name=newname,format=get_fits_type(names[i]),array=data_tab[names[i]]))
                prev_colname = colname

        if "-mof" in data_file or "-sof" in data_file:
            for b,band in enumerate(bands):
                psf_mag_err = 1.0857362*data_tab['psf_flux_err'][:,b]/data_tab['psf_flux'][:,b]
                psf_mag_err[psf_mag_err == 1.0857362] = -9999
                cm_mag_err = 1.0857362*np.sqrt(data_tab['cm_flux_cov'][:,b,b])/data_tab['cm_flux'][:,b]
                cm_mag_err[cm_mag_err == 1.0857362] = -9999
                cm_mag_err = [-9999 if np.isnan(x) else x for x in cm_mag_err]
    	        cols.append(fits.Column(name=mofsofstr+'psf_mag_err_'+band,format='D',array=psf_mag_err))
                cols.append(fits.Column(name=mofsofstr+'cm_mag_err_'+band,format='D',array=cm_mag_err))

        if fid == 0 and tilename != '':
            cols.append(fits.Column(name='TILENAME',format='27A',array=data_tab.size*[tilename]))

        new_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
        new_tbdata = new_hdu.data

        # Original Nacho code
        if save_all is True:
            new_hdu = fits.BinTableHDU(data=new_tbdata)
            out_file = os.path.join(out_dir,os.path.splitext(data_file)[0]+'_flat.fits')
            new_hdu.writeto(out_file,clobber='True')
            print 'Wrote',out_file

        # Not saving individual files, less IO intensive:
        else:
            if fid == 0:
                merged = Table(new_tbdata)
            else:
                merged = join(Table(new_tbdata), merged)
            print 'Merged',data_file

    if save_all is False:
        return merged

def main():

    # Parse command line
    args = parser.parse_args()

    if not os.path.isdir(args.data_dir):
        print 'Path with data not found at',args.data_dir
        sys.exit(1)

    if args.out_dir is None:
        args.out_dir = args.data_dir
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    print 'Reading data from',args.data_dir
    print 'Writing files to upload to',args.out_dir

    print 'Getting list of files...'
    data_files = [f for f in os.listdir(args.data_dir) if os.path.isfile(os.path.join(args.data_dir, f))]
    check_string = ['g_cat','r_cat','i_cat','z_cat','-mof','-sof']

    chk = 0
    select_files = []
    tilename = ''
    for check in check_string:
        for data_file in data_files:
            if tilename == '':
                idx = data_file.find('DES')
                tilename = data_file[idx:idx+12]
            if check in data_file:
                chk += 1
                select_files.append(data_file)
                print 'Appending',data_file
                break
    if chk != len(check_string):
        print 'Missing file of one of these types:',check_string
        sys.exit()

    # New setup flattens during merge
    merge(args, tilename, select_files)

if __name__ == '__main__':
    sys.exit(main())


