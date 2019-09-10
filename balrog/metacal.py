import numpy as np
import fitsio
import h5py
import glob
import sys

def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

def convert_mcal_to_h5(catdir, outfile, bands, version_tag=None, match_type='default', vb=False):
    """
    Converts metacal fits files into a single h5 file with separate tables for each of the unsheared and four sheared values.
    This form is much faster to access, but doesn't preserve the complete recarray table structure. Instead, each column is a different data array in the h5 file. The output of this function is sorted by coadd id.

    This is a modified version of Troxel's convert_mcal_to_h5 to work with Balrog outputs
    See original at https://github.com/des-science/y3_matched_cat/blob/master/make.py
    """

    # Some mcal col definitions
    # TODO: Wrap this up in config classes
    gold_cols_default = {'FLAGS_GOLD':'flags_gold',
                        'EXTENDED_CLASS_MOF':'extended_class_mof',
                        'EXTENDED_CLASS_SOF':'extended_class_sof',
                        'TILENAME':'tilename'
                        }

    gold_cols_mof = {'FLAGS_GOLD_MOF_ONLY':'flags_gold',
                    'EXTENDED_CLASS_MOF':'extended_class_mof',
                    'TILENAME':'tilename'
                    }

    gold_cols_sof= {'FLAGS_GOLD_SOF_ONLY':'flags_gold',
                    'EXTENDED_CLASS_SOF':'extended_class_sof',
                    'TILENAME':'tilename'
                    }

    mcal_flat = {'id':'coadd_object_id',
                'bal_id':'bal_id',
                'flags':'flags',
                'mask_frac':'mask_frac',
                'ra':'ra',
                'dec':'dec',
                'psfrec_T':'psf_T',
                'mcal_Tpsf':'mcal_psf_T'
                }

    mcal_flats = {'mcal_T_r':'T',
                'mcal_T_err':'T_err',
                'mcal_s2n_r':'snr'
                }

    mcal_vec2 = {'psfrec_g':'psf_e',
                'mcal_gpsf':'mcal_psf_e'
                }
    mcal_vec2_ext = ['1','2']

    mcal_vec3 = {'nimage_tot':'nimage_tot_',
                'nimage_use':'nimage_use_'
                }
    mcal_vec3_ext = ['r','i','z']

    mcal_vec4 = {'nimage_tot':'nimage_tot_',
                'nimage_use':'nimage_use_'
                }
    mcal_vec4_ext = ['g','r','i','z']

    # Determine if using vec3 or vec4 (riz or griz)
    if bands == 'riz':
        mcal_vecb = mcal_vec3
        mcal_vecb_ext = mcal_vec3_ext
        mcal_bindx = dict(zip('riz', range(3)))
        mcal_pars_bindx = dict(zip('riz', [5,6,7]))
    elif bands == 'griz':
        mcal_vecb = mcal_vec4
        mcal_vecb_ext = mcal_vec4_ext
        mcal_bindx = dict(zip('griz', range(4)))
        mcal_pars_bindx = dict(zip('griz', [5,6,7,8]))
    else:
        raise ValueError('Only `riz` and `griz` currently allowed for arg `bands`!')

    if match_type == 'default':
        mcal_flat = merge_two_dicts(mcal_flat, gold_cols_default)
    elif match_type == 'mof_only':
        mcal_flat = merge_two_dicts(mcal_flat, gold_cols_gold)
    elif match_type == 'sof_only':
        mcal_flat = merge_two_dicts(mcal_flat, gold_cols_sof)
    else:
        raise ValueError('Not a valid match_type!')

    if version_tag is None:
        version_tag = ''

    catfiles = glob.glob('{}/*{}*.fits'.format(catdir, version_tag))

    total_length = 0
    for catfile in catfiles:
        try:
            h = fitsio.read_header(catfile, ext=1)
            total_length += h['NAXIS2']
        except IOError as e:
            print('{}\nSkipping file.'.format(e))

    print('Total stack length will be {}'.format(total_length))
    chunk_size = 1000000
    if total_length < chunk_size:
        # Can't have a chunk size greater than length
        chunk_size = total_length

    print('Writing stack to {}'.format(outfile))
    f = h5py.File(outfile, 'w')
    iter_end = 0 # Total length of data at each iteration step.
    for i, fname in enumerate(sorted(catfiles)): # Loop over fits files

        if version_tag not in fname: # Ditch extra files in dir that don't match version_tag
            continue

        print('Converting {} ({} of {})...'.format(fname, i+1, len(catfiles)))
        cat  = fitsio.read(fname) # Read fits file into np recarray
        # mask = np.searchsorted(gold,cat['id']) # Remove things that aren't in gold
        # cat = cat[ (mask>0) & (mask<len(gold)) ]
        lencat   = len(cat) # Get length after culling of cat array

        if vb is True:
            print('Writing unsheared cols...')

        for name in mcal_flat: # Iterate over columns
            if i==0: # If first fits file, create the data sets in the h5 file
                f.create_dataset( 'catalog/unsheared/'+mcal_flat[name], maxshape=(total_length,), shape=(total_length,), dtype=cat[name].dtype, chunks=(chunk_size,) )
            # Write the cut cat array into h5 file
            f['catalog/unsheared/'+mcal_flat[name]][iter_end:iter_end+lencat] = cat[name]

        for name in mcal_vec2: # Iterate over columns
            for j in range(2):
                if i==0: # If first fits file, create the data sets in the h5 file
                    f.create_dataset( 'catalog/unsheared/'+mcal_vec2[name]+mcal_vec2_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat[name].dtype, chunks=(chunk_size,) )
                # Write the cut cat array into h5 file
                f['catalog/unsheared/'+mcal_vec2[name]+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat[name][:,j]

        for name in mcal_vecb: # Iterate over columns
            for j in range(len(bands)):
                if i==0: # If first fits file, create the data sets in the h5 file
                    f.create_dataset( 'catalog/unsheared/'+mcal_vecb[name]+mcal_vecb_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat[name].dtype, chunks=(chunk_size,) )
                # Write the cut cat array into h5 file
                f['catalog/unsheared/'+mcal_vecb[name]+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat[name][:,j]

        # Create new column size_ratio = size/psf_size
        if i==0:
            f.create_dataset( 'catalog/unsheared/size_ratio', maxshape=(total_length,), shape=(total_length,), dtype=cat['mcal_T_r'].dtype, chunks=(chunk_size,) )
        f['catalog/unsheared/size_ratio'][iter_end:iter_end+lencat] = (cat['mcal_T_r']/cat['psfrec_T'])

        for name in mcal_flats: # Iterate over columns
            if i==0: # If first fits file, create the data sets in the h5 file
                f.create_dataset( 'catalog/unsheared/'+mcal_flats[name], maxshape=(total_length,), shape=(total_length,), dtype=cat[name].dtype, chunks=(chunk_size,) )
            # Write the cut cat array into h5 file
            f['catalog/unsheared/'+mcal_flats[name]][iter_end:iter_end+lencat] = cat[name]

        # Do some more complicated cases explicitly
        for j in range(2):
            if i==0:
                f.create_dataset( 'catalog/unsheared/e_'+mcal_vec2_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat['mcal_g'].dtype, chunks=(chunk_size,) )
                f.create_dataset( 'catalog/unsheared/covmat_'+mcal_vec2_ext[j]+'_'+mcal_vec2_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat['mcal_g_cov'].dtype, chunks=(chunk_size,) )
                if j==0:
                    f.create_dataset( 'catalog/unsheared/covmat_0_1', maxshape=(total_length,), shape=(total_length,), dtype=cat['mcal_g_cov'].dtype, chunks=(chunk_size,) )
            f['catalog/unsheared/e_'+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat['mcal_g'][:,j]
            f['catalog/unsheared/covmat_'+mcal_vec2_ext[j]+'_'+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat['mcal_g_cov'][:,j,j]
            if j==0:
                f['catalog/unsheared/covmat_0_1'][iter_end:iter_end+lencat] = cat['mcal_g_cov'][:,0,1]

        # Do some more complicated cases explicitly
        for j in range(len(bands)):
            if i==0:
                # mcal -> gauss for dtype as columns are not natively created for mcal (but share dtype)
                f.create_dataset( 'catalog/unsheared/flux_'+mcal_vecb_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat['gauss_flux'].dtype, chunks=(chunk_size,) )
                f.create_dataset( 'catalog/unsheared/flux_err_'+mcal_vecb_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat['gauss_flux_cov'].dtype, chunks=(chunk_size,) )
            f['catalog/unsheared/flux_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars'][:,5+j]
            f['catalog/unsheared/flux_err_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars_cov'][:,5+j,5+j]

        for sheared in ['_1p','_1m','_2p','_2m']:
            if vb is True:
                print('Writing sheared{} cols...'.format(sheared))

            if i==0:
                f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/size_ratio', maxshape=(total_length,), shape=(total_length,), dtype=cat['mcal_T_r'].dtype, chunks=(chunk_size,) )
            f['catalog/sheared_'+sheared[1:]+'/size_ratio'][iter_end:iter_end+lencat] = (cat['mcal_T_r'+sheared]/cat['psfrec_T'])

            for name in mcal_flats: # Iterate over columns
                if i==0: # If first fits file, create the data sets in the h5 file
                    f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/'+mcal_flats[name], maxshape=(total_length,), shape=(total_length,), dtype=cat[name].dtype, chunks=(chunk_size,) )
                # Write the cut cat array into h5 file
                f['catalog/sheared_'+sheared[1:]+'/'+mcal_flats[name]][iter_end:iter_end+lencat] = cat[name+sheared]

            # Do some more complicated cases explicitly
            for j in range(2):
                if i==0:
                    f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/e'+mcal_vec2_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat['mcal_g'].dtype, chunks=(chunk_size,) )
                    f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/covmat_'+mcal_vec2_ext[j]+'_'+mcal_vec2_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat['mcal_g_cov'].dtype, chunks=(chunk_size,) )
                    if j==0:
                        f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/covmat_0_1', maxshape=(total_length,), shape=(total_length,), dtype=cat['mcal_g_cov'].dtype, chunks=(chunk_size,) )
                f['catalog/sheared_'+sheared[1:]+'/e'+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat['mcal_g'+sheared][:,j]
                f['catalog/sheared_'+sheared[1:]+'/covmat_'+mcal_vec2_ext[j]+'_'+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat['mcal_g_cov'+sheared][:,j,j]
                if j==0:
                    f['catalog/sheared_'+sheared[1:]+'/covmat_0_1'][iter_end:iter_end+lencat] = cat['mcal_g_cov'+sheared][:,0,1]

            # Do some more complicated cases explicitly
            for j in range(len(bands)):
                if i==0:
                    f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/flux_'+mcal_vecb_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat['gauss_flux'].dtype, chunks=(chunk_size,) )
                    f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/flux_err_'+mcal_vecb_ext[j], maxshape=(total_length,), shape=(total_length,), dtype=cat['gauss_flux_cov'].dtype, chunks=(chunk_size,) )
                f['catalog/sheared_'+sheared[1:]+'/flux_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars'+sheared][:,5+j]
                f['catalog/sheared_'+sheared[1:]+'/flux_err_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars_cov'+sheared][:,5+j,5+j]

        iter_end += lencat

    # Close h5 file to dump cache and clear memory used in cat
    # f.close()
    cat = None

    # f2           = h5py.File(goldfile, 'r')
    # gold         = f2['catalog']['gold']['coadd_object_id'][:]
    # gold         = np.sort(gold)
    # f2.close()

    # # Reopen h5 file to do some sorting
    # f = h5py.File(outfile, 'r+')
    # # kill missing objects relative to gold
    # coadd = f['catalog']['unsheared']['coadd_object_id'][:]
    # mask = np.in1d(coadd,gold,assume_unique=False)
    # coadd = coadd[mask]

    # sort against gold
    # s0 = np.argsort(coadd)
    # coadd = coadd[s0]
    # s = np.where(np.in1d(gold,coadd))[0]
    # # Loop over tables and data sets to reorder h5 columns to match gold coadd id ordering. At this point gold (and now this file) are sorted by coadd id
    # for table in ['unsheared','sheared_1p','sheared_1m','sheared_2p','sheared_2m']:
    #     for col in f['catalog'][table].keys():
    #         print table,col
    #         x = f['catalog'][table][col][:][mask][s0]
    #         c = np.zeros(len(gold),dtype=x.dtype)
    #         c[s] = x # Insert objects into dummy array at correct position to match gold
    #         del(f['catalog'][table][col])
    #         f.create_dataset( 'catalog/'+table+'/'+col, maxshape=(len(gold),), shape=(len(gold),), dtype=c.dtype, chunks=(chunk_size,) )
    #         f['catalog'][table][col][:] = c
    #         # This is complicated, but necessary because the shape catalog has fewer objects than gold in total, but also has objects that gold doesn't (and duplicates in those extra objects). 

    # # Add bit to flag to represent objects not in gold
    # c = f['catalog']['unsheared']['flags'][:]
    # c[f['catalog']['unsheared']['coadd_object_id'][:]==0] += 2**28
    # f['catalog']['unsheared']['flags'][:] = c

    # NOTE: We don't apply calibrations to the fluxes in Balrog as we can't compute them.
    # But this is where flux_err is properly square rooted in the data, so we will do it
    # here as well
    # Calibrate fluxes
    # f2 = h5py.File(goldfile, 'r')
    # s  = np.argsort(f2['catalog']['gold']['coadd_object_id'][:])
    for band in bands:
        # corr = 10**(-0.4*(f2['catalog/gold/delta_mag_chrom_'+band][:]+f2['catalog/gold/delta_mag_y4_'+band][:]-f2['catalog/gold/a_sed_sfd98_'+band][:]))[s]
        for table in ['unsheared','sheared_1p','sheared_1m','sheared_2p','sheared_2m']:
            # c = f['catalog/'+table+'/flux_'+band][:]
            # f['catalog/'+table+'/flux_'+band][:]  = c*corr
            c = f['catalog/'+table+'/flux_err_'+band][:]
            f['catalog/'+table+'/flux_err_'+band][:]  = np.sqrt(c)

    # f2.close()
    # Create Rxx columns
    if vb is True:
        print('Calculating shear responses...')
    f.create_dataset( 'catalog/unsheared/R11', maxshape=(total_length,), shape=(total_length,), dtype=float, chunks=(chunk_size,) )
    f.create_dataset( 'catalog/unsheared/R12', maxshape=(total_length,), shape=(total_length,), dtype=float, chunks=(chunk_size,) )
    f.create_dataset( 'catalog/unsheared/R21', maxshape=(total_length,), shape=(total_length,), dtype=float, chunks=(chunk_size,) )
    f.create_dataset( 'catalog/unsheared/R22', maxshape=(total_length,), shape=(total_length,), dtype=float, chunks=(chunk_size,) )
    f['catalog/unsheared/R11'][:] = (f['catalog/sheared_1p/e1'][:] - f['catalog/sheared_1m/e1'][:])/0.02
    f['catalog/unsheared/R12'][:] = (f['catalog/sheared_2p/e1'][:] - f['catalog/sheared_2m/e1'][:])/0.02
    f['catalog/unsheared/R21'][:] = (f['catalog/sheared_1p/e2'][:] - f['catalog/sheared_1m/e2'][:])/0.02
    f['catalog/unsheared/R22'][:] = (f['catalog/sheared_2p/e2'][:] - f['catalog/sheared_2m/e2'][:])/0.02

    # Close h5 file and dump cache
    f.close()

    return
