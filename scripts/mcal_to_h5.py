import h5py
import fitsio

mcal_flat = {'id':'coadd_object_id',
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

def mcal_to_h5(mcal_file, h5_filename, bands, max_shape=450000000, vb=False):

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

    # Create the shape h5 file
    f = h5py.File(h5_filename, 'w')
    iter_end = 0 # Total length of data at each iteration step.

    cat  = fitsio.read(mcal_file) # Read fits file into np recarray
    # mask = np.searchsorted(gold,cat['id']) # Remove things that aren't in gold
    # cat = cat[ (mask>0) & (mask<len(gold)) ]
    lencat   = len(cat) # Get length after culling of cat array

    if vb is True:
        print('Writing unsheared cols...')
    for name in mcal_flat: # Iterate over columns
        f.create_dataset( 'catalog/unsheared/'+mcal_flat[name], maxshape=(max_shape,), shape=(lencat,), dtype=cat[name].dtype, chunks=(1000000,) )
        # Write the cut cat array into h5 file
        f['catalog/unsheared/'+mcal_flat[name]][iter_end:iter_end+lencat] = cat[name]

    for name in mcal_vec2: # Iterate over columns
        for j in range(2):
            f.create_dataset( 'catalog/unsheared/'+mcal_vec2[name]+mcal_vec2_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat[name].dtype, chunks=(1000000,) )
            # Write the cut cat array into h5 file
            f['catalog/unsheared/'+mcal_vec2[name]+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat[name][:,j]

    for name in mcal_vecb: # Iterate over columns
        for j in range(len(bands)):
            f.create_dataset( 'catalog/unsheared/'+mcal_vecb[name]+mcal_vecb_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat[name].dtype, chunks=(1000000,) )
            # Write the cut cat array into h5 file
            f['catalog/unsheared/'+mcal_vecb[name]+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat[name][:,j]

    # Create new column size_ratio = size/psf_size
    f.create_dataset( 'catalog/unsheared/size_ratio', maxshape=(max_shape,), shape=(lencat,), dtype=cat['mcal_T_r'].dtype, chunks=(1000000,) )
    f['catalog/unsheared/size_ratio'][iter_end:iter_end+lencat] = (cat['mcal_T_r']/cat['psfrec_T'])

    for name in mcal_flats: # Iterate over columns
        f.create_dataset( 'catalog/unsheared/'+mcal_flats[name], maxshape=(max_shape,), shape=(lencat,), dtype=cat[name].dtype, chunks=(1000000,) )
        # Write the cut cat array into h5 file
        f['catalog/unsheared/'+mcal_flats[name]][iter_end:iter_end+lencat] = cat[name]

    # Do some more complicated cases explicitly
    for j in range(2):
        f.create_dataset( 'catalog/unsheared/e_'+mcal_vec2_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat['mcal_g'].dtype, chunks=(1000000,) )
        f.create_dataset( 'catalog/unsheared/covmat_'+mcal_vec2_ext[j]+'_'+mcal_vec2_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat['mcal_g_cov'].dtype, chunks=(1000000,) )
        if j==0:
            f.create_dataset( 'catalog/unsheared/covmat_0_1', maxshape=(max_shape,), shape=(lencat,), dtype=cat['mcal_g_cov'].dtype, chunks=(1000000,) )
        f['catalog/unsheared/e_'+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat['mcal_g'][:,j]
        f['catalog/unsheared/covmat_'+mcal_vec2_ext[j]+'_'+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat['mcal_g_cov'][:,j,j]
        if j==0:
            f['catalog/unsheared/covmat_0_1'][iter_end:iter_end+lencat] = cat['mcal_g_cov'][:,0,1]

    # Do some more complicated cases explicitly
    for j in range(len(bands)):
        # mcal -> gauss for dtype as columns are not natively created for mcal
        f.create_dataset( 'catalog/unsheared/flux_'+mcal_vecb_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat['gauss_flux'].dtype, chunks=(1000000,) )
        f.create_dataset( 'catalog/unsheared/flux_err_'+mcal_vecb_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat['gauss_flux'].dtype, chunks=(1000000,) )
        # Map j to band to pars index
        # pars_indx = mcal_pars_bindx[mcal_vecb_ext[j]]
        # f['catalog/unsheared/flux_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars'][:,pars_indx]
        # f['catalog/unsheared/flux_err_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars_cov'][:,pars_indx,pars_indx]
        f['catalog/unsheared/flux_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars'][:,5+j]
        f['catalog/unsheared/flux_err_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars_cov'][:,5+j,5+j]

    for sheared in ['_1p','_1m','_2p','_2m']:
        if vb is True:
            print('Writing sheared{} cols...'.format(sheared))

        f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/size_ratio', maxshape=(max_shape,), shape=(lencat,), dtype=cat['mcal_T_r'].dtype, chunks=(1000000,) )
        f['catalog/sheared_'+sheared[1:]+'/size_ratio'][iter_end:iter_end+lencat] = (cat['mcal_T_r'+sheared]/cat['psfrec_T'])

        for name in mcal_flats: # Iterate over columns
            f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/'+mcal_flats[name], maxshape=(max_shape,), shape=(lencat,), dtype=cat[name].dtype, chunks=(1000000,) )
            # Write the cut cat array into h5 file
            f['catalog/sheared_'+sheared[1:]+'/'+mcal_flats[name]][iter_end:iter_end+lencat] = cat[name+sheared]
        # Do some more complicated cases explicitly

        for j in range(2):
            f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/e'+mcal_vec2_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat['mcal_g'].dtype, chunks=(1000000,) )
            f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/covmat_'+mcal_vec2_ext[j]+'_'+mcal_vec2_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat['mcal_g_cov'].dtype, chunks=(1000000,) )
            if j==0:
                f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/covmat_0_1', maxshape=(max_shape,), shape=(lencat,), dtype=cat['mcal_g_cov'].dtype, chunks=(1000000,) )
            f['catalog/sheared_'+sheared[1:]+'/e'+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat['mcal_g'+sheared][:,j]
            f['catalog/sheared_'+sheared[1:]+'/covmat_'+mcal_vec2_ext[j]+'_'+mcal_vec2_ext[j]][iter_end:iter_end+lencat] = cat['mcal_g_cov'+sheared][:,j,j]
            if j==0:
                f['catalog/sheared_'+sheared[1:]+'/covmat_0_1'][iter_end:iter_end+lencat] = cat['mcal_g_cov'+sheared][:,0,1]

        # Do some more complicated cases explicitly
        for j in range(len(bands)):
            f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/flux_'+mcal_vecb_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat['gauss_flux'].dtype, chunks=(1000000,) )
            f.create_dataset( 'catalog/sheared_'+sheared[1:]+'/flux_err_'+mcal_vecb_ext[j], maxshape=(max_shape,), shape=(lencat,), dtype=cat['gauss_flux_cov'].dtype, chunks=(1000000,) )
            f['catalog/sheared_'+sheared[1:]+'/flux_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars'+sheared][:,5+j]
            f['catalog/sheared_'+sheared[1:]+'/flux_err_'+mcal_vecb_ext[j]][iter_end:iter_end+lencat] = cat['mcal_pars_cov'+sheared][:,5+j,5+j]

        # Only used for multiple files! Will bug here
        # iter_end += lencat

    if vb is True:
        print('Calculating shear responses...')
    # Create Rxx columns
    f.create_dataset( 'catalog/unsheared/R11', maxshape=(max_shape,), shape=(lencat,), dtype=float, chunks=(1000000,) )
    f.create_dataset( 'catalog/unsheared/R12', maxshape=(max_shape,), shape=(lencat,), dtype=float, chunks=(1000000,) )
    f.create_dataset( 'catalog/unsheared/R21', maxshape=(max_shape,), shape=(lencat,), dtype=float, chunks=(1000000,) )
    f.create_dataset( 'catalog/unsheared/R22', maxshape=(max_shape,), shape=(lencat,), dtype=float, chunks=(1000000,) )
    f['catalog/unsheared/R11'][:] = (f['catalog/sheared_1p/e1'][:] - f['catalog/sheared_1m/e1'][:])/0.02
    f['catalog/unsheared/R12'][:] = (f['catalog/sheared_2p/e1'][:] - f['catalog/sheared_2m/e1'][:])/0.02
    f['catalog/unsheared/R21'][:] = (f['catalog/sheared_1p/e2'][:] - f['catalog/sheared_1m/e2'][:])/0.02
    f['catalog/unsheared/R22'][:] = (f['catalog/sheared_2p/e2'][:] - f['catalog/sheared_2m/e2'][:])/0.02

    # Close h5 file to dump cache and clear memory used in cat
    f.close()

    return
