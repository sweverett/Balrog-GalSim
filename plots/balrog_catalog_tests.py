#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import fitsio
import treecorr
import ngmix
import os
import errno
from esutil import htm
# import lfunc

# Can use for debugging
import pudb

def mag_flux_tick_function(flux):
    return ["%.1f" % mag for mag in  30 - 2.5*np.log10(flux)]

def flux_mag_convert(mag, zeropt = 30.0):
    return 10**((30 - mag)/2.5)

def match_catalogs(truth_catalog = None, meas_catalog = None,ratag_truth = 'ra', dectag_truth = 'dec',\
                       ratag_meas = 'ra',dectag_meas = 'dec', match_radius = 1./3600):
    matcher = htm.Matcher(depth=14,ra=truth_catalog[ratag_truth], dec = truth_catalog[dectag_truth])
    id_m, id_t, dist = matcher.match(ra = meas_catalog[ratag_meas], dec = meas_catalog[dectag_meas],radius=match_radius)
    truth_cut = truth_catalog[id_t]
    meas_cut = meas_catalog[id_m]
    return truth_cut, meas_cut, dist

def remove_close_pairs(catalog,ratag='ra',dectag='dec',radius=5./3600):
    matcher = htm.Matcher(depth=14,ra=catalog[ratag],dec=catalog[dectag])
    ind1,ind2,dist = matcher.match(ra=catalog[ratag],dec=catalog[dectag],radius=radius,maxmatch=0)
    nonself = dist>0
    ind1 = ind1[nonself]
    ind2= ind2[nonself]
    all_inds = np.arange(catalog.size)
    keep = np.in1d(all_inds,ind1,invert=True)
    return catalog[keep]

def quality_cuts(catalog, band = 'i'):
    bands = ['g','r','i','z','Y']
    bandind = np.arange(len(bands))[np.in1d(bands,band)][0]

    keep = (catalog['flags'] == 0) & (catalog['cm_s2n_r'] > 10) & (catalog['cm_T']/catalog['cm_T_err'] > .5) & (catalog['cm_T']/catalog['psfrec_T'] > 0.5)
    catalog = catalog[keep]
    return catalog

def get_catalogs(path = '.', tilename = 'DES0347-5540', re_id = '0', minsep = 0., stars=False, ratag = 'ra', dectag = 'dec'):
    if type(re_id) is not type(''):
        re_id = str(re_id)
    if stars is False:
        cat = tilename+'_'+re_id+'_balrog_truth_cat_gals.fits'
        filename = os.path.join(path, re_id, tilename, cat)
        truth_catalog = fitsio.read(filename)
    else:
        cat = tilename+'_'+re_id+'_balrog_truth_cat_stars.fits'
        filename = os.path.join(path, re_id, tilename, cat)
        truth_catalog = fitsio.read(filename)
    if minsep > 0.:
        truth_catalog = remove_close_pairs(truth_catalog,radius=minsep,ratag=ratag,dectag=dectag)

    cat = tilename + '_mof.fits'
    filename = os.path.join(path, re_id, tilename, 'mof', cat)
    meas_catalog = fitsio.read(filename)
    meas_catalog = quality_cuts(meas_catalog)
    truth_matched, meas_gal_matched, dist = match_catalogs(truth_catalog = truth_catalog, meas_catalog =meas_catalog,ratag_truth = ratag, dectag_truth = dectag)

    return truth_catalog, truth_matched, meas_gal_matched, dist

def make_plots(truth = None, meas_matched = None, truth_matched = None, sep = None, filetag = '',
               bandind = 1, outdir=None):

    # TODO: fontsize as input
    matplotlib.rcParams.update({'font.size': 16})

    if outdir is None:
        outdir = os.getcwd()
    else:
        try:
            os.makedirs(outdir)
        except OSError as e:
            # Ignore error if dir already exists 
            if e.errno != errno.EEXIST:
                raise e

    # flux, size:
    tag = 'cm_flux'
    errtag = 'cm_flux_cov'

    bands = ['g','r','i','z']
    mag_tick_locations = np.array([0,5,10,15,20,25,30,35,40])
    flux_tick_locations = flux_mag_convert(mag_tick_locations)

    err = np.sqrt(meas_matched['cm_flux_cov'][:,bandind,bandind])
    magerr = err / meas_matched['cm_flux'][:,bandind]

    dy = [-4, 2]
    a = 0.5
    fig,(ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(21,7))
    ax1.plot(truth_matched['cm_mag'][:,bandind],meas_matched['cm_mag'][:,bandind] - truth_matched['cm_mag'][:,bandind],'.k',alpha=a)
    ax1.axhline(0,color='red',linestyle='--',alpha=0.5)
    ax1.set_ylim(dy[0],dy[1])
    #ax1.set_yscale('symlog')
    ax1.set_xlim(16,26)
    ax1b = ax1.twiny()
    ax1b.set_xlim(ax1.get_xlim())
    ax1b.set_xticks(mag_tick_locations)
    ax1b.set_xticklabels(np.log10(flux_tick_locations))
    ax1b.set_xlabel('log_10 flux')
    ax1.set_xlabel('magnitude')
    ax1.set_xlabel('cm_mag ('+bands[bandind]+')')
    ax1.set_ylabel('magnitude difference ('+bands[bandind]+')')


    ax2.plot(truth_matched['cm_mag'][:,bandind],(meas_matched[tag][:,bandind] - truth_matched[tag][:,bandind])/np.sqrt(meas_matched[errtag][:,bandind,bandind]),'.k',alpha=a)
    #ax2.set_ylim(-250,250)
    ax2.set_xlim(16,26)
    ax2.axhline(0,color='red',linestyle='--',alpha=0.5)
    ax2.set_xlabel('cm_mag ('+bands[bandind]+')')
    ax2.set_ylabel('flux chi ('+bands[bandind]+')')
    ax3.hist((meas_matched[tag][:,bandind] - truth_matched[tag][:,bandind])/np.sqrt(meas_matched[errtag][:,bandind,bandind]),bins=np.linspace(-25,25,250),color='k')
    ax3.axvline(0,color='red',linestyle='--',alpha=a)
    ax3.set_xlabel('flux chi ('+bands[bandind]+')')

    outfile = os.path.join(outdir, 'cm_flux-'+filetag+bands[bandind])
    fig.savefig(outfile)
    plt.close(fig)

    # Now make a detection efficiency plot.
    size_bin_edges = np.linspace(0,10,26)
    size_bin_centers = (size_bin_edges[1:] + size_bin_edges[0:-1])/2.

    flux_bin_edges = np.linspace(15,28,11)
    flux_bin_centers = (flux_bin_edges[1:] + flux_bin_edges[0:-1])/2.

    nobj_truth_meas_mag,_ = np.histogram(truth_matched['cm_mag'],bins=flux_bin_edges)
    nobj_truth_input_mag,_= np.histogram(truth['cm_mag'],bins=flux_bin_edges)

    nobj_truth_meas_size,_ = np.histogram(np.log10(truth_matched['cm_T']),bins=size_bin_edges)
    nobj_truth_input_size,_= np.histogram(np.log10(truth['cm_T']),bins=size_bin_edges)

    fig,(ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(21,7))
    ax1.plot(flux_bin_centers, nobj_truth_meas_mag*1./nobj_truth_input_mag)
    ax1.set_ylabel('magnitude completeness')
    ax1.set_xlabel('magnitude')
    ax2.plot(size_bin_centers, nobj_truth_meas_size*1./nobj_truth_input_size)
    ax2.set_ylabel('size completeness')
    ax2.set_xlabel('size (T)')
    ax3.plot(truth['cm_mag'][:,1],np.log10(truth['cm_T']),'.')
    ax3.plot(truth_matched['cm_mag'][:,1],np.log10(truth_matched['cm_T']),'.',alpha=0.5)
    ax3.set_xlabel('mag')
    ax3.set_ylabel(r'log_10 cm_T')

    outfile = os.path.join(outdir, 'completeness-'+filetag+bands[bandind])
    fig.savefig(outfile)
    plt.close(fig)

    # Now do this for colors.

    gr_meas = meas_matched['cm_mag'][:,0] - meas_matched['cm_mag'][:,1]
    ri_meas = meas_matched['cm_mag'][:,1] - meas_matched['cm_mag'][:,2]
    iz_meas = meas_matched['cm_mag'][:,2] - meas_matched['cm_mag'][:,3]

    gr_truth = truth_matched['cm_mag'][:,0] - truth_matched['cm_mag'][:,1]
    ri_truth = truth_matched['cm_mag'][:,1] - truth_matched['cm_mag'][:,2]
    iz_truth = truth_matched['cm_mag'][:,2] - truth_matched['cm_mag'][:,3]

    dy = 2.0
    fig,(ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(21,7))
    ax1.plot(gr_truth, gr_meas - gr_truth,'.')
    ax1.axhline(0,color='red',linestyle='--')
    ax1.set_ylim(-dy,dy)
    ax1.set_xlim(-0.5,2.5)
    ax1.set_xlabel('g-r (truth)')
    ax1.set_ylabel('g-r (meas) - g-r (truth)')

    ax2.plot(ri_truth,ri_meas - ri_truth,'.')
    ax2.axhline(0,color='red',linestyle='--')
    ax2.set_ylim(-dy,dy)
    ax2.set_xlabel('r-i (truth)')
    ax2.set_ylabel('r-i (meas) - r-i (truth)')

    ax3.plot(iz_truth,iz_meas - iz_truth,'.')
    ax3.axhline(0,color='red',linestyle='--')
    ax3.set_ylim(-dy,dy)
    ax3.set_xlabel('i-z (truth)')
    ax3.set_ylabel('i-z (meas) - i-z (truth)')

    outfile = os.path.join(outdir, 'colors_vs_colors')
    fig.savefig(outfile)
    plt.close(fig)


    # Same for colors, but vs mag.
    dy = 1.6
    fig,(ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(21,7))
    ax1.plot(truth_matched['cm_mag'][:,1], gr_meas - gr_truth,'.')
    ax1.axhline(0,color='red',linestyle='--')
    ax1.set_ylim(-dy,dy)
    ax1.set_xlabel('r (truth)')
    ax1.set_ylabel('g-r (meas) - g-r (truth)')

    ax2.plot(truth_matched['cm_mag'][:,1],ri_meas - ri_truth,'.')
    ax2.axhline(0,color='red',linestyle='--')
    ax2.set_ylim(-dy,dy)
    ax2.set_xlabel('r (truth)')
    ax2.set_ylabel('r-i (meas) - r-i (truth)')

    ax3.plot(truth_matched['cm_mag'][:,1],iz_meas - iz_truth,'.')
    ax3.axhline(0,color='red',linestyle='--')
    ax3.set_ylim(-dy,dy)
    ax3.set_xlabel('r (truth)')
    ax3.set_ylabel('i-z (meas) - i-z (truth)')

    outfile = os.path.join(outdir, 'colors_vs_rmag')
    fig.savefig(outfile)
    plt.close(fig)

    # And finally, vs T

    fig,(ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(21,7))
    ax1.plot(truth_matched['cm_T'], gr_meas - gr_truth,'.')
    ax1.axhline(0,color='red',linestyle='--')
    ax1.set_ylim(-3,3)
    ax1.set_xscale('log')
    ax1.set_xlabel(' cm_T (truth)')
    ax1.set_ylabel('g-r (meas) - g-r (truth)')

    ax2.plot(truth_matched['cm_T'],ri_meas - ri_truth,'.')
    ax2.axhline(0,color='red',linestyle='--')
    ax2.set_ylim(-3,3)
    ax2.set_xscale('log')
    ax2.set_xlabel('cm_T (truth)')
    ax2.set_ylabel('r-i (meas) - r-i (truth)')

    ax3.plot(truth_matched['cm_T'],iz_meas - iz_truth,'.')
    ax3.axhline(0,color='red',linestyle='--')
    ax3.set_ylim(-3,3)
    ax3.set_xscale('log')
    ax3.set_xlabel('cm_T (truth)')
    ax3.set_ylabel('i-z (meas) - i-z (truth)')

    outfile = os.path.join(outdir, 'colors_vs_T')
    fig.savefig(outfile)
    plt.close(fig)

    # Variant of the below plot
    fig,ax1 = plt.subplots(nrows=1,ncols=1,figsize=(7,7))
    color = ['red','blue','black','cyan','orange']
    delta_mag =  meas_matched['cm_mag'][:,1] - truth_matched['cm_mag'][:,1]
    n_delta_bins = 5
    delta_mag_bds = np.arange(n_delta_bins+1)*0.5
    delta_mag_bds[-1] = 100
    for i in xrange(n_delta_bins):
        these = (np.abs(delta_mag) > delta_mag_bds[i]) & (np.abs(delta_mag) <= delta_mag_bds[i+1])
        ax1.plot(meas_matched['cm_T'][these],delta_mag[these],'.',alpha=0.5,color='k')
    ax1.set_xscale('log')
    ax1.set_xlabel('cm_T (r)')
    ax1.set_ylabel('cm_mag Measured - True (r)')

    outfile = os.path.join(outdir, 't_vs_mag_dif'+bands[bandind])
    fig.savefig(outfile)
    plt.close(fig)

    # Fun color-coded magnitude plot.
    fig,ax1 = plt.subplots(nrows=1,ncols=1,figsize=(7,7))
    color = ['red','blue','black','cyan','orange']
    delta_mag =  meas_matched['cm_mag'][:,1] - truth_matched['cm_mag'][:,1]
    n_delta_bins = 5
    delta_mag_bds = np.arange(n_delta_bins+1)*0.5
    delta_mag_bds[-1] = 100
    for i in xrange(n_delta_bins):
        these = (np.abs(delta_mag) > delta_mag_bds[i]) & (np.abs(delta_mag) <= delta_mag_bds[i+1])
        ax1.plot(meas_matched['cm_T'][these],meas_matched['cm_mag'][:,1][these],'.',alpha=0.5,color=color[i],label='%.1f <'%delta_mag_bds[i]+'|delta mag| < '+'%.1f'%delta_mag_bds[i+1])
    ax1.set_xscale('log')
    ax1.set_xlabel('cm_T (r)')
    ax1.set_ylabel('cm_mag (r)')
    ax1.legend(loc='best')

    outfile = os.path.join(outdir, 't_vs_mag_color-coded-by-delta-mag'+bands[bandind])
    fig.savefig(outfile)
    plt.close(fig)

    # Do matching errors matter?
    fig,ax1 = plt.subplots(figsize=(7,7))
    color = ['red','blue','black','cyan','orange']
    delta_mag =  meas_matched['cm_mag'][:,1] - truth_matched['cm_mag'][:,1]
    n_delta_bins = 5
    delta_mag_bds = np.arange(n_delta_bins+1)*0.5
    delta_mag_bds[-1] = 100
    for i in xrange(n_delta_bins):
        these = (np.abs(delta_mag) > delta_mag_bds[i]) & (np.abs(delta_mag) <= delta_mag_bds[i+1])
        ax1.plot(sep[these]*3600, meas_matched['cm_T'][these] - truth_matched['cm_T'][these],'.',alpha=0.5,color=color[i],label='%.1f <'%delta_mag_bds[i]+'|delta mag| < '+'%.1f'%delta_mag_bds[i+1])
    ax1.set_xlabel('match separation (arcsec)')
    ax1.set_ylabel('error in T')
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.legend(loc='best')

    outfile = os.path.join(outdir, 'sep_vs_t')
    fig.savefig(outfile)
    plt.close(fig)

    fig,ax1 = plt.subplots(figsize=(7,7))
    #small = (truth_matched['cm_T'] / truth_matched['cm_T_err'] ) > 10.
    small = truth_matched['cm_T']  < 10.
    ax1.plot(meas_matched['cm_T_err'],meas_matched['cm_T'] - truth_matched['cm_T'],'.',markersize=5)
    ax1.plot(meas_matched['cm_T_err'][small],meas_matched['cm_T'][small] - truth_matched['cm_T'][small],'.',markersize=4)

    ax1.axhline(0,color='red',linestyle='--')
    ax1.set_xlabel('reported error on cm_T')
    ax1.set_ylabel('actual error in cm_T')
    ax1.set_yscale('symlog')
    ax1.set_xscale('log')
    xlo,xhi = ax1.get_xlim()
    ax1.plot(np.logspace(np.log10(xlo),np.log10(xhi),1000),np.logspace(np.log10(xlo),np.log10(xhi),1000),color='red',linestyle='--')
    ax1.plot(np.linspace(np.log10(xlo),np.log10(xhi),1000),-np.linspace(np.log10(xlo),np.log10(xhi),1000),color='red',linestyle='--')

    outfile = os.path.join(outdir, 'Terr')
    fig.savefig(outfile)
    plt.close(fig)

    fig,ax1 = plt.subplots(figsize=(7,7))
    for i in xrange(n_delta_bins):
        these = (np.abs(delta_mag) >= delta_mag_bds[i]) & (np.abs(delta_mag) <= delta_mag_bds[i+1])
        ax1.plot(truth_matched['cm_logsb'][these,1],np.abs(meas_matched['cm_T'][these] - truth_matched['cm_T'][these])/meas_matched['cm_T_err'][these],'.',color=color[i],label='%.1f <'%delta_mag_bds[i]+'|delta mag| < '+'%.1f'%delta_mag_bds[i+1],alpha=0.5)

    ax1.set_yscale('symlog')
    #ax1.set_ylim(0,10)
    ax1.set_xlabel('cm_logsb')
    ax1.set_ylabel('| measured - input| cm_T / reported error')
    ax1.legend(loc='best')

    outfile = os.path.join(outdir, 'logsb_vs_cm_T')
    fig.savefig(outfile)
    plt.close(fig)

    # Shape and Correlation plots
    # fig, ax = ...
    # pudb.set_trace()

    if bandind == 1:
        # only do once!

        # min/max separation and bin size
        mins, maxs = 0.1, 10
        bs = 0.075

        # First do the truth catalog
        g1, g2 = truth_matched['cm_g'][:,0], truth_matched['cm_g'][:,1]
        e1, e2 = ngmix.shape.g1g2_to_e1e2(g1, g2)

        # Easiest to use treecorr by making a new temporary catalog
        truth_outfile = os.path.join(outdir,'treecorr_temp_file_truth.fits')
        delete_file(truth_outfile)
        fits = fitsio.FITS(truth_outfile,'rw')
        data = [truth_matched['ra'], truth_matched['dec'], g1, g2, e1, e2]
        names = ['ra', 'dec', 'g1', 'g2', 'e1', 'e2']
        fits.write(data, names=names)
        fits.close()

        cat = treecorr.Catalog(truth_outfile, ra_col='ra', dec_col='dec',
                            ra_units='degrees', dec_units='degrees',
                            g1_col='e1', g2_col='e2')
        gg = treecorr.GGCorrelation(min_sep=mins, max_sep=maxs, bin_size=bs,
                                    sep_units='arcmin')
        gg.process(cat)
        fig = plot_gg_corr(gg, plt_type='Truth')
        outfile = os.path.join(outdir, 'gg_corr_truth.png')
        fig.savefig(outfile)
        plt.close(fig)

        # Now for measured catalog
        g1, g2 = meas_matched['cm_g'][:,0], meas_matched['cm_g'][:,1]
        e1, e2 = ngmix.shape.g1g2_to_e1e2(g1, g2)

        # Easiest to use treecorr by making a new temporary catalog
        meas_outfile = os.path.join(outdir,'treecorr_temp_file_meas.fits')
        delete_file(meas_outfile)
        fits = fitsio.FITS(meas_outfile,'rw', clobber=True)
        data = [meas_matched['ra'], meas_matched['dec'], g1, g2, e1, e2]
        names = ['ra', 'dec', 'g1', 'g2', 'e1', 'e2']
        fits.write(data, names=names)
        fits.close()

        cat = treecorr.Catalog(meas_outfile, ra_col='ra', dec_col='dec',
                            ra_units='degrees', dec_units='degrees',
                            g1_col='e1', g2_col='e2')
        gg = treecorr.GGCorrelation(min_sep=mins, max_sep=maxs, bin_size=bs,
                                    sep_units='arcmin')
        gg.process(cat)
        fig = plot_gg_corr(gg, plt_type='Measured')
        outfile = os.path.join(outdir, 'gg_corr_meas.png')
        fig.savefig(outfile)
        plt.close(fig)

        # Now for differences
        g1t, g2t = truth_matched['cm_g'][:,0], truth_matched['cm_g'][:,1]
        e1t, e2t = ngmix.shape.g1g2_to_e1e2(g1t, g2t)
        g1m, g2m = meas_matched['cm_g'][:,0], meas_matched['cm_g'][:,1]
        e1m, e2m = ngmix.shape.g1g2_to_e1e2(g1m, g2m)
        g1d, g2d = g1m-g1t, g2m-g2t
        e1d, e2d = e1m-e1t, e2m-e2t

        # Easiest to use treecorr by making a new temporary catalog
        diff_outfile = os.path.join(outdir,'treecorr_temp_file_diff.fits')
        delete_file(diff_outfile)
        fits = fitsio.FITS(diff_outfile,'rw', clobber=True)
        data = [truth_matched['ra'], truth_matched['dec'], g1d, g2d, e1d, e2d]
        names = ['ra', 'dec', 'g1', 'g2', 'e1', 'e2']
        fits.write(data, names=names)
        fits.close()

        cat = treecorr.Catalog(diff_outfile, ra_col='ra', dec_col='dec',
                            ra_units='degrees', dec_units='degrees',
                            g1_col='e1', g2_col='e2')
        gg = treecorr.GGCorrelation(min_sep=mins, max_sep=maxs, bin_size=bs,
                                    sep_units='arcmin')
        gg.process(cat)
        fig = plot_gg_corr(gg, plt_type='Measured-True')
        outfile = os.path.join(outdir, 'gg_corr_diff.png')
        fig.savefig(outfile)
        plt.close(fig)

def delete_file(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise
def plot_gg_corr(gg, plt_type=None):
    r = np.exp(gg.meanlogr)
    xip = gg.xip
    xim = gg.xim
    sig = np.sqrt(gg.varxi)

    plt.plot(r, xip, color='blue')
    # plt.plot(r, -xip, color='blue', ls=':')
    plt.errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='blue', lw=1, ls='')
    # plt.errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='blue', lw=0.1, ls='')
    lp = plt.errorbar(-r, xip, yerr=sig, color='blue')

    plt.plot(r, xim, color='green')
    # plt.plot(r, -xim, color='green', ls=':')
    plt.errorbar(r[xim>0], xim[xim>0], yerr=sig[xim>0], color='green', lw=1, ls='')
    # plt.errorbar(r[xim<0], -xim[xim<0], yerr=sig[xim<0], color='green', lw=0.1, ls='')
    lm = plt.errorbar(-r, xim, yerr=sig, color='green')

    # Reference line
    plt.axhline(0, linestyle='--', c='k')

    plt.xscale('log')
    # if plt_type == 'Measured' or plt_type=='Truth':
    #     plt.yscale('log', nonposy='clip')
    # plt.yscale('log', nonposy='clip')
    plt.xlabel(r'$\theta$ (arcmin)')

    plt.legend([lp, lm], [r'$\xi_+(\theta)$', r'$\xi_-(\theta)$'])
    # plt.xlim( [1,100] )
    plt.ylabel(r'$\xi_{+,-}$')
    plt.title(plt_type, fontsize=16)

    plt.gcf().set_size_inches(7,7)

    return plt.gcf()

def make_star_plots(truth = None, meas_matched = None, truth_matched = None, filetag = ''):
    tag = 'cm_flux'
    errtag = 'cm_flux_cov'
    bandind = 2
    bands = ['g','r','i','z']

    pass

def main(argv):
    run_name = 'shear_test_ps'
    # run_name = 'sof_stars'
    # run_name = 'grid_with_noise'
    # run_name = 'grid_test_ngmix'
    # run_name = 'grid_test_shear_sof'
    path = os.path.join('/home/spencer/research/balrog/outputs/' + run_name)
    tilename = 'DES0347-5540'
    re_id = '0'
    min_sep = 0./3600 # minimum input object separation, degrees
    #truth_gal_catalog = fitsio.read(path+tilename+'_0/'+tilename+'_'+re_id+'_balrog_truth_cat_gals.fits')
    #truth_star_catalog = fitsio.read(path+tilename+'_0/'+tilename+'_'+re_id+'_balrog_truth_cat_stars.fits')
    #meas_catalog = fitsio.read(path+tilename+'_0/'+tilename+'_mof.fits')
    #truth_gal_matched, meas_gal_matched = match_catalogs(truth_catalog = truth_gal_catalog, meas_catalog =meas_catalog)

    truth_gal_catalog, truth_gal_matched, meas_gal_matched, sep = get_catalogs(path = path, tilename = tilename, re_id = re_id, minsep = min_sep,\
                                                                                ratag = 'ra', dectag = 'dec',stars=False)

    for i in xrange(4):
        make_plots(truth=truth_gal_catalog, meas_matched = meas_gal_matched, truth_matched=truth_gal_matched,
                   sep = sep,filetag='',bandind=i, outdir=run_name)
    #truth_star_catalog, truth_star_matched, meas_star_matched = get_catalogs(path = path, tilename = tilename, re_id = re_id, minsep = min_sep,\
    #                                                                            ratag = 'RA_new', dectag = 'DEC_new',stars=True)
    #make_star_plots(truth=truth_star_catalog,meas_matched = meas_star_matched, truth_matched = truth_star_matched)

if __name__ == "__main__":
    import pdb, traceback,sys
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
