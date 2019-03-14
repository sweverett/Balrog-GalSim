import fitsio
from astropy.table import Table, Column, vstack, hstack
from astropy.io import fits
import corner
import numpy as np
from scipy.stats import norm
import scipy.stats as stats
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as Colors
import matplotlib
matplotlib.rcParams.update({'font.size': 24})

import seaborn as sb
import os
from glob import glob
import esutil.htm as htm
plt.style.use('seaborn')

# import pudb

class BalrogMatchedCatalog(object):

    def __init__(self, filename, cat_type='mof', prefix='cm_', use_deredden=True, real=0):

        # TODO: Implement full type-checking below
        self.filename = filename
        self.cat_type = cat_type
        self.prefix = prefix
        self.use_deredden = use_deredden
        self.b_indx = {'g':0, 'r':1, 'i':2, 'z':3}
        self.cov_indx = {'g':0, 'r':5, 'i':10, 'z':15}

        # Could be optional parameters in the future, but enforced for now
        self.true_prefix = 'true_'
        self.meas_prefix = 'meas_'

        p = prefix
        if self.use_deredden is True:
            self.true_flux_colname = self.true_prefix + p + 'flux_deredden'
            self.true_mag_colname  = self.true_prefix + p + 'mag_deredden'
            self.meas_flux_colname = self.meas_prefix + p + 'flux_deredden'
            self.meas_mag_colname  = self.meas_prefix + p + 'mag_deredden'
        else:
            self.true_flux_colname = self.true_prefix + p + 'flux'
            self.true_mag_colname  = self.true_prefix + p + 'mag'
            self.meas_flux_colname = self.meas_prefix + p + 'flux'
            self.meas_mag_colname  = self.meas_prefix + p + 'mag'

        if not isinstance(real, int):
            raise TypeError('real must be an int!')
        self.real = real

        return

    def plot_flux_chis(self, bands='griz', xlim=[-10, 10.0],
                       S=16, title=None, bins=40, show=False):
        p = self.prefix
        qt = self.true_flux_colname
        qm = p + 'flux'

        for band in bands:
            bi = self.b_indx[band]
            ci = self.cov_indx[band]
            true = self.true[qt][:,bi]
            meas = self.meas[qm][:,bi] / self.ext_flux[bi]

            diff = meas - true
            err = np.sqrt(self.meas[qm+'_cov'][:,bi,bi])#self.cov_indx[band]])
            chi = diff / err
            cuts = np.where( (chi > xlim[0]) & (chi < xlim[1]) )
            plt.subplot(2, 2, bi+1)
            n, bins, patches = plt.hist(chi[cuts], bins=bins, ec='k', normed=True)

            # add a 'best fit' line
            (mu, sigma) = norm.fit(chi[cuts])
            yy = mlab.normpdf(bins, 0.0, 1.0)
            y = mlab.normpdf(bins, mu, sigma)
            plt.plot(bins, yy, 'k:', linewidth=4, label=r'$\mu=%.2f,\ \sigma=%.2f$' %(0.0, 1.0))
            plt.plot(bins, y, 'r:', linewidth=4, label=r'$\mu=%.2f,\ \sigma=%.2f$' %(mu, sigma))

            ax = plt.gca()
            ax.axvline(0, ls='--', c='k', lw=3)
            plt.xlabel('Flux Chi (sigma)')
            plt.ylabel('Normed Counts')
#             plt.title(r'$\mathrm{Best Fit line:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))
            plt.legend(fontsize=12)
            if title: plt.suptitle(title)
        plt.gcf().set_size_inches(S, S)

        if show is True:
            plt.show()

        return

    def plot_flux_vs_err(self, bands='griz', title=None, xlim=[1.0, 5e5], ylim=[1e-5, 5e5],
                         logx=True, logy=True, cmap='inferno', alpha=0.15, S=16, show=True):
        p = self.prefix
        qt = self.true_flux_colname
        qm = p + 'flux'
        qe = qm + '_cov'

        for band in bands:
            bi = self.b_indx[band]
            true = self.true[qt][:,bi]
            meas = self.meas[qm][:,bi] / self.ext_flux[bi]
            true_err = self.true[qe][:,bi,bi]
            meas_err = self.meas[qe][:,bi,bi]

            cuts = np.where( (true>xlim[0]) & (true<xlim[1]) & (true_err>ylim[0]) & (true_err<ylim[1]) &
                             (meas>xlim[0]) & (meas<xlim[1]) & (meas_err>ylim[0]) & (meas_err<ylim[1]) )

            plt.subplot(2, 2, bi+1)
            ax = plt.gca()
            hb = plt.scatter(true[cuts], true_err[cuts], alpha=alpha, label='True')
            hb = plt.scatter(meas[cuts], meas_err[cuts], alpha=0.5*alpha, label='SOF')
#             sb.kdeplot(true[cuts], true_err[cuts], cmap='Greens',shade=True, shade_lowest=False, cut=5, ax=ax, gridsize=200)
#             sb.kdeplot(meas[cuts], meas_err[cuts], cmap='Blues',shade=True, shade_lowest=False, cut=5, ax=ax, gridsize=200)

            # Add contour lines
#             X, Y, Z = self.density_estimation(np.log10(true[cuts]), np.log10(true_err[cuts]),
#                                               np.log10(xlim), np.log10(ylim), dz=1.0)
#             plt.contour(X, Y, Z)
#             X, Y, Z = self.density_estimation(np.log10(meas[cuts]), np.log10(meas_err[cuts]),
#                                               np.log10(xlim), np.log10(ylim), dz=1.0)
#             print 'starting 1'
#             X, Y, Z = self.density_estimation(true[cuts], true_err[cuts],
#                                               xlim, ylim, dz=20000.0)
#             print 'starting 2'
#             plt.contour(X, Y, Z)
#             print 'starting 3'
#             X, Y, Z = self.density_estimation(meas[cuts], meas_err[cuts],
#                                               xlim, ylim, dz=20000.0)
#             print 'starting 4'
#             plt.contour(X, Y, Z) 
            plt.xlabel('True %s-band %s' % (band, q) )
            plt.ylabel('Meas %s-band %s err' % (band, q) )
            plt.ylim([1e-2, ylim[1]])
            plt.legend()
            if title: plt.suptitle(title)
            if logx: plt.xscale('log')
            if logy: plt.yscale('log')
        plt.gcf().set_size_inches(S, S)

        if show is True:
            plt.show()

        return

    def plot_over_bands(self, val, bands='griz', xlim=[0.0, 1.0], ylim=[0.0, 1.0],
                        S=16, title=None, cmap='inferno', dim=2, logx=False, gs=100,
                        show=True):
        p = self.prefix

        if val == 'flux':
            qt = self.true_flux_colname
            qm = self.meas_flux_colname
        elif (val == 'mag') or (val == 'color'):
            qt = self.true_mag_colname
            qm = self.meas_mag_colname
        else:
            # Not guaranteed to work, but try something sensible:
            qt = self.true_prefix + p + val
            qm = self.meas_prefix + p + val

        columns = [qt, qm]
        cat = fitsio.read(self.filename, columns=columns)

        bk = 0
        for band in bands:
            if len(band) != 1:
                if val != 'color':
                    raise ValueError('Can only pass single bands unless plotting over colors!')
                if len(band) != 3:
                    raise ValueError('If passing a color, must pass each band as '
                                     '\'{band1}-{band2}\'')

                b1 = self.b_indx[band[0]]
                b2 = self.b_indx[band[2]]
                # true = cat[qt][:,b1] - cat[qt][:,b2]
                true = cat[qt][:,b1]
                meas1 = cat[qm][:,b1]
                meas2 = cat[qm][:,b2]
                diff = meas1 - meas2

                # Keep track of which color we're on
                bk += 1

            else:
                bi = self.b_indx[band]
                true = cat[qt][:,bi]
                meas = cat[qm][:,bi]
                diff = meas - true

            cuts = np.where( (true > xlim[0]) & (true < xlim[1]) &
                             (diff > ylim[0]) & (diff < ylim[1]) )

            if logx:
                x = np.log10(true[cuts])
                lx = 'log10'
            else:
                x = true[cuts]
                lx = ''

            # fig, ax = plt.figure()

            if val == 'color':
                ax = plt.subplot(3, 1, bk)
            else:
                ax = plt.subplot(2, 2, bi+1)

            hb = plt.hexbin(x, diff[cuts], gridsize=gs, cmap=cmap, norm=Colors.LogNorm())
            ax.axhline(0, ls='--', c='k', lw=4)
            med = np.median(diff[cuts])
            ax.axhline(med, ls=':', c='w', lw=3, label='Median={:.3f}'.format(med))
            cb = plt.colorbar(hb, ax=ax)
            legend = plt.legend(bbox_to_anchor=(0.6, 0.925), bbox_transform=ax.transAxes, fontsize=18)
            plt.setp(legend.get_texts(), color='w')

            label = val
            if self.use_deredden is True:
                label += ' (dereddened)'
            ax.set_xlabel('%s True %s-band %s' % (lx, band, label))
            ax.set_ylabel('Meas-True %s-band %s' % (band, label))

            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(20)
            cb.ax.tick_params(labelsize=16)

            if title: plt.suptitle(title)
#             if logx: plt.xscale('log')

        if val == 'color':
            plt.gcf().set_size_inches(S/2., 2.*S)
        else:
            plt.gcf().set_size_inches(S, S)

        if show is True:
            plt.show()

        return

    # TODO: Try to incorporate this into `plot_over_bands()`. Almost works, but some details
    # to figure out
    def plot_colors(self, colors=['g-r','r-i','i-z'], xlim=[16., 26.], ylim=[-5., 5.],
                        S=16, title=None, cmap='inferno', dim=2, logx=False, gs=200,
                        show=True):
        p = self.prefix
        qt = self.true_mag_colname
        qm = self.meas_mag_colname

        columns = [qt, qm]
        cat = fitsio.read(self.filename, columns=columns)

        bk = 0
        for color in colors:
            if len(color) != 3:
                raise ValueError('If passing a color, must pass each band as '
                                    '\'{band1}-{band2}\'')

            b1 = self.b_indx[color[0]]
            b2 = self.b_indx[color[2]]
            # true = cat[qt][:,b1] - cat[qt][:,b2]
            true = cat[qt][:,b1]
            diff = cat[qm][:,b1] - cat[qm][:,b2]

            # Keep track of which color we're on
            bk += 1

            cuts = np.where( (true > xlim[0]) & (true < xlim[1]) &
                             (diff > ylim[0]) & (diff < ylim[1]) )

            plt.subplot(3, 1, bk)
            hb = plt.hexbin(true[cuts], diff[cuts], gridsize=gs, cmap=cmap, norm=Colors.LogNorm())
            ax = plt.gca()
            ax.axhline(0, ls='--', c='k', lw=4)
            med = np.median(diff[cuts])
            ax.axhline(med, ls=':', c='w', lw=3, label='Median={:.3f}'.format(med))
            cb = plt.colorbar(hb, ax=ax)
            legend = plt.legend(bbox_to_anchor=(0.6, 0.925), bbox_transform=ax.transAxes)
            plt.setp(legend.get_texts(), color='w')
            plt.xlabel('True %s-band mag (dereddened)' % (color[0]) )
            plt.ylabel('Meas-True %s-band mags (dereddened)' % (color) )
            if title:
                plt.suptitle(title)

        plt.gcf().set_size_inches(S/2., 1.5*S)

        if show is True:
            plt.show()

        return

    def plot_mags(self, bands='griz', xlim=[0.0, 1.0], ylim=[-1.0, 1.0],
                  S=16, title=None, cmap='inferno', logx=False, gs=100, show=True):
        self.plot_over_bands('mag', bands=bands, xlim=xlim, ylim=ylim,
                             S=S, title=title, cmap=cmap, logx=logx, gs=gs, show=show)
        return

    # TODO: Use this if we can incorporate it into `plot_over_bands()`
    # def plot_colors(self, bands=['g-r','r-i','i-z'], xlim=[0.0, 1.0], ylim=[-1.0, 1.0],
    #               S=16, title=None, cmap='inferno', logx=False, gs=100, show=True):
    #     self.plot_over_bands('color', bands=bands, xlim=xlim, ylim=ylim,
    #                          S=S, title=title, cmap=cmap, logx=logx, gs=gs, show=show)
    #     return

    def plot_fluxes(self, bands='griz', xlim=[0.0, 1e5], ylim=[-2e3, 2e3],
                    S=16, title=None, cmap='inferno', logx=True, show=True):
        self.plot_over_bands('flux', bands=bands, xlim=xlim, ylim=ylim,
                             S=S, title=title, cmap=cmap, logx=logx, show=show)
        return

    def plot_dist(self, val, minD=0.0, maxD=1.0, log=False, alpha=0.5, title=None, S=6, dD=0.25,
                  show=True):
        qt = self.true_prefix + self.prefix + val
        qm = self.meas_prefix + self.prefix + val

        columns = [qt, qm]
        cat = fitsio.read(self.filename, columns=columns)

        cuts = np.where( (cat[qt]<=maxD) &
                         (cat[qm]<=maxD) &
                         (cat[qm]>=minD) &
                         (cat[qt]>=minD) )
        if log:
            bins = np.power(10.0, np.arange(np.log10(minD), np.log10(maxD)+dD, dD))
        else:
            bins = np.arange(minD, maxD+dD, dD)

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        plt.hist(cat[qt][cuts], ec='k', bins=bins, label='True');
        plt.hist(cat[qm][cuts], ec='k', bins=bins, label='Meas', alpha=alpha);
        plt.legend()
        if title:
            plt.title(title)
        ax.set_xlabel(qt)
        ax.set_ylabel('Counts')
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)
        fig.set_size_inches(S,S)
        if log:
            plt.xscale('log', nonposx='clip')

        if show is True:
            plt.show()

        return

    def plot_T(self, minT=1e-4, maxT=100.0, alpha=0.5, title=None, S=6, dT=0.25, show=True):
        self.plot_dist('T', minD=minT, maxD=maxT, log=True, alpha=alpha, title=title, S=S, dD=dT,
                       show=show)
        return

    def plot_fracdev(self, minF=0.0, maxF=1.0, alpha=0.5, title=None, S=6, dF=0.05, show=True):
        self.plot_dist('fracdev', minD=minF, maxD=maxF, log=False, alpha=alpha, title=title, S=S, dD=dF,
                       show=show)
        return

    def plot_TdByTe(self, minTT=0.1, maxTT=1e2, alpha=0.5, title=None, S=6, dTT=0.25, show=True):
        self.plot_dist('TdByTe', minD=minTT, maxD=maxTT, log=True, alpha=alpha, title=title, S=S, dD=dTT,
                       show=show)
        return

    def density_estimation(self, m1, m2, xlim, ylim, dz=100.0):
        X, Y = np.mgrid[xlim[0]:xlim[1]:dz, ylim[0]:ylim[1]:dz]
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([m1, m2])
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)
        return X, Y, Z
