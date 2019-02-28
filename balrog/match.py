import fitsio
from astropy.table import Table, Column, vstack, hstack
from astropy.io import fits
import corner
import numpy as np
from scipy.stats import norm
import scipy.stats as stats
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib

import seaborn as sb
import os
from glob import glob
import esutil.htm as htm
rc = matplotlib.rcParams.update({'font.size': 20})
plt.style.use('seaborn')

import pudb

# TODO: When constructing MatchedCatalogs, grab the extinction factors, etc. from the truth catalogs!
class MatchedCatalogs(object):

    # TODO: Add stars into all of this!
    def __init__(self, basedir, meds_conf='y3v02', real=0, tile_list=None, inj_type='gals', **kwargs):
        if not isinstance(basedir, str):
            raise TypeError('basedir must be a string!')
        if not os.path.isdir(basedir):
            raise ValueError('{} is not a directory!'.format(basedir))
        if not os.path.exists(basedir):
            raise ValueError('{} does not exist!'.format(basedir))
        self.basedir = os.path.abspath(basedir)

        if tile_list is not None:
            if not isinstance(tile_list, list):
                raise TypeError('tile_list must be a list if passed!')
        self.tile_list = tile_list

        if tile_list is not None:
            self.tiles = tile_list
            self.ntiles = len(self.tile_list)
        else:
            self.tiles = []
            self.ntiles = 0

        if not isinstance(meds_conf, str):
            raise TypeError('meds_conf must be a string!')
        self.meds_conf = meds_conf

        if not isinstance(real, int):
            raise TypeError('The realization `real` must be an int!')
        self.real = real

        if not isinstance(inj_type, str):
            raise TypeError('inj_type must be a string!')
        if inj_type == 'gals':
            self.inj_type = 'ngmixGalaxy'
        elif inj_type == 'stars':
            self.inj_type = 'desStar'
        # TODO: Implement!
        # elif inj_type.lower() == 'both':
        else:
            raise ValueError('Must pass inj_type as `gals` or `stars`!')

        self.true_stack = None
        self.meas_stack = None
        self.full_true_stack = None

        self._has_matched = False

        # TODO: Can add things from the kwargs here as well!

        self.tiledir = {}
        self.cats = {}

        self._added_base_tiles = False
        self.add_catalogs(tile_list, **kwargs)

        return

    def add_catalogs(self, tiles=None, **kwargs):
        if tiles is None:
            assert self._added_base_tiles is False
        if tiles is not None:
            if not isinstance(tiles, (str, list)):
                raise TypeError('tiles must be passed as a str (for 1) or in a list!')

        if self._added_base_tiles is False:
            if tiles is None:
                # Grab all tiles from base directory
                tilepaths = glob(self.basedir+'/*/')
                tiles = [os.path.basename(os.path.normpath(tilepath)) for tilepath in tilepaths
                         if 'DES' in tilepath]

        real = self.real
        self.nobjects = 0
        for tile in tiles:
            tdir = os.path.join(self.basedir, tile)
            self.tiledir[tile] = tdir
            true_name = '{}_{}_balrog_truth_cat_{}.fits'.format(tile, real, self.inj_type)
            meas_name = 'real_{}_{}-{}-mof.fits'.format(real, tile, self.meds_conf)
            true_file = os.path.join(tdir, true_name)
            meas_file = os.path.join(tdir, meas_name)

            self.cats[tile] = MatchedCatalog(true_file, meas_file, **kwargs)
            assert len(self.cats[tile].meas) == len(self.cats[tile].true)
            self.nobjects += len(self.cats[tile].meas)

        self.ntiles += len(tiles)
        self.tiles.append(tiles)

        return

    def get_full_true_stack(self):
        if self.full_true_stack is not None:
            return self.full_true_stack

        full_true_stack = None

        for tile, cat in self.cats.items():
            if full_true_stack is None:
                full_true_stack = cat.true_cat
            else:
                full_true_stack = vstack([full_true_stack, cat.true_cat])

        self.full_true_stack = full_true_stack

        return full_true_stack

    def get_matched_stack(self):

        if self._has_matched is True:
            return self.true_stack, self.meas_stack

        matched_stack = None
        true_stack = None
        meas_stack = None

        for tile, cat in self.cats.items():
            if matched_stack is None:
                true_stack = cat.true
                meas_stack = cat.meas
                matched_stack = True
            else:
                true_stack = vstack([true_stack, cat.true])
                meas = cat.meas
                meas['separation'] = cat.dist
                meas_stack = vstack([meas_stack, meas])

        self.true_stack = true_stack
        self.meas_stack = meas_stack

        self._has_matched = True

        return true_stack, meas_stack

    def write_stacks(self, outdir=None):
        write_full_truth_stack(outdir=outdir)
        write_matched_stacks(outdir=outdir)

        return

    def _write_truth_base(self, outfile, outdir):
        if self.full_true_stack is None:
            stack = self.get_full_true_stack()
        else:
            # Don't re-stack unless specifically asked to do so
            stack = self.full_true_stack

        if outdir is None:
            outdir = self.basedir
        outfile = os.path.join(outdir, outfile)

        if os.path.exists(outfile):
            os.remove(outfile)

        return stack, outfile

    def write_full_truth_stack(self, outfile='full_truth_stack.fits', outdir=None):
        stack, outfile = self._write_truth_base(outfile=outfile, outdir=outdir)

        stack.write(outfile)

        # TODO: Add some metadata information to PHU!

        return

    def write_truth_det_stack(self, outfile='truth_det_stack.fits', outdir=None):
        stack, outfile = self._write_truth_base(outfile=outfile, outdir=outdir)

        # Detection stack only needs limited information
        det_stack = Table()
        det_stack['id'] = stack['id']
        det_stack['id_number'] = stack['number'] # 'number' is a reserved name in DB
        det_stack['ra'] = stack['ra']
        det_stack['dec'] = stack['dec']
        det_stack['detected'] = stack['detected']

        det_stack.write(outfile)

        # TODO: Add some metadata information to PHU!

        return

    def write_matched_stacks(self, outfile_base=None, outdir=None):
        if self.meas_stack is not None:
            true_stack, meas_stack = self.get_matched_stack()
        else:
            # Don't re-stack unless specifically asked to do so
            true_stack, meas_stack = self.true_stack, self.meas_stack

        outfiles = {'truth_stack.fits': true_stack,
                    'meas_stack.fits' : meas_stack }

        if outfile_base is None:
            outfile_base = ''
        if outdir is None:
            outdir = self.basedir

        for outfile, stack in outfiles.items():
            outfile = os.path.join(outdir, outfile_base, outfile)

            if os.path.exists(outfile):
                os.remove(outfile)

            stack.write(outfile)

        # TODO: Add some metadata information to PHU!

        return

    def write_combined_stack(self, outdir=None, outfile='balrog_matched_catalog.fits'):
        if self._has_matched is False:
            true_stack, meas_stack = self.get_matched_stack()
        else:
            # Don't re-stack unless specifically asked to do so
            true_stack, meas_stack = self.true_stack, self.meas_stack

        # Add column prefixes (avoids table copying)
        for col in true_stack.colnames:
            true_stack.rename_column(col, 'true_'+col)
        for col in meas_stack.colnames:
            if col == 'separation':
                continue
            meas_stack.rename_column(col, 'meas_'+col)

        self.combined_stack = hstack([true_stack, meas_stack],
                                      join_type='exact',
                                      table_names=['true', 'meas'],
                                      uniq_col_name='{table_name}_{col_name}')

        # Remove column prefixes (avoids table copying)
        for col in true_stack.colnames:
            clist = col.split('_')[1:]
            true_stack.rename_column(col, '_'.join(clist))
        for col in meas_stack.colnames:
            if col == 'separation':
                continue
            clist = col.split('_')[1:]
            meas_stack.rename_column(col, '_'.join(clist))

        if outdir is None:
            outdir = self.basedir

        outfile = os.path.join(outdir, outfile)

        if os.path.exists(outfile):
            os.remove(outfile)

        self.combined_stack.write(outfile)

        # TODO: Add some metadata information to PHU!

        return

    # TODO: make single function that handles any MatchedCatalog plot func
    # def ...

    def plot_flux_chis(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_flux_chis(show=show, **kwargs)

    def plot_flux_vs_error(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_flux_vs_error(show=show, **kwargs)

    def plot_over_bands(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_over_bands(show=show, **kwargs)

    def plot_mags(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_mags(show=show, **kwargs)

        return

    def plot_fluxes(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_fluxes(show=show, **kwargs)

    def plot_dist(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_dist(show=show, **kwargs)

    def plot_T(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_T(show=show, **kwargs)

    def plot_fracdev(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_fracdev(show=show, **kwargs)

    def plot_TdByTe(self, **kwargs):
        count = 0
        for tile, cat in self.cats.items():
            count +=1
            if count < self.ntiles:
                show = False
            else:
                show = True
            cat.plot_TdByTe(show=show, **kwargs)

    def plot_density_estimation(self, **kwargs):
        pass

class MatchedCatalog(object):

    def __init__(self, true_file, meas_file, prefix='cm_', ratag='ra', dectag='dec',
                 match_radius=1.0/3600, depth=14, de_reddened=False, ext_flux=None,
                 ext_mag=None):

        self.true_file = true_file
        self.meas_file = meas_file
        self.prefix = prefix
        self.ratag = ratag
        self.dectag = dectag
        self.match_radius = match_radius
        self.depth = depth
        self.de_reddened = de_reddened
        self.b_indx = {'g':0, 'r':1, 'i':2, 'z':3}
        self.cov_indx = {'g':0, 'r':5, 'i':10, 'z':15}

        p = prefix
        if self.de_reddened is True:
            self.true_flux_colname = p + 'flux_deredden'
            self.true_mag_colname = p + 'mag_deredden'
        else:
            self.true_flux_colname = p + 'flux'
            self.true_mag_colname = p + 'mag'

        if ext_flux is None:
            self.ext_flux = [1., 1., 1., 1.]
        else:
            assert ext_mag is not None
            self.ext_flux = ext_flux
        if ext_mag is None:
            self.ext_mag = [0., 0., 0., 0.]
        else:
            assert ext_flux is not None
            self.ext_mag = ext_mag

        self._match()

        return

    def _match(self):
        true_cat, meas_cat = self._load_cats()
        self.matcher = htm.Matcher(depth=self.depth, ra=true_cat[self.ratag], dec=true_cat[self.dectag])
        id_m, id_t, dist = self.matcher.match(ra=meas_cat[self.ratag], dec=meas_cat[self.dectag],
                                              radius=self.match_radius)
        self.true = true_cat[id_t]
        self.meas = meas_cat[id_m]
        self.dist = dist

        assert len(self.true) == len(self.meas)

        # Need for detection efficiency plots
        self.true_cat = true_cat
        col = Column(name='detected', data=np.zeros(len(self.true_cat)), dtype=int)
        self.true_cat.add_column(col)
        self.true_cat['detected'][id_t] = 1

        self._make_cuts()

        return

    def _load_cats(self):
        true_cat = Table(fits.getdata(self.true_file, ext=1))
        meas_cat = Table(fits.getdata(self.meas_file, ext=1))

        return true_cat, meas_cat

    def _make_cuts(self):
        # TODO: Add more flexibility in future!
        cuts = np.where(self.meas['flags']==0)
        bad = np.where(self.meas['flags']!=0)

        # Don't want to count these in efficiency plots
        ids = self.true['id'][bad]
        true_cat_ids = np.where(self.true_cat['id'] in ids)
        self.true_cat['detected'][true_cat_ids] = -1

#         cuts = np.ones(len(self.meas['flags']), ctype=bool)
        self.true = self.true[cuts]
        self.meas = self.meas[cuts]
        self.dist = self.dist[cuts]

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
        qm = p + val

        if val == 'flux':
            qt = self.true_flux_colname
        elif val == 'mag':
            qt = self.true_mag_colname
        else:
            qt = p + val

        for band in bands:
            bi = self.b_indx[band]
            true = self.true[qt][:,bi]

            if val == 'flux':
                meas = self.meas[qm][:,bi] / self.ext_flux[bi]
            elif val == 'mag':
                meas = self.meas[qm][:,bi] - self.ext_mag[bi]
            else:
                meas = self.meas[qm][:,bi]

            diff = meas - true
            cuts = np.where( (true > xlim[0]) & (true < xlim[1]) &
                             (diff > ylim[0]) & (diff < ylim[1]) )
            plt.subplot(2, 2, bi+1)
            if logx:
                x = np.log10(true[cuts])
                lx = 'log10'
            else:
                x = true[cuts]
                lx = ''
            hb = plt.hexbin(x, diff[cuts], gridsize=gs, cmap=cmap, norm=colors.LogNorm())
            ax = plt.gca()
            ax.axhline(0, ls='--', c='k', lw=4)
            med = np.median(diff[cuts])
#             if show is True:
            ax.axhline(med, ls=':', c='w', lw=3, label='Median={:.3f}'.format(med))
            cb = plt.colorbar(hb, ax=ax)
            legend = plt.legend(bbox_to_anchor=(0.6, 0.925), bbox_transform=ax.transAxes)
            plt.setp(legend.get_texts(), color='w')
            plt.xlabel('%s True %s-band %s' % (lx, band, qt) )
            plt.ylabel('Meas-True %s-band (%s-%s)' % (band, qm, qt) )
            if title: plt.suptitle(title)
#             if logx: plt.xscale('log')
        plt.gcf().set_size_inches(S, S)

        if show is True:
            plt.show()

        return

    def plot_mags(self, bands='griz', xlim=[0.0, 1.0], ylim=[-1.0, 1.0],
                  S=16, title=None, cmap='inferno', logx=False, gs=100, show=True):
        self.plot_over_bands('mag', bands=bands, xlim=xlim, ylim=ylim,
                             S=S, title=title, cmap=cmap, logx=logx, gs=gs, show=show)
        return

    def plot_fluxes(self, bands='griz', xlim=[0.0, 1e5], ylim=[-2e3, 2e3],
                    S=16, title=None, cmap='inferno', logx=True, show=True):
        self.plot_over_bands('flux', bands=bands, xlim=xlim, ylim=ylim,
                             S=S, title=title, cmap=cmap, logx=logx, show=show)
        return

    def plot_dist(self, val, minD=0.0, maxD=1.0, log=False, alpha=0.5, title=None, S=6, dD=0.25,
                  show=True):
        q = self.prefix + val
        cuts = np.where( (self.true[q]<maxD) &
                         (self.meas[q]<maxD) &
                         (self.meas[q]>minD) &
                         (self.true[q]>minD) )
        if log:
            bins = np.power(10.0, np.arange(np.log10(minD), np.log10(maxD)+dD, dD))
        else:
            bins = np.arange(minD, maxD+dD, dD)
        plt.hist(self.true[q][cuts], ec='k', bins=bins, label='True');
        plt.hist(self.meas[q][cuts], ec='k', bins=bins, label='SOF', alpha=alpha);
        plt.legend()
        if title:
            plt.title(title)
        plt.xlabel(q)
        plt.ylabel('Counts')
        plt.gcf().set_size_inches(S,S)
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
