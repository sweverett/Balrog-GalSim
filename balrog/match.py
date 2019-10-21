import fitsio
from astropy.table import Table, Column, vstack, hstack, join
from astropy.io import fits
import numpy as np

# TODO: These should be moved to a plotting library, rather than
# part of the class
#import corner
#from scipy.stats import norm
#import scipy.stats as stats
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
#import matplotlib.colors as colors
#import matplotlib
#rc = matplotlib.rcParams.update({'font.size': 20})

import os
from glob import glob
import esutil.htm as htm

from blacklisted import _blacklisted_tiles

#-------------------------------------------------------------------------------
# Photometry catalogs
class MatchedCatalogs(object):
    _extra_basename = 'TILENAME_merged.fits'

    # TODO: Add stars into all of this!
    def __init__(self, basedir, meds_conf='y3v02', real=0, tile_list=None, inj_type='gals',
                 match_type='default', ngmix_type='mof', vb=False, **kwargs):
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

        if ngmix_type is not None:
            if ngmix_type not in ['mof', 'sof']:
                raise ValueError('`ngmix_type` must be mof or sof!')
        self.ngmix_type = ngmix_type

        # This is checked in `build_matched_catalog()`
        self.match_type = match_type

        if not isinstance(vb, bool):
            raise TypeError('vb must be a bool!')
        self.vb = vb

        self.true_stack = None
        self.meas_stack = None
        self.full_true_stack = None

        self._has_matched = False

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
        Nt = len(tiles)
        removed = []
        k = 0
        for tile in tiles:
            k += 1
            if self.vb:
                print('Matching tile {} ({} of {})'.format(tile, k, Nt))

            if tile in _blacklisted_tiles:
                removed.append(tile)
                print('Tile is blacklisted - skipping tile')
                continue

            tdir = os.path.join(self.basedir, tile)
            self.tiledir[tile] = tdir
            true_name = '{}_{}_balrog_truth_cat_{}.fits'.format(tile, real, self.inj_type)
            meas_name = 'real_{}_{}-{}-{}.fits'.format(real, tile, self.meds_conf, self.ngmix_type)
            true_file = os.path.join(tdir, true_name)
            meas_file = os.path.join(tdir, meas_name)

            # The kwargs for each tile will be a little different, and need to
            # clean up a few
            tile_kwargs = kwargs.copy()
            del tile_kwargs['extra_base']
            del tile_kwargs['extra_subdir']

            if kwargs['extra_base'] is not None:
                extra_tbase = os.path.join(kwargs['extra_base'], tile)
                extra_filename = self._extra_basename.replace('TILENAME', tile)
                if kwargs['extra_subdir'] is not None:
                    extra_subdir = kwargs['extra_subdir']
                else:
                    extra_subdir = ''
                tile_kwargs['extra_catfile'] = os.path.join(extra_tbase,
                                                            extra_subdir,
                                                            extra_filename)

            tile_kwargs['tilename'] = tile
            tile_kwargs['real'] = real
            try:
                self.cats[tile] = build_matched_catalog(self.match_type, true_file, meas_file, **tile_kwargs)
            except IOError as e:
                print('Following IO error occured:\n{}\nSkipping tile.'.format(e))
                removed.append(tile)
                continue
            # TODO: Can turn on once we understand what's causing current assertion error
            # (i.e. len(gold) != len(sof))
            #except AssertionError as e:
            #    print('Following assertion error occured:\n{}\nSkipping tile.'.format(e))
            #    removed.append(tile)
            #    continue

            assert len(self.cats[tile].meas) == len(self.cats[tile].true)
            self.nobjects += len(self.cats[tile].meas)

        for t in removed:
            tiles.remove(t)

        self.ntiles += len(tiles)
        self.tiles += tiles

        return

    def get_full_true_stack(self):
        if self.full_true_stack is not None:
            return self.full_true_stack

        k = 0
        Nt = len(self.tiles)
        cats = []
        for tile, cat in self.cats.items():
            k += 1
            if self.vb:
                print('Stacking det tile {} ({} of {})'.format(tile, k, Nt))
            # if full_true_stack is None:
            #     full_true_stack = cat.det_cat
            # else:
            #     full_true_stack = vstack([full_true_stack, cat.det_cat])
            cats.append(cat.det_cat)

        if self.vb:
            print('Stacking all...')
        self.full_true_stack = vstack(cats)

        return self.full_true_stack

    def get_matched_stack(self):

        if self._has_matched is True:
            return self.true_stack, self.meas_stack

        true_stack = None
        meas_stack = None

        k = 0
        Nt = len(self.tiles)
        true_cats = []
        meas_cats = []
        for tile, cat in self.cats.items():
            k += 1
            if self.vb:
                print('Stacking matched tile {} ({} of {})'.format(tile, k, Nt))

            true_cats.append(cat.true)
            meas_cats.append(cat.meas)
            # if matched_stack is None:
            #     true_stack = cat.true
            #     meas_stack = cat.meas
            #     matched_stack = True
            # else:
            #     true_stack = vstack([true_stack, cat.true])
            #     meas_stack = vstack([meas_stack, cat.meas])

        if self.vb:
            print('Stacking all...')
        self.true_stack = vstack(true_cats)
        self.meas_stack = vstack(meas_cats)

        self._has_matched = True

        return self.true_stack, self.meas_stack

    def write_stacks(self, outdir=None, clobber=False):
        write_full_truth_stack(outdir=outdir, clobber=clobber)
        write_matched_stacks(outdir=outdir, clobber=clobber)

        return

    def _write_truth_base(self, outfile, outdir, clobber=False):
        if self.full_true_stack is None:
            stack = self.get_full_true_stack()
        else:
            # Don't re-stack unless specifically asked to do so
            stack = self.full_true_stack

        if outdir is None:
            outdir = self.basedir
        outfile = os.path.join(outdir, outfile)

        if os.path.exists(outfile) and (clobber is True):
            os.remove(outfile)

        return stack, outfile

    def write_full_truth_stack(self, outfile='full_truth_stack.fits', outdir=None, clobber=False):
        stack, outfile = self._write_truth_base(outfile=outfile, outdir=outdir, clobber=clobber)

        stack.write(outfile)

        # TODO: Add some metadata information to PHU!

        return

    def write_truth_det_stack(self, outfile='truth_det_stack.fits', outdir=None, clobber=False,
                              save_mags=True, save_gap_flux=False):
        stack, outfile = self._write_truth_base(outfile=outfile, outdir=outdir, clobber=clobber)

        # Detection stack only needs limited information
        det_stack = self._setup_det_cat(stack, save_mags=save_mags,
                                        save_gap_flux=save_gap_flux)
        det_stack.write(outfile)

        # TODO: Add some metadata information to PHU!

        return outfile

    def write_matched_stacks(self, outfile_base=None, outdir=None, clobber=False):
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

            if os.path.exists(outfile) and (clobber is True):
                os.remove(outfile)

            stack.write(outfile)

        # TODO: Add some metadata information to PHU!

        return

    def write_combined_stack(self, outdir=None, outfile='balrog_matched_catalog.fits',
                             cache=False, clobber=False, table_format='fits'):
        if self._has_matched is False:
            true_stack, meas_stack = self.get_matched_stack()
            if self.vb:
                print('Combined stack built')
        else:
            # Don't re-stack unless specifically asked to do so
            true_stack, meas_stack = self.true_stack, self.meas_stack
            if self.vb:
                print('Loaded combined stack')

        combined_stack = self._merge_matches(true_stack, meas_stack)

        if outdir is None:
            outdir = self.basedir

        outfile = os.path.join(outdir, outfile)

        if os.path.exists(outfile) and (clobber is True):
            os.remove(outfile)

        if self.vb:
            print('Combined stack writing...')
        combined_stack.write(outfile)

        # TODO: Add some metadata information to PHU!

        return

    def write_combined_cats(self, outdir=None, outbase='balrog_matched_cat', clobber=False):
        if outdir is None:
            outdir = self.basedir

        Nt = len(self.tiles)
        i = 0
        for tile, cat in self.cats.items():
            i += 1
            if self.vb is True:
                print('Writing tile {} ({} of {})'.format(tile, i, Nt))

            outfile = os.path.join(outdir, '{}_{}.fits'.format(tile, outbase))
            if os.path.exists(outfile) and (clobber is True):
                os.remove(outfile)

            true_stack = cat.true
            meas_stack = cat.meas
            combined_stack = self._merge_matches(true_stack, meas_stack)

            combined_stack.write(outfile)

        return

    def _setup_det_cat(self, cat, save_mags=True, save_gap_flux=False):
        # Detection stack only needs limited information
        det_cat= Table()
        det_cat['bal_id'] = cat['bal_id']
        det_cat['true_id'] = cat['id']
        det_cat['meas_tilename'] = cat['meas_tilename']
        # NOTE: Not present in prerun2
        #det_cat['true_number'] = cat['number'] # 'number' is a reserved name in DB
        det_cat['true_ra'] = cat['ra']
        det_cat['true_dec'] = cat['dec']
        det_cat['meas_id'] = cat['meas_id']
        det_cat['detected'] = cat['detected']

        if save_mags is True:
            det_cat['true_'+self.true_mag_colname] = cat[self.true_mag_colname]

        if save_gap_brightness is True:
            det_cat['true_gap_riz_flux_deredden'] = cat['gap_riz_flux_deredden']

        return det_cat

    def write_det_cats(self, outdir=None, outbase='balrog_det_cat', save_mags=save_mags,
                       clobber=False, save_gap_flux=False):
        if outdir is None:
            outdir = self.basedir

        Nt = len(self.tiles)
        i = 0
        for tile, cat in self.cats.items():
            i += 1
            if self.vb is True:
                print('Writing det tile {} ({} of {})'.format(tile, i, Nt))

            outfile = os.path.join(outdir, '{}_{}.fits'.format(tile, outbase))
            if os.path.exists(outfile) and (clobber is True):
                os.remove(outfile)

            det = cat.det_cat
            det_cat = self._setup_det_cat(det, save_mags=save_mags,
                                          save_gap_flux=save_gap_flux)
            det_cat.write(outfile)

        return

    def _merge_matches(self, true_stack, meas_stack):
        skip_cols = ['separation', 'bal_id']

        # Add column prefixes (avoids table copying)
        for col in true_stack.colnames:
            if col in skip_cols:
                continue
            true_stack.rename_column(col, 'true_'+col)
        for col in meas_stack.colnames:
            if col in skip_cols:
                continue
            meas_stack.rename_column(col, 'meas_'+col)

        combined_stack = hstack([true_stack, meas_stack],
                                join_type='exact',
                                table_names=['true', 'meas'],
                                uniq_col_name='{table_name}_{col_name}')

        # Remove column prefixes (avoids table copying)
        for col in true_stack.colnames:
            if col in skip_cols:
                continue
            clist = col.split('_')[1:]
            true_stack.rename_column(col, '_'.join(clist))
        for col in meas_stack.colnames:
            if col in skip_cols:
                continue
            clist = col.split('_')[1:]
            meas_stack.rename_column(col, '_'.join(clist))

        return combined_stack

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

    # In principle we'd match gold catalogs directly. Until then, we may
    # want to include certain value-added cols from Gold-like cats in the
    # matched ngmix catalogs
    _allowed_extra_cols = ['COADD_OBJECT_ID',
                           'FLAGS_GOLD',
                           'EXTENDED_CLASS_MOF',
                           'EXTENDED_CLASS_SOF']

    b_indx   = {'g':0, 'r':1, 'i':2,  'z':3}
    cov_indx = {'g':0, 'r':5, 'i':10, 'z':15}

    def __init__(self, true_file, meas_file, prefix='cm_', ratag='ra', dectag='dec',
                 match_radius=1.0/3600, depth=14, de_reddened=False, ext_flux=None,
                 ext_mag=None, tilename=None, real=None, make_cuts=False,
                 extra_catfile=None, save_gap_flux=False):

        self.true_file = true_file
        self.meas_file = meas_file
        self.prefix = prefix
        self.ratag = ratag
        self.dectag = dectag
        self.match_radius = match_radius
        self.depth = depth
        self.de_reddened = de_reddened
        self.save_gap_flux = save_gap_flux

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

        # For when the matches are entirely contained in a tile
        if (tilename is not None) and (not isinstance(tilename, str)):
            raise TypeError('tilename must be a string!')
        self.tilename = tilename

        if (real is not None) and (not isinstance(real, int)):
            raise TypeError('real must be an int!')
        self.real = real

        if not isinstance(make_cuts, bool):
            raise TypeError('make_cuts must be a bool!')
        self.make_cuts = make_cuts

        if extra_catfile is not None:
            if not isinstance(extra_catfile, str):
                raise TypeError('extra_catfile must be a string!')
        self.extra_catfile = extra_catfile

        self._match()

        return

    def _match(self):
        true_cat, meas_cat = self._load_cats()
        h = htm.HTM(self.depth)
        self.matcher = htm.Matcher(depth=self.depth, ra=true_cat[self.ratag], dec=true_cat[self.dectag])
        id_m, id_t, dist = self.matcher.match(ra=meas_cat[self.ratag], dec=meas_cat[self.dectag],
                                              radius=self.match_radius)
        self.true = true_cat[id_t]
        self.meas = meas_cat[id_m]
        self.meas['separation'] = dist
        self.dist = dist

        if self.tilename is not None:
            tilenames = np.array(len(dist) * [self.tilename])
            self.meas['tilename'] = tilenames

        assert len(self.true) == len(self.meas)

        # Need for detection efficiency plots
        self.det_cat = true_cat
        dcol = Column(name='detected', data=np.zeros(len(self.det_cat)), dtype=int)
        icol = Column(name='meas_id', data=-1.*np.ones(len(self.det_cat)), dtype=int)
        self.det_cat.add_column(dcol)
        self.det_cat.add_column(icol)
        self.det_cat['detected'][id_t] = 1
        self.det_cat['meas_id'][id_t] = self.meas['id']

        # Sometimes want to compute an avg dereddened riz Gaussian
        # aperture flux for brightness comparisons
        if self.save_gap_flux is True:
            self._compute_gap_flux()

        if self.tilename is not None:
            L = len(self.det_cat)
            tilenames = np.array(L * [self.tilename])
            col = Column(name='meas_tilename', data=np.zeros(L), dtype=str)
            self.det_cat.add_column(col)
            self.det_cat['meas_tilename'] = tilenames

        if self.make_cuts is True:
            self._make_cuts()

        return

    def _load_cats(self):
        true_cat = Table(fits.getdata(self.true_file, ext=1))
        meas_cat = Table(fits.getdata(self.meas_file, ext=1))

        # Sometimes might want a few extra columns not in MOF/SOF
        # (e.g. FLAGS_GOLD) appended to the matched catalog
        if self.extra_catfile is not None:
            extra_cat = Table(fitsio.read(self.extra_catfile,
                                          columns=self._allowed_extra_cols))

            assert len(meas_cat) == len(extra_cat)
            # Need the ID colnames to match
            extra_cat.rename_column('COADD_OBJECT_ID', 'id')
            meas_cat = join(meas_cat, extra_cat, keys='id', join_type='left')

        self._assign_bal_id(true_cat)

        return true_cat, meas_cat

    def _compute_gap_flux(self):
        gap_riz_deredden = Column(name='gap_riz_flux_deredden',
                                  data=-1.*np.ones(len(self.det_cat)),
                                  dtype=float)

        flux_factor = self.det_cat['bdf_flux_deredden'] /
                      self.det_cat['bdf_flux']

        gap_deredden = self.det_cat['gap_flux'] * flux_factor

        # Average over riz gap fluxes
        gap_riz_deredden[:] = np.mean(gap_deredden[:,1:], axis=1)

        self.det_cat.add_column(gap_riz_deredden)

        return

    def _make_cuts(self):
        # TODO: Add more flexibility in future!
        cuts = np.where(self.meas['flags']==0)
        bad = np.where(self.meas['flags']!=0)

        # Don't want to count these in efficiency plots
        ids = self.true['id'][bad]
        det_cat_ids = np.where(self.det_cat['id'] in ids)
        self.det_cat['detected'][det_cat_ids] = -1

        self.true = self.true[cuts]
        self.meas = self.meas[cuts]
        self.dist = self.dist[cuts]

        return

    def _assign_bal_id(self, true_cat):
        if self.tilename is not None:
            # base = self.tilename + '-'
            base = self.tilename.replace('DES','')
            # Replace +/- with 1/0
            base = base.replace('+', '1').replace('-', '0')
        else:
            base = ''

        if self.real is not None:
            #real = str(self.real) + '-'
            real = str(self.real)
        else:
            real = ''

        # bal_id = [int(base + real + str(i)) for i in range(len(true_cat))]
        bal_id = [int(str(1) + real + base + str(i)) for i in range(len(true_cat))]
        # bal_id = [int(real + base + str(meas_cat['id'][i])) for i in range(len(meas_cat))]

        # Save to self.true as it's also needed in det stack
        true_cat['bal_id'] = bal_id

        return

    def plot_flux_chis(self, bands='griz', xlim=[-10, 10.0],
                       S=16, title=None, bins=40, show=False):
        p = self.prefix
        qt = self.true_flux_colname
        qm = p + 'flux'

        for band in bands:
            bi = b_indx[band]
            ci = cov_indx[band]
            true = self.true[qt][:,bi]
            meas = self.meas[qm][:,bi] / self.ext_flux[bi]

            diff = meas - true
            err = np.sqrt(self.meas[qm+'_cov'][:,bi,bi])#cov_indx[band]])
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
            bi = b_indx[band]
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
            bi = b_indx[band]
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
            # legend = plt.legend(bbox_to_anchor=(0.6, 0.925), bbox_transform=ax.transAxes)
            legend = plt.legend()
            # plt.setp(legend.get_texts(), color='w')
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

class MOFOnlyMatchedCatalog(MatchedCatalog):
    '''
    Use for MOF-only Balrog runs that won't have SOF columns for extra cols
    '''

    _allowed_extra_cols = ['COADD_OBJECT_ID',
                           'FLAGS_GOLD_MOF_ONLY',
                           'EXTENDED_CLASS_MOF']

class SOFOnlyMatchedCatalog(MatchedCatalog):
    '''
    Use for SOF-only Balrog runs that won't have MOF columns for extra cols
    '''

    _allowed_extra_cols = ['COADD_OBJECT_ID',
                           'FLAGS_GOLD_SOF_ONLY',
                           'EXTENDED_CLASS_SOF']

def build_matched_catalog(match_type, *args, **kwargs):
    if match_type in MATCHED_CATALOG_TYPES:
        return MATCHED_CATALOG_TYPES[match_type](*args, **kwargs)
    else:
        raise ValueError('{} is not a valid matched catalog type. Allowed tyes are: {}'.format(
            match_type, MATCHED_CATALOG_TYPES.keys()))

MATCHED_CATALOG_TYPES = {
    'default' : MatchedCatalog,
    'matched_catalog' : MatchedCatalog,
    'mof_only' : MOFOnlyMatchedCatalog,
    'sof_only' : SOFOnlyMatchedCatalog
    }
