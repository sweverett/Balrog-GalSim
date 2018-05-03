#####################################################################
#
# This file contains code relevant to the construction and use of
# an injection grid for Balrog objects (for now, likely just
# galaxies). Only rectangular and hexagonal grids currently
# supported.
#
#
# Spencer Everett
# UCSC
# 5/2/18
#####################################################################

import numpy as np
import random as rand
import os, sys, errno
import warnings
from astropy.wcs import WCS

class Grid(object):

    def __init__(self, grid_spacing, wcs, Npix_x=10000, Npix_y=10000, pixscale=0.2631):
        # stuff...
        self.grid_spacing = grid_spacing # arcsec
        self.Npix_x, self.Npix_y = Npix_x, Npix_y
        self.pixscale = pixscale # arcsec
        self.wcs = wcs

        return

class RectGrid(Grid):
    def __init__(self, grid_spacing, wcs, Npix_x=10000, Npix_y=10000, pixscale=0.2631):
        super(RectGrid, self).__init__(grid_spacing, wcs, Npix_x=Npix_x, Npix_y=Npix_y,
                                       pixscale=pixscale)
        self._create_grid()

        return

    def _create_grid(self):
        # gs = 40.0 / 3600.0 # grid spacing in deg
        # pixscale = 7.306e-5 # taken from tile header
        im_gs = self.grid_spacing * (1.0 / self.pixscale) # pixels
        self.im_grid_spacing = im_gs

        self.im_ra  = np.arange(im_gs, self.Npix_x-im_gs, im_gs)
        self.im_dec = np.arange(im_gs, self.Npix_y-im_gs, im_gs)

        # Get all image coordinate pairs
        self.im_pos = np.array(np.meshgrid(self.im_ra, self.im_dec)).T.reshape(-1, 2)
        self.pos = self.wcs.wcs_pix2world(self.im_pos, 1)
        self.ra  = self.pos[:, 0]
        self.dec = self.pos[:, 1]

        return

# TODO: Complete HexGrid class
class HexGrid(Grid):
    def __init__(self, grid_spacing, wcs, Npix_x=10000, Npix_y=10000, pixscale=0.2631):
        super(RectGrid, self).__init__(grid_spacing, wcs, Npix_x=Npix_x, Npix_y=Npix_y,
                                       pixscale=pixscale)
        self._create_grid()

        return

    def _create_grid(self):
        im_gs = self.grid_spacing * (1.0 / self.pixscale) # pixels
        self.im_grid_spacing = im_gs

        # TODO: Change to make a hexagonal grid
        # self.im_ra  = np.arange(im_gs, self.Npix_x-im_gs, im_gs)
        # self.im_dec = np.arange(im_gs, self.Npix_y-im_gs, im_gs)
        # # Get all image coordinate pairs
        # self.im_pos = np.array(np.meshgrid(self.im_ra, self.im_dec)).T.reshape(-1, 2)
        # self.pos = self.wcs.wcs_pix2world(self.im_pos, 1)

        return
