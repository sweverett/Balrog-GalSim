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
import matplotlib.pyplot as plt
# import pudb

class Grid(object):

    def __init__(self, grid_spacing, wcs, Npix_x=10000, Npix_y=10000, pixscale=0.2631):
        self.grid_spacing = grid_spacing  # arcsec
        self.Npix_x, self.Npix_y = Npix_x, Npix_y
        self.pixscale = pixscale  # arcsec
        self.wcs = wcs

        return

class RectGrid(Grid):
    def __init__(self, grid_spacing, wcs, Npix_x=10000, Npix_y=10000,
                 pixscale=0.2631):
        super(RectGrid, self).__init__(grid_spacing, wcs, Npix_x=Npix_x,
                                       Npix_y=Npix_y, pixscale=pixscale)
        self._create_grid()

        return

    def _create_grid(self):
        # gs = 40.0 / 3600.0 # grid spacing in deg
        # pixscale = 7.306e-5 # taken from tile header
        im_gs = self.grid_spacing * (1.0 / self.pixscale)  # pixels
        self.im_grid_spacing = im_gs

        self.im_ra  = np.arange(im_gs, self.Npix_x-im_gs, im_gs)
        self.im_dec = np.arange(im_gs, self.Npix_y-im_gs, im_gs)

        # Get all image coordinate pairs
        self.im_pos = np.array(np.meshgrid(self.im_ra, self.im_dec)).T.reshape(
                -1, 2)
        self.pos = self.wcs.wcs_pix2world(self.im_pos, 1)
        self.ra = self.pos[:, 0]
        self.dec = self.pos[:, 1]

        return

class HexGrid(Grid):
    def __init__(self, grid_spacing, wcs, Npix_x=10000, Npix_y=10000,
                 pixscale=0.2631):
        super(HexGrid, self).__init__(grid_spacing, wcs, Npix_x=Npix_x, Npix_y=Npix_y,
                                       pixscale=pixscale)
        self._create_grid()

        return

    def _create_grid(self):
        im_gs = self.grid_spacing * (1.0 / self.pixscale) # pixels
        self.im_grid_spacing = im_gs

        # p = HexGrid.calc_polygons(im_gs, im_gs, self.Npix_x-im_gs, self.Npix_y-im_gs, im_gs)
        # p = HexGrid.calc_polygons(0, 0, self.Npix_x, self.Npix_y, im_gs)
        # self.im_pos = HexGrid.polygons2coords(p)
        self.im_pos = HexGrid.calc_hex_coords(im_gs, im_gs, self.Npix_x-im_gs, self.Npix_y-im_gs, im_gs)

        self.im_ra  = self.im_pos[:, 0]
        self.im_dec = self.im_pos[:, 1]

        # Get all image coordinate pairs
        self.pos = self.wcs.wcs_pix2world(self.im_pos, 1)

        return

    @classmethod
    def calc_hex_coords(cls, HexGrid, startx, starty, endx, endy, radius):
        sl = (2 * radius) * np.tan(np.pi / 6)

        # Geoemtric factors of given hexagon
        p = sl * 0.5
        b = radius
        h = 2 * sl
        dx = b * 2
        dy = 3 * p

        row = 1

        xs = []
        ys = []

        while startx < endx:
            x = [startx, startx, startx + b, startx + dx, startx + dx, startx + b, startx]
            xs.append(x)
            startx += dx

        while starty < endy:
            if row % 2:
                y = [starty + p, starty + dy, starty + h, starty + dy, starty + p, starty, starty + p]
                ys.append(y)
            starty += dy
            row += 1

        polygons = [zip(x, y) for x in xs for y in ys]
        hexgrid = HexGrid.polygons2coords(polygons)

        # Some hexagonal elements go beyond boundary; cut these out
        indx = np.where( (hexgrid[:,0]<endx) & (hexgrid[:,1]<endy) )
        return hexgrid[indx]

    @classmethod
    def calc_polygons(HexGrid, startx, starty, endx, endy, radius):
        sl = (2 * radius) * np.tan(np.pi / 6)

        # calculate coordinates of the hexagon points
        p = sl * 0.5
        b = sl * np.cos(np.radians(30))
        w = b * 2
        h = 2 * sl

        # offsets for moving along and up rows
        xoffset = b
        yoffset = 3 * p

        row = 1

        shifted_xs = []
        straight_xs = []
        shifted_ys = []
        straight_ys = []

        while startx < endx:
            xs = [startx, startx, startx + b, startx + w, startx + w, startx + b, startx]
            straight_xs.append(xs)
            shifted_xs.append([xoffset + x for x in xs])
            startx += w

        while starty < endy:
            ys = [starty + p, starty + (3 * p), starty + h, starty + (3 * p), starty + p, starty, starty + p]
            (straight_ys if row % 2 else shifted_ys).append(ys)
            starty += yoffset
            row += 1

        polygons = [zip(xs, ys) for xs in shifted_xs for ys in shifted_ys] + [zip(xs, ys) for xs in straight_xs for ys in straight_ys]

        return polygons

    @classmethod
    def polygons2coords(HexGrid, p):
        s = np.shape(p)
        L = s[0]*s[1]
        pp = np.array(p).reshape(L,2)
	return np.vstack({tuple(row) for row in pp})
        # NOTE: Requires numpy 1.3.3
        # return np.unique(pp, axis=0)

    def rotate_polygons():
        return

class FibonacciGrid(Grid):
    def __init__(self, wcs, N=1000000):
        pass

    def make_fib_grid(self, N=1000000):
        self.golden_angle = np.pi * (3 - np.sqrt(5))
        theta = self.golden_angle * np.arange(N)
        z = np.linspace(1 - 1.0 / N, 1.0 / N - 1, N)
        radius = np.sqrt(1 - z * z)

        points = np.zeros((N, 3))
        points[:, 0] = radius * np.cos(theta)  # rad
        points[:, 1] = radius * np.sin(theta)  # rad
        points[:, 2] = z

        return points

def plot_fib_grid(points):
    fig = plt.figure()
    fig.set_size_inches(10,10)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:,0], points[:,1], zs=points[:,2])
    return

def plot_sphere(points):
    fig = plt.figure()
    fig.set_size_inches(10, 10)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:,0], points[:,1], zs=points[:,2])
    return
