# Adapted from code by Daniel Gruen
import numpy as np
import h5py
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument(
    'mcal_file',
    type=str,
    help='Balrog mcal file to add weights to'
    )
# parser.add_argument(
#     '--max_percentile',
#     default=97.5,
#     type=float,
#     help='Maximum percentile before pileup'
#     )
parser.add_argument(
    '--make_plots',
    action='store_true',
    default=False,
    help='Set to make plots'
    )
parser.add_argument(
    '--save_weight_table',
    action='store_true',
    default=False,
    help='Set to save weight table as a txt file'
    )
parser.add_argument(
    '--vb',
    action='store_true',
    default=False,
    help='Set to print out more information'
    )

def assign_loggrid(x, y, xmin=snmin, xmax=snmax, xsteps=steps, ymin=sizemin, ymax=sizemax, ysteps=steps):
    # return x and y indices of data (x,y) on a log-spaced grid that runs from [xy]min to [xy]max in [xy]steps
    x = np.maximum(x, xmin)
    x = np.minimum(x, xmax)

    y = np.maximum(y, ymin)
    y = np.minimum(y, ymax)

    logstepx = np.log10(xmax/xmin)/xsteps
    logstepy = np.log10(ymax/ymin)/ysteps

    indexx = (np.log10(x/xmin)/logstepx).astype(int)
    indexy = (np.log10(y/ymin)/logstepy).astype(int)

    indexx = np.minimum(indexx, xsteps-1)
    indexy = np.minimum(indexy, ysteps-1)

    return indexx, indexy

def mesh_average(quantity, indexx, indexy, steps, count):
    m = np.zeros((steps,steps))
    np.add.at(m,(indexx,indexy),quantity)
    m /= count
    return m

def apply_loggrid(x, y, grid, xmin=snmin, xmax=snmax, xsteps=steps, ymin=sizemin, ymax=sizemax, ysteps=steps):
    indexx,indexy = assign_loggrid(x, y, xmin, xmax, xsteps, ymin, ymax, ysteps)
    res = np.zeros(len(x))
    res = grid[indexx,indexy]
    return res

def main():
    args = parser.parse_args()
    mcal_file = args.mcal_file
    make_plots = args.make_plots
    vb = args.vb
    size_ratio=np.array(f['catalog/metacal/unsheared/size_ratio'])[:]
    if vb is True:
        print('read R')
    R11=np.array(f['catalog/metacal/unsheared/R11'])[:]
    R22=np.array(f['catalog/metacal/unsheared/R22'])[:]
    if vb is True:
        print('read e')
    e1=np.array(f['catalog/metacal/unsheared/e_1'])[:]
    e2=np.array(f['catalog/metacal/unsheared/e_2'])[:]
    f.close()

    # NOTE: The below is taken from Y3 values that we need to reproduce.
    # Can generalize in the future if needed
    max_percentile = 97.5
    np.percentile(snr, max_percentile)
    np.percentile(size_ratio, max_percentile)
    # definitions of grid - span inner 95% of range (i.e., upper limit is 97.5th percentile)
    snmin = 10
    snmax = 275
    sizemin = 0.5
    sizemax = 5
    steps = 20
    save_weight_table = True
    indexx, indexy = assign_loggrid(snr, size_ratio, snmin, snmax, steps, sizemin, sizemax, steps)
    count = np.zeros((steps,steps))
    np.add.at(count,(indexx,indexy), 1)
    meanes = mesh_average(np.sqrt(e1**2+e2**2)/2, indexx, indexy, steps, count)
    response = mesh_average((R11+R22)/2, indexx, indexy, steps, count)
    w = 1/(meanes/response)**2
    if save_weight_table:
            np.savetxt("balrog_y3_shape_w_grid.txt")

    if make_plots is True:
        H, xedges, yedges = np.histogram2d(snr,
                                           size_ratio,
                                           bins=[np.logspace(np.log10(snmin),np.log10(snmax),steps+1),
                                                 np.logspace(np.log10(sizemin),np.log10(sizemax),steps+1)]
                                           )
        X, Y = np.meshgrid(xedges, yedges)
        plt.pcolormesh(Y, X, H.transpose())
        plt.xscale('log')
        plt.yscale('log')
        plt.colorbar(label="number of galaxies in Y3")
        plt.xlabel("metacal size_ratio")
        plt.ylabel("metacal SNR")
        plt.tight_layout()
        plt.savefig("y3_sn_size.png", dpi=300)

        plt.figure()
        logmeshplot(count, xedges, yedges, r"count in Y3")
        plt.tight_layout()
        plt.savefig("y3_count.png", dpi=200)

        plt.figure()
        logmeshplot(response, xedges, yedges, r"$(R11+R22)/2$ in Y3")
        plt.tight_layout()
        plt.savefig("y3_response.png", dpi=200)

        plt.figure()
        logmeshplot(meanes, xedges, yedges, r"$\sqrt{e_1^2+e_2^2}$ in Y3")
        plt.tight_layout()
        plt.savefig("y3_ellipticity_noise.png", dpi=200)

        plt.figure()
        logmeshplot(meanes/response, xedges, yedges, r"$\sqrt{e_1^2+e_2^2}/\langle R\rangle$ in Y3")
        plt.tight_layout()
        plt.savefig("y3_shape_noise.png", dpi=200)

        plt.figure()
        logmeshplot(w, xedges, yedges, r"$w=1/[(e_1^2+e_2^2)/2/\langle R\rangle]$")
        plt.tight_layout()
        plt.savefig("y3_weight.png", dpi=200)

    # Step 2 - assign weight to each galaxy
    w = np.genfromtxt("y3_shape_w_grid.txt")
    weights = apply_loggrid(snr, size_ratio, w)

    # TODO: Save to cols...
    f = h5py.File(mcal_file, 'r+')
    f[wgt_colname] = weights

    return

if __name__ == '__main__':
    main()
