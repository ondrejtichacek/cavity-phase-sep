import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from bp import bp_main
from common import timer
from numpy.lib.function_base import append, interp
from optim import Optimizer

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--molecule", default="all", help="")
parser.add_argument("--maxiter", type=int, default=100, help="")
parser.add_argument("--nrepeat", type=int, default=1, help="")
parser.add_argument("--plot", action="store_true", default=False, help="")
parser.add_argument("--test", action="store_true", default=False, help="")
parser.add_argument("--mock", action="store_true", default=False, help="")
parser.add_argument("--mock-data-file", type=str, default=None, help="")
parser.add_argument("--new", action="store_true", default=False, help="")
parser.add_argument("--exp-data-file", type=str, default=None, help="")

args = parser.parse_args()

mpl.rcParams['figure.figsize'] = (7,6)

scale = {
    'amp': 232,
    'adp': 50,
    'atp': 20,
    'poly_arg': 20,
    'Lys10_ADP': 1,
}

scale = {
    'mock': 1,
    'amp': 1,
    'adp': 1,
    'atp': 1,
    'Lys10_ADP': 1,
}

exp_limit = {
    'atp': {
        'x': [0, 1],
        'y': [0, 1],
    },
    'adp': {
        'x': [0, 0.001],
        'y': [0.0005, 1],
    },
    'amp': {
        'x': [0, 0.0012],
        'y': [0, 1],
    },
}

exp_limit_plot = {
    'atp': {
        'x': [0, 1],
        'y': [0, 1],
    },
    'adp': {
        'x': [0, 1],
        'y': [0, 1],
    },
    'amp': {
        'x': [0, 1],
        'y': [0, 1],
    },
}

def main():    

    print(f"We will perform {args.nrepeat} optimizations")

    for i in np.arange(args.nrepeat):
        sel = args.molecule

        opt = Optimizer(
            lattice_connectivity=6,
            sel=sel,
            vol_frac_scaling_x=scale[sel],
            # vol_frac_scaling_y=scale['poly_arg'],
            vol_frac_scaling_y=scale[sel],
        )

        opt.load_exp_data(exp_limit_plot[sel])

        # p ... positive ... ARG
        # m ... negative ... ATP/ADP/AMP
        # 0 ... solvent ... water

        J = 10

        bounds = {
            'lp': (0, 10), # rel. to l0
            'lm': (0, 10), # rel. to l0
            'l0': (1, 1),  # reference
            'Jp': (-J, 0),
            'Jm': (-J, 0),
            'Jpm': (0, J),
            'J0': (-1, 1),
            'J0p': (-1, 1),
            'J0m': (-1, 1),
            'scale_x': (100, 500),
            'scale_y': (100, 500),
            # 'rel_scale_y': (0.1, 10),
        }

        
        opt.optimize_cy(bounds, maxiter=args.maxiter)

def main_new():

    print(f"We will perform {args.nrepeat} optimizations")

    for i in np.arange(args.nrepeat):
        sel = args.molecule

        opt = Optimizer(
            use_cont_err=True,
            lattice_connectivity=6,
            sel=sel,
            vol_frac_scaling_x=scale[sel],
            # vol_frac_scaling_y=scale['poly_arg'],
            vol_frac_scaling_y=scale[sel],
        )

        opt.load_exp_data_new(args.exp_data_file)

        # p ... positive ... ARG
        # m ... negative ... ATP/ADP/AMP
        # 0 ... solvent ... water

        J = 50

        bounds = {
            'lp': (0, 10), # rel. to l0
            'lm': (0, 10), # rel. to l0
            'l0': (1, 1),  # reference
            # ASSOCIATIVE
            # 'Jp': (-J, 0),
            # 'Jm': (-J, 0),
            # 'Jpm': (0, J),
            # SEGREGATIVE
            'Jp': (0, J),
            'Jm': (0, J),
            'Jpm': (-J, 0),
            # FREE
            # 'Jp': (-J, J),
            # 'Jm': (-J, J),
            # 'Jpm': (-J, J),
            #
            # 'J0': (-1, 1),
            # 'J0p': (-1, 1),
            # 'J0m': (-1, 1),
            'J0': (-J, J),
            'J0p': (-J, J),
            'J0m': (-J, J),
            'scale_x': (0, 0.3),
            'scale_y': (0, 0.3),
            # 'rel_scale_y': (0.1, 10),
        }

        
        opt.optimize_cy(bounds, maxiter=args.maxiter)


def main_mock():
    
    print(f"We will perform {args.nrepeat} optimizations")

    for i in np.arange(args.nrepeat):

        sel = args.molecule

        opt = Optimizer(
            lattice_connectivity=4,
            sel=sel,
            vol_frac_scaling_x=1,
            vol_frac_scaling_y=1,
            use_cont_err=True,
        )

        print(args.mock_data_file)

        opt.load_mock_data(args.mock_data_file)

        J1 = 4
        eps = 0.4
        J0 = 0

        bounds = {
            'lp': (1, 1), # rel. to l0
            'lm': (1, 1), # rel. to l0
            'l0': (1, 1),  # reference
            'Jp': (-10, 10),
            'Jm': (-10, 10),
            'Jpm': (-10, 10),
            'J0': (0, 0),
            'J0p': (0, 0),
            'J0m': (0, 0),
            'scale_x': (1, 1),
            'scale_y': (1, 1),
        }

        opt.optimize_cy(bounds, maxiter=args.maxiter)

def test(
    sel, 
    lp = 1,
    lm = 1,
    l0 = 1,
    Jp = -1.3,
    Jm = -2,
    Jpm = 2.2,
    J0 =  0.1,
    J0p = 0.2,
    J0m = 0.6,
    opt_scale_x = 1,
    opt_scale_y = 1,
):

    if args.mock_data_file is None:
        lattice_connectivity = 6

        opt = Optimizer(
            lattice_connectivity=lattice_connectivity,
            sel=sel,
            vol_frac_scaling_x=scale[sel] * opt_scale_x,
            vol_frac_scaling_y=scale[sel] * opt_scale_y,
            use_cont_err=True,
        )

    else:
        lattice_connectivity = 4

        opt = Optimizer(
            lattice_connectivity=lattice_connectivity,
            sel=sel,
            vol_frac_scaling_x=opt_scale_x,
            vol_frac_scaling_y=opt_scale_y,
            use_cont_err=False,
        )

    if args.mock_data_file is not None:
        opt.load_mock_data(args.mock_data_file)
    elif args.exp_data_file is not None:
        opt.load_exp_data_new(args.exp_data_file)
    else:
        opt.load_exp_data(exp_limit_plot[sel], extended=False)
        
    if args.exp_data_file is not None:
        opt.x *= opt_scale_x
        opt.y *= opt_scale_y

    x = opt.x
    y = opt.y
    sep_exp = opt.sep
    sep_exp_cont = opt.sep_cont

    err, sep = opt.bp_wrapper(x, y,
        lp, lm, l0, Jp, Jm, Jpm, J0, J0p, J0m, 
        sep_exp, sep_exp_cont, opt.is_on_boundary)

    print(f"Err: {err}")

    i1 = sep != 0
    i2 = sep == 0

    i3 = sep_exp != 0
    i4 = sep_exp == 0

    plt.figure()
    plt.plot(x[i1 & i3], y[i1 & i3], 's', color='b', alpha=0.3)
    plt.plot(x[i2 & i4], y[i2 & i4], 's', color='#999', alpha=0.3)
    plt.plot(x[i1 & i4], y[i1 & i4], 's', color='r', alpha=0.3)
    plt.plot(x[i2 & i3], y[i2 & i3], 's', color='#222', alpha=0.3)
    plt.xlabel(sel)
    plt.ylabel('polyARG')
    # plt.show()

    plt.figure(figsize=(3.3,3))
    ax = plt.gca()
    z = opt.sep_cont
    X, Y = np.meshgrid(np.unique(x), np.unique(y))
    Z = np.zeros_like(X)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            n = (X[i,j] - x)**2 + (Y[i,j] - y)**2
            k = np.argmin(n)
            Z[i,j] = z[k]
    # ax.pcolormesh(X*100, Y*100, Z, alpha=0.3, shading='gouraud')
    
    plot_vf = True
    plot_vf = False
    if plot_vf:
        plt.scatter(x*100, y*100, c=opt.sep_cont, cmap='viridis', marker='s')
        plt.xlabel('Lys10 (% vf)')
        plt.ylabel('ADP (% vf)')
    else:
        plt.scatter(opt.orig_x, opt.orig_y, c=opt.sep_cont, cmap='viridis', marker='s')
        plt.xlabel('Lys10 (mM)')
        plt.ylabel('ADP (mM)')

    cbar = plt.colorbar()
    cbar.set_ticks([])
    cbar.set_label('sep. intensity (au)')
    plt.tight_layout(pad=2)
    # plt.show()

    max_x = x.max()
    max_y = y.max()

    xp = np.linspace(0.001, 0.75, 128)
    yp = np.linspace(0.001, 0.75, 128)
    X, Y = np.meshgrid(xp, yp)
    X = X.flatten()
    Y = Y.flatten()
    ind = ((X + opt.dx + Y) < 1) & ((X + Y + opt.dy) < 1) & (X > 0)  & (Y > 0)
    print(ind.shape)
    x = np.ascontiguousarray(X[ind])
    y = np.ascontiguousarray(Y[ind])

    n = x.size

    print(n)

    err, sep = opt.bp_wrapper(x, y, lp, lm, l0, Jp, Jm, Jpm, J0, J0p, J0m)

    # print(f"Err: {err}")

    print(min(sep), max(sep))

    Z = -np.ones_like(X)
    Z[ind] = sep
    X, Y = np.meshgrid(xp, yp)
    Z = Z.reshape(X.shape)

    print(X.shape)

    plt.figure()
    plt.pcolormesh(X, Y, Z, shading='auto', cmap='gray_r', vmax=3)

    print(max_x)
    print(max_y)
    xp = np.linspace(0.001, max_x, 128)
    yp = np.linspace(0.001, max_y, 128)
    X, Y = np.meshgrid(xp, yp)
    X = X.flatten()
    Y = Y.flatten()
    ind = ((X + 10*opt.dx + Y) < 1) & ((X + Y + 10*opt.dy) < 1) & (X > 0)  & (Y > 0)
    print(ind.shape)
    x = np.ascontiguousarray(X[ind])
    y = np.ascontiguousarray(Y[ind])

    print(x)

    n = x.size

    err, sep = opt.bp_wrapper(x, y, lp, lm, l0, Jp, Jm, Jpm, J0, J0p, J0m)

    print(sep)

    Z = -np.ones_like(X)
    Z[ind] = sep
    X, Y = np.meshgrid(xp, yp)
    Z = Z.reshape(X.shape)

    plt.figure()
    plt.pcolormesh(X, Y, Z, shading='auto', cmap='gray_r', vmax=3)

    opt.plot_phase_diagram(X, Y, Z)

if __name__ == '__main__':
    if args.test:
        test('adp')

    elif args.mock:
        main_mock()

    elif args.new:
        main_new()

    else:

        if not args.plot:
            main()
        else:
            params = {}
            if args.molecule == "mock_0":
                pass
                par = []
                par.append([1, 1, 1, -1, -1, 3, 0, 0, 0, 1, 1])
                params['mock_0'] = par

            if args.molecule == "mock_1":
                pass
                par = []
                par.append([1, 1, 1, -1, -1, 3, 0, 0, 0, 1, 1])
                params['mock_1'] = par

            if args.molecule == "mock_2":
                pass
                par = []
                par.append([1, 1, 1, 1, 1, -3, 0, 0, 0, 1, 1])
                params['mock_2'] = par

            if args.molecule == "mock_3":
                pass
                par = []
                par.append([1, 1, 1, 0, 2, 0.5, 0, 0, 0, 1, 1])
                params['mock_3'] = par


            if args.molecule == "Lys10_ADP" or args.molecule == "all":
                pass
                par = []
                par.append([4.41604, 2.19046,   1, 0.546494, -3.21525, 3.52499, 0.680959, 0.971192, -0.492456, 0.0669402, 0.0326574])
                params['Lys10_ADP'] = par

            for mol in params:
                par = params[mol]
                for p in par:
                    test(mol, *p)

plt.show()