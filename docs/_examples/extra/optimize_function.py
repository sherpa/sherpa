#!/usr/bin/env python

# PNG files to go to
#   ../../_static/extra/
#
# an example of optimizing a function

import numpy as np

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import colors

from astropy.table import Table

from sherpa.optmethods import optfcts


def savefig(outfile):
    plt.savefig(outfile)
    print(f'Created: {outfile}')

a = 2
b = 10

def rosenbrock(x, y):
    return (a - x) * (a - x) + b * (y - x * x) * (y - x * x)


y, x = np.mgrid[-2:5.1:0.1, -3:3.1:0.1]
surface = rosenbrock(x, y)
print(surface.shape)
print(surface.min())
print(surface.max())

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#surf = ax.plot_surface(x, y, surface)
surf = ax.plot_surface(x, y, surface,
                       # cmap=cm.YlOrRd_r,
                       alpha=0.8,
                       cmap=cm.plasma_r,
                       vmin=0, vmax=300
                       )

# ax.scatter([2], [4], zs=0, c='k', label='Minimum')
# ax.legend()

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.view_init(30, -110)

savefig('rosenbrock-surface.png')

# How do we optimize this function?
#
def to_optimize(args):
    x = args[0]
    y = args[1]
    pa = (a - x) * (a - b)
    pb = b * (y - x * x) * (y - x * x)
    stat = pa * pa + pb * pb
    return stat, None


start = [-1.2, 1]
lo = [-100, -100]
hi = [100, 100]
res = optfcts.minim(to_optimize, start, lo, hi)

print(f"Success: {res[0]}")
print(f"Message: {res[3]}")
print(f"extra:   {res[4]}")

print(f"best-fit location: {res[1]}")
print(f"          minimum: {res[2]}")


tbl = Table(names=['method', 'stat0', 'x', 'y'],
            dtype=[str, float, float, float])

for method in [optfcts.minim, optfcts.neldermead, optfcts.lmdif, optfcts.montecarlo]:
    res = method(to_optimize, start, lo, hi)
    if res[0]:
        tbl.add_row([method.__name__, res[2], res[1][0], res[1][1]])
    else:
        print(f"Failed {method.__name__}: {res[3]}")

print(tbl)

def to_optimize2(args):
    x = args[0]
    y = args[1]
    pa = (a - x) * (a - b)
    pb = b * (y - x * x) * (y - x * x)
    stat = pa * pa + pb * pb
    return stat, [pa, pb]

res2 = optfcts.lmdif(to_optimize2, start, lo, hi)

print(res2[0])
print(res2[1])

print(res2[4])


# Now try to optimize a more-realistic case, using the data from
#
# https://sherpa.readthedocs.io/en/latest/examples/fit_peaked_data.html


d = np.loadtxt('../examples/fit_peaked_data/test_peak.dat')
x = d[:, 0]
y = d[:, 1]

def ngauss(x, ampl, pos, fwhm):
    term = 4 * np.log(2)
    numerator = ampl * np.exp(-term * (x - pos) * (x - pos) / (fwhm * fwhm))
    denominator = np.sqrt(np.pi / term) * fwhm
    return numerator / denominator

def cb(pars):
    model = ngauss(x, pars[0], pars[1], pars[2])
    delta = model - y
    statval = (delta * delta).sum()
    # Keep a record of the parameters we've visited
    store.add_row([statval] + list(pars))
    return statval, delta


start = [y.max(), (x[0] + x[-1]) / 2, (x[-1] - x[0]) / 10]
lows = [0, x[0], 0]
his = [1e4, x[-1], x[-1] - x[0]]

store = Table(names=['stat', 'ampl', 'pos', 'fwhm'])
flag, bestfit, statval, msg, opts = optfcts.neldermead(cb, start, lows, his, ftol=1e-4)
print(flag)
print(bestfit)
print(statval)
print(opts)
print(len(store))

plt.clf()
plt.plot(x, y, label='Data', alpha=0.5)
plt.plot(x, ngauss(x, *start), label='Start')
plt.plot(x, ngauss(x, *bestfit), label='Best fit', c='k')
plt.legend()

savefig('normgauss1d-example.png')

print(store)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

vmin = store['stat'].min()
vmax = store['stat'].max()
norm = colors.LogNorm(vmin=vmin, vmax=vmax)

ax.plot(store['ampl'], store['pos'], store['fwhm'], alpha=0.4)
scatter = ax.scatter(store['ampl'], store['pos'], store['fwhm'],
                     c=store['stat'], norm=norm)

ax.set_xlabel('ampl')
ax.set_ylabel('pos')
ax.set_zlabel('fwhm')

cbar = fig.colorbar(scatter, shrink=0.5, orientation='horizontal')
cbar.set_label('least-squares statistic')

savefig('normgauss1d-trail.png')
