# for fit/index.rst

import numpy as np
import matplotlib.pyplot as plt
from sherpa.data import Data1D
d = Data1D('fit example', [-13, -5, -3, 2, 7, 12],
           [102.3, 16.7, -0.6, -6.7, -9.9, 33.2],
           np.ones(6) * 5)


def report(name):
    print("# print({})".format(name))
    print(eval(name))
    print("----------------------------------------")


def dump(name):
    print("# dump")
    print("{}".format(name))
    print(repr(eval(name)))
    print("----------------------------------------")


def savefig(name):
    plt.savefig(name)
    print("# Created: {}".format(name))


from sherpa.models.basic import Polynom1D
mdl = Polynom1D()

mdl.c2.thaw()
report("mdl")


from sherpa.fit import Fit
f = Fit(d, mdl)
report("f")
report("f.data")
report("f.model")

dump("f.model.c2.val")
mdl.c2 = 1
dump("f.model.c2.val")

print("Starting statistic: {:.3f}".format(f.calc_stat()))
sinfo1 = f.calc_stat_info()
report("sinfo1")

d.ignore(0, 5)
sinfo2 = f.calc_stat_info()
d.notice()
dump("sinfo1.numpoints")
dump("sinfo2.numpoints")

res = f.fit()
if res.succeeded: print("Fit succeeded")
if not res.succeeded: print("**** ERRRR, the fit failed folks")

report("res.format()")
report("res")

from sherpa.plot import DataPlot, ModelPlot
dplot = DataPlot()
dplot.prepare(f.data)
mplot = ModelPlot()
mplot.prepare(f.data, f.model)
dplot.plot()
mplot.overplot()

savefig("data_model_c0_c2.png")

dump("f.method.name")
original_method = f.method

from sherpa.optmethods import NelderMead
f.method = NelderMead()
resn = f.fit()
print("Change in statistic: {}".format(resn.dstatval))

fit2 = Fit(d, mdl, method=NelderMead())
fit2.fit()

mdl.c1.thaw()
f.method = original_method

res2 = f.fit()
report("res2.format()")


from sherpa.plot import DelchiPlot, FitPlot, SplitPlot
fplot = FitPlot()
rplot = DelchiPlot()
splot = SplitPlot()
mplot.prepare(f.data, f.model)
fplot.prepare(dplot, mplot)
splot.addplot(fplot)
rplot.prepare(f.data, f.model, f.stat)
splot.addplot(rplot)

savefig("fit_delchi_c0_c1_c2.png")

mdl.reset()

report("[(p.name, p.val, p.frozen) for p in mdl.pars[:3]]")

for p in mdl.pars[:3]:
    p.val = 10

report("mdl.pars[1]")

res3 = f.fit(outfile='fitpath.txt', clobber=True)
report("res3.format()")

print("*** cat fitpath.txt")
with open('fitpath.txt', 'r') as fh:
    for l in fh.readlines():
        print(l.strip())

print("----------------------------------------")

report("res3.statval == res2.statval")
report("res3.statval - res2.statval")

report("res2.parvals")
report("res3.parvals")

for p2, p3 in zip(res2.parvals, res3.parvals):
    print("{:+.2e}".format(p3 - p2))

from sherpa.utils import calc_mlr
report("calc_mlr(res.dof - res2.dof, res.statval - res2.statval)")

report("f.estmethod.name")

coverrs = f.est_errors()
report("coverrs.format()")

report("coverrs")

dump("f.estmethod.sigma")
f.estmethod.sigma = 1.6
coverrs90 = f.est_errors()
report("coverrs90.format()")

dump("coverrs90.percent")

report("coverrs.extra_output")

print([p.split('.')[1] for p in coverrs.parnames])


from sherpa.estmethods import Confidence
f.estmethod = Confidence()
print("*** about to start confidence run")
conferrs = f.est_errors()
print("*** finished confidence run")

report("conferrs.format()")

dump("f.estmethod.name")
dump("f.estmethod.parallel")


print("*** about to start confidence run")
c1errs = f.est_errors(parlist=(mdl.c1, ))
print("*** finished confidence run")
report("c1errs")

from sherpa.plot import IntervalProjection
iproj = IntervalProjection()
iproj.calc(f, mdl.c1)
iproj.plot()
savefig("iproj_c1_auto.png")

iproj.prepare(fac=5, nloop=51)
iproj.calc(f, mdl.c1)
iproj.plot()
pmin = c1errs.parvals[0] + c1errs.parmins[0]
pmax = c1errs.parvals[0] + c1errs.parmaxes[0]
iproj.vline(pmin, overplot=True, linestyle='dot')
iproj.vline(pmax, overplot=True, linestyle='dot')
savefig("iproj_c1_manual.png")


from sherpa.plot import RegionProjection
rproj = RegionProjection()
rproj.calc(f, mdl.c0, mdl.c2)
rproj.contour()
savefig("rproj_c0_c2_auto.png")

rproj.prepare(min=(-22, 0.35), max=(3, 0.6), nloop=(41, 41))
rproj.calc(f, mdl.c0, mdl.c2)
rproj.contour()
xlo, xhi = plt.xlim()
ylo, yhi = plt.ylim()

def get_limits(i):
    return conferrs.parvals[i] + \
           np.asarray([conferrs.parmins[i],
                       conferrs.parmaxes[i]])


plt.vlines(get_limits(0), ymin=ylo, ymax=yhi)
plt.hlines(get_limits(2), xmin=xlo, xmax=xhi)

# TODO: is this better with
#   plt.axhline(val) ?

savefig("rproj_c0_c2_manual.png")


xmin, xmax = rproj.x0.min(), rproj.x0.max()
ymin, ymax = rproj.x1.min(), rproj.x1.max()
nx, ny = rproj.nloop
hx = 0.5 * (xmax - xmin) / (nx - 1)
hy = 0.5 * (ymax - ymin) / (ny - 1)
extent = (xmin - hx, xmax + hx, ymin - hy, ymax + hy)

y = rproj.y.reshape((ny, nx))

plt.clf()
# plt.imshow(y, origin='lower', extent=extent, aspect='auto', cmap='cubehelix_r')
plt.imshow(y, origin='lower', extent=extent, aspect='auto', cmap='viridis_r')
plt.colorbar()
plt.xlabel(rproj.xlabel)
plt.ylabel(rproj.ylabel)
rproj.contour(overplot=True)

savefig("rproj_c0_c2_image.png")
