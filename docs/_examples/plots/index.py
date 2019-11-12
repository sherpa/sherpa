# for plot/index.rst

import numpy as np
import matplotlib.pyplot as plt

def savefig(name):
    plt.savefig(name)
    print("# Created: {}".format(name))

def report(name):
    print("# print({})".format(name))
    print(eval(name))
    print("----------------------------------------")


def dump(name):
    print("# dump")
    print("{}".format(name))
    print(repr(eval(name)))
    print("----------------------------------------")


edges = np.asarray([-10, -5, 5, 12, 17, 20, 30, 56, 60])
y = np.asarray([28, 62, 17, 4, 2, 4, 125, 55])

from sherpa.data import Data1DInt
d = Data1DInt('example histogram', edges[:-1], edges[1:], y)

from sherpa.plot import DataPlot
dplot = DataPlot()
dplot.prepare(d)
dplot.plot()

savefig('dataplot_histogram.png')

from sherpa.plot import Histogram
hplot = Histogram()
hplot.overplot(d.xlo, d.xhi, d.y)

savefig('dataplot_histogram_overplot.png')

from sherpa.models.basic import Const1D, Gauss1D
mdl = Const1D('base') - Gauss1D('line')
mdl.pars[0].val = 10
mdl.pars[1].val = 25
mdl.pars[2].val = 22
mdl.pars[3].val = 10
report("mdl")

from sherpa.plot import ModelPlot
mplot = ModelPlot()
mplot.prepare(d, mdl)
mplot.plot()
dplot.overplot()

savefig('modelplot_histogram_overplot.png')

from sherpa.plot import FitPlot
fplot = FitPlot()
fplot.prepare(dplot, mplot)
fplot.plot()

savefig('fitplot_histogram.png')

# settings

report("dplot.plot_prefs")

dplot.plot_prefs['ylog'] = True
dplot.plot(marker='s', linestyle='dashed')

savefig('settings_dataplot_combined.png')

dplot.plot()
savefig('settings_dataplot_ylog.png')

dplot.plot_prefs['ylog'] = False


from sherpa.optmethods import NelderMead
from sherpa.stats import Cash
from sherpa.fit import Fit
f = Fit(d, mdl, stat=Cash(), method=NelderMead())
out = f.fit()

# How bad a fit is this?
print(out)
# f.est_errors()

mplot.prepare(d, mdl)
fplot.plot()

savefig('fitplot_histogram_after.png')

from sherpa.plot import IntervalProjection
iproj = IntervalProjection()
iproj.calc(f, mdl.pars[2])

iproj.plot()

savefig('intproj_histogram_pos.png')
