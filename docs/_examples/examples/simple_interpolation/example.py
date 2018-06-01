import numpy as np
import matplotlib.pyplot as plt

def savefig(name):
    plt.savefig(name)
    print("Created: {}".format(name))

def report(name):
    print("# print({})".format(name))
    print(eval(name))

def dump(name):
    print("# dump")
    print(name)
    print(repr(eval(name)))


from sherpa.data import Data1D
x = [1, 1.5, 2, 4, 8, 17]
y = [1, 1.5, 1.75, 3.25, 6, 16]
d = Data1D('interpolation', x, y)

report("print(d)")

from sherpa.plot import DataPlot
dplot = DataPlot()
dplot.prepare(d)
dplot.plot()
savefig('data.png')

# Note: can not print(dplot) as there is a problem with the fact
#       the input to the data object is a list, not ndarray
#       Sherpa 4.10.0

from sherpa.models.basic import Polynom1D
mdl = Polynom1D()
report("print(mdl)")

mdl.c2.thaw()

from sherpa.plot import ModelPlot
mplot = ModelPlot()
mplot.prepare(d, mdl)
dplot.plot()
mplot.overplot()
savefig("data_model_initial.png")

from sherpa.stats import LeastSq
from sherpa.optmethods import NelderMead
from sherpa.fit import Fit
f = Fit(d, mdl, stat=LeastSq(), method=NelderMead())
report("print(f)")

res = f.fit()
dump("res.succeeded")

report("res.format()")
report("res")

report("mdl")

stat2 = f.calc_stat()
print("Statistic = {:.4f}".format(stat2))


from sherpa.plot import FitPlot, ResidPlot, SplitPlot
fplot = FitPlot()
mplot.prepare(d, mdl)
fplot.prepare(dplot, mplot)
splot = SplitPlot()
splot.addplot(fplot)
rplot = ResidPlot()
print("### should get a WARNING about missing errors")
rplot.prepare(d, mdl, stat=LeastSq())

rplot.plot_prefs['yerrorbars'] = False
splot.addplot(rplot)

savefig("data_model_resid.png")


report("mdl([2, 5, 10])")
report("mdl([-100])")
report("mdl([234.56])")


mdl.c1.thaw()
mdl.c2 = 0
mdl.c2.freeze()
dump("f.fit()")

report("mdl")

stat1 = f.calc_stat()
print("Statistic: order 1 = {:.3f} order 2 = {:.3f}".format(stat1, stat2))


mplot2 = ModelPlot()
mplot2.prepare(d, mdl)
mplot.plot()
mplot2.overplot()
savefig("model_comparison.png")

xgrid = np.linspace(0, 20, 21)
y1 = mdl(xgrid)
mdl.c0 = res.parvals[0]
mdl.c1 = 0
mdl.c2 = res.parvals[1]
y2 = mdl(xgrid)
plt.clf()
plt.plot(xgrid, y2, label='order=2')
plt.plot(xgrid, y1, label='order=1')
plt.legend();
plt.title("Manual evaluation of the models")
savefig("model_comparison_manual.png")
