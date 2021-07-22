# for data/index.rst

import numpy as np
import matplotlib.pyplot as plt
from sherpa.stats import LeastSq
from sherpa.optmethods import LevMar
from sherpa import data

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


np.random.seed(0)
x = np.arange(20, 40, 0.5)
y = x**2 + np.random.normal(0, 10, size=x.size)
d1 = data.Data1D('test', x, y)
report("d1")

plt.plot(d1.x, d1.y, 'o')
savefig("data1d.png")


z = np.random.gamma(20, scale=0.5, size=1000)
(y, edges) = np.histogram(z)
d2 = data.Data1DInt('gamma', edges[:-1], edges[1:], y)
report("d2")

plt.clf()
plt.bar(d2.xlo, d2.y, d2.xhi - d2.xlo, align='edge')
savefig("data1dint.png")

d1.ignore(21.2, 22.8)
dump("d1.x[np.invert(d1.mask)]")

report('d1.get_filter()')

from sherpa.models import Polynom1D
from sherpa.fit import Fit
mdl = Polynom1D()
mdl.c2.thaw()
fit = Fit(d1, mdl, stat=LeastSq(), method=LevMar())
res1 = fit.fit()

d1.notice()
res2 = fit.fit()

report('"Degrees of freedom: {} vs {}".format(res1.dof, res2.dof)')

d1.notice()
report("d1.get_filter(format='%.1f')")
d1.notice(25, 27)
report("d1.get_filter(format='%.1f')")
d1.notice(30, 35)
report("d1.get_filter(format='%.1f')")

d1.notice()
report('d1.get_x()')
report('d1.get_y()')

d1.notice(21.1, 23.5)

report('d1.get_x(filter=True)')
report('d1.get_y(filter=True)')

from sherpa.plot import DataHistogramPlot
pdata = DataHistogramPlot()
pdata.prepare(d2)
pdata.plot()

savefig('data_int_to_plot.png')

pdata.plot(linestyle='solid', marker='')
savefig('data_int_line_to_plot.png')

from sherpa.plot import DataPlot
pdata = DataPlot()
d1.notice()
d1.ignore(25, 30)
d1.notice(26, 27)
pdata.prepare(d1)
pdata.plot()

savefig('data_to_plot.png')

d1.notice()

d1.notice(22, 25)
y1 = d1.eval_model(mdl)
y2 = d1.eval_model_to_fit(mdl)
x2 = d1.x[d1.mask]
plt.plot(d1.x, d1.y, 'ko', label='Data')
plt.plot(d1.x, y1, '--', label='Model (all points)')
plt.plot(x2, y2, linewidth=2, label='Model (filtered)')
plt.legend(loc=2)

savefig('data_eval_model.png')
