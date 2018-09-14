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


from openpyxl import load_workbook
wb = load_workbook('pone.0171996.s001.xlsx')
fig4 = wb['Fig4data']
t = []; y = []; dy = []
for r in list(fig4.values)[2:]:
    t.append(r[0])
    y.append(r[3])
    dy.append(r[4])

from sherpa.data import Data1D
d = Data1D('NaNO_3', t, y, dy)

from sherpa.plot import DataPlot
dplot = DataPlot()
dplot.prepare(d)
dplot.plot()
savefig("data.png")

report("d")
dump("d.get_filter()")

d.ignore(None, 1)
dump("d.get_filter()")

dump("d.get_filter(format='%d')")

dplot.prepare(d)

from sherpa.models.basic import Const1D, Exp
plateau = Const1D('plateau')
rise = Exp('rise')
mdl = plateau / (1 + rise)
report("mdl")

rise.ampl.freeze()
report("mdl")

from sherpa.plot import ModelPlot
mplot = ModelPlot()
mplot.prepare(d, mdl)
plt.subplot(2, 1, 1)
mplot.plot(clearwindow=False)
plt.subplot(2, 1, 2)
dplot.plot(clearwindow=False)
plt.title('')

savefig("model_data_before_fit.png")

from sherpa.stats import Chi2
from sherpa.fit import Fit
f = Fit(d, mdl, stat=Chi2())
report("f")

print("Starting statistic: {}".format(f.calc_stat()))

fitres = f.fit()
report("fitres.format()")

print("Reduced chi square = {:.2f}".format(fitres.rstat))

mplot.prepare(d, mdl)
dplot.plot()
mplot.overplot()
savefig("model_data_fit1.png")

from sherpa.optmethods import NelderMead
f.method = NelderMead()
fitres2 = f.fit()
report("mdl")

dump("fitres2.dstatval")

mdl.reset()
report("mdl")

plateau.c0 = np.max(d.y)
mplot.prepare(d, mdl)
dplot.plot()
mplot.overplot()
savefig("model_data_reset.png")

fitres3 = f.fit()
report("fitres3.format()")

mplot.prepare(d, mdl)
dplot.plot()
mplot.overplot()
savefig("model_data_fit2.png")

from sherpa.plot import DelchiPlot
residplot = DelchiPlot()
residplot.prepare(d, mdl, f.stat)
residplot.plot()
savefig("model_data_delchi.png")

d.notice()
dump("d.get_filter(format='%d')")

from sherpa.plot import FitPlot
fitplot = FitPlot()
dplot.prepare(d)
mplot.prepare(d, mdl)
fitplot.prepare(dplot, mplot)
fitplot.plot()
savefig("model_data_fit_all.png")

# do we get an error? Actually, it looks to not be the divide-by-zero
# being the problem but list/list instead:
#
"""
residplot.prepare(d, mdl, f.stat)

/home/djburke/miniconda2/envs/sherpa410-py35/lib/python3.5/site-packages/sherpa-4.10.0-py3.5-linux-x86_64.egg/sherpa/plot/__init__.py:1128: RuntimeWarning: divide by zero encountered in true_divide
  return (ylist[0] - ylist[1]) / staterr
Traceback (most recent call last):
  File "example.py", line 125, in <module>
    residplot.prepare(d, mdl, f.stat)
  File "/home/djburke/miniconda2/envs/sherpa410-py35/lib/python3.5/site-packages/sherpa-4.10.0-py3.5-linux-x86_64.egg/sherpa/plot/__init__.py", line 1140, in prepare
    self.yerr = staterr / staterr
TypeError: unsupported operand type(s) for /: 'list' and 'list'
"""

d.ignore(None, 1)

statinfo = f.calc_stat_info()
report("statinfo")

dump("statinfo.rstat == fitres3.rstat")
dump("f.estmethod.name")

coverrs = f.est_errors()
report("coverrs.format()")

dump("f.estmethod.sigma")

f.estmethod.sigma = 1.6
coverrs90 = f.est_errors()
report("coverrs90.format()")

from sherpa.estmethods import Confidence
f.estmethod = Confidence()
print("*** start confidence errors")
conferrs = f.est_errors()
print("*** end   confidence errors")

report("conferrs.format()")

print("*** start confidence errors")
offseterrs = f.est_errors(parlist=(mdl.pars[1], ))
print("*** end   confidence errors")
report("offseterrs")

fmt = "{:13s} covar=±{:4.2f}  conf={:+5.2f} {:+5.2f}"
for i in range(len(conferrs.parnames)):
    print(fmt.format(conferrs.parnames[i], coverrs.parmaxes[i],
                     conferrs.parmins[i], conferrs.parmaxes[i]))

from sherpa.plot import IntervalProjection
intproj = IntervalProjection()
intproj.calc(f, plateau.c0)
intproj.plot()
savefig("intproj_c0_auto.png")

intproj.prepare(min=12.5, max=20, nloop=51)
intproj.calc(f, plateau.c0)
intproj.plot()
s0 = f.calc_stat()
for ds in [1, 4, 9]:
    intproj.hline(s0 + ds, overplot=True, linestyle='dot', linecolor='gray')

savefig("intproj_c0_manual.png")

from sherpa.plot import RegionProjection
regproj = RegionProjection()
regproj.calc(f, rise.offset, rise.coeff)
regproj.contour()
savefig("regproj_offset_coeff_auto.png")

regproj.prepare(min=(2, -1.2), max=(8, -0.1), nloop=(21, 21))
regproj.calc(f, rise.offset, rise.coeff)
regproj.contour()
savefig("regproj_offset_coeff_manual.png")

from sherpa.models.basic import ArithmeticModel
from sherpa.models.parameter import Parameter


class MyExp(ArithmeticModel):
    """A simpler form of the Exp model.

    The model is f(x) = exp(a + b * x).
    """

    def __init__(self, name='myexp'):

        self.a = Parameter(name, 'a', 0)
        self.b = Parameter(name, 'b', -1)

        # The _exp instance is used to perform the model calculation,
        # as shown in the calc method.
        self._exp = Exp('hidden')

        return ArithmeticModel.__init__(self, name, (self.a, self.b))

    def calc(self, pars, *args, **kwargs):
        """Calculate the model"""

        # Tell the exp model to evaluate the model, after converting
        # the parameter values to the required form, and order, of:
        # offset, coeff, ampl.
        #
        coeff = pars[1]
        offset = -1 * pars[0] / coeff
        ampl = 1.0
        return self._exp.calc([offset, coeff, ampl], *args, **kwargs)


plateau2 = Const1D('plateau2')
rise2 = MyExp('rise2')
mdl2 = plateau2 / (1 + rise2)
report("mdl2")

fit2 = Fit(d, mdl2, stat=Chi2())
res2 = fit2.fit()
report("res2.format()")

dplot.prepare(d)
mplot2 = ModelPlot()
mplot2.prepare(d, mdl2)
dplot.plot()
mplot2.overplot()
savefig("model_data_myexp.png")


fit2.estmethod = Confidence()
print("*** start confidence errors")
conferrs2 = fit2.est_errors()
print("*** end   confidence errors")
report("conferrs2.format()")


regproj2 = RegionProjection()
regproj2.prepare(min=(0.5, -1.2), max=(5, -0.1), nloop=(21, 21))
regproj2.calc(fit2, rise2.a, rise2.b)
regproj2.contour()
plt.plot(1.941, -0.453, 'ko', label='NaNO$_3$ Table 5')
plt.legend(loc=1)
savefig("regproj_a_b_manual.png")
