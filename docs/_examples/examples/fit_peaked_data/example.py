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

from sherpa.io import read_data

d = read_data('test_peak.dat')
report("print(d)")

from sherpa.plot import DataPlot
dplot = DataPlot()
dplot.prepare(d)
dplot.plot()

savefig('data.png')

from sherpa.models.basic import NormGauss1D
from sherpa.astro.models import Lorentz1D, PseudoVoigt1D, Voigt1D
from sherpa.stats import LeastSq
from sherpa.optmethods import NelderMead
from sherpa.fit import Fit
models = [NormGauss1D, Lorentz1D, PseudoVoigt1D, Voigt1D]
stat = LeastSq()
method = NelderMead()
fits = {}
for model in models:
    mdl = model()
    # how to call guess? for now skipping
    fit = Fit(d, mdl, stat, method)
    fits[mdl.name] = fit

results = fits['normgauss1d'].fit()
dump("print(results.format())")

results = fits['lorentz1d'].fit()
dump("print(results.format())")

from sherpa.plot import ModelPlot
mplot1 = ModelPlot()
mplot2 = ModelPlot()
mplot1.prepare(d, fits['normgauss1d'].model)
mplot2.prepare(d, fits['lorentz1d'].model)

dplot.plot(alpha=0.5)
ax = plt.gca()
ax.lines[-1].set_label('data')
mplot1.overplot()
ax.lines[-1].set_label('NormGauss1D')
mplot2.overplot()
ax.lines[-1].set_label('Lorentz1D')
plt.legend()

savefig('models-normgauss1d-lorentz1d.png')

results = fits['pseudovoigt1d'].fit()
dump("print(results.format())")

results = fits['voigt1d'].fit()
dump("print(results.format())")

mplot3 = ModelPlot()
mplot4 = ModelPlot()
mplot3.prepare(d, fits['voigt1d'].model)
mplot4.prepare(d, fits['pseudovoigt1d'].model)

from sherpa.plot import SplitPlot

splot = SplitPlot()
splot.addplot(dplot, alpha=0.5)
splot.overlayplot(mplot3)
ax = plt.gca()
ax.set_title('Voigt1D')
ax.set_xlabel('')
splot.addplot(dplot, alpha=0.5)
splot.overlayplot(mplot4)
ax = plt.gca()
ax.set_title('PseudoVoigt1D')

savefig('models-voigts.png')

mplot1.plot()
ax = plt.gca()
ax.lines[-1].set_label('NormGauss1D')
mplot2.overplot()
ax.lines[-1].set_label('Lorentz1D')
mplot3.overplot()
ax.lines[-1].set_label('Voigt1D')
plt.legend()

savefig('models.png')
