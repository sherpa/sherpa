# for evaluation/examples#model_evaluate_example_pha_directly

import os

import numpy as np
import matplotlib.pyplot as plt


# as we jump to the data directory need to know the working directory so
# can write the PNG files there
#
savedir = os.getcwd()


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
    outfile = os.path.join(savedir, name)
    plt.savefig(outfile)
    print("# Created: {}".format(name))


os.chdir('../../../../sherpa-test-data/sherpatest/')

from sherpa.astro.io import read_pha

pha = read_pha('9774.pi')

dump("pha")
dump("pha.get_background()")
dump("pha.get_arf()")
dump("pha.get_rmf()")

dump("pha.header['INSTRUME']")
dump("pha.header['DETNAM']")
dump("pha.channel.size")

pha.set_analysis('energy')
pha.notice(0.3, 7)
tabs = ~pha.mask
pha.group_counts(20, tabStops=tabs)

from sherpa.plot import DataPlot
dplot = DataPlot()
dplot.prepare(pha)
dplot.plot(xlog=True, ylog=True)
savefig('pha_data.png')

chans, = pha.get_indep(filter=True)
counts = pha.get_dep(filter=True)
dump("chans.size, counts.size")

gchans = pha.apply_filter(chans, pha._middle)
dump("gchans.size")

plt.clf()
plt.plot(gchans, counts, 'o')
plt.xlabel('Channel')
plt.ylabel('Counts')
savefig('pha_data_manual.png')

x = pha.get_x()
dump("x.min(), x.max()")
x = pha.apply_filter(x, pha._middle)
y = pha.get_y(filter=True)
dplot.plot(xlog=True, ylog=True)
plt.plot(x, y)
savefig('pha_data_compare.png')

pha.set_analysis('wave')
dump("pha.get_x().max()")
wplot = DataPlot()
wplot.prepare(pha)
wplot.plot()
savefig('pha_data_wave.png')
pha.set_analysis('energy')

from sherpa.models.basic import PowLaw1D
from sherpa.astro.xspec import XSphabs
pl = PowLaw1D()
gal = XSphabs()
mdl = gal * pl
pl.gamma = 1.7
gal.nh = 0.2

report("mdl")

egrid = np.arange(0.1, 10, 0.01)
elo, ehi = egrid[:-1], egrid[1:]
emid = (elo + ehi) / 2

plt.clf()
plt.plot(emid, mdl(elo, ehi), label='Absorbed')
plt.plot(emid, pl(elo, ehi), ':', label='Unabsorbed')
plt.xscale('log')
plt.ylim(0, 0.01)
plt.legend()
savefig('pha_model_energy.png')

from sherpa.astro.instrument import Response1D
rsp = Response1D(pha)
full = rsp(mdl)

report("full")

dump("elo.size")
dump("full(elo, ehi).size")
dump("full([1, 2, 3]).size")
dump("np.all(full(elo, ehi) == full([1, 2, 3]))")

plt.clf()
plt.plot(pha.channel, full(pha.channel))
plt.xlabel('Channel')
plt.ylabel('Counts')

savefig('pha_fullmodel_manual.png')

y1 = pha.eval_model(full)
y2 = pha.eval_model_to_fit(full)
dump("y1.size, y2.size")

rmf = pha.get_rmf()
dump("rmf.e_min.size, rmf.e_max.size")

xlo = pha.apply_filter(rmf.e_min, pha._min)
xhi = pha.apply_filter(rmf.e_max, pha._max)

x2 = pha.get_x()
xmid = pha.apply_filter(x2, pha._middle)
plt.clf()
plt.plot(xmid, y2 / (xhi - xlo) / pha.exposure)

plt.xlabel('Energy (keV)')
plt.ylabel('Counts/sec/keV')

savefig('pha_eval_model_to_fit.png')

from sherpa.astro.plot import ModelHistogram
mplot = ModelHistogram()
mplot.prepare(pha, full)
mplot.plot()

savefig('pha_fullmodel_model.png')

from sherpa.fit import Fit
fit = Fit(pha, full)
res = fit.fit()

report('res.format()')

from sherpa.plot import ModelPlot

dplot.prepare(pha)
dplot.plot(xlog=True)

mplot2 = ModelPlot()
mplot2.prepare(pha, full)
mplot2.overplot()

savefig('pha_fullmodel_fit.png')
