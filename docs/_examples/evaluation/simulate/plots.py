from pathlib import Path
import sys

from matplotlib import pyplot as plt

import numpy as np

from sherpa.utils.random import poisson_noise

# Assume run from docs directory.
#
indir = Path("../sherpa-test-data/sherpatest")
outdir = Path("_static/evaluation/")

for wanted in [indir, outdir]:
    if not wanted.is_dir():
        sys.stderr.write(f"ERROR: no directory {wanted}\n")
        sys.exit(1)

def savefig(name):
    plt.savefig(f"_static/evaluation/{name}.png")
    print(f"** Created: {name}.png")

data_dir = str(indir)
if not data_dir.endswith("/"):
    data_dir += "/"

#

rng = np.random.default_rng(235)  # a "repeatable" set of random values

from sherpa.models.basic import Polynom1D
mdl = Polynom1D('mdl')
mdl.offset = 35
mdl.c1 = 0.5
mdl.c2 = 0.12
print(mdl)

x = np.arange(10, 100, 12)
ymdl = mdl(x)
plt.plot(x, ymdl, label="model")
plt.legend()

savefig("model")

ypoisson = poisson_noise(ymdl, rng=rng)
plt.clf()
plt.plot(x, ymdl, label="model")
plt.plot(x, ypoisson, ".", label="with noise")
plt.legend()

savefig("model_with_noise")


from sherpa.astro.data import DataPHA
from sherpa.astro.io import read_arf, read_rmf, read_pha
data = DataPHA(name='simulation', channel=None, counts=None, exposure=10000.)
data.set_arf(read_arf(data_dir + '9774.arf'))
data.set_rmf(read_rmf(data_dir + '9774.rmf'))
data.channel = np.arange(1, 1025, dtype=np.int16)

from sherpa.models.basic import PowLaw1D, Gauss1D
pl = PowLaw1D()
line = Gauss1D()
pl.gamma = 1.8
pl.ampl = 2e-05
line.pos = 6.7
line.ampl = .0003
line.fwhm = .1
srcmdl = pl + line
print(srcmdl)

resp = data.get_full_response()
instmdl = resp(srcmdl)
print(instmdl)
from sherpa.astro.fake import fake_pha
from sherpa.astro.plot import DataPHAPlot
fake_pha(data, instmdl)  # NOTE: include the instrument response
dplot = DataPHAPlot()
dplot.prepare(data)
dplot.plot()

savefig("simulated_pha")

data.set_analysis('energy')
data.notice(0.3, 8)
data.group_counts(5)
dplot.prepare(data)
dplot.plot(xlog=True, ylog=True)

savefig("simulated_pha_grouped")

from sherpa.astro.plot import ModelPHAHistogram, ModelHistogram
dplot.plot(xlog=True, ylog=True)
m1plot = ModelPHAHistogram()
m1plot.prepare(data, instmdl)
m1plot.overplot(alpha=0.5, label="Grouped model")
m2plot = ModelHistogram()
m2plot.prepare(data, instmdl)
m2plot.overplot(alpha=0.5, label="Ungrouped model")

savefig("compare_data_and_model")


bgdata = read_pha(data_dir + '9774_bg.pi')
data.set_background(bgdata)
data.backscal = 9.6e-06
fake_pha(data, instmdl, add_bkgs=True)
bplot = DataPHAPlot()
bplot.prepare(data)
dplot.plot(linestyle="solid", ylog=True, label="data")
bplot.overplot(linestyle="solid", alpha=0.5, label="with background")

savefig("simulated_pha_with_background")


# from sherpa.models.basic import Const1D
# bkgmdl = Const1D('bmdl')
# bkgmdl.c0 = 2
# fake_pha(data, mdl, add_bkgs=True, bkg_models={'1': bkgmdl})  # doctest: +SKIP
