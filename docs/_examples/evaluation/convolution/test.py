from matplotlib import pyplot as plt
import numpy as np

from sherpa import models, instrument
from sherpa.data import Data1D

# Copy all the objects just to make sure
g1 = models.Gauss1D('g1')
g1.fwhm = 3.0

g2 = models.Gauss1D('g2')
g2.fwhm = 3.0

psf1 = instrument.PSFModel('psf1', g1)
psf2 = instrument.PSFModel('psf2', g2)

mdl_orig = models.Box1D('box')
mdl_orig.xlow = 10.0
mdl_orig.xhi = 20.0

mdl_conv1 = psf1(mdl_orig)
mdl_conv2 = psf2(mdl_orig)

x1 = np.arange(-5.0, 30, 0.5)
d1 = Data1D('grid1', x1, x1 * 0)

x2 = np.arange(1.0, 29.0, 0.2)
d2 = Data1D('grid2', x2, x2 * 0)

y_orig = mdl_orig(x1)

psf1.fold(d1)
y_conv1 = mdl_conv1(x1)

psf2.fold(d2)
y_conv2 = mdl_conv2(x2)

plt.plot(x1, y_conv1, '.-', label='x1')
plt.plot(x2, y_conv2, label='x2')
plt.plot(x1, y_orig, '.', label='original', alpha=0.6)
plt.legend(loc='upper left')
plt.xlim(-6, 30)
plt.ylim(-0.2, 1.2)
