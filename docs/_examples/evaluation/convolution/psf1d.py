# This isn't currently used

from matplotlib import pyplot as plt
import numpy as np

from sherpa import models, instrument
from sherpa.data import Data1D
from sherpa.utils.err import PSFErr

# For now make sure the convolved models are 0 at both edges of the
# grid, to avoid having to explain the current behavior, which I
# don't understand.
#
g3 = models.Gauss1D('g3')
g3.fwhm = 3.0

# g6 = models.Gauss1D('g6')
# g6.fwhm = 6.0

psf3 = instrument.PSFModel('psf3', g3)
# psf6 = instrument.PSFModel('psf6', g6)

mdl_orig = models.Box1D('box')
mdl_orig.xlow = 10.0
mdl_orig.xhi = 20.0

mdl_conv3 = psf3(mdl_orig)
# mdl_conv6 = psf6(mdl_orig)

print("# Model: PSF=3")
print(mdl_conv3)

# print("# Model: PSF=6")
# print(mdl_conv6)

print("# Kernel: PSF=3")
print(psf3)

# print("# Kernel: PSF=6")
# print(psf6)

print("----------------")

# x = np.arange(2.0, 23, 0.5)
x = np.arange(-5.0, 30, 0.5)
# x = np.arange(1.0, 29.0, 0.2)
d = Data1D('grid', x, x * 0)

y_orig = mdl_orig(x)

# This will fail
try:
    y_conv3 = mdl_conv3(x)
except PSFErr as exc:
    print("** failed with message: {}".format(exc))

# Need this
psf3.fold(d)
# psf6.fold(d)

# Qus: does the grid have to match d.x here?
y_conv3 = mdl_conv3(x)
# y_conv6 = mdl_conv6(x)

plt.plot(x, y_orig, label='original')
plt.plot(x, y_conv3, label='FWHM=3')
# plt.plot(x, y_conv6, label='FWHM=6')
plt.legend(loc='upper left')
plt.xlim(0, 30)
plt.ylim(-0.2, 1.2)
