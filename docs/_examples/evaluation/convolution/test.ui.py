from matplotlib import pyplot as plt
import numpy as np

from sherpa import ui

x1 = np.arange(-5.0, 30, 0.5)

x2 = np.arange(1.0, 29.0, 0.2)

ui.load_arrays(1, x1, x1 * 0, ui.Data1D)
ui.load_arrays(2, x2, x2 * 0, ui.Data1D)

ui.set_source(1, ui.box1d.box)
box = ui.get_model_component('box')
ui.set_source(2, box)

box.xlow = 10.0
box.xhi = 20.0

# Copy all the objects just to make sure
g1 = ui.gauss1d('g1')
g1.fwhm = 3.0

g2 = ui.gauss1d('g2')
g2.fwhm = 3.0

ui.load_psf('psf1', g1)
ui.load_psf('psf2', g2)

ui.set_psf(1, 'psf1')
ui.set_psf(2, 'psf2')

ui.plot_model(1)
ui.plot_model(2, overplot=True)
