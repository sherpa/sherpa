from matplotlib import pyplot as plt

from sherpa.ui.utils import Session
from sherpa.data import Data1DInt
from sherpa.models.basic import Polynom1D

def savefig(name):
    plt.savefig(name)
    print("# Created: {}".format(name))


s = Session()
xlo = [2, 3, 5, 7, 8]
xhi = [3, 5, 6, 8, 9]
y = [10, 27, 14, 10, 14]
s.load_arrays(1, xlo, xhi, y, Data1DInt)
mdl = Polynom1D('mdl')
mdl.c0 = 6
mdl.c1 = 1
s.set_source(mdl)
s.plot_fit()

savefig('ui_plot_fit_basic.png')

s.plot_data(color='black')
p = s.get_model_plot_prefs()
p['marker'] = '*'
p['markerfacecolor'] = 'green'
p['markersize'] = 12
s.plot_model(linestyle=':', alpha=0.7, overplot=True)

savefig('ui_plot_fit_manual.png')

print("**** s.get_model_plot()")
plot = s.get_model_plot(recalc=False)
print(type(plot))
print(plot)
