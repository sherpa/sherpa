import numpy as np
import matplotlib.pyplot as plt
from sherpa import data, models, stats, optmethods, fit, plot

np.random.seed(42)

s1 = models.Gauss1D('sim1')
s2 = models.Gauss1D('sim2')

print(s1)
print([p.name for p in s1.pars])

sim_model = s1 + s2

print(sim_model)
print([p.fullname for p in sim_model.pars])

s1.ampl = 1.0
s1.pos = 0.0
s1.fwhm = 0.5
s2.ampl = 2.5
s2.pos = 0.5
s2.fwhm = 0.25

x = np.linspace(-1, 1, 200)
y = sim_model(x) + np.random.normal(0., 0.2, x.shape)

d = data.Data1D('simulated', x, y)
dplot = plot.DataPlot()
dplot.prepare(d)
print(">>> dplot.plot()")
print(">>> save to _static/models/combine/model_combine_data.png")

print(sim_model)
print(sim_model.op)
print(repr(sim_model.lhs))
print(repr(sim_model.rhs))
print(sim_model.parts)
for cpt in sim_model.parts:
    print(cpt)

print(sim_model([-1.0, 0, 1]))

g1 = models.Gauss1D('g1')
g2 = models.Gauss1D('g2')
mdl = g1 + g2

g2.pos = g1.pos + 0.5
g1.fwhm = 0.1
g2.fwhm = 0.1

print(mdl)

ystart = mdl(x)

mplot = plot.ModelPlot()
mplot.prepare(d, mdl)

print(">>> dplot.plot()")
print(">>> mplot.plot(overplot=True)")
print(">>> save to _static/models/combine/model_combine_start.png")

f = fit.Fit(d, mdl, stats.LeastSq())
res = f.fit()
print(res.succeeded)

fplot = plot.FitPlot()
mplot.prepare(d, mdl)
fplot.prepare(dplot, mplot)
print(">>> fplot.plot()")
print(">>> plt.plot(x, ystart, label='Start')")
print(">>> plt.legend(loc=2)")
print(">>> save as docs/_static/models/combine/model_combine.png")

print(mdl)

for p in mdl.pars:
    if p.link is None:
        print("{:10s} -> {:.3f}".format(p.fullname, p.val))
    else:
        print("{:10s} -> link to {}".format(p.fullname, p.link.name))

print(repr(g2.pos))
print(g2.pos.link)
print(repr(g2.pos.link))
print(g2.pos.link)
