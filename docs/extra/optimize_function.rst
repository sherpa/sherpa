********************
Optimizing functions
********************

The :py:mod:`sherpa.optmethods` module provides classes that let you
optimize :doc:`models <../models/index>` when applyed to :doc:`data
<../data/index>` instances, and this is handled by the :doc:`fit
<../fit/index>` module.  You can apply the optimizers directly to
functions if you do not need to take advantage of Sherpa's model and
data classes.

The optimizers in :py:mod:`sherpa.optmethods.optfcts` are functions
which follow the same interface::

  optimizer(cb, start, lolims, hilims, **kwargs)

where `cb` is a callable that given a set or parameters returns the
current statistic, `start` is the starting positions of the
parameters, `lolims` and `hilims` are the bounds of the parameters,
and the remaining keyword arguments are specific to the optimizer.

Examples
========

The following imports have been made::

  >>> import numpy as np
  >>> from matplotlib import pyplot as plt
  >>> from matplotlib import colors
  >>> from mpl_toolkits.mplot3d import Axes3D
  >>> from astropy.table import Table
  >>> from sherpa.optmethods import optfcts

A simple function
-----------------

The function we want to optimize is a version of the
`Rosenbrock function <https://en.wikipedia.org/wiki/Rosenbrock_function>`_,
where `a=2` and `b=10`, so the best-fit location is at::

  x = 2
  y = 2 * 2

.. note::
   The Sherpa test suite contains a number of test cases based on
   functions like this, taken from
   `More, J. J., Garbow, B. S., and Hillstrom, K. E.. Testing unconstrained optimization software. United States: N. p., 1978. Web. doi:10.2172/6650344. <https://www.osti.gov/biblio/6650344>`_.

The function is::

  def rosenbrock(x, y):
      a, b = 2, 10
      return (a - x) * (a - x) + b * (y - x * x) * (y - x * x)

We can look at the surface close to the best-fit location
What does the function look like?

  >>> y, x = np.mgrid[-2:5.1:0.1, -3:3.1:0.1]
  >>> surface = rosenbrock(x, y)
  >>> print(surface.shape)
  (71, 61)
  >>> print(surface.min())
  0.0
  >>> print(surface.max())
  1235.0

This can be visualized as::

  >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
  >>> surf = ax.plot_surface(x, y, surface,
  ...                        cmap=cm.plasma_r, alpha=0.8,
  ...                        vlim=0, vmax=300)
  >>> ax.set_xlabel('x')
  >>> ax.set_ylabel('y')
  >>> ax.view_init(30, -100)

.. image:: ../_static/extra/rosenbrock-surface.png

In order to optimize this function we need to determine how we
want to define the "search surface"; that is the value that the
optimizer is going to try and minimize. For this dataset, where
the minimum is 0, we can use the square of the function (i.e. a
least-squares model). Now, the function we pass to the minimizer
requires two values - the "statistic" and a per-parameter breakdown
of the parameter - but for now we are going to ignore the latter,
as it's only needed for some optimizers. This gives us::

  def to_optimize1(args):
      x = args[0]
      y = args[1]
      pa = (2 - x) * (2 - x)
      pb = 10 * (y - x * x) * (y - x * x)
      stat = pa * pa + pb * pb
      return stat, None

The function can tehn be monimized by passing the routine, the
starting point, and the bounds, to the optimizer- in this case the
:py:func:`~sherpa.optmethods.optfcts.minim` routine:

  >>> start = [-1.2, 1]
  >>> lo = [-100, -100]
  >>> hi = [100, 100]
  >>> res = optfcts.minim(to_optimize, start, lo, hi)

The return value is a tuple which contains a sucess flag,
the best-fit parameters, the value at this location, a
message, and a dictionary with optimizer-specific information:

  >>> print(f"Success: {res[0]}")
  Success: True
  >>> print(f"Message: {res[3]}")
  Message: successful termination
  >>> print(f"extra:   {res[4]}")
  extra:   {'info': 0, 'nfev': 88}
  >>> print(f"best-fit location: {res[1]}")
  best-fit location: [1.99999996 4.00064474]
  >>> print(f"          minimum: {res[2]}")
            minimum: 1.73832000295664e-11

So, the correct location is (2, 4) and we can see we got close.

As the different optimizers use the same interface we can easily try
different optimizers:

  >>> tbl = Table(names=['method', 'stat0', 'x', 'y'],
  ...             dtype=[str, float, float, float])
  >>> for method in [optfcts.minim, optfcts.neldermead, optfcts.lmdif, optfcts.montecarlo]:
  ...     res = method(to_optimize, start, lo, hi)
  ...     if res[0]:
  ...         tbl.add_row([method.__name__, res[2], res[1][0], res[1][1]])
  ...     else:
  ...         print(f"Failed {method.__name__}: {res[3]}")
  ...
  Failed lmdif: improper input parameters
  >>> print(tbl)
      method           stat0                 x                  y
    ---------- --------------------- ------------------ ------------------
         minim  1.73832000295664e-11 1.9999999630648686  4.000644742396361
    neldermead 4.932797434506576e-11 2.0000007686008545 3.9994204823962214
    montecarlo 6.629393582604866e-12 1.9999997816942297  3.999564165982678

In this particular case :py:func:`~sherpa.optmethods.optfcts.lmdif`
failed, and this is because the ``to_optimize`` function returned None
as its second argument. For the other cases we can see that they all
found the minimum location.

In order to use :py:func:`~sherpa.optmethods.optfcts.lmdif` we need
the per-parameter statistic, which we can get with a small tweak::

  def to_optimize2(args):
      x = args[0]
      y = args[1]
      pa = (a - x) * (a - b)
      pb = b * (y - x * x) * (y - x * x)
      stat = pa * pa + pb * pb
      return stat, [pa, pb]

This lets us use `lmdif`::

  >>> res2 = optfcts.lmdif(to_optimize2, start, lo, hi)
  >>> res2[0]
  True
  >>> res2[1]
  [2.        4.0000846]

The `lmdif` optimizer is one of those that returns a number of
elements in the "extra" dictionary::

  >>> res2[4]
  {'info': 0, 'nfev': 108, 'covar': array([[ 1.56250000e-02, -1.04361202e-01],
         [-1.04361202e-01,  4.15098457e+03]]), 'num_parallel_map': 0}

Optimizing a model
------------------

We can re-create the :ref:`Fitting the NormGauss1D and Lorentz1D models <fit_peaked_data_normgauss1d_lorentz1d>`
section - at least for the :py:class:`~sherpa.models.basic.NormGauss1D` fit.

The normalized gaussian model can be expressed as

::

  def ngauss(x, ampl, pos, fwhm):
      term = 4 * np.log(2)
      numerator = ampl * np.exp(-term * (x - pos) * (x - pos) / (fwhm * fwhm))
      denominator = np.sqrt(np.pi / term) * fwhm
      return numerator / denominator

and the data we are going to fit is read from a file:

  >>> d = np.loadtxt('test_peak.dat')
  >>> x = d[:, 0]
  >>> y = d[:, 1]

In this case the "callback" routine we want to give the minimizer a
routine which will apply a least-squares minimizer to compare the
model to the data (the `x` and `y` values for the independent and
dependent axes respectively).

::

  def cb(pars):
      model = ngauss(x, pars[0], pars[1], pars[2])
      delta = model - y
      statval = (delta * delta).sum()
      # Keep a record of the parameters we've visited
      store.add_row([statval] + list(pars))
      return statval, delta

.. note::
   This function would normally be written as a closure to ensure
   that the use of global terms like `x`, `y`, `ngauss`, and `store` do not
   cause problems.

For the starting point and bounds I'm using to use a "guess" based on
the data:

  >>> start = [y.max(), (x[0] + x[-1]) / 2, (x[-1] - x[0]) / 10]
  >>> lows = [0, x[0], 0]
  >>> his = [1e4, x[-1], x[-1] - x[0]]

which can be used with the :py:func:`~sherpa.optmethods.optfcts.neldermead` optimizer:

  >>> store = Table(names=['stat', 'ampl', 'pos', 'fwhm'])
  >>> flag, bestfit, statval, msg, opts = optfcts.neldermead(cb, start, lows, his)
  >>> flag
  True
  >>> bestfits
  array([30.31354379,  9.24277056,  2.90156713])
  >>> statval
  29.994315740119312
  >>> opts
  {'info': True, 'nfev': 421}
  >>> len(store)

and displayed:

  >>> plt.plot(x, y, label='Data', alpha=0.5)
  >>> plt.plot(x, ngauss(x, *start), label='Start')
  >>> plt.plot(x, ngauss(x, *bestfit), label='Best fit', c='k')
  >>> plt.legend()

.. image:: ../_static/extra/normgauss1d-example.png

This result matches the
:ref:`Fitting the NormGauss1D and Lorentz1D models <fit_peaked_data_normgauss1d_lorentz1d>`
result. Note that at this point you are close to replicating the
main parts of Sherpa (but with a lot-less functionality)!

One tweak that was added to the `cb` routine, compared to `to_optimize`,
was the ability to store the parameter values at each iteration:

  >>> print(store)
         stat               ampl               pos                fwhm
  ------------------ ------------------ ----------------- -------------------
  1995.1782885013076          10.452393              10.0                 2.0
  1830.6327671752736          11.652393              10.0                 2.0
  3522.2622560146397          10.452393              11.2                 2.0
    2156.39741128647          10.452393              10.0                 3.2
   1764.261448064872 11.252393000000001               8.8                 2.8
   2715.012206450446 11.652393000000004 7.600000000000001  3.1999999999999993
   1366.678033289349 11.785726333333333               9.2   1.333333333333333
    4770.95772056462          12.452393 8.799999999999997 0.39999999999999947
  1504.4139142274962 12.674615222222219 8.666666666666668  2.0888888888888886
  2731.7112487850172 12.156096703703703 7.777777777777779   2.148148148148148
                 ...                ...               ...                 ...
  29.994315975508734 30.313386240009088 9.242774794233913  2.9015358619056464
  29.994316369122366  30.31358595155565 9.242759714112402  2.9016069773705953
  29.994318894587913 30.313745720792898 9.242747650015195  2.9016638697425545
  29.994315764507792 30.313466124627716 9.242768762185309  2.9015643080916256
  29.994315898963315 30.313608124884944 9.242779534765045  2.9015614069320215
  29.994316025654843 30.313467246158012 9.242779279880374  2.9015835332886226
  29.994316138869017  30.31320541897285 9.242777555373543   2.901543765169257
   29.99431599926865 30.313726830282576 9.242759968997074   2.901584851013994
   29.99431606407549 30.313465003097413 9.242758244490243  2.9015450828946285
  29.994316180878897  30.31332412437048 9.242757989605572  2.9015672092512297
  29.994315740119312 30.313543791452684 9.242770559884388   2.901567125484238
  Length = 421 rows

We can use this to see how the optimizer worked, color-coding each point by
the statistic (using a log-normalized scale as we go from
:math:`\gt 2000` to :math:`\sim 30`):

  >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
  >>> vmin = store['stat'].min()
  >>> vmax = store['stat'].max()
  >>> norm = colors.LogNorm(vmin=vmin, vmax=vmax)
  >>> ax.plot(store['ampl'], store['pos'], store['fwhm'], alpha=0.4)
  >>> scatter = ax.scatter(store['ampl'], store['pos'], store['fwhm'],
  ...                      c=store['stat'], norm=norm)
  >>> ax.set_xlabel('ampl')
  >>> ax.set_ylabel('pos')
  >>> ax.set_zlabel('fwhm')
  >>> cbar = fig.colorbar(scatter, shrink=0.5. orientation='horizontal')
  >>> cbar.set_label('least-squares statistic')

.. image:: ../_static/extra/normgauss1d-trail.png
