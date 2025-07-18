#
#  Copyright (C) 2007, 2015, 2018, 2020, 2021, 2023 - 2025
#  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""Optimization classes.

The `OptMethod` class provides an interface to a number of optimisers.
When creating an optimizer an optional name can be added; this name is
only used in string representations of the class:

>>> from sherpa.optmethods import NelderMead
>>> opt = NelderMead()
>>> print(opt)
name         = simplex
ftol         = 1.1920928955078125e-07
maxfev       = None
initsimplex  = 0
finalsimplex = 9
step         = None
iquad        = 1
verbose      = 0
reflect      = True

A model is fit by providing the ``fit`` method a callback, the starting
point (parameter values), and parameter ranges. The callback should
match::

    callback(pars, *statargs, **statkwargs)

and return the statistic value to minimise along with the per-bin
statistic values (note that per-bin here refers to the independent
axis of the data being fit, and not the number of parameters being
fit).

Notes
-----

Each optimizer has certain classes of problem where it is more, or
less, successful. For instance, the `NelderMead` class should
only be used with chi-square based statistics.

Examples
--------

Using Sherpa classes for data, models, and statistics we can create a
callback, in this case using a least-squared statistic to fit a
constant model to a 1D dataset (we do not need to send any extra
arguments to the callback other than the parameter values in this
case):

>>> import numpy as np
>>> from sherpa.data import Data1D
>>> from sherpa.models.basic import Const1D
>>> from sherpa.stats import LeastSq
>>> x = np.asarray([1, 2, 5])
>>> y = np.asarray([3, 2, 7])
>>> d = Data1D('data', x, y)
>>> mdl = Const1D()
>>> stat = LeastSq()
>>> def cb(pars):
...     mdl.thawedpars = pars
...     return stat.calc_stat(d, mdl)

We can check the model before the optimisaton run:

>>> print(mdl)
const1d
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   const1d.c0   thawed            1 -3.40282e+38  3.40282e+38

The model can be fit using the ``fit`` method:

>>> from sherpa.optmethods import NelderMead
>>> opt = NelderMead()
>>> res = opt.fit(cb, mdl.thawedpars, mdl.thawedparmins, mdl.thawedparmaxes)

The return from ``fit`` is a tuple where the first element indicates
whether the fit was successful, then the best-fit parameters, the
best-fit statistic, a string message, along with a dictionary
depending on the optimiser:

>>> print(res)
(True, array([4.]), 14.0, 'Optimization terminated successfully', {'info': True, 'nfev': 98})
>>> print(f"Best-fit value: {res[1][0]}")
Best-fit value: 4.0

We can see that the model has been updated thanks to the use of
``cb``, which sets the ``thawedpars`` attribute of ``mdl``:

>>> print(mdl)
const1d
   Param        Type          Value          Min          Max      Units
   -----        ----          -----          ---          ---      -----
   const1d.c0   thawed            4 -3.40282e+38  3.40282e+38

"""
from sherpa.optmethods import buildin
from sherpa.optmethods import optscipy
from sherpa.optmethods import optoptimagic
from sherpa.optmethods.buildin import *
from sherpa.optmethods.optscipy import *
from sherpa.optmethods.optoptimagic import *

__all__ = buildin.__all__ + optscipy.__all__ + optoptimagic.__all__