.. _usermodel:

**********************
Writing your own model
**********************

A model class can be created to fit any function, or interface with
external code.

.. todo::

   There should be some description of what needs to be done, as well
   as examples.
   Does the 1D example need to be cleaned up to separate out unnescessary code,
   perhaps just hiding the setup code (and it would be nice if
   this could be shared with the setup).
   Should the 2D example add commentary to point out the following
   note I added at the time:
   "Hmmm, this looks similar to the Sherpa results. In particular
   the 0,0 value is -80 not 1. Aha, is it a normalization at
   (0,0) vs (1,1) sort of thing?"
   

Overview
========

A single model is defined by the parameters of the model - stored
as `sherpa.models.parameter.Parameter` instances - and the function that
takes the parameter values along with an array of grid values to calculate
values of the model. This is packaged together in a class for each model,
and instances of this class can then be used to fit data.

Sherpa provides several base classes in `sherpa.models.model`; the most
commonly used ones are:

* `~sherpa.models.model.ArithmeticModel` is the main base class for deriving user models since
  it supports combining models (e.g. by addition or multiplication) and
  a cache to reduce evaluation time at the expense of memory use.
  (`~sherpa.models.model.ArithmeticConstantModel` and `ArithmeticFunctionModel` are less general
  than `~sherpa.models.model.ArithmeticModel` and can be used to represent
 a constant value or a function.)

* `~sherpa.models.model.RegriddableModel1D` and
  `~sherpa.models.model.RegriddableModel2D` which build on
  `~sherpa.models.model.ArithmeticModel` to allow a model to be
  evaluated on a different grid to that requested.

To define a new model, a new class inherits from one of these base classes and
implements an ``__init__`` method that defines the model parameters
and a ``calc`` method that, given the parameter values and the coordinates
at which the model is to be evaluated, returns the model values.

Calculating the model values
============================

Many Sherpa models are set up to work with both integrated and non-integrated data,
some might even work with data of different dimensionality. Because of that, the input
to the ``calc`` method can take several forms:

 - For `~sherpa.data.Data1D` data, the input is ``calc(pars, x, **kwargs)``
 - For `~sherpa.data.Data1DInt` data, the input is ``calc(pars, xlo, xhi, **kwargs)``, where
   ``xlo`` and ``xhi`` are the lower and upper bounds of the bins in the independent coordinate.
 - For `~sherpa.data.Data2D` data, the input is ``calc(pars, x0, x1, **kwargs)``. ``x0`` and ``x1``
   are flattened arrays for the two independent coordinates, so they are one-dimensional.
 - For `~sherpa.data.Data2DInt` data, the input is ``calc(pars, x0lo, x1lo, x0hi, x1hi, **kwargs)``.
   Here, ``x0lo`` and ``x0hi`` are the lower and upper bounds of the bins in the first coordinate of
   the independent variable, and ``x1lo`` and ``x1hi`` are the lower and upper bounds of the bins in
   the second coordinate.
 - For `~sherpa.data.astro.DataPHA` data, the input is ``calc(pars, xlo, xhi, **kwargs)``,
   where x could be given in either energy or wavelength units, depending on the setting
   of ``set_analysis``.

Sherpa provides other data classes and users can implement their own data classes, that might follow
different conventions, but the cases listed above cover the typical use of Sherpa.

In all cases, ``**kwargs`` are keywords that the model might accept (this feature is not used
by most Sherpa models) and ``pars`` is a list of parameter values, in the same order as they are
listed in the model instances::

    >>> from sherpa.models.basic import Gauss1D
    >>> mdl = Gauss1D(name='mygauss')
    >>> print([p.fullname for p in mdl.pars])
    ['mygauss.fwhm', 'mygauss.pos', 'mygauss.ampl']

Often, in Sherpa, we want models to work with both integrated and non-integrated datasets, and thus
``calc`` is defined to accept a flexible number of arguments with
``def calc(self, pars, x, *args, **kwargs)`` and the code in ``calc`` can branch::

    >>> from sherpa.models import model
    >>> class LinearModel(model.RegriddableModel1D):
    ...     def __init__(self, name='linemodel'):
    ...         self.slope = model.Parameter(name, 'slope', 1, min=-10, hard_min=-1e20, max=10, hard_max=1e20)
    ...         super().__init__(name, (self.slope, ))
    ...
    ...     def calc(self, pars, x, *args, **kwargs):
    ...         if len(args) == 0:
    ...             return x * pars[0]
    ...         elif len(args) == 1:
    ...             xlo, xhi = x, args[0]
    ...             return (xlo + xhi) / 2 * pars[0]
    ...         else:
    ...             raise ValueError('This model only works in 1 dimension')

For a linear model, the integrated value over the bin is just the linear model evaluated in
the middle of the bin, but in general, the implementation can be more complex.

We can now evaluate this model for points or for bins::

    >>> import numpy as np
    >>> linear = LinearModel()
    >>> linear.calc([2], np.array([1.3, 3.4]))
    array([2.6, 6.8])
    >>> linear.calc([2.], np.array([0, 2, 3]), np.array([2, 3, 5]))
    array([2., 5., 8.])

In the examples below, we will set up full data classes and fits and not just pass
the numbers directly into the ``calc`` method.

A one-dimensional model
=======================

An example is the
`AstroPy trapezoidal model <https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Trapezoid1D.html>`_,
which has four parameters: the amplitude of the central region, the center
and width of this region, and the slope. The following model class,
which was not written for efficiency or robustness, implements this
interface:

.. literalinclude:: ../code/trap.py

This can be used in the same manner as the
:py:class:`~sherpa.models.basic.Gauss1D` model
in the :ref:`quick guide to Sherpa<quick-gauss1d>`.

First, create the data to fit::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> np.random.seed(0)
    >>> x = np.linspace(-5., 5., 200)
    >>> ampl_true = 3
    >>> pos_true = 1.3
    >>> sigma_true = 0.8
    >>> err_true = 0.2
    >>> y = ampl_true * np.exp(-0.5 * (x - pos_true)**2 / sigma_true**2)
    >>> y += np.random.normal(0., err_true, x.shape)

Now create a Sherpa data object::
  
    >>> from sherpa.data import Data1D
    >>> d = Data1D('example', x, y)

Set up the user model::
  
    >>> from trap import Trap1D
    >>> t = Trap1D()
    >>> print(t)
    trap1d
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       trap1d.ampl  thawed            1            0  3.40282e+38           
       trap1d.center thawed            1 -3.40282e+38  3.40282e+38           
       trap1d.width thawed            1            0  3.40282e+38           
       trap1d.slope thawed            1            0  3.40282e+38           

Finally, perform the fit::
  
    >>> from sherpa.fit import Fit
    >>> from sherpa.stats import LeastSq
    >>> from sherpa.optmethods import LevMar
    >>> tfit = Fit(d, t, stat=LeastSq(), method=LevMar())
    >>> tres = tfit.fit()
    >>> if not tres.succeeded: print(tres.message)

Rather than use a :py:class:`~sherpa.plot.ModelPlot` object,
the ``overplot`` argument can be set to allow multiple values
in the same plot::

    >>> from sherpa import plot
    >>> dplot = plot.DataPlot()
    >>> dplot.prepare(d)
    >>> dplot.plot()
    >>> mplot = plot.ModelPlot()
    >>> mplot.prepare(d, t)
    >>> mplot.plot(overplot=True)

.. image:: ../_static/model_classes/usermodel/data1d_trap_fit.png

.. _example-usermodel-2d:

A two-dimensional model
=======================

The two-dimensional case is similar to the one-dimensional case,
with the major difference being the number of independent axes to
deal with. In the following example the model is assumed to only be
applied to non-integrated data sets, as it simplifies the implementation
of the ``calc`` method.

It also shows one way of embedding models from a different system,
in this case the
`two-dimemensional polynomial model 
<https://docs.astropy.org/en/stable/api/astropy.modeling.polynomial.Polynomial2D.html>`_
from the AstroPy package.

.. literalinclude:: ../code/poly.py

Repeating the 2D fit by first setting up the data to fit::

    >>> np.random.seed(0)
    >>> y2, x2 = np.mgrid[:128, :128]
    >>> z = 2. * x2 ** 2 - 0.5 * y2 ** 2 + 1.5 * x2 * y2 - 1.
    >>> z += np.random.normal(0., 0.1, z.shape) * 50000.

Put this data into a Sherpa data object::
    
    >>> from sherpa.data import Data2D
    >>> x0axis = x2.ravel()
    >>> x1axis = y2.ravel()
    >>> d2 = Data2D('img', x0axis, x1axis, z.ravel(), shape=(128,128))

Create an instance of the user model::
  
    >>> from poly import WrapPoly2D
    >>> wp2 = WrapPoly2D('wp2')
    >>> wp2.c1_0.frozen = True
    >>> wp2.c0_1.frozen = True

Finally, perform the fit::
  
    >>> f2 = Fit(d2, wp2, stat=LeastSq(), method=LevMar())
    >>> res2 = f2.fit()
    >>> if not res2.succeeded: print(res2.message)
    >>> print(res2)
    datasets       = None
    itermethodname = none
    methodname     = levmar
    statname       = leastsq
    succeeded      = True
    parnames       = ('wp2.c0_0', 'wp2.c2_0', 'wp2.c0_2', 'wp2.c1_1')
    parvals        = (-80.289475553599914, 1.9894112623565667, -0.4817452191363118, 1.5022711710873158)
    statval        = 400658883390.6685
    istatval       = 6571934382318.328
    dstatval       = 6.17127549893e+12
    numpoints      = 16384
    dof            = 16380
    qval           = None
    rstat          = None
    message        = successful termination
    nfev           = 80
    >>> print(wp2)
    wp2
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       wp2.c0_0     thawed     -80.2895 -3.40282e+38  3.40282e+38           
       wp2.c1_0     frozen            0 -3.40282e+38  3.40282e+38           
       wp2.c2_0     thawed      1.98941 -3.40282e+38  3.40282e+38           
       wp2.c0_1     frozen            0 -3.40282e+38  3.40282e+38           
       wp2.c0_2     thawed    -0.481745 -3.40282e+38  3.40282e+38           
       wp2.c1_1     thawed      1.50227 -3.40282e+38  3.40282e+38
   
