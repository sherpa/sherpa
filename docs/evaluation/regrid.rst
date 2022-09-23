****************************************
Evaluating the model on a different grid
****************************************

Sherpa now provides *experimental* support for evaluating a model
on a different grid to the independent axis. This can be used to
better model the underlying distribution (by use of a finer grid),
or to include features that lie outside the data range but, due
to the use of a convolution model, can affect the values within
the data range.

Existing Sherpa models take advantage of this support by inheriting
from the
:py:class:`~sherpa.models.model.RegriddableModel1D`
or
:py:class:`~sherpa.models.model.RegriddableModel2D` classes.
The idea is to create a copy of the model which evaluates the
model on the new grid, and then resamples it onto the output grid.
This is handled by the `regrid` method from these two
classes, which is given the grid to evaluate the model on, and it
returns a copy of the original model.

The examples assume the following:

    >>> import numpy as np

One dimensional models
======================

Non-integrated models
---------------------

For this example we shall evaluate the model on the
grid x=4450.25, 4450.75, ... 4549.75.

::

   >>> xgrid = np.arange(4450.25, 4550, 0.5)

The "regridded" model is created by passing the grid
to the :py:meth:`~sherpa.models.model.RegriddableModel1D.regrid`
method::

   >>> from sherpa.models.basic import Gauss1D
   >>> mdl = Gauss1D('orig')
   >>> rmdl = mdl.regrid(xgrid)

The "new" model is really just a wrapper around the original model,
which means that it shares the parameter values and settings. This
is shown below where the ``pos`` and ``fwhm`` parameters of the
``rmdl`` are labeled by the label given to ``mdl`` (namely
"orig")::

   >>> mdl.pos = 4500
   >>> mdl.fwhm = 20
   >>> print(mdl)
   orig
      Param        Type          Value          Min          Max      Units
      -----        ----          -----          ---          ---      -----
      orig.fwhm    thawed           20  1.17549e-38  3.40282e+38
      orig.pos     thawed         4500 -3.40282e+38  3.40282e+38
      orig.ampl    thawed            1 -3.40282e+38  3.40282e+38
   >>> print(rmdl)
   regrid1d(orig)
      Param        Type          Value          Min          Max      Units
      -----        ----          -----          ---          ---      -----
      orig.fwhm    thawed           20  1.17549e-38  3.40282e+38
      orig.pos     thawed         4500 -3.40282e+38  3.40282e+38
      orig.ampl    thawed            1 -3.40282e+38  3.40282e+38

The regridded model can be evaluated as any other model, but behind the
scenes the data is calculated on the original grid (in this case ``xgrid``)
and then converted to match the evaluation grid.

::

   >>> x = np.arange(4460, 4540, 2)
   >>> y1 = mdl(x)
   >>> y2 = rmdl(x)
   >>> np.testing.assert_allclose(y1, y2)

Integrated models
-----------------

The approach is similar, except that the
:py:meth:`~sherpa.models.model.RegriddableModel1D.regrid` method is
given two arrays, for the low and high edges of each bin.

::

   >>> imdl = mdl.regrid(xgrid[:-1], xgrid[1:])
   >>> y1 = mdl(x[:-1], x[1:])
   >>> y2 = imdl(x[:-1], x[1:])
   >>> np.testing.assert_allclose(y1, y2, rtol=1e-2)

The fact that the models are integrated across the bins means that the
results do not agree that closely in this case, which is why the
relative tolerance is only 0.01 in the check above. The results would
be closer if the two grids had a common origin (e.g. if ``xgrid`` were
defined over the bins 0-0.5, 0.5-1 rather than 0.25-0.75, 0.75-1.25).

Mixing the two
--------------

It is an error to try to evaluate a regridded model with the wrong
style of grid, and will raise a :py:class:`~sherpa.utils.err.ModelErr`
error::

   >>> rmdl(x[:-1], x[1:])
   Traceback (most recent call last):
   ...
   sherpa.utils.err.ModelErr: A non-integrated grid is required for model evaluation


   >>> imdl(x)
   Traceback (most recent call last):
   ...
   sherpa.utils.err.ModelErr: A non-overlapping integrated grid is required for model evaluation,
   e.g. [0.1,0.2],[0.2,0.3]


Two dimensional models
======================

The two-dimensional non-integrated case is created by passing the
``x0`` and ``x1`` arrays to the
:py:meth:`~sherpa.models.model.RegriddableModel1D.regrid` method::

   >>> from sherpa.models.basic import Gauss2D
   >>> g2 = Gauss2D('g2')
   >>> x1grid, x0grid = np.mgrid[100:200:0.5, 50:150:0.5]
   >>> rg2 = g2.regrid(x0grid.flatten(), x1grid.flatten())

As with the one-dmensional case, the regridded model uses the
parameters of the original model::

   >>> g2.xpos = 100
   >>> g2.ypos = 150
   >>> g2.fwhm = 25
   >>> print(rg2)
   regrid2d(g2)
      Param        Type          Value          Min          Max      Units
      -----        ----          -----          ---          ---      -----
      g2.fwhm      thawed           25  1.17549e-38  3.40282e+38
      g2.xpos      thawed          100 -3.40282e+38  3.40282e+38
      g2.ypos      thawed          150 -3.40282e+38  3.40282e+38
      g2.ellip     frozen            0            0        0.999
      g2.theta     frozen            0     -6.28319      6.28319    radians
      g2.ampl      thawed            1 -3.40282e+38  3.40282e+38

.. note::

   Evaluation of the 2D model is complicated by the current implementation.
   Please see `issue 840 <https://github.com/sherpa/sherpa/issues/840>`_
   for more information.
