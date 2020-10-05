************************
Creating model instances
************************

The :mod:`sherpa.models` and :mod:`sherpa.astro.models` namespaces
provides a collection of one- and two-dimensional models. There
are also more specialised models, such as those in
:mod:`sherpa.astro.optical`, :mod:`sherpa.astro.xspec`,
:mod:`sherpa.instrument`, and :mod:`sherpa.astro.instrument`.

The following modules are assumed to have been imported for this
section::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from sherpa import models

Creating a model instance
=========================

Models must be created before there parameter values can
be set. In this case a one-dimensional gaussian using the
:py:class:`~sherpa.models.basic.Gauss1D` class::

    >>> g = models.Gauss1D()
    >>> print(g)
    gauss1d
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       gauss1d.fwhm thawed           10  1.17549e-38  3.40282e+38
       gauss1d.pos  thawed            0 -3.40282e+38  3.40282e+38
       gauss1d.ampl thawed            1 -3.40282e+38  3.40282e+38

A description of the model is provided by ``help(g)``.

The parameter values have a current value, a valid range
(as given by the the minimum and maximum columns in the table above),
and a units field. The units field is a string, describing the
expected units for the parameter; there is currently *no support* for
using `astropy.units
<https://docs.astropy.org/en/stable/units/index.html>`_ to set a
parameter value.  The "Type" column refers to whether the parameter is
fixed, (``frozen``) or can be varied during a fit (``thawed``),
as described below, in the :ref:`params-freeze` section.

Models can be given a name, to help distinguish multiple versions
of the same model type. The default value is the lower-case version
of the class name.

::

    >>> g.name
    'gauss1d'
    >>> h = models.Gauss1D('other')
    >>> print(h)
    other
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       other.fwhm   thawed           10  1.17549e-38  3.40282e+38
       other.pos    thawed            0 -3.40282e+38  3.40282e+38
       other.ampl   thawed            1 -3.40282e+38  3.40282e+38
    >>> h.name
    'other'

The model classes are expected to derive from the
:py:class:`~sherpa.models.model.ArithmeticModel` class, although
more-complicated cases, such as :doc:`convolution models
<../evaluation/convolution>`, may extend other classes.

.. _model-combine:

Combining models
================

Models can be combined and shared by using the standard Python
numerical operators. For instance, a one-dimensional gaussian
plus a flat background - using the
:py:class:`~sherpa.models.basic.Const1D` class - would be
represented by the following model::

    >>> src1 = models.Gauss1D('src1')
    >>> back = models.Const1D('back')
    >>> mdl1 = src1 + back
    >>> print(mdl1)
    (src1 + back)
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       src1.fwhm    thawed           10  1.17549e-38  3.40282e+38
       src1.pos     thawed            0 -3.40282e+38  3.40282e+38
       src1.ampl    thawed            1 -3.40282e+38  3.40282e+38
       back.c0      thawed            1 -3.40282e+38  3.40282e+38

Now consider fitting a second dataset where it is known that the background
is two times higher than the first::

    >>> src2 = models.Gauss1D('src2')
    >>> mdl2 = src2 + 2 * back
    >>> print(mdl2)
    (src2 + (2 * back))
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       src2.fwhm    thawed           10  1.17549e-38  3.40282e+38
       src2.pos     thawed            0 -3.40282e+38  3.40282e+38
       src2.ampl    thawed            1 -3.40282e+38  3.40282e+38
       back.c0      thawed            1 -3.40282e+38  3.40282e+38

The two models can then be fit separately or simultaneously. In this
example the two source models (the Gaussian component) were completely
separate, but they could have been identical - in which case
``mdl2 = src1 + 2 * back`` would have been used instead - or
:ref:`parameter linking <params-link>` could be used to constrain the
models. An example of the use of linking would be to force the two
FWHM (full-width half-maximum)
parameters to be the same but to let the position and amplitude
values vary independently.

More information is available in the
:doc:`combining models <../evaluation/combine>`
and
:doc:`convolution <../evaluation/convolution>`
documentation.

Changing a parameter
====================

The parameters of a model - those numeric variables that control the
shape of the model, and that can be varied during a fit -
can be accesed as attributes, both to read or change
the current settings. The
:py:attr:`~sherpa.models.parameter.Parameter.val` attribute
contains the current value::

    >>> print(h.fwhm)
    val         = 10.0
    min         = 1.17549435082e-38
    max         = 3.40282346639e+38
    units       =
    frozen      = False
    link        = None
    default_val = 10.0
    default_min = 1.17549435082e-38
    default_max = 3.40282346639e+38
    >>> h.fwhm.val
    10.0
    >>> h.fwhm.min
    1.1754943508222875e-38
    >>> h.fwhm.val = 15
    >>> print(h.fwhm)
    val         = 15.0
    min         = 1.17549435082e-38
    max         = 3.40282346639e+38
    units       =
    frozen      = False
    link        = None
    default_val = 15.0
    default_min = 1.17549435082e-38
    default_max = 3.40282346639e+38

Assigning a value to a parameter directly (i.e. without using the
``val`` attribute) also works::

    >>> h.fwhm = 12
    >>> print(h.fwhm)
    val         = 12.0
    min         = 1.17549435082e-38
    max         = 3.40282346639e+38
    units       =
    frozen      = False
    link        = None
    default_val = 12.0
    default_min = 1.17549435082e-38
    default_max = 3.40282346639e+38

.. _params-limits:

The soft and hard limits of a parameter
=======================================

Each parameter has two sets of limits, which are referred to as
"soft" and "hard". The soft limits are shown when the model
is displayed, and refer to the
:py:attr:`~sherpa.models.parameter.Parameter.min`
and
:py:attr:`~sherpa.models.parameter.Parameter.max`
attributes for the parameter, whereas the hard limits are
given by the
:py:attr:`~sherpa.models.parameter.Parameter.hard_min`
and
:py:attr:`~sherpa.models.parameter.Parameter.hard_max`
(which are not displayed, and can not be changed).

    >>> print(h)
    other
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       other.fwhm   thawed           12  1.17549e-38  3.40282e+38
       other.pos    thawed            0 -3.40282e+38  3.40282e+38
       other.ampl   thawed            1 -3.40282e+38  3.40282e+38
    >>> print(h.fwhm)
    val         = 12.0
    min         = 1.17549435082e-38
    max         = 3.40282346639e+38
    units       =
    frozen      = False
    link        = None
    default_val = 12.0
    default_min = 1.17549435082e-38
    default_max = 3.40282346639e+38

These limits act to bound the acceptable parameter range; this
is often because certain values are physically impossible, such
as having a negative value for the full-width-half-maxium value
of a Gaussian, but can also be used to ensure that the fit is
restricted to a meaningful part of the search space. The hard
limits are set by the model class, and represent the full
valid range of the parameter, whereas the soft limits can be
changed by the user, although they often default to the same
values as the hard limits.

Setting a parameter to a value outside its soft limits will
raise a :py:exc:`~sherpa.utils.err.ParameterErr` exception.

During a fit the paramater values are bound by the soft limits,
and a screen message will be displayed if an attempt to move
outside this range was made. During error analysis the parameter
values are allowed outside the soft limits, as long as they remain
inside the hard limits.

.. _params-guess:

Guessing a parameter's value from the data
==========================================

Sherpa models have a
:py:meth:`~sherpa.models.model.Model.guess`
method which is used to seed the paramters (or
parameter) with values and
:ref:`soft-limit ranges <params-limits>`
which match the data.
The idea is to move the parameters to values appropriate
for the data, which can avoid un-needed computation by
the optimiser.

The existing ``guess`` routines are very basic - such as
picking the index of the largest value in the data for
the peak location - and do not always account for the
full complexity of the model expression, so care should
be taken when using this functionality.

The arguments depend on the model type, since both the
independent and dependent axes may be used, but the
:py:meth:`~sherpa.data.Data.to_guess` method of
a data object will return the correct data (assuming the
dimensionality and type match)::

    >>> mdl.guess(*data.to_guess())

Note that the soft limits can be changed, as in this example
which ensures the position of the gaussian falls within the
grid of points (since this is the common situation; if the source
is meant to lie outside the data range then the limits will
need to be increased manually)::

    >>> yg, xg = np.mgrid[4000:4050:10, 3000:3070:10]
    >>> r2 = (xg - 3024.2)**2 + (yg - 4011.7)**2
    >>> zg = 2400 * np.exp(-r2 / 1978.2)
    >>> d2d = Data2D('example', xg.flatten(), yg.flatten(), zg.flatten(),
                     shape=zg.shape)
    >>> mdl = Gauss2D('mdl')
    >>> print(mdl)
    mdl
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       mdl.fwhm     thawed           10  1.17549e-38  3.40282e+38
       mdl.xpos     thawed            0 -3.40282e+38  3.40282e+38
       mdl.ypos     thawed            0 -3.40282e+38  3.40282e+38
       mdl.ellip    frozen            0            0        0.999
       mdl.theta    frozen            0     -6.28319      6.28319    radians
       mdl.ampl     thawed            1 -3.40282e+38  3.40282e+38
    >>> mdl.guess(*d2d.to_guess())
    >>> print(mdl)
    mdl
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       mdl.fwhm     thawed           10  1.17549e-38  3.40282e+38
       mdl.xpos     thawed         3020         3000         3060
       mdl.ypos     thawed         4010         4000         4040
       mdl.ellip    frozen            0            0        0.999
       mdl.theta    frozen            0     -6.28319      6.28319    radians
       mdl.ampl     thawed      2375.22      2.37522  2.37522e+06

.. _params-freeze:

Freezing and Thawing parameters
===============================

Not all model parameters should be varied during a fit: perhaps
the data quality is not sufficient to constrain all the parameters,
it is already known, the parameter is highly correlated with
another, or perhaps the parameter value controls a behavior of the
model that should not vary during a fit (such as the interpolation
scheme to use). The :py:attr:`~sherpa.models.parameter.Parameter.frozen`
attribute controls whether a fit
should vary that parameter or not; it can be changed directly,
as shown below::

    >>> h.fwhm.frozen
    False
    >>> h.fwhm.frozen = True

or via the :py:meth:`~sherpa.models.parameter.Parameter.freeze`
and :py:meth:`~sherpa.models.parameter.Parameter.thaw`
methods for the parameter.

::

    >>> h.fwhm.thaw()
    >>> h.fwhm.frozen
    False

There are times when a model parameter should *never* be varied
during a fit. In this case the
:py:attr:`~sherpa.models.parameter.Parameter.alwaysfrozen`
attribute will be set to ``True`` (this particular
parameter is read-only).

.. _params-link:

Linking parameters
==================

There are times when it is useful for one parameter to be
related to another: this can be equality, such as saying that
the width of two model components are the same, or a functional
form, such as saying that the position of one component is a
certain distance away from another component. This concept
is refererred to as linking parameter values. The second case
incudes the first - where the functional relationship is equality -
but it is treated separately here as it is a common operation.
Lnking parameters also reduces the number of free parameters in a fit.

The following examples use the same two model components::

    >>> g1 = models.Gauss1D('g1')
    >>> g2 = models.Gauss1D('g2')

Linking parameter values requires referring to the parameter, rather
than via the :py:attr:`~sherpa.models.parameter.Parameter.val` attribute.
The :py:attr:`~sherpa.models.parameter.Parameter.link` attribute
is set to the link value (and is ``None`` for parameters that are
not linked).

Equality
--------

After the following, the two gaussian components have the same
width::

    >>> g2.fwhm.val
    10.0
    >>> g2.fwhm = g1.fwhm
    >>> g1.fwhm = 1024
    >>> g2.fwhm.val
    1024.0
    >>> g1.fwhm.link is None
    True
    >>> g2.fwhm.link
    <Parameter 'fwhm' of model 'g1'>

When displaying the model, the value and link expression are included::

    >>> print(g2)
    g2
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       g2.fwhm      linked         1024            expr: g1.fwhm
       g2.pos       thawed            0 -3.40282e+38  3.40282e+38
       g2.ampl      thawed            1 -3.40282e+38  3.40282e+38

Functional relationship
-----------------------

The link can accept anything that evaluates to a value,
such as adding a constant.

::

    >>> g2.pos = g1.pos + 8234
    >>> g1.pos = 1200
    >>> g2.pos.val
    9434.0

The :py:class:`~sherpa.models.parameter.CompositeParameter` class
controls how parameters are combined. In this case the result
is a :py:class:`~sherpa.models.parameter.BinaryOpParameter` object.

Including another parameter
---------------------------

It is possible to include other parameters in a link expression,
which can lead to further constraints on the fit. For instance,
rather than using a fixed separation, a range can be used. One
way to do this is to use a :py:class:`~sherpa.models.basic.Const1D`
model, restricting the value its one parameter can vary.

::

    >>> sep = models.Const1D('sep')
    >>> print(sep)
    sep
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       sep.c0       thawed            1 -3.40282e+38  3.40282e+38
    >>> g2.fwhm = g1.fwhm + sep.c0
    >>> sep.c0 = 1200
    >>> sep.c0.min = 800
    >>> sep.c0.max = 1600

In this example, the separation of the two components is restricted
to lie in the range 800 to 1600.

In order for the optimiser to recognize that it needs to vary the
new parameter (``sep.c0``), the component *must* be included in the
model expression. As it does not contribute to the model output
directly, it should be multiplied by zero. So, for this example
the model to be fit would be given by an expression like::

   >>> mdl = g1 + g2 + 0 * sep

.. _parameter_reset:

Resetting parameter values
==========================

.. todo::

   Needs work, including discussing the
   :py:attr:`~sherpa.models.parameter.Parameter.default_val` attribute?

The
:py:meth:`~sherpa.models.parameter.Parameter.reset`
method of a parameter will change the parameter settings (which
includes the status of the thawed flag and allowed ranges,
as well as the value) to the values they had the last time
the parameter was *explicitly* set. That is, it does not restore
the initial values used when the model was created, but the
last values the user set.

The model class has its own
:py:meth:`~sherpa.models.model.Model.reset`
method which calls reset on the thawed parameters. This can be used to
:ref:`change the starting point of a fit <change_fit_starting_point>`
to see how robust the optimiser is by:

- explicitly setting parameter values (or using the default values)
- fit the data
- call reset
- change one or more parameters
- refit


Inspecting models and parameters
================================

Models, whether a single component or composite, contain a
``pars`` attribute which is a tuple of all the parameters
for that model. This can be used to programatically query
or change the parameter values.
There are several attributes that return arrays of values
for the thawed parameters of the model expression: the most
useful is :py:attr:`~sherpa.models.model.Model.thawedpars`,
which gives the current values.

Composite models can be queried to find the individual
components using the ``parts`` attribute, which contains
a tuple of the components (these components can themselves
be composite objects).
