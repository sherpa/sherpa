=========================
Caching model evaluations
=========================

Sherpa contains a rudimentary system for caching the results
of 1D model evaluations in order to speed up the time to evaluate
models, at the expense of using more memory.
The :py:func:`~sherpa.models.model.modelCacher1d`
function decorator is applied to the
:py:meth:`~sherpa.models.model.ArithmeticModel.calc` method of
:py:class:`~sherpa.models.model.ArithmeticModel` models, and this then
uses the parameter values, evaluation grid, and integrate setting to
look for a value from that model's cache. If found the value is returned,
otherwise the model is evaluated and the result is added to the cache.

Unfortunately it is not always obvious if a model uses caching, or how
effective it is.

What models are cached?
=======================

There is unfortunately no easy way to determine whether a model
uses the cache without either viewing the model definition - looking
for the application of ``@modelCacher1d`` to the ``calc`` method - or
by running a test as shown below,
:ref:`in the example section <example-modelcacher1d>`.

When is the cache useful?
=========================

At present most 1D models use the cache by default.
It is intended to improve fit performance, but the actual
time saved depends on the model and the data being fit.
Compared to most built-in sherpa models, models in the optional XSPEC model
library (:py:mod:`sherpa.astro.xspec`) tend to be more complex and
thus benefit more from caching.

Can I turn off this behavior?
=============================

The size of the cache for a specific model component called ``mdl`` can
be set to zero (``mdl.cache=0``) to turn off the cache behavior.
This may be useful if you are evaluating models over a large grid,
to save memory. For a composite model (e.g. a sum of models) you need
to set the cache for each component.

When we use the fit method of the UI or the fit method of an optimizer, it is possible to
turn off caching for all models with the ``cache`` parameter, e.g.
``ui.fit(..., cache=False)``.
This works by iterating over all the model components and setting the
``cache`` attribute to zero. Note that the opposite is not true: setting
``cache=True`` does not turn on the cache for all model components, it simply
leaves it at the previous setting. This is because some models may not work with
caching at all and need to stay at ``cache=0`` at all times.
The cache has to be manually set to a positive number for all models that should use the cache
to allow caching again.

How does the cache work?
========================

The parameter values, integrate setting, and grid values are used to
create a unique token - the SHA256 hash of the values - which is used
to look up a value in the `_cache` dictionary. If it exists then the
stored value is returned, otherwise the model is evaluated and added
to the `_cache` dictionary. In order to keep the cache size small, the
oldest element in the cache is removed when the number of entries becomes
larger than :py:attr:`~sherpa.models.model.ArithmeticModel.cache` elements (the
default value for this attribute is 5).


Examples
========

.. _example-modelcacher1d:

Checking the cache
------------------

In the following example we evaluate a model and check the `_cache`
attribute, and see that it has been updated by the model evaluation.

>>> from sherpa.models.basic import Box1D
>>> m = Box1D()
>>> m.xlow = 1.5
>>> m.xhi = 4.5
>>> print(m._cache)
{}
>>> print(m([1, 2, 3, 4, 5, 6]))
[0. 1. 1. 1. 0. 0.]
>>> print(m._cache)  # doctest: +SKIP
{b'<random byte string>': array([0., 1., 1., 1., 0., 0.])}


Fit and the startup method
--------------------------

The fit method can also be seen to use the cache (although in this
case it isn't worth it!). First we set up the data::

    >>> import numpy as np
    >>> from sherpa.data import Data1D
    >>> x = np.arange(1, 4)
    >>> y = [4, 5, 2]
    >>> data = Data1D('example', x, y)

A simple model is used::

    >>> from sherpa.models.basic import Const1D
    >>> mdl = Const1D()
    >>> print(mdl.c0.val)
    1.0
    >>> print(mdl._cache)
    {}

The fit only takes 4 iterations, so the cache doesn't help here! Note that
the `startup` and `teardown` methods are called automatically by
:py:meth:`~sherpa.fit.Fit.fit`:

    >>> from sherpa.fit import Fit
    >>> f = Fit(data, mdl)
    >>> result = f.fit()
    >>> print(result.format())
    Method                = levmar
    Statistic             = chi2gehrels
    Initial fit statistic = 2.4176
    Final fit statistic   = 0.534697 at function evaluation 4
    Data points           = 3
    Degrees of freedom    = 2
    Probability [Q-value] = 0.765406
    Reduced statistic     = 0.267349
    Change in statistic   = 1.8829
       const1d.c0     3.39944      +/- 1.74862

The cache contains 4 elements which we can display::

    >>> print(mdl.c0.val)
    3.399441714533379
    >>> print(len(mdl._cache))
    4
    >>> for v in mdl._cache.values():
    ...     print(v)
    ...
    [1. 1. 1.]
    [1.00034527 1.00034527 1.00034527]
    [3.39944171 3.39944171 3.39944171]
    [3.40061543 3.40061543 3.40061543]

