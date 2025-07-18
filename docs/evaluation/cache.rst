.. _cache:

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

It is hard to predict how effective caching is so
the defaults are set for typical Sherpa use cases, based on some benchmarking,
but your performance might be improved with different settings.

What models are cached?
=======================

A model uses the cache if ``@modelCacher1d`` is applied to the ``calc`` method
**and** `model.cache` is set to a positive integer.
Unfortunately it is not easy to check weather ``@modelCacher1d`` is applied without
looking at the source codeb- or by running a test as shown below,
:ref:`in the example section <example-modelcacher1d>`.

When is the cache useful?
=========================

At present most 1D models use the cache by default.
It is intended to improve fit performance, but the actual
time saved depends on the model and the data being fit.
Compared to most built-in sherpa models, models in the optional XSPEC model
library (:py:mod:`sherpa.astro.xspec`) tend to be more complex and
thus benefit more from caching.

By default, the cache is switched off (`mdl.cache=0`) for simple models where the model
evaluation is fast, sometimes even faster than the hashing that is needed to look up
a value in the cache. An obvious example is a scale or constant model
(`sherpa.models.basic.Scale1D` or `sherpa.models.basic.Constant1D`),
but this also applies to some fast analytical
models such as an XSPEC black body (`sherpa.astro.xspec.XSbbody`).

Can I turn off this behavior for other models?
==============================================

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

How do I set the default in my own models to use or not use the cache?
======================================================================

In order to be cacheable at all, a model must have the ``@modelCacher1d`` decorator
applied to the ``calc`` method. The default for an `~sherpa.models.model.ArithmeticModel`
is a cache size of 5, but this can be changed by setting the
``cache`` attribute in the model's
constructor. For example, here we deactivate the cache by default by setting the value to 0,
but we still decorate the ``calc`` method so that a user can switch if back on for
individual instances of the model::

    >>> from sherpa.models.model import ArithmeticModel, modelCacher1d, Parameter
    >>> class MyModel(ArithmeticModel):
    ...     def __init__(self, name='mymodel'):
    ...         self.xpos = Parameter(name, 'offset', 0)
    ...         super().__init__(name, (self.offset,))
    ...         self.cache = 0
    ...
    ...     @modelCacher1d
    ...     def calc(self, p, *args, **kwargs):
    ...         # do something
    ...         return p[0] + args[0]
    >>> m = MyModel()
    >>> m.cache = 3  # use the cache in models instance m

Do **not** set `cache = 0` as a class attribute. If you use `~sherpa.models.model.modelCacher1d`
``cache`` is actually a property that does other things when the cache is set (e.g. reset the
cache content). Setting a class attribute will lead to errors when the decorated ``calc`` method
is called, i.e. the following will not work::

    >>> class MyModel(ArithmeticModel):
    ...     cache = 0  # DO NOT DO THIS!


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
    >>> x = np.arange(0, 3)
    >>> y = [2, 0.3, 0.02]
    >>> data = Data1D('example', x, y)

A simple model is used::

    >>> from sherpa.models.basic import Exp10
    >>> mdl = Exp10()
    >>> mdl.offset.frozen = True
    >>> mdl.offset = 1.0
    >>> mdl.coeff.frozen = True
    >>> mdl.coeff = -1.0
    >>> print(mdl.ampl.val)
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
    Initial fit statistic = 9.178
    Final fit statistic   = 0.00239806 at function evaluation 4
    Data points           = 3
    Degrees of freedom    = 2
    Probability [Q-value] = 0.998802
    Reduced statistic     = 0.00119903
    Change in statistic   = 9.1756
       exp10.ampl     0.201694     +/- 0.263543

The cache contains 4 elements which we can display::

    >>> print(len(mdl._cache))
    4
    >>> for v in mdl._cache.values():
    ...     print(v)
    ...
    [10.   1.   0.1]
    [10.00345267  1.00034527  0.10003453]
    [2.01694277 0.20169428 0.02016943]
    [2.01763916 0.20176392 0.02017639]
