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

At present most 1D models use the cache by default when evaluated
normally, but not during a fit. It is intended to improve fit
performance - that is, reduce the time taken to fit a dataset - but
there has been limited effort to evaluate its efficiency.

Can I turn off this behavior?
=============================

The `_use_caching` attribute of a model can be set to `False` to stop
the cache behavior. This may be useful if you are evaluating models
over a large grid, to save memory, or the model calculation is not
expensive, and so the extra time used to store the result is not
beneficial.

Note that the :ref:`the startup method <startup-modelcacher1d>` can
change this value, but it depends if you are calling `startup`
directly or indirectly, via the ``fit`` and ``est_errors`` methods of
a fit object.

How does the cache work?
========================

The parameter values, integrate setting, and grid values are used to
create an unique token - the SHA256 hash of the values - which is used
to look up a value in the `_cache` dictionary. If it exists then the
stored value is returned, otherwise the model is evaluated and added
to the `_cache` dictionary. In order to keep the cache size small, the
`_queue` array is used to remove an existing value from the store when
a new value is added. The default size for the `_queue` array is
a single value, but it can be changed by
:ref:`the startup method <startup-modelcacher1d>`.

.. _startup-modelcacher1d:

The startup method
==================

The model :py:meth:`~sherpa.models.model.ArithmeticModel.startup`
method is automatically called by the :py:meth:`~sherpa.fit.Fit.fit`
method, but can also be called manually. It sets the `_use_caching`
attribute and sets the `_queue` array to have
:py:attr:`~sherpa.models.model.ArithmeticModel.cache` elements (the
default value for this attribute is 5).

Although the default value for the `cache` argument to `startup` is
set to `False`, the `sherpa.fit.evaluates_model` decorator - which
is used to wrap the :py:meth:`~sherpa.fit.Fit.fit`,
:py:meth:`~sherpa.fit.Fit.simulfit`, and
:py:meth:`~sherpa.fit.Fit.est_errors` methods - over-rides this value
and uses a value of `True`. Therefore, to turn off the cache you
have to explicitly pass ``cache=False`` to the fit method::

    f.fit(cache=False)

The teardown method
===================

The model :py:meth:`~sherpa.models.model.ArithmeticModel.teardown`
method is run after the fit is done - to match `startup` - and
currently sets the `_use_caching` setting to `False`.

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
>>> print(m._use_caching)
True
>>> print(m._cache)
{}
>>> print(m([1, 2, 3, 4, 5, 6]))
[0. 1. 1. 1. 0. 0.]
>>> print(m._cache)  # doctest: +SKIP
{b'<random byte string>': array([0., 1., 1., 1., 0., 0.])}
>>> print(m._queue)  # doctest: +SKIP
[b'<random byte string>']

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

Note that if we had called::

    f.fit(cache=False)

then the cache would not have been used (e.g. `mdl._cache` would
have remained empty).
