#
#  Copyright (C) 2010, 2016 - 2024
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

"""Allow models to be defined and combined.

A single model is defined by the parameters of the model - stored
as `sherpa.models.parameter.Parameter` instances - and the function that
takes the parameter values along with an array of grid values. The
main classes are:

* `Model` which is the base class and defines most of the interfaces.

* `ArithmeticConstantModel` and `ArithmeticFunctionModel` for representing
  a constant value or a function.

* `ArithmeticModel` is the main base class for deriving user models since
  it supports combining models (e.g. by addition or multiplication) and
  a cache to reduce evaluation time at the expense of memory use.

* `RegriddableModel` builds on ArithmeticModel to allow a model to be
  evaluated on a different grid to that requested: most model classes
  are derived from the 1D (`RegriddableModel1D`) and 2D
  (`RegriddableModel2D`) variants of RegriddableModel.

* `CompositeModel` which is used to represent a model expression, that
  is combined models, such as `m1 * (m2 + m3)`

  * `UnaryOpModel` for model expressions such as `-m1`.

  * `BinaryOpModel` for model expressions such as `m1 + m2`.

  * `NestedModel` for applying one model to another.

* `SimulFitModel` for fitting multiple models and datasets.

Creating a model
================

Models can be created with an optional name, which is useful for
identifying a component in an expression:

    >>> from sherpa.models.basic import Gauss1D
    >>> m1 = Gauss1D()
    >>> m2 = Gauss1D('gmdl')
    >>> print(m1)
    gauss1d
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       gauss1d.fwhm thawed           10  1.17549e-38  3.40282e+38
       gauss1d.pos  thawed            0 -3.40282e+38  3.40282e+38
       gauss1d.ampl thawed            1 -3.40282e+38  3.40282e+38

    >>> print(m2)
    gmdl
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       gmdl.fwhm    thawed           10  1.17549e-38  3.40282e+38
       gmdl.pos     thawed            0 -3.40282e+38  3.40282e+38
       gmdl.ampl    thawed            1 -3.40282e+38  3.40282e+38

Changing parameters
===================

The parameters are the model values that control the output of the
model. A particular model has a fixed set of parameters that can
be inspected with print or the pars attribute:

    >>> print(m2)
    gmdl
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       gmdl.fwhm    thawed           10  1.17549e-38  3.40282e+38
       gmdl.pos     thawed            0 -3.40282e+38  3.40282e+38
       gmdl.ampl    thawed            1 -3.40282e+38  3.40282e+38

    >>> print(m2.pars)
    (<Parameter 'fwhm' of model 'gmdl'>, <Parameter 'pos' of model 'gmdl'>, <Parameter 'ampl' of model 'gmdl'>)

The parameters are instances of the `sherpa.models.parameter.Parameter`
class:

    >>> print(m2.fwhm)
    val         = 10.0
    min         = 1.1754943508222875e-38
    max         = 3.4028234663852886e+38
    units       =
    frozen      = False
    link        = None
    default_val = 10.0
    default_min = 1.1754943508222875e-38
    default_max = 3.4028234663852886e+38

    >>> print(m2.fwhm.val)
    10.0

Setting the model parameter does not require going through the val
attribute as you can say:

    >>> m2.fwhm = 20

Accessing parameter values
--------------------------

The model class is set up so that any attribute access is case
insensitive, so the following are all ways to change the ``fwhm``
parameter:

    >>> m2.fwhm = 10
    >>> m2.FWHM = 10
    >>> m2.FwHm = 10

Linking parameters
------------------

One parameter can be made to reference one or more other parameters, a
process called "linking". The linked parameter is no-longer considered
a free parameter in a fit since its value is derived from the other
parameters. This link can be a simple one-to-one case, such as
ensuring the fwhm parameter of one model is the same as the other:

    >>> m2.fwhm = m1.fwhm

It can be more complex, such as ensuring the position of one line
is a fixed distance from another:

    >>> m2.pos = m1.pos + 23.4

It can even include multiple parameters:

    >>> m3 = Gauss1D("m3")
    >>> m3.ampl = (m1.ampl + m2.ampl) / 2

Requesting the parameter value will return the evaluated expression,
and the expression is stored in the link attribute:

    >>> m1.ampl = 10
    >>> m2.ampl = 12
    >>> m3.ampl.val
    11.0
    >>> m3.ampl.link
    <BinaryOpParameter '(gauss1d.ampl + gmdl.ampl) / 2'>

The string representation of the model changes for linked parameters
to indicate the expression:

    >>> print(m3)
    m3
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       m3.fwhm      thawed           10  1.17549e-38  3.40282e+38
       m3.pos       thawed            0 -3.40282e+38  3.40282e+38
       m3.ampl      linked           11 expr: (gauss1d.ampl + gmdl.ampl) / 2

Model evaluation
================

With a `sherpa.data.Data` instance a model can be evaluated with the
eval_model method of the object. For example:

    >>> import numpy as np
    >>> from sherpa.data import Data1D
    >>> from sherpa.models.basic import Gauss1D
    >>> x = np.asarray([4000, 4100, 4250, 4300, 4400])
    >>> y = np.asarray([10, 20, 50, 40, 30])
    >>> d = Data1D('example', x, y)
    >>> mdl = Gauss1D()
    >>> mdl.pos = 4200
    >>> mdl.fwhm = 200
    >>> mdl.ampl = 50
    >>> ymdl1 = d.eval_model(mdl)
    >>> print(ymdl1)
    [ 3.125      25.         42.04482076 25.          3.125     ]

The model can also be evaluated directly with the independent axis
values:

    >>> ymdl2 = mdl(x)
    >>> print(ymdl2)
    [ 3.125      25.         42.04482076 25.          3.125     ]

Integrated bins
---------------

If given the low and high edges of the bins then the model will - if
supported - evaluate the integral of the model across the bins:

    >>> xlo = np.asarray([4180, 4190, 4195, 4200, 4210])
    >>> xhi = np.asarray([4190, 4194, 4200, 4210, 4220])
    >>> y = mdl(xlo, xhi)
    >>> print(y)
    [491.98725233 199.0964993  249.85566938 498.847153   491.98725233]

Note that the bins are expected to be in ascending order and do not
overlap, but they do not need to be consecutive.

The behavior of a model when given low and high edges depends on
whether the model is written to support this mode - that is,
integrating the model across the bin - and the setting of the
integrate flag of the model. For example, the
`sherpa.models.basic.Gauss1D` model will, by default, integrate the
model across each bin when given the bin edges, but if the flag is set
to `False` then just the first array (here ``xlo``) is used:

    >>> print(mdl.integrate)
    True
    >>> mdl.integrate = False
    >>> y2 = mdl(xlo, xhi)
    >>> print(y2)
    [48.63274737 49.65462477 49.91343163 50.         49.65462477]
    >>> y3 = mdl(xlo)
    >>> y2 == y3
    array([ True,  True,  True,  True,  True])

Direct access
-------------

The calc method of a model can also be used to evaluate the model, and
this requires a list of the parameters and the independent axes:

    >>> pars = [200, 4200, 50]
    >>> y4 = mdl.calc(pars, x)
    >>> y5 = mdl.calc(pars, xlo, xhi)

The parameter order matches the pars attribute of the model:

    >>> print([p.fullname for p in mdl.pars])
    ['gauss1d.fwhm', 'gauss1d.pos', 'gauss1d.ampl']

Model expressions
=================

The `CompositeModel` class is the base class for creating model
expressions - that is the overall model that is combined of one or
more model objects along with possible numeric terms, such as a
model containing two gaussians and a polynomial:

    >>> from sherpa.models.basic import Gauss1D, Polynom1D
    >>> l1 = Gauss1D('l1')
    >>> l2 = Gauss1D('l2')
    >>> l1.pos = 5
    >>> l2.pos = 20
    >>> l2.ampl = l1.ampl
    >>> c = Polynom1D('c')
    >>> mdl = l1 + (0.5 * l2) + c

The resulting model can be evaluated just like an individual
component:

    >>> x = np.arange(-10, 40, 2)
    >>> y = mdl(x)

This model is written so that the amplitude of the ``l2`` component is
half the ``l1`` component by linking the two ``ampl`` parameters and then
including a scaling factor in the model expression for ``l2``. An
alternative would have been to include this scaling factor in the link
expression:

    >>> l2.ampl = l1.ampl / 2

Model cache
===========

The `ArithmeticModel` class and `modelCacher1d` decorator provide basic
support for caching one-dimensional model evaluations - that is, to
avoid re-calculating the model. The idea is to save the results of the
latest calls to a model and return the values from the cache,
hopefully saving time at the expense of using more memory. This is
most effective when the same model is used with multiple datasets
which all have the same grid.

The `_use_caching` attribute of the model is used to determine whether
the cache is used, but this setting can be over-ridden by the startup
method, which is automatically called by the fit and est_errors
methods of a `sherpa.fit.Fit` object.

The `cache_clear` and `cache_status` methods of the `ArithmeticModel`
and `CompositeModel` classes allow you to clear the cache and display
to the standard output the cache status of each model component.

Example
=======

The following class implements a simple scale model which has a single
parameter (``scale``) which defaults to 1. It can be used for both
non-integrated and integrated datasets of any dimensionality (see
`sherpa.models.basic.Scale1D` and `sherpa.models.basic.Scale2D`)::

    class ScaleND(ArithmeticModel):
        '''A constant value per element.'''

        def __init__(self, name='scalend'):
            self.scale = Parameter(name, 'scale', 1)
            self.integrate = False
            pars = (self.scale, )
            ArithmeticModel.__init__(self, name, pars)

        def calc(self, p, *args, **kwargs):
            return p[0] * np.ones_like(args[0])

"""

from __future__ import annotations

import functools
import itertools
import logging
from typing import TYPE_CHECKING, Any, Callable, Iterator, Optional, \
    Sequence, SupportsFloat, SupportsIndex, Type, TypeVar, Union
import warnings

import numpy as np

from sherpa.models.regrid import EvaluationSpace1D, ModelDomainRegridder1D, EvaluationSpace2D, ModelDomainRegridder2D
from sherpa.utils import NoNewAttributesAfterInit, formatting
from sherpa.utils.err import ModelErr, ParameterErr
from sherpa.utils.numeric_types import SherpaFloat

from .op import get_precedences_op, get_precedence_expr, \
    get_precedence_lhs, get_precedence_rhs
from .parameter import Parameter, expand_par

# What routine do we use for the hash in modelCacher1d?  As we do not
# need cryptographic security go for a "quick" algorithm, but md5 is
# not guaranteed to always be present.  There has been no attempt to
# check the run times of these routines for the expected data sizes
# they will be used with.
#
try:
    from hashlib import md5 as hashfunc
except ImportError:
    from hashlib import sha256 as hashfunc

info = logging.getLogger(__name__).info
warning = logging.getLogger(__name__).warning


__all__ = ('Model', 'CompositeModel', 'SimulFitModel',
           'ArithmeticConstantModel', 'ArithmeticModel', 'RegriddableModel1D', 'RegriddableModel2D',
           'UnaryOpModel', 'BinaryOpModel', 'FilterModel', 'modelCacher1d',
           'ArithmeticFunctionModel', 'NestedModel', 'MultigridSumModel')


# These tests refer to other variables which are just not worth
# setting up (or would require too much work to create any useful
# state, such as the cache status lines).
#
__doctest_skip__ = ['ArithmeticModel.cache_status',
                    'CompositeModel.cache_status',
                    'SimulFitModel']


def boolean_to_byte(boolean_value: bool) -> bytes:
    """Convert a boolean to a byte value.

    Parameters
    ----------
    boolean_value : bool
        The value to convert. If not a boolean then it is
        treated as `False`.

    Returns
    -------
    val : bytes
        b'1' if `True` otherwise b'0'.
    """

    bmap = {True: b'1', False: b'0'}
    return bmap.get(boolean_value, b'0')


def modelCacher1d(func: Callable) -> Callable:
    """A decorator to cache 1D ArithmeticModel evaluations.

    Apply to the `calc` method of a 1D model to allow the model
    evaluation to be cached. The decision is based on the
    `_use_caching` attribute of the cache along with the `integrate`
    setting, the evaluation grid, parameter values, and the keywords
    sent to the model.

    Notes
    -----
    The keywords are included in the hash calculation even if they are
    not relevant for the model (as there's no easy way to find this
    out).

    Example
    -------

    Allow `MyModel` model evaluations to be cached::

        def MyModel(ArithmeticModel):
            ...
            @modelCacher1d
            def calc(self, p, *args, **kwargs):
                ...

    """

    @functools.wraps(func)
    def cache_model(cls, pars, xlo, *args, **kwargs):
        # Counts all accesses, even those that do not use the cache.
        cache_ctr = cls._cache_ctr
        cache_ctr['check'] += 1

        # Short-cut if the cache is not being used.
        #
        if not cls._use_caching:
            return func(cls, pars, xlo, *args, **kwargs)

        # Up until Sherpa 4.12.2 we used the kwargs to define the
        # integrate setting, with
        # boolean_to_byte(kwargs.get('integrate', False)) but
        # unfortunately this is used in code like
        #
        #    @modelCacher1d
        #    def calc(..):
        #        kwargs['integrate'] = self.integrate
        #        return somefunc(... **kwargs)
        #
        # and the decorator is applied to calc, which is not
        # called with a integrate kwarg, rather than the call to
        # somefunc, which was sent an integrate setting.
        #
        try:
            integrate = cls.integrate
        except AttributeError:
            # Rely on the integrate kwarg as there's no
            # model setting.
            #
            integrate = kwargs.get('integrate', False)

        data = [np.array(pars).tobytes(),
                boolean_to_byte(integrate),
                np.asarray(xlo).tobytes()]
        if args:
            data.append(np.asarray(args[0]).tobytes())

        # Add any keyword arguments to the list. This will
        # include the xhi named argument if given. Can the
        # value field fail here?
        #
        for k, v in kwargs.items():
            data.extend([k.encode(), np.asarray(v).tobytes()])

        # Is the value cached?
        #
        token = b''.join(data)
        digest = hashfunc(token).digest()
        cache = cls._cache
        if digest in cache:
            cache_ctr['hits'] += 1
            return cache[digest].copy()

        # Evaluate the model.
        #
        vals = func(cls, pars, xlo, *args, **kwargs)

        # remove first item in queue and remove from cache
        queue = cls._queue
        key = queue.pop(0)
        cache.pop(key, None)

        # append newest model values to queue
        queue.append(digest)
        cache[digest] = vals.copy()

        cache_ctr['misses'] += 1

        return vals

    return cache_model


# It is tempting to convert the explicit class names below into calls
# to super(), but this is problematic since it ends up breaking a
# number of invariants the classes rely on. An example is that
# instances of Unary/BinaryOpModel classes should not cache-related
# attributes, but they can do if we change to using super. There is
# more discussion of this in
# https://www.artima.com/weblogs/viewpost.jsp?thread=237121 which
# points out that you should either always use super or never do (or,
# that multiple inheritance is tricky in Python).
#

class Model(NoNewAttributesAfterInit):
    """The base class for Sherpa models.

    A model contains zero or more parameters that control the
    predictions of the model when given one or more coordinates.
    These parameters may represent some variable that describes
    the model, such as the temperature of a black body, or the
    computation, such as what form of interpolation to use.

    Parameters
    ----------
    name : str
        A label for the model instance.
    pars : sequence of sherpa.parameter.Parameter objects
        The parameters of the model.

    Notes
    -----
    Parameters can be accessed via the ``pars`` attribute, but it is
    expected that they will generally be accessed directly, as the
    class provides case-insensitive access to the parameter names
    as object attributes. That is, if the model contains parameters
    called ``breakFreq`` and ``norm``, and the instance is stored in
    the variable ``mdl``, then the following can be used to access
    the parameters::

        print(f"Break frequency = {mdl.breakfreq.val}")
        mdl.norm = 1.2e-3

    """

    ndim: Optional[int] = None
    "The dimensionality of the model, if defined, or None."

    def __init__(self,
                 name: str,
                 pars: Sequence[Parameter] = ()) -> None:
        self.name = name
        self.type = self.__class__.__name__.lower()
        self._pars = tuple(pars)
        self.is_discrete = False
        NoNewAttributesAfterInit.__init__(self)

    @property
    def pars(self) -> tuple[Parameter, ...]:
        """Return the parameters of the model.

        This does not include any linked parameters.

        .. versionchanged:: 4.16.1
           The pars field can no-longer be set directly. Individual
           elements can still be changed.

        See Also
        --------
        lpars

        """

        return tuple(par for par in self._pars)

    @property
    def lpars(self) -> tuple[Parameter, ...]:
        """Return any linked parameters.

        This only returns linked parameters that are not related
        to the model, and each parameter is not repeated.

        .. versionadded:: 4.16.1

        See Also
        --------
        pars

        Examples
        --------

        By default there are no linked parameters:

        >>> from sherpa.models.basic import Gauss2D
        >>> mdl = Gauss2D("mdl")
        >>> len(mdl.pars)
        6
        >>> mdl.lpars
        ()

        Force the model to have identical xpos and ypos parameters.
        Since the linked parameter value (mdl.xpos) is part of the
        model it is not included in `lpars`:

        >>> mdl.ypos = mdl.xpos
        >>> len(mdl.pars)
        6
        >>> mdl.lpars
        ()

        Add a link to allow the sigma term to be fit rather than
        FWHM. Since the linked parameter - here from the Const1D
        model - is not a part of the model it is included in
        `lpars`:

        >>> import numpy as np
        >>> from sherpa.models.basic import Const1D
        >>> sigma = Const1D("sigma")
        >>> mdl.fwhm = 2 * np.sqrt(2 * np.log(2)) * sigma.c0
        >>> len(mdl.pars)
        6
        >>> mdl.lpars
        (<Parameter 'c0' of model 'sigma'>,)

        """

        # Find all the linked parameters, but only report the first
        # occurrence.
        #
        out = []
        for par in self._pars:
            if not par.link:
                continue

            for lpar in expand_par(par.link):
                # This could be a parameter we've already seen: e.g.
                #    mdl.x2 = mdl.x1 + 5
                #
                if lpar in self._pars:
                    continue

                # The parameter could be used in several expressions: e.g.
                #    mdl.x1 = other.c0 + 5
                #    mdl.x2 = other.c0 + 7
                #
                if lpar in out:
                    continue

                out.append(lpar)

        return tuple(out)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} model instance '{self.name}'>"

    def __str__(self) -> str:
        s = self.name
        sep5 = '-' * 5
        sep4 = '-' * 4
        sep3 = '-' * 3
        hfmt = '\n   %-12s %-6s %12s %12s %12s %10s'
        s += hfmt % ('Param', 'Type', 'Value', 'Min', 'Max', 'Units')
        s += hfmt % (sep5, sep4, sep5, sep3, sep3, sep5)
        for p in self.pars:
            if p.hidden:
                continue

            if p.link is not None:
                tp = 'linked'
            elif p.frozen:
                tp = 'frozen'
            else:
                tp = 'thawed'

            s += f'\n   {p.fullname:<12s} {tp:<6s} {p.val:12g} '
            if p.link is not None:
                linkstr = f'expr: {p.link.fullname}'
                s += f'{linkstr:>24s}'
            else:
                s += f'{p.min:12g} {p.max:12g}'

            s += f" {p.units:>10s}"

        return s

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the model
        """
        return html_model(self)

    # This allows all models to be used in iteration contexts, whether or
    # not they're composite
    def __iter__(self) -> Iterator[Model]:
        return iter([self])

    def __getattr__(self, name: str) -> Any:
        """Access to parameters is case insensitive. Other fields are exact."""

        # Any model that defines a parameter will end up looking for
        # the _par_index field, and this is done before Model.__init__
        # is reached (i.e. setting it there does not help), so create
        # it if needed.
        #
        if "_par_index" == name:
            if self.__dict__.get('_par_index') is None:
                self.__dict__['_par_index'] = {}
            return self.__dict__['_par_index']

        lowered_name = name.lower()
        parameter = self._par_index.get(lowered_name)
        if parameter is not None:
            if lowered_name in parameter.aliases:
                wmsg = f'Parameter name {name} is deprecated for model ' + \
                    f'{type(self).__name__}, use {parameter.name} instead'
                warnings.warn(wmsg, DeprecationWarning)

            return parameter

        NoNewAttributesAfterInit.__getattribute__(self, name)

    def __setattr__(self, name: str, val) -> None:
        # case sensitivity is handled by getattr
        par = getattr(self, name, None)
        if isinstance(par, Parameter):
            # When setting an attribute that is a Parameter, set the
            # parameter's value instead.
            par.val = val
            return

        NoNewAttributesAfterInit.__setattr__(self, name, val)

        # If adding a parameter, add the entry to _par_index and
        # handle any aliases.
        #
        if not isinstance(val, Parameter):
            return

        # We use lower-case for parameter access in _par_index.
        vname = val.name.lower()
        self._par_index[vname] = val
        if not val.aliases:
            return

        # Update index of aliases, if necessary. The aliases array
        # is assumed to be in lower case form.
        #
        for alias in val.aliases:
            self._par_index[alias] = val

    def startup(self, cache: bool = False) -> None:
        """Called before a model may be evaluated multiple times.

        Parameters
        ----------
        cache : bool, optional
            Should a cache be used when evaluating the models.

        See Also
        --------
        teardown
        """
        raise NotImplementedError

    def calc(self,
             p: Sequence[SupportsFloat],
             *args,
             **kwargs) -> np.ndarray:
        """Evaluate the model on a grid.

        Parameters
        ----------
        p : sequence of numbers
            The parameter values to use. The order matches the
            ``pars`` field.
        *args
            The model grid. The values can be scalar or arrays,
            and the number depends on the dimensionality of the
            model and whether it is being evaluated over an
            integrated grid or at a point (or points).
        **kwargs
            Any model-specific values that are not parameters.
        """
        raise NotImplementedError

    def teardown(self) -> None:
        """Called after a model may be evaluated multiple times.

        See Also
        --------
        startup
        """
        raise NotImplementedError

    def guess(self, dep, *args, **kwargs):
        """Set an initial guess for the parameter values.

        Attempt to set the parameter values, and ranges, for
        the model to match the data values. This is intended
        as a rough guess, so it is expected that the model
        is only evaluated a small number of times, if at all.
        """
        raise NotImplementedError

    def get_center(self):
        raise NotImplementedError

    def set_center(self, *args, **kwargs):
        raise NotImplementedError

    def __call__(self, *args, **kwargs) -> Union[Model, np.ndarray]:
        # A bit of trickery, to make model creation
        # in IPython happen without raising errors, when
        # model is made automatically callable
        if len(args) == 0 and len(kwargs) == 0:
            return self

        # By accessing the val field of each parameter in self.pars we
        # process all the linked parameters (if there are any) and so
        # do not need to worry about the lpars field.
        #
        return self.calc([p.val for p in self.pars], *args, **kwargs)

    def get_thawed_pars(self) -> list[Parameter]:
        """Return the thawed parameter objects.

        This includes linked parameters, which complicates the min/max
        settings, since the range on the components of a linked
        parameter does not match that of the original parameter, which
        is an issue when the limits are exceeded.

        .. versionadded:: 4.16.1

        """

        pars = [p for p in self.pars if not p.frozen]
        pars.extend(p for p in self.lpars if not p.frozen)
        return pars

    def _get_thawed_par_vals(self) -> list[SupportsFloat]:
        return [p.val for p in self.get_thawed_pars()]

    def _set_thawed_par_vals(self, vals: Sequence[SupportsFloat]) -> None:
        tpars = self.get_thawed_pars()

        ngot = len(vals)
        nneed = len(tpars)
        if ngot != nneed:
            raise ModelErr('numthawed', nneed, ngot)

        # Note that this check ignores the soft limits. However,
        # it sets the limits to min/max and not hard_min/max,
        # which is issue #1980.
        #
        for p, v in zip(tpars, vals):
            v = SherpaFloat(v)
            if v < p.hard_min:
                p.val = p.min
                warning('value of parameter %s is below minimum; '
                        'setting to minimum', p.fullname)
            elif v > p.hard_max:
                p.val = p.max
                warning('value of parameter %s is above maximum; '
                        'setting to maximum', p.fullname)
            else:
                # We do not want to set val directly because we do not
                # want the default field to change.
                #
                p._val = v

        # Check that each linked parameter lies within it's limits
        # (since we can not guarantee it from the limits of the
        # linking parameters).
        #
        # Unfortunately, as noted in #1981, it's not obvious what we
        # should do if a linked parameter is now out of range. For now
        # we just evaluate the parameter value, which will trigger a
        # ParameterErr.  This triggers if the *soft* limits are
        # exceeded, rather than the hard limits, which is slightly
        # different to the above check.
        #
        for par in self.pars:
            if par.link is None:
                continue

            # This relies on the parameter validation logic and we do
            # not care about the return value.
            #
            # We could change par._val but we can not "feed" that
            # value back to the system to know what the linked
            # parameter should be.
            #
            _ = par.val

    thawedpars = property(_get_thawed_par_vals, _set_thawed_par_vals,
                          doc="""The thawed parameters of the model.

Get or set the thawed parameters of the model as a list of
numbers. If there are no thawed parameters then [] is used.
The ordering matches that of the pars attribute.

See Also
--------
thawedparmaxes, thawedparmins
""")

    def _get_thawed_par_mins(self) -> list[SupportsFloat]:
        return [p.min for p in self.get_thawed_pars()]

    def _set_thawed_pars_mins(self, vals: Sequence[SupportsFloat]) -> None:
        tpars = self.get_thawed_pars()

        ngot = len(vals)
        nneed = len(tpars)
        if ngot != nneed:
            raise ModelErr('numthawed', nneed, ngot)

        for p, v in zip(tpars, vals):
            v = SherpaFloat(v)
            if v < p.hard_min:
                p.min = p.hard_min
                warning('value of parameter %s minimum is below '
                        'hard minimum; setting to hard minimum',
                        p.fullname)
            elif v > p.hard_max:
                p.min = p.hard_max
                warning('value of parameter %s minimum is above '
                        'hard maximum; setting to hard maximum',
                        p.fullname)
            else:
                p._min = v

    thawedparmins = property(_get_thawed_par_mins, _set_thawed_pars_mins,
                             doc="""The minimum limits of the thawed parameters.

Get or set the minimum limits of the thawed parameters
of the model as a list of numbers. If there are no
thawed parameters then [] is used. The ordering matches
that of the pars attribute.

See Also
--------
thawedpars, thawedarhardmins, thawedparmaxes
""")

    def _get_thawed_par_maxes(self) -> list[SupportsFloat]:
        return [p.max for p in self.get_thawed_pars()]

    def _set_thawed_pars_maxes(self, vals: Sequence[SupportsFloat]) -> None:
        tpars = self.get_thawed_pars()

        ngot = len(vals)
        nneed = len(tpars)
        if ngot != nneed:
            raise ModelErr('numthawed', nneed, ngot)

        for p, v in zip(tpars, vals):
            v = SherpaFloat(v)
            if v < p.hard_min:
                p.max = p.hard_min
                warning('value of parameter %s maximum is below '
                        'hard minimum; setting to hard minimum',
                        p.fullname)
            elif v > p.hard_max:
                p.max = p.hard_max
                warning('value of parameter %s maximum is above '
                        'hard maximum; setting to hard maximum',
                        p.fullname)
            else:
                p._max = v

    thawedparmaxes = property(_get_thawed_par_maxes, _set_thawed_pars_maxes,
                              doc="""The maximum limits of the thawed parameters.

Get or set the maximum limits of the thawed parameters
of the model as a list of numbers. If there are no
thawed parameters then [] is used. The ordering matches
that of the pars attribute.

See Also
--------
thawedpars, thawedarhardmaxes, thawedparmins
""")

    def _get_thawed_par_hardmins(self) -> list[SupportsFloat]:
        return [p.hard_min for p in self.get_thawed_pars()]

    thawedparhardmins = property(_get_thawed_par_hardmins,
                                 doc="""The hard minimum values for the thawed parameters.

The minimum and maximum range of the parameters can be
changed with thawedparmins and thawedparmaxes but only
within the range given by thawedparhardmins
to thawparhardmaxes.

See Also
--------
thawedparhardmaxes, thawedparmins
""")

    def _get_thawed_par_hardmaxes(self) -> list[SupportsFloat]:
        return [p.hard_max for p in self.get_thawed_pars()]

    thawedparhardmaxes = property(_get_thawed_par_hardmaxes,
                                  doc="""The hard maximum values for the thawed parameters.

The minimum and maximum range of the parameters can be
changed with thawedparmins and thawedparmaxes but only
within the range given by thawedparhardmins
to thawparhardmaxes.

See Also
--------
thawedparhardmins, thawedparmaxes
""")

    # TODO: should this reset linked parameters? Or does a reset clear
    # the link?
    #
    def reset(self) -> None:
        """Reset the parameter values.

        Restores each parameter to the last value it was set to.
        This allows the parameters to be easily reset after a
        fit.
        """

        for p in self.pars:
            p.reset()

    # TODO: should this freeze linked parameters?
    def freeze(self) -> None:
        """Freeze any thawed parameters of the model."""

        for p in self.pars:
            p.freeze()

    # TODO: should this thaw linked parameters?
    def thaw(self) -> None:
        """Thaw any frozen parameters of the model.

        Those parameters that are marked as "always frozen" are
        skipped.

        """

        # Note that we have to handle "always frozen" cases, but rather
        # than asking for permission we just handle the failure case.
        #
        for p in self.pars:
            try:
                p.thaw()
            except ParameterErr:
                continue


class CompositeModel(Model):
    """Represent a model with composite parts.

    This is the base class for representing expressions that combine
    multiple models and values.

    Parameters
    ----------
    name : str
        The name for the collection of models.
    parts : sequence of Model objects
        The models.

    Attributes
    ----------
    parts : sequence of Model

    Notes
    -----
    Composite models can be iterated through to find their
    components:

    >>> from sherpa.models.basic import Gauss1D, Polynom1D
    >>> l1 = Gauss1D('l1')
    >>> l2 = Gauss1D('l2')
    >>> b = Polynom1D('b')
    >>> mdl = l1 + (0.5 * l2) + b
    >>> mdl
    <BinaryOpModel model instance 'l1 + 0.5 * l2 + b'>
    >>> for cpt in mdl:
    ...     print(type(cpt))
    ...
    <class 'sherpa.models.model.BinaryOpModel'>
    <class 'sherpa.models.basic.Gauss1D'>
    <class 'sherpa.models.model.BinaryOpModel'>
    <class 'sherpa.models.model.ArithmeticConstantModel'>
    <class 'sherpa.models.basic.Gauss1D'>
    <class 'sherpa.models.basic.Polynom1D'>

    """

    def __init__(self, name: str, parts: Sequence[Model]) -> None:
        self.parts = tuple(parts)
        allpars = []
        model_with_dim = None
        for part in self.parts:

            ndim = part.ndim
            if ndim is not None:
                if self.ndim is None:
                    self.ndim = ndim
                    model_with_dim = part
                elif self.ndim != ndim:
                    if TYPE_CHECKING:
                        # help the type checker out
                        assert model_with_dim is not None

                    raise ModelErr('Models do not match: ' +
                                   f'{self.ndim}D ({model_with_dim.name}) and ' +
                                   f'{ndim}D ({part.name})')

            for p in part.pars:
                if p in allpars:
                    # If we already have a reference to this
                    # parameter, store a hidden, linked proxy
                    # instead. This is presumably to ensure that we
                    # have the correct number of degrees of freedom
                    # (as pnew is frozen) while still sending the
                    # correct parameters to the different components.
                    #
                    pnew = Parameter(p.modelname, p.name, 0.0, hidden=True)
                    pnew.link = p
                    p = pnew

                allpars.append(p)

        Model.__init__(self, name, allpars)

        # This check should probably be removed since it refers to
        # loading a very-old session, and that is unlikely to work due
        # to other changes in the code base.
        #
        for part in self.parts:
            try:
                self.is_discrete = self.is_discrete or part.is_discrete
            except AttributeError:
                warning("Could not determine whether the model is discrete.\n" +
                        "This probably means that you have restored a session saved with a previous version of Sherpa.\n" +
                        "Falling back to assuming that the model is continuous.\n")
                self.is_discrete = False

    def __iter__(self) -> Iterator[Model]:
        return iter(self._get_parts())

    def _get_parts(self) -> list[Model]:
        parts = []

        for p in self.parts:
            # A CompositeModel should not hold a reference to itself
            assert (p is not self), f"'{type(self).__name__}' " + \
                "object holds a reference to itself"

            # Including itself seems a bit strange if it's a CompositeModel
            # but is used by sherpa.astro.instrument.has_pha_instance (and
            # possibly elsewhere).
            #
            parts.append(p)
            if isinstance(p, CompositeModel):
                parts.extend(p._get_parts())

        # FIXME: do we want to remove duplicate components from parts?

        return parts

    def startup(self, cache: bool = False) -> None:
        pass

    def teardown(self) -> None:
        pass

    def cache_clear(self) -> None:
        """Clear the cache for each component."""
        for p in self.parts:
            try:
                p.cache_clear()
            except AttributeError:
                pass

    def cache_status(self) -> None:
        """Display the cache status of each component.

        Information on the cache - the number of "hits", "misses", and
        "requests" - is displayed at the INFO logging level.

        Example
        -------

        >>> mdl.cache_status()
         xsphabs.gal                size:    5  hits:   715  misses:   158  check=  873
         powlaw1d.pl                size:    5  hits:   633  misses:   240  check=  873

        """
        for p in self.parts:
            try:
                p.cache_status()
            except AttributeError:
                pass

    def guess(self, dep, *args, **kwargs):
        """Call guess on each component.

        At the moment there is no recognition of the full
        model expression - e.g. cpt1 * cpt2 and cpt1 + cpt2
        would ideally have different scalings applied here.

        .. versionchanged:: 4.17.0
           Prior to 4.17.0 the guess method could not be called on
           composite models.

        """

        # It would be good to apply the various "guess" routines
        # to the combined model expression, but that would
        # significantly complicate the analysis.
        #
        seen = set()
        for p in self.parts:
            # Only call guess the first time a model is seen.
            #
            if p in seen:
                continue

            try:
                p.guess(dep, *args, **kwargs)
            except (AttributeError, NotImplementedError):
                # Skip those models without a guess. This could just
                # be a call to pass, but follow the behaviour of
                # sherpa.ui.utils.guess. This is not ideal as it is
                # better decided at the ui layer than here.
                #
                warning('No guess found for %s', p.name)

            seen.add(p)


class SimulFitModel(CompositeModel):
    """Store multiple models.

    This class is for use with sherpa.data.DataSimulFit.

    Parameters
    ----------
    name : str
        The name for the collection of models.
    parts : sequence of Model objects
        The models.

    Attributes
    ----------
    parts : sequence of Model

    See Also
    --------
    sherpa.data.DataSimulFit

    Examples
    --------

    >>> from sherpa.models.basic import Gauss1D, Polynom1D
    >>> m1 = Polynom1D('m1')
    >>> m2 = Gauss1D('g1')
    >>> mall = SimulFitModel('comp', (m1, m1 + m2))

    If dall is a DataSimulFit object then the model components
    can be evaluated for the composite object using:

    >>> ymdl = dall.eval_model_to_fit(mall)

    """

    def __iter__(self) -> Iterator[Model]:
        return iter(self.parts)

    # Why is this not defined in CompositeModel?
    #
    def startup(self, cache: bool = False) -> None:
        for part in self:
            part.startup(cache)
        CompositeModel.startup(self, cache)

    def teardown(self) -> None:
        for part in self:
            part.teardown()
        CompositeModel.teardown(self)


# TODO: what benefit does this provide versus just using the number?
# I guess it does simplify any attempt to parse the components of
# a model expression.
#
class ArithmeticConstantModel(Model):
    """Represent a constant value, or values.

    Parameters
    ----------
    val : number or sequence
        The number, or numbers, to store.
    name : str or None, optional
        The display name. If not set the value is used when the value
        is a scalar, otherwise it indicates an array of elements.
    """

    def __init__(self,
                 # Use SupportsIndex rather than Sequence[SupportsFloat] to
                 # avoid mypy warnings.
                 val: Union[SupportsFloat, SupportsIndex],
                 name: Optional[str] = None) -> None:

        store = SherpaFloat(val)
        if name is None:
            if np.isscalar(store):
                name = str(store)
            else:
                if TYPE_CHECKING:
                    assert isinstance(val, np.ndarray)

                # For some reason mypy didn't like
                # '[str(s) for s in store.shape]`
                dims = map(str, store.shape)
                nstr = ','.join(dims)
                name = f'{store.dtype.name}[{nstr}]'

        self.name = name
        self.val = store

        # store has to be a scalar or 1D array, even if used with a 2D
        # model, due to the way model evaluation works, so as we
        # can't easily define the dimensionality of this model, we
        # remove any dimensionality checking for this class.
        #
        self.ndim = None

        Model.__init__(self, self.name)

    def _get_val(self) -> Union[SherpaFloat, np.ndarray]:
        return self._val

    def _set_val(self, val:Union[SupportsFloat, SupportsIndex]) -> None:
        val = SherpaFloat(val)
        if val.ndim > 1:
            raise ModelErr('The constant must be a scalar or 1D, not 2D')

        self._val = val

    val = property(_get_val, _set_val,
                   doc='The constant value (scalar or 1D).')

    def startup(self, cache: bool = False) -> None:
        pass

    # This doesn't match superclass, as we can return a scalar here
    # and the superclass assumes it returns a ndarray, so we do not
    # type this routine.
    #
    def calc(self, p, *args, **kwargs):
        return self.val

    def teardown(self) -> None:
        pass

    def guess(self, dep, *args, **kwargs):
        """There is nothing to guess here."""
        pass


def _make_unop(op: Callable, opstr: str) -> Callable:
    def func(self):
        return UnaryOpModel(self, op, opstr)
    return func


def _make_binop(op: Callable, opstr: str) -> tuple[Callable, Callable]:
    def func(self, rhs):
        return BinaryOpModel(self, rhs, op, opstr)

    def rfunc(self, lhs):
        return BinaryOpModel(lhs, self, op, opstr)

    return (func, rfunc)


class ArithmeticModel(Model):
    """Support combining model expressions and caching results."""

    cache = 5
    """The maximum size of the cache."""

    def __init__(self,
                 name: str,
                 pars: Sequence[Parameter] = ()) -> None:
        self.integrate = True

        # Model caching ability
        self.cache = 5  # repeat the class definition
        self._use_caching = True  # FIXME: reduce number of variables?
        self.cache_clear()
        Model.__init__(self, name, pars)

    def cache_clear(self) -> None:
        """Clear the cache."""
        # It is not obvious what to set the queue length to
        self._queue = ['']

        self._cache: dict[bytes, np.ndarray] = {}
        self._cache_ctr = {'hits': 0, 'misses': 0, 'check': 0}

    def cache_status(self) -> None:
        """Display the cache status.

        Information on the cache - the number of "hits", "misses", and
        "requests" - is displayed at the INFO logging level.

        Example
        -------

        >>> pl.cache_status()
         powlaw1d.pl                size:    5  hits:   633  misses:   240  check=  873

        """
        c = self._cache_ctr
        info(f" {self.name:25s}  size: {len(self._queue):4d}  " +
             f"hits: {c['hits']:5d}  misses: {c['misses']:5d}  " +
             f"check: {c['check']:5d}")

    # Unary operations
    __neg__ = _make_unop(np.negative, '-')
    __pos__ = _make_unop(np.positive, '+')
    __abs__ = _make_unop(np.absolute, 'abs')

    # Binary operations
    __add__, __radd__ = _make_binop(np.add, '+')
    __sub__, __rsub__ = _make_binop(np.subtract, '-')
    __mul__, __rmul__ = _make_binop(np.multiply, '*')
    __div__, __rdiv__ = _make_binop(np.divide, '/')
    __floordiv__, __rfloordiv__ = _make_binop(np.floor_divide, '//')
    __truediv__, __rtruediv__ = _make_binop(np.true_divide, '/')
    __mod__, __rmod__ = _make_binop(np.remainder, '%')
    __pow__, __rpow__ = _make_binop(np.power, '**')

    def __setstate__(self, state):
        self.__dict__.update(state)

        if '_use_caching' not in state:
            self.__dict__['_use_caching'] = True

        if '_queue' not in state:
            self.__dict__['_queue'] = ['']

        if '_cache' not in state:
            self.__dict__['_cache'] = {}
            self.__dict__['_cache_ctr'] = {'hits': 0, 'misses': 0, 'check': 0}

        if 'cache' not in state:
            self.__dict__['cache'] = 5

    def __getitem__(self, filter):
        return FilterModel(self, filter)

    def startup(self, cache: bool = False) -> None:
        self.cache_clear()
        self._use_caching = cache
        if int(self.cache) <= 0:
            return

        self._queue = [''] * int(self.cache)
        frozen = np.array([par.frozen for par in self.pars], dtype=bool)
        if len(frozen) > 0 and frozen.all():
            self._use_caching = cache

    def teardown(self) -> None:
        self._use_caching = False

    def apply(self, outer, *otherargs, **otherkwargs):
        return NestedModel(outer, self, *otherargs, **otherkwargs)


class RegriddableModel(ArithmeticModel):
    """Support models that can be evaluated on a different grid.

    """

    def regrid(self, *args, **kwargs):
        """Allow a model to be evaluated on a different grid than requested.

        The return value is a new instance of the model, set up to
        evaluate the model on the supplied axes which will be
        regridded onto the requested grid.

        """
        raise NotImplementedError


class RegriddableModel1D(RegriddableModel):
    """Allow 1D models to be regridded."""

    ndim = 1
    "A one-dimensional model."

    def regrid(self, *args, **kwargs):
        """
        The class RegriddableModel1D allows the user to evaluate in the
        requested space then interpolate onto the data space. An optional
        argument 'interp' enables the user to change the interpolation method.

        Examples
        --------
        >>> import numpy as np
        >>> from sherpa.models.basic import Box1D
        >>> from sherpa.utils import linear_interp
        >>> mybox = Box1D()
        >>> request_space = np.arange(1, 10, 0.1)
        >>> regrid_model = mybox.regrid(request_space, interp=linear_interp)
        """
        valid_keys = ('interp',)
        for key in kwargs.keys():
            if key not in valid_keys:
                raise TypeError(f"unknown keyword argument: '{key}'")
        eval_space = EvaluationSpace1D(*args)
        regridder = ModelDomainRegridder1D(eval_space, **kwargs)
        regridder._make_and_validate_grid(args)
        return regridder.apply_to(self)


class RegriddableModel2D(RegriddableModel):
    """Allow 2D models to be regridded."""

    ndim = 2
    "A two-dimensional model."

    def regrid(self, *args, **kwargs):
        eval_space = EvaluationSpace2D(*args)
        regridder = ModelDomainRegridder2D(eval_space)
        return regridder.apply_to(self)


class UnaryOpModel(CompositeModel, ArithmeticModel):
    """Apply an operator to a model expression.

    Parameters
    ----------
    arg : Model instance
        The expression.
    op : function reference
        The ufunc to apply to the model values.
    opstr : str
        The symbol used to represent the operator.

    Attributes
    ----------
    arg : Model instance
        The model.
    op : function reference
        The ufunc to apply to the model values.
    opstr : str
        The symbol used to represent the operator.

    See Also
    --------
    BinaryOpModel

    Examples
    --------

    >>> import numpy
    >>> from sherpa.models.basic import Gauss1D
    >>> m1 = Gauss1D()
    >>> m2 = UnaryOpModel(m1, numpy.negative, '-')

    """

    @staticmethod
    def wrapobj(obj) -> Model:
        return _wrapobj(obj, ArithmeticConstantModel)

    def __init__(self,
                 arg: Any,
                 op: Callable,
                 opstr: str) -> None:
        self.arg = self.wrapobj(arg)
        self.op = op
        self.opstr = opstr
        self.opprec = get_precedences_op(op)[0]

        # We do not simplify this (e.g. remove brackets if self.arg
        # is not a composite model).
        #
        name = f'{opstr}({self.arg.name})'
        CompositeModel.__init__(self, name, (self.arg,))

    def calc(self, p: Sequence[SupportsFloat],
             *args, **kwargs) -> np.ndarray:
        return self.op(self.arg.calc(p, *args, **kwargs))


class BinaryOpModel(CompositeModel, RegriddableModel):

    """Combine two model expressions.

    Parameters
    ----------
    lhs : Model instance
        The left-hand side of the expression.
    rhs : Model instance
        The right-hand side of the expression.
    op : function reference
        The ufunc which combines two array values.
    opstr : str
        The symbol used to represent the operator.

    Attributes
    ----------
    lhs : Model instance
        The left-hand sides of the expression.
    rhs : Model instance
        The right-hand sides of the expression.
    op : function reference
        The ufunc which combines two array values.
    opstr : str
        The symbol used to represent the operator.

    See Also
    --------
    UnaryOpModel

    Examples
    --------

    >>> import numpy
    >>> from sherpa.models.basic import Gauss1D, Polynom1D
    >>> m1 = Gauss1D()
    >>> m2 = Polynom1D()
    >>> m = BinaryOpModel(m1, m2, numpy.add, '+')

    """

    @staticmethod
    def wrapobj(obj) -> Model:
        return _wrapobj(obj, ArithmeticConstantModel)

    def __init__(self,
                 lhs: Any,
                 rhs: Any,
                 op: Callable,
                 opstr: str) -> None:
        self.lhs = self.wrapobj(lhs)
        self.rhs = self.wrapobj(rhs)
        self.op = op
        self.opstr = opstr

        p, a =  get_precedences_op(op)
        self.opprec = p

        # Simplify the expression if possible.
        #
        lp = get_precedence_expr(self.lhs)
        rp = get_precedence_expr(self.rhs)

        lstr = get_precedence_lhs(self.lhs.name, lp, p, a)
        rstr = get_precedence_rhs(self.rhs.name, opstr, rp, p)

        name = f'{lstr} {opstr} {rstr}'
        CompositeModel.__init__(self, name, (self.lhs, self.rhs))

    def regrid(self, *args, **kwargs):
        for part in self.parts:
            # ArithmeticConstantModel does not support regrid by design
            if not hasattr(part, 'regrid'):
                continue
            # The full model expression must be used
            return part.__class__.regrid(self, *args, **kwargs)
        raise ModelErr('Neither component supports regrid method')

    def startup(self, cache: bool = False) -> None:
        self.lhs.startup(cache)
        self.rhs.startup(cache)
        CompositeModel.startup(self, cache)

    def teardown(self) -> None:
        self.lhs.teardown()
        self.rhs.teardown()
        CompositeModel.teardown(self)

    def calc(self, p: Sequence[SupportsFloat],
             *args, **kwargs) -> np.ndarray:
        # Note that the kwargs are sent to both model components.
        #
        nlhs = len(self.lhs.pars)
        lhs = self.lhs.calc(p[:nlhs], *args, **kwargs)
        rhs = self.rhs.calc(p[nlhs:], *args, **kwargs)
        try:
            val = self.op(lhs, rhs)
        except ValueError as ve:
            raise ValueError("shape mismatch between " +
                             f"'{type(self.lhs).__name__}: {len(lhs)}' and " +
                             f"'{type(self.rhs).__name__}: {len(rhs)}'") from ve
        return val


# TODO: do we actually make use of this functionality anywhere?
# We only have 1 test that checks this class, and it is an existence
# test (check that it works), not that it is used anywhere.
#
class FilterModel(CompositeModel, ArithmeticModel):

    def __init__(self, model, filter):
        self.model = model
        self.filter = filter

        if isinstance(filter, tuple):
            filter_str = ','.join([self._make_filter_str(f) for f in filter])
        else:
            filter_str = self._make_filter_str(filter)

        CompositeModel.__init__(self,
                                f'({self.model.name})[{filter_str}]',
                                (self.model,))

    @staticmethod
    def _make_filter_str(filter):
        if not isinstance(filter, slice):
            if filter is Ellipsis:
                return '...'
            return str(filter)

        s = ''
        if filter.start is not None:
            s += str(filter.start)
        s += ':'
        if filter.stop is not None:
            s += str(filter.stop)
        if filter.step is not None:
            s += f':{filter.step}'

        return s

    def calc(self, p, *args, **kwargs):
        return self.model.calc(p, *args, **kwargs)[self.filter]


class ArithmeticFunctionModel(Model):
    """Represent a callable function.

    Parameters
    ----------
    func : function reference
        A callable function. It is called with the grid arguments and
        any keyword arguments sent to calc(), but not the model
        parameter values.

    Attributes
    ----------
    func : function reference

    """

    def __init__(self, func: Callable) -> None:
        if isinstance(func, Model):
            raise ModelErr('badinstance', type(self).__name__)
        if not callable(func):
            raise ModelErr('noncall', type(self).__name__, type(func).__name__)
        self.func = func
        Model.__init__(self, func.__name__)

    def calc(self, p: Sequence[SupportsFloat],
             *args, **kwargs) -> np.ndarray:
        return self.func(*args, **kwargs)

    def startup(self, cache: bool = False) -> None:
        pass

    def teardown(self) -> None:
        pass

    def guess(self, dep, *args, **kwargs):
        """There is nothing to guess here."""
        pass


class NestedModel(CompositeModel, ArithmeticModel):
    """Apply a model to the results of a model.

    Parameters
    ----------
    outer : Model instance
        The model to apply second.
    inner : Model instance
        The model to apply first.
    *otherargs
        Arguments that are to be applied to the outer model.
    **otherkwargs
        Keyword arguments that are to be applied to the outer model.

    Attributes
    ----------
    outer : Model instance
        The outer model.
    inner : Model instance
        The inner model.
    otherargs
        Arguments that are to be applied to the outer model.
    otherkwargs
        Keyword arguments that are to be applied to the outer model.

    """

    @staticmethod
    def wrapobj(obj) -> Model:
        return _wrapobj(obj, ArithmeticFunctionModel)

    def __init__(self, outer, inner, *otherargs, **otherkwargs) -> None:
        self.outer = self.wrapobj(outer)
        self.inner = self.wrapobj(inner)
        self.otherargs = otherargs
        self.otherkwargs = otherkwargs
        CompositeModel.__init__(self, f'{self.outer.name}({self.inner.name})',
                                (self.outer, self.inner))

    def startup(self, cache: bool = False) -> None:
        self.inner.startup(cache)
        self.outer.startup(cache)
        CompositeModel.startup(self, cache)

    def teardown(self) -> None:
        self.inner.teardown()
        self.outer.teardown()
        CompositeModel.teardown(self)

    def calc(self, p: Sequence[SupportsFloat],
             *args, **kwargs) -> np.ndarray:
        nouter = len(self.outer.pars)
        return self.outer.calc(p[:nouter],
                               self.inner.calc(p[nouter:], *args, **kwargs),
                               *self.otherargs, **self.otherkwargs)


# TODO: can we remove this as it is now unused?
#
class MultigridSumModel(CompositeModel, ArithmeticModel):

    def __init__(self, models: Sequence[Model]) -> None:
        self.models = tuple(models)
        arg = ','.join([m.name for m in models])
        name = f'{type(self).__name__}({arg})'
        CompositeModel.__init__(self, name, self.models)

    # This does not match the superclass so do not type it
    def calc(self, p, arglist):
        vals = []
        for model, args in zip(self.models, arglist):
            # FIXME: we're not using p here (and therefore assuming that the
            # parameter values have already been updated to match the contents
            # of p)
            vals.append(model(*args))
        return sum(vals)


class RegridWrappedModel(CompositeModel, ArithmeticModel):

    def __init__(self, model, wrapper: Model) -> None:
        self.model = self.wrapobj(model)
        self.wrapper = wrapper

        if hasattr(model, 'integrate'):
            self.wrapper.integrate = model.integrate

        CompositeModel.__init__(self,
                                f"{self.wrapper.name}({self.model.name})",
                                (self.model, ))

    def calc(self, p: Sequence[SupportsFloat],
             *args, **kwargs) -> np.ndarray:
        return self.wrapper.calc(p, self.model.calc, *args, **kwargs)

    def get_center(self):
        return self.model.get_center()

    def set_center(self, *args, **kwargs):
        return self.model.set_center(*args, **kwargs)

    def guess(self, dep, *args, **kwargs):
        return self.model.guess(dep, *args, **kwargs)

    @property
    def grid(self):
        return self.wrapper.grid

    @grid.setter
    def grid(self, value):
        self.wrapper.grid = value

    @property
    def evaluation_space(self):
        return self.wrapper.evaluation_space

    @staticmethod
    def wrapobj(obj) -> Model:
        # TODO: choice of ArithmeticConstandModel or
        #       ArithmeticFunctionModel?
        return _wrapobj(obj, ArithmeticFunctionModel)


def _wrapobj(obj, wrapper: Callable[..., Model]) -> Model:
    """Wrap an object with the wrapper if needed.

    Parameters
    ----------
    obj
        The input object.
    wrapper
        The wrapper class which is applied to obj if needed.

    Returns
    -------
    newobj
        This is either obj or wrapper(obj).

    """

    # This has been placed outside the classes as
    # the full list of classes are needed to be accessible
    # when called.
    #
    if isinstance(obj, (ArithmeticModel,
                        ArithmeticConstantModel,
                        ArithmeticFunctionModel)):
        return obj

    return wrapper(obj)


# Simple model "deconstruction" - that is, turn a model expression
# like "mdl1 * (mdl2 + mdl3)" into "mdl1 * mdl2" and "mdl1 * mdl3".
#
# The code could be re-worked to use continuation-passing style to
# avoid the recursion error, but that would significantly complicate
# the code. For example
# https://www.tweag.io/blog/2023-01-19-fp2-dial-m-for-monoid/#what-the-thunk
# For now we try to catch any recursion errors and treat the model
# as a "singleton" for that component.
#
def model_deconstruct(model: Model) -> list[Model]:
    """Separate models into additive components.

    Identify the separate "additive" components in the model
    expression, that is the terms separated by addition or
    subtraction. The resulting list can be summed together to recreate
    the original model expression. This is not guaranteed to create
    the expand the expression completely, and the resulting terms are
    not guaranteed to be "simplified".

    .. versionadded:: 4.16.1

    Parameters
    ----------
    model : Model instance
       The model expression to separate.

    Returns
    -------
    terms : list of Model instances
       The separated terms (this list may contain a single
       element).

    Notes
    -----

    If the model includes a subtracted component, such as:

        a - b

    then the b component will be negated in the output. This negation
    is applied directly rather than applied to any term combined with
    it (such as "a * (b - c)", which will create terms

        a * b
        a * -(c)

    While an expression like

        -(a + b * c + d)

    could be split up, at present it is not. However, when written as
    part of a binary expression it is split up, so

        x - (a + b * c + d)

    is split into

        x, -(a), -(b * x), -(d)

    Examples
    --------

    >>> from sherpa.models.basic import Box1D, Gauss1D
    >>> b1 = Box1D("b1")
    >>> g1 = Gauss1D("g1")
    >>> g2 = Gauss1D("g2")
    >>> model_deconstruct(b1)
    [<Box1D model instance 'b1'>]
    >>> model_deconstruct(b1 * (g1 + g2))
    [<BinaryOpModel model instance 'b1 * g1'>, <BinaryOpModel model instance 'b1 * g2'>]
    >>> model_deconstruct(b1 * (g1 - g2))
    [<BinaryOpModel model instance 'b1 * g1'>, <BinaryOpModel model instance 'b1 * -(g2)'>]

    Unary operators are not expanded by this routine, which makes the
    behaviour a bit different than the binary-operator case, as shown
    below:

    >>> model_deconstruct(-(g1 + b1 * g2))
    [<UnaryOpModel model instance '-(g1 + b1 * g2)'>]

    >>> model_deconstruct(b1 - (g1 + b1 * g2))
    [<Box1D model instance 'b1'>, <UnaryOpModel model instance '-(g1)'>, <UnaryOpModel model instance '-(b1 * g2)'>]

    >>> model_deconstruct(-g1 - g2)
    [<UnaryOpModel model instance '-(g1)'>, <UnaryOpModel model instance '-(g2)'>]

    """

    # The return is a list of components that add together to match
    # the input model. So, if we are not a binary operation the result
    # is simple.
    #
    # Note that there is no attempt to simplify an expression, that is,
    #
    #     -a * -b
    #
    # is not changed to
    #
    #     a * b
    #
    if not isinstance(model, BinaryOpModel):
        return [model]

    # The idea is to separate terms for model expressions that users
    # are expected to write, not to separate all possible
    # expressions. This means that the code forgoes the possibility of
    # certain cases to simplify things. Also, the binary operator
    # could technically be anything, and all that is checked here are
    # the NumPy versions.
    #
    # Note that
    #     (a + b) / c -> a / b, a / c
    # but
    #     (a + b) // c -> unchanged
    #
    if model.op not in [np.add, np.multiply,
                        np.subtract, np.divide]:
        return [model]

    try:
        lhs = model_deconstruct(model.lhs)
    except RecursionError:
        # If we can not recurse treat this as a single term.
        lhs = [model.lhs]

    # Special case division, since we do only want to deconstruct the
    # lhs for this case.
    #
    if model.op == np.divide:
        return [BinaryOpModel(lterm, model.rhs, np.divide, '/')
                for lterm in lhs]

    try:
        rhs = model_deconstruct(model.rhs)
    except RecursionError:
        # If we can not recurse treat this as a single term.
        rhs = [model.rhs]

    if model.op == np.multiply:
        # Expand the term, so
        #
        #     (a + b) * (c - d)
        #
        # will go to
        #
        #     (*, a, c)
        #     (*, a, -d)
        #     (*, b, c)
        #     (*, b, -d)
        #
        return [BinaryOpModel(lterm, rterm, np.multiply, '*')
                for lterm, rterm in itertools.product(lhs, rhs)]

    # The terms are intended to be summed together, so for subtraction
    # the RHS terms must be negated.
    #
    if model.op == np.subtract:
        # Since we claim that model is a Model and not ArithmeticModel
        # we explicitly create the UnaryOpModel term (rather than just
        # say "-term").
        #
        rhs = [UnaryOpModel(term, np.negative, '-')
               for term in rhs]

    # This is either addition or subtraction, so we return the
    # combined list.
    #
    lhs.extend(rhs)
    return lhs


# Notebook representation
#
def modelcomponents_to_list(model: Model) -> list[Model]:
    if hasattr(model, 'parts'):
        modellist = []
        for p in model.parts:
            modellist.extend(modelcomponents_to_list(p))
        return modellist

    return [model]


def html_model(mdl: Model) -> str:
    """Construct the HTML to display the model."""

    # Note that as this is a specialized table we do not use
    # formatting.html_table but create everything directly.
    #
    complist = []
    nrows = []
    for comp in modelcomponents_to_list(mdl):
        this_comp_nrows = 0
        for par in comp.pars:
            if par.hidden:
                continue

            this_comp_nrows +=1

        if this_comp_nrows > 0:
            complist.append(comp)
            nrows.append(this_comp_nrows)

    out = '<table class="model">'
    expr = formatting.clean_bracket(mdl.name)
    out += f'<caption>Expression: {expr}</caption>'
    out += '<thead><tr>'
    cols = ['Component', 'Parameter', 'Thawed', 'Value',
            'Min', 'Max', 'Units']
    for col in cols:
        out += f'<th>{col}</th>'

    out += '</tr></thead><tbody>'

    for mcount, (n, comp) in enumerate(zip(nrows, complist)):
        for i, par in enumerate(comp.pars):
            style = '' if ((i > 0) or (mcount == 0)) else ' class="block"'

            if par.hidden:
                continue

            def addtd(val):
                "Use the parameter to convert to HTML"
                return f'<td>{par._val_to_html(val)}</td>'

            out += f'<tr{style}>'
            if i == 0:
                cls = "model-" + ("even" if mcount % 2 == 1 else "odd")
                out += f'<th class="{cls}" scope="rowgroup" '
                out += f'rowspan={n}>{comp.name}</th>'

            out += f'<td>{par.name}</td>'

            if par.link is not None:
                out += "<td>linked</td>"
            else:
                out += '<td><input disabled type="checkbox"'
                if not par.frozen:
                    out += ' checked'
                out += '></input></td>'

            out += addtd(par.val)

            if par.link is not None:
                # 8592 is single left arrow
                # 8656 is double left arrow
                #
                linkstr = formatting.clean_bracket(par.link.fullname)
                out += f'<td colspan=2>&#8656; {linkstr}</td>'
            else:
                out += addtd(par.min)
                out += addtd(par.max)

            out += f'<td>{par._units_to_html()}</td>'
            out += '</tr>'

    out += '</tbody></table>'

    ls = ['<details open><summary>Model</summary>' + out + '</details>']
    return formatting.html_from_sections(mdl, ls)
