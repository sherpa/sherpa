#
#  Copyright (C) 2010, 2016, 2017, 2018, 2019, 2020, 2021
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

A single model is defined by the parameters of the model - represented
as sherpa.models.model.Parameter instances - and the function that
takes the parameter values along with an array of grid values. The
main classes are:

* Model which is the base class and defines most of the interfaces.

* ArithmeticConstantModel and ArithmeticFunctionModel for representing
  a constant value or a function.

* ArithmeticModel is the main base class for deriving user models since
  it supports combining models (e.g. by addition or multiplication) and
  a cache to reduce evaluation time at the expense of memory use.

* RegriddableModel builds on ArithmeticModel to allow a model to be
  evaluated on a different grid to that requested: most model classes
  are derived from the 1D and 2D variants of RegriddableModel.

* CompositeModel which is used to represent a model expression, that
  is combined models, such as `m1 * (m2 + m3)`

  * UnaryOpModel for model expressions such as `- m1`.

  * BinaryOpModel for model expressions such as `m1 + m2`.

  * NestedModel for applying one model to another.

* SimulFitModel for fitting multiple models and datasets.

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

The parameters are instances of the sherpa.models.parameter.Parameter
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
process called "linking". The lniked is no-longer considered a free
parameter in a fit since it's value is derived from the other
parameters. This link can be a simple one-to-one case, such as
ensuring the fwhm parameter of one model is the same as the other:

    >>> m2.fwhm = m1.fwhm

It can be more complex, such as ensuring the position of one line
is a fixed distance from another:

    >>> l2.pos = l1.pos + 23.4

It can even include multiple parameters:

    >>> l3.ampl = (l1.ampl + l2.ampl) / 2

Requesting the parameter value will return the evaluated expression,
and the expression is stored in the link attribute:

    >>> l1.ampl = 10
    >>> l2.ampl = 12
    >>> l3.ampl.val
    11.0
    >>> l3.ampl.link
    <BinaryOpParameter '((l1.ampl + l2.ampl) / 2)'>

The string representation of the model changes for linked parameters
to indicate the expression:

    >>> print(l3)
    l3
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       l3.fwhm      thawed           10  1.17549e-38  3.40282e+38
       l3.pos       thawed            0 -3.40282e+38  3.40282e+38
       l3.ampl      linked           11 expr: ((l1.ampl + l2.ampl) / 2)

Model evaluation
================

With a sherpa.data.Data instance a model can be evaluated with the
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
integrate flag of the model. The Gauss1D model does support an
integrated mode, so switching the integrate flag will change the model
output:

    >>> print(mdl.integrate)
    True
    >>> mdl.integrate = False
    >>> y2 = mdl(xlo, xhi)
    >>> print(y2)
    [48.63274737 49.65462477 49.91343163 50.         49.65462477]

The behavior when the integrate flag is False depends on the model but
it normally just uses the low edge, as shown for the Gauss1D case:

    >>> y3 = mdl(xlo)
    >>> print(y2 == y3)
    [ True  True  True  True  True]

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

The CompositeModel class is the base class for creating model
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

This model is written so that the amplitude of the `l2` component is
half the `l1` component by linking the two `ampl` parameters and then
including a scaling factor in the model expression for `l2`. An
alternative would have been to include this scaling factor in the link
expression:

    >>> l2.ampl = l1.ampl / 2

Model cache
===========

The ArithmeticModel class and modelCacher1d decorator provide basic
support for caching one-dimensional model evaluations - that is, to
avoid re-calculating the model. The idea is to save the results of the
latest calls to a model and return the values from the cache,
hopefully saving time at the expense of using more memory. This is
most effective when the same model is used with multiple datasets
which all have the same grid.

The _use_caching attribute of the model is used to determine whether
the cache is used, but this setting can be over-ridden by the startup
method, which is automatically called by the fit and est_errors
methods of a sherpa.fit.Fit object.

The cache_clear and cache_status methods of ArithmeticModel and
CompositeModel allow you to clear the cache and display to the
standard output the cache status of each model component.

Example
=======

The following class implements a simple scale model which has a single
parameter (`scale`) which defaults to 1. It can be used for both
non-integrated and integrated datasets of any dimensionality (see
sherpa.models.basic.Scale1D and sherpa.models.basic.Scale2D)::

    class ScaleND(ArithmeticModel):
        '''A constant value per element.'''

        def __init__(self, name='scalend'):
            self.scale = Parameter(name, 'scale', 1)
            self.integrate = False
            pars = (self.scale, )
            ArithmeticModel.__init__(self, name, pars)

        def calc(self, *args, **kwargs):
            return self.scale.val * np.ones(len(args[0]))

"""


import functools
import logging
import warnings

import numpy

from sherpa.models.regrid import EvaluationSpace1D, ModelDomainRegridder1D, EvaluationSpace2D, ModelDomainRegridder2D
from sherpa.utils import SherpaFloat, NoNewAttributesAfterInit
from sherpa.utils.err import ModelErr
from sherpa.utils import formatting

from .parameter import Parameter

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


def boolean_to_byte(boolean_value):
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


def modelCacher1d(func):
    """A decorater to cache 1D ArithmeticModel evalutions.

    Apply to the `calc` method of a 1D model to allow the model
    evaluation to be cached. The decision is based on the
    `_use_caching` attribute of the cache along with the `integrate`
    setting, the evaluation grid, and parameter values.

    Example
    -------

    Allow `MyModel` model evaluations to be cached::

        def MyModel(ArithmeticModel):
            ...
            @modelCacher1d
            def calc(self, *args, **kwargs):
                ...

    """

    @functools.wraps(func)
    def cache_model(cls, pars, xlo, *args, **kwargs):
        use_caching = cls._use_caching
        cache = cls._cache
        cache_ctr = cls._cache_ctr
        queue = cls._queue

        # Counts all accesses, even those that do not use the cache.
        cache_ctr['check'] += 1

        digest = ''
        if use_caching:

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

            data = [numpy.array(pars).tobytes(),
                    boolean_to_byte(integrate),
                    numpy.asarray(xlo).tobytes()]
            if args:
                data.append(numpy.asarray(args[0]).tobytes())

            token = b''.join(data)
            digest = hashfunc(token).digest()
            if digest in cache:
                cache_ctr['hits'] += 1
                return cache[digest].copy()

        vals = func(cls, pars, xlo, *args, **kwargs)

        if use_caching:
            # remove first item in queue and remove from cache
            key = queue.pop(0)
            cache.pop(key, None)

            # append newest model values to queue
            queue.append(digest)
            cache[digest] = vals.copy()

            cache_ctr['misses'] += 1

        return vals

    return cache_model


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

    ndim = None
    "The dimensionality of the model, if defined, or None."

    def __init__(self, name, pars=()):
        self.name = name
        self.type = self.__class__.__name__.lower()
        self.pars = tuple(pars)
        self.is_discrete = False
        NoNewAttributesAfterInit.__init__(self)

    def __repr__(self):
        return "<%s model instance '%s'>" % (type(self).__name__, self.name)

    def __str__(self):
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

            if tp == 'linked':
                linkstr = 'expr: %s' % p.link.fullname
                s += ('\n   %-12s %-6s %12g %24s %10s' %
                      (p.fullname, tp, p.val, linkstr, p.units))
            else:
                s += ('\n   %-12s %-6s %12g %12g %12g %10s' %
                      (p.fullname, tp, p.val, p.min, p.max, p.units))
        return s

    def _repr_html_(self):
        """Return a HTML (string) representation of the model
        """
        return html_model(self)

    # This allows all models to be used in iteration contexts, whether or
    # not they're composite
    def __iter__(self):
        return iter([self])

    def __getattr__(self, name):
        """Access to parameters is case insensitive."""

        if "_par_index" == name:
            if self.__dict__.get('_par_index') is None:
                self.__dict__['_par_index'] = {}
            return self.__dict__['_par_index']

        lowered_name = name.lower()

        def warn(oname, nname):
            wmsg = 'Parameter name {} is deprecated'.format(oname) + \
                ' for model {}, '.format(type(self).__name__) + \
                'use {} instead'.format(nname)
            warnings.warn(wmsg, DeprecationWarning)

        parameter = self._par_index.get(lowered_name)

        if parameter is not None:
            if lowered_name in parameter.aliases:
                warn(lowered_name, parameter.name)
            return parameter

        NoNewAttributesAfterInit.__getattribute__(self, name)

    def __setattr__(self, name, val):
        lname = name.lower()
        par = getattr(self, lname, None)
        if (par is not None) and isinstance(par, Parameter):
            # When setting an attribute that is a Parameter, set the parameter's
            # value instead.
            par.val = val
        else:
            NoNewAttributesAfterInit.__setattr__(self, name, val)
            if isinstance(val, Parameter):
                vname = val.name.lower()

                # Check the parameter names match - as this is a
                # 'development' error then just make this an assert.
                # Ideally it should be exact but support lower case
                # comparison only
                assert lname == vname, (name, val.name, self.name)

                # Update parameter index
                self._par_index[vname] = val
                if val.aliases:
                    # Update index of aliases, if necessary
                    for alias in val.aliases:
                        self._par_index[alias] = val

    def startup(self, cache=False):
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

    def calc(self, p, *args, **kwargs):
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
        """
        raise NotImplementedError

    def teardown(self):
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

    def __call__(self, *args, **kwargs):
        # A bit of trickery, to make model creation
        # in IPython happen without raising errors, when
        # model is made automatically callable
        if (len(args) == 0 and len(kwargs) == 0):
            return self
        return self.calc([p.val for p in self.pars], *args, **kwargs)

    def _get_thawed_pars(self):
        return [p.val for p in self.pars if not p.frozen]

    def _set_thawed_pars(self, vals):
        tpars = [p for p in self.pars if not p.frozen]

        ngot = len(vals)
        nneed = len(tpars)
        if ngot != nneed:
            raise ModelErr('numthawed', nneed, ngot)

        for p, v in zip(tpars, vals):
            v = SherpaFloat(v)
            if v < p.hard_min:
                p.val = p.min
                warning(('value of parameter %s is below minimum; ' +
                         'setting to minimum') % p.fullname)
            elif v > p.hard_max:
                p.val = p.max
                warning(('value of parameter %s is above maximum; ' +
                         'setting to maximum') % p.fullname)
            else:
                p._val = v

    thawedpars = property(_get_thawed_pars, _set_thawed_pars,
                          doc='The thawed parameters of the model.\n\n' +
                          'Get or set the thawed parameters of the model as a list of\n' +
                          'numbers. If there are no thawed parameters then [] is used.\n' +
                          'The ordering matches that of the pars attribute.\n\n' +
                          'See Also\n' +
                          '--------\n' +
                          'thawedparmaxes, thawedparmins\n')

    def _get_thawed_par_mins(self):
        return [p.min for p in self.pars if not p.frozen]

    def _set_thawed_pars_mins(self, vals):
        tpars = [p for p in self.pars if not p.frozen]

        ngot = len(vals)
        nneed = len(tpars)
        if ngot != nneed:
            raise ModelErr('numthawed', nneed, ngot)

        for p, v in zip(tpars, vals):
            v = SherpaFloat(v)
            if v < p.hard_min:
                p.min = p.hard_min
                warning(('value of parameter %s minimum is below ' +
                         'hard minimum; ' +
                         'setting to hard minimum') % p.fullname)
            elif v > p.hard_max:
                p.min = p.hard_max
                warning(('value of parameter %s minimum is above ' +
                         'hard maximum; ' +
                         'setting to hard maximum') % p.fullname)
            else:
                p._min = v

    thawedparmins = property(_get_thawed_par_mins, _set_thawed_pars_mins,
                             doc='The minimum limits of the thawed parameters.\n\n' +
                             'Get or set the minimum limits of the thawed parameters\n' +
                             'of the model as a list of numbers. If there are no\n' +
                             'thawed parameters then [] is used. The ordering matches\n' +
                             'that of the pars attribute.\n\n' +
                             'See Also\n' +
                             '--------\n' +
                             'thawedpars, thawedarhardmins, thawedparmaxes\n')

    def _get_thawed_par_maxes(self):
        return [p.max for p in self.pars if not p.frozen]

    def _set_thawed_pars_maxes(self, vals):
        tpars = [p for p in self.pars if not p.frozen]

        ngot = len(vals)
        nneed = len(tpars)
        if ngot != nneed:
            raise ModelErr('numthawed', nneed, ngot)

        for p, v in zip(tpars, vals):
            v = SherpaFloat(v)
            if v < p.hard_min:
                p.max = p.hard_min
                warning(('value of parameter %s maximum is below ' +
                         'hard minimum; ' +
                         'setting to hard minimum') % p.fullname)
            elif v > p.hard_max:
                p.max = p.hard_max
                warning(('value of parameter %s maximum is above ' +
                         'hard maximum; ' +
                         'setting to hard maximum') % p.fullname)
            else:
                p._max = v

    thawedparmaxes = property(_get_thawed_par_maxes, _set_thawed_pars_maxes,
                              doc='The maximum limits of the thawed parameters.\n\n' +
                              'Get or set the maximum limits of the thawed parameters\n' +
                              'of the model as a list of numbers. If there are no\n' +
                              'thawed parameters then [] is used. The ordering matches\n' +
                              'that of the pars attribute.\n\n' +
                              'See Also\n' +
                              '--------\n' +
                              'thawedpars, thawedarhardmaxes, thawedparmins\n')

    def _get_thawed_par_hardmins(self):
        return [p.hard_min for p in self.pars if not p.frozen]

    thawedparhardmins = property(_get_thawed_par_hardmins,
                                 doc='The hard minimum values for the thawed parameters.\n\n' +
                                 'The minimum and maximum range of the parameters can be\n' +
                                 'changed with thawedparmins and thawedparmaxes but only\n' +
                                 'within the range given by thawedparhardmins\n' +
                                 'to thawparhardmaxes.\n\n' +
                                 'See Also\n' +
                                 '--------\n' +
                                 'thawedparhardmaxes, thawedparmins\n')

    def _get_thawed_par_hardmaxes(self):
        return [p.hard_max for p in self.pars if not p.frozen]

    thawedparhardmaxes = property(_get_thawed_par_hardmaxes,
                                  doc='The hard maximum values for the thawed parameters.\n\n' +
                                 'The minimum and maximum range of the parameters can be\n' +
                                 'changed with thawedparmins and thawedparmaxes but only\n' +
                                 'within the range given by thawedparhardmins\n' +
                                 'to thawparhardmaxes.\n\n' +
                                  'See Also\n' +
                                  '--------\n' +
                                  'thawedparhardmins, thawedparmaxes\n')

    def reset(self):
        """Reset the parameter values.

        Restores each parameter to the last value it was set to.
        This allows the parameters to be easily reset after a
        fit.
        """

        for p in self.pars:
            p.reset()


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

       >>> l1 = Gauss1D('l1')
       >>> l2 = Gauss1D('l2')
       >>> b = Polynom1D('b')
       >>> mdl = l1 + (0.5 * l2) + b
       >>> mdl
       <BinaryOpModel model instance '((l1 + (0.5 * l2)) + polynom1d)'>
       >>> for cpt in mdl:
       ...     print(type(c))
       ...
       <class 'BinaryOpModel'>
       <class 'sherpa.models.basic.Gauss1D'>
       <class 'BinaryOpModel'>
       <class 'ArithmeticConstantModel'>
       <class 'sherpa.models.basic.Gauss1D'>
       <class 'sherpa.models.basic.Polynom1D'>

    """

    def __init__(self, name, parts):
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
                    raise ModelErr('Models do not match: ' +
                                   '{}D ({}) and '.format(self.ndim, model_with_dim.name) +
                                   '{}D ({})'.format(ndim, part.name))

            for p in part.pars:
                if p in allpars:
                    # If we already have a reference to this parameter, store
                    # a hidden, linked proxy instead
                    pnew = Parameter(p.modelname, p.name, 0.0, hidden=True)
                    pnew.link = p
                    p = pnew
                allpars.append(p)

        Model.__init__(self, name, allpars)

        for part in self.parts:
            try:
                self.is_discrete = self.is_discrete or part.is_discrete
            except AttributeError:
                warning("Could not determine whether the model is discrete.\n" +
                        "This probably means that you have restored a session saved with a previous version of Sherpa.\n" +
                        "Falling back to assuming that the model is continuous.\n")
                self.is_discrete = False

    def __iter__(self):
        return iter(self._get_parts())

    def _get_parts(self):
        parts = []

        for p in self.parts:
            # A CompositeModel should not hold a reference to itself
            assert (p is not self), (("'%s' object holds a reference to " +
                                      "itself") % type(self).__name__)

            # Including itself seems a bit strange if it's a CompositeModel
            # but is used by sherpa.astro.instrument.has_pha_instance (and
            # possibly elsewhere).
            #
            parts.append(p)
            if isinstance(p, CompositeModel):
                parts.extend(p._get_parts())

        # FIXME: do we want to remove duplicate components from parts?

        return parts

    def startup(self, cache=False):
        pass

    def teardown(self):
        pass

    def cache_clear(self):
        """Clear the cache for each component."""
        for p in self.parts:
            try:
                p.cache_clear()
            except AttributeError:
                pass

    def cache_status(self):
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

    def create(self, *args):
        """Create a copy of the object.

        Note that it will re-use model components rather than
        re-creating them.
        """

        raise NotImplementedError()

    def expand(self):
        """Expand the model expression to separate out "additive" terms.

        Note that it will re-use model components rather than
        re-creating them.
        """

        raise NotImplementedError()


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

    >>> m1 = Polynom1D('m1')
    >>> m2 = Gauss1D('g1')
    >>> mall = SimulFitModel('comp', (m1, m1 + m2))

    If dall is a DataSimulFit object then the model components
    can be evaluated for the composite object using:

    >>> ymdl = dall.eval_model_to_fit(mall)

    """

    def __iter__(self):
        return iter(self.parts)

    def startup(self, cache=False):
        for part in self:
            part.startup(cache)
        CompositeModel.startup(self, cache)

    def teardown(self):
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

    Attributes
    ----------
    val : number

    """

    def __init__(self, val, name=None):
        val = SherpaFloat(val)
        if name is None:
            if numpy.isscalar(val):
                name = str(val)
            else:
                name = '{}[{}]'.format(val.dtype.name,
                                       ','.join([str(s) for s in val.shape]))

        self.name = name
        self.val = val

        # val has to be a scalar or 1D array, even if used with a 2D
        # model, due to the way model evaluation works, so as we
        # can't easily define the dimensionality of this model, we
        # remove any dimensionality checking for this class.
        #
        self.ndim = None

        Model.__init__(self, self.name)

    def _get_val(self):
        return self._val

    def _set_val(self, val):
        val = SherpaFloat(val)
        if val.ndim > 1:
            raise ModelErr('The constant must be a scalar or 1D, not 2D')

        self._val = val

    val = property(_get_val, _set_val,
                   doc='The constant value (scalar or 1D).')

    def startup(self, cache=False):
        pass

    def calc(self, p, *args, **kwargs):
        return self.val

    def teardown(self):
        pass


def _make_unop(op, opstr):
    def func(self):
        return UnaryOpModel(self, op, opstr)
    return func


def _make_binop(op, opstr):
    def func(self, rhs):
        return BinaryOpModel(self, rhs, op, opstr)

    def rfunc(self, lhs):
        return BinaryOpModel(lhs, self, op, opstr)

    return (func, rfunc)


class ArithmeticModel(Model):
    """Support combining model expressions and caching results."""

    cache = 5
    """The maximum size of the cache."""

    def __init__(self, name, pars=()):
        self.integrate = True

        # Model caching ability
        self.cache = 5  # repeat the class definition
        self._use_caching = True  # FIXME: reduce number of variables?
        self.cache_clear()
        Model.__init__(self, name, pars)

    def cache_clear(self):
        """Clear the cache."""
        # It is not obvious what to set the queue length to
        self._queue = ['']
        self._cache = {}
        self._cache_ctr = {'hits': 0, 'misses': 0, 'check': 0}

    def cache_status(self):
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
    __neg__ = _make_unop(numpy.negative, '-')
    __abs__ = _make_unop(numpy.absolute, 'abs')

    # Binary operations
    __add__, __radd__ = _make_binop(numpy.add, '+')
    __sub__, __rsub__ = _make_binop(numpy.subtract, '-')
    __mul__, __rmul__ = _make_binop(numpy.multiply, '*')
    __div__, __rdiv__ = _make_binop(numpy.divide, '/')
    __floordiv__, __rfloordiv__ = _make_binop(numpy.floor_divide, '//')
    __truediv__, __rtruediv__ = _make_binop(numpy.true_divide, '/')
    __mod__, __rmod__ = _make_binop(numpy.remainder, '%')
    __pow__, __rpow__ = _make_binop(numpy.power, '**')

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

    def startup(self, cache=False):
        self.cache_clear()
        self._use_caching = cache
        if int(self.cache) > 0:
            self._queue = [''] * int(self.cache)
            frozen = numpy.array([par.frozen for par in self.pars], dtype=bool)
            if len(frozen) > 0 and frozen.all():
                self._use_caching = cache

    def teardown(self):
        self._use_caching = False

    def apply(self, outer, *otherargs, **otherkwargs):
        return NestedModel(outer, self, *otherargs, **otherkwargs)


class RegriddableModel(ArithmeticModel):
    def regrid(self, *args, **kwargs):
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
        >>> mybox = Box1D()
        >>> request_space = np.arange(1, 10, 0.1)
        >>> regrid_model = mybox.regrid(request_space, interp=linear_interp)
        """
        valid_keys = ('interp',)
        for key in kwargs.keys():
            if key not in valid_keys:
                raise TypeError("unknown keyword argument: '%s'" % key)
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

    >>> m1 = Gauss1d()
    >>> m2 = UnaryOpModel(m1, numpy.negative, '-')

    """

    @staticmethod
    def wrapobj(obj):
        return _wrapobj(obj, ArithmeticConstantModel)

    def __init__(self, arg, op, opstr):
        self.arg = self.wrapobj(arg)
        self.op = op
        self.opstr = opstr
        CompositeModel.__init__(self, f'{opstr}({self.arg.name})',
                                (self.arg,))

    def calc(self, p, *args, **kwargs):
        return self.op(self.arg.calc(p, *args, **kwargs))

    def create(self, *args):
        """Create a version of unary-op model applied to model.

        Parameters
        ----------
        model : Model instance

        Returns
        -------
        applied : UnaryOpModel instance
        """

        if len(args) != 1:
            raise TypeError("create() missing 1 required positional argument: model")

        return self.__class__(args[0], self.op, self.opstr)

    def expand(self):
        """Expand the model expression to separate out "additive" terms.

        Returns
        -------
        expanded : UnaryOpModel instance
        """

        try:
            expanded = self.arg.expand()
        except AttributeError:
            return self

        return self.create(expanded)

    def separate(self):
        """Separate out the 'additive' terms in the model expression.

        Returns
        -------
        models : list of Model instances
            The additive terms of the model. Summing the elements will
            return the equivalent model expression to that stored in the
            object.

        Notes
        -----
        We separate out the components and then apply the operator to
        each term.

        """

        try:
            out = self.arg.separate()
        except AttributeError:
            return [self]

        return [self.create(m) for m in out]


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

    >>> m1 = Gauss1d()
    >>> m2 = Polynom1D()
    >>> m = BinaryOpModel(m1, m2, numpy.add, '+')

    """

    @staticmethod
    def wrapobj(obj):
        return _wrapobj(obj, ArithmeticConstantModel)

    def __init__(self, lhs, rhs, op, opstr):
        self.lhs = self.wrapobj(lhs)
        self.rhs = self.wrapobj(rhs)
        self.op = op
        self.opstr = opstr

        CompositeModel.__init__(self,
                                f'({self.lhs.name} {opstr} {self.rhs.name})',
                                (self.lhs, self.rhs))

    def regrid(self, *args, **kwargs):
        for part in self.parts:
            # ArithmeticConstantModel does not support regrid by design
            if not hasattr(part, 'regrid'):
                continue
            # The full model expression must be used
            return part.__class__.regrid(self, *args, **kwargs)
        raise ModelErr('Neither component supports regrid method')

    def startup(self, cache=False):
        self.lhs.startup(cache)
        self.rhs.startup(cache)
        CompositeModel.startup(self, cache)

    def teardown(self):
        self.lhs.teardown()
        self.rhs.teardown()
        CompositeModel.teardown(self)

    def calc(self, p, *args, **kwargs):
        nlhs = len(self.lhs.pars)
        lhs = self.lhs.calc(p[:nlhs], *args, **kwargs)
        rhs = self.rhs.calc(p[nlhs:], *args, **kwargs)
        try:
            val = self.op(lhs, rhs)
        except ValueError:
            raise ValueError("shape mismatch between '%s: %i' and '%s: %i'" %
                             (type(self.lhs).__name__, len(lhs),
                              type(self.rhs).__name__, len(rhs)))
        return val

    def create(self, *args):
        """Create a version of binary-op model applied to the models.

        Parameters
        ----------
        lhs, rhs : Model instance

        Returns
        -------
        applied : BinaryOpModel instance
        """

        if len(args) != 2:
            raise TypeError("create() missing 2 required positional arguments: lhs, rhs")

        return self.__class__(args[0], args[1], self.op, self.opstr)

    def expand(self):
        """Expand the model expression.

        Returns
        -------
        expanded : BinaryOpModel instance
        """

        if is_add(self):
            try:
                lhs = self.lhs.expand()
            except AttributeError:
                lhs = self.lhs

            try:
                rhs = self.rhs.expand()
            except AttributeError:
                rhs = self.rhs

            return self.create(lhs, rhs)

        if not is_mul(self):
            # This is technically possible but requires using BinaryOpModel
            # in a way that Sherpa doesn't do.
            #
            return self

        # Try to avoid isinstance checks and be more python-like.
        #
        try:
            is_lhs_binop = len(self.lhs.parts) == 2
            lhs_add = is_add(self.lhs)
            lhs_mul = is_mul(self.lhs)
        except AttributeError:
            is_lhs_binop = False
            lhs_add = False
            lhs_mul = False

        try:
            is_rhs_binop = len(self.rhs.parts) == 2
            rhs_add = is_add(self.rhs)
            rhs_mul = is_mul(self.rhs)
        except AttributeError:
            is_rhs_binop = False
            rhs_add = False
            rhs_mul = False

        def mcomb(a, b):
            """Multiply/Divide/... the two components

            Note we expand the components and then expand the
            result - ie expand(expand(a) * expand(b))
            for a multiplicative model.
            """
            lhs = expand(a)
            rhs = expand(b)
            mul = self.create(lhs, rhs)
            return expand(mul)

        # The idea is that we expand the components and then, if
        # we create a (a * b) term, we expand that. We don't expand
        # addition terms.
        #
        # In the following I use + and * but they are meant to be
        # generic, since it could be - or /. We therefore use
        # mcomb and acomb to combine the components.
        #

        # (l1 + l2) * notbinop =
        #    expand (expand l1 * notbinop) +
        #    expand (expand l2 * notbinop)
        #
        if lhs_add and not is_rhs_binop:
            a = mcomb(self.lhs.lhs, self.rhs)
            b = mcomb(self.lhs.rhs, self.rhs)
            return self.lhs.create(a, b)

        # notbinop * (r1 + r2) =
        #   expand (notbinop * expand r1) +
        #   expand (notbinop * expand r2)
        #
        if not is_lhs_binop and rhs_add:
            a = mcomb(self.lhs, self.rhs.lhs)
            b = mcomb(self.lhs, self.rhs.rhs)
            return self.rhs.create(a, b)

        # (l1 + l2) * (r1 + r2) =
        #   expand (expand l1 * expand r1) +
        #   expand (expand l1 * expand r2) +
        #   expand (expand l2 * expand r1) +
        #   expand (expand l2 * expand r2) +
        #
        # Need to handle potentially different addition terms for the
        # two sides, so
        #
        # (l1 `f` l2) * (r1 `g` r2) =
        #   (l1 * r1) `f` (l2 * r1) `g` (l1 * r2) `g` (l2 * r2)
        #
        if lhs_add and rhs_add:
            # This duplicates the expanson of the components
            # but it is unlikely to be expensive and I think
            # the simplification in the code is worth it.
            #
            l1 = self.lhs.lhs
            l2 = self.lhs.rhs
            r1 = self.rhs.lhs
            r2 = self.rhs.rhs
            t1 = self.lhs.create(mcomb(l1, r1), mcomb(l2, r1))
            t2 = self.rhs.create(t1, mcomb(l1, r2))
            t3 = self.rhs.create(t2, mcomb(l2, r2))
            return t3

        # Can this be handled by (l1 + l2) * notbinop?
        #
        # (l1 + l2) * (r1 * r2) =
        #   expand (expand l1 * expand (expand r1 * expand r2)) +
        #   expand (expand l2 * expand (expand r1 * expand r2)) +
        #
        if lhs_add and rhs_mul:
            # er gets deconstructed and reconstructed in mcomb
            er = mcomb(self.rhs.lhs, self.rhs.rhs)
            return self.lhs.create(mcomb(self.lhs.lhs, er),
                                   mcomb(self.lhs.rhs, er))

        # Can this be handled by notbinop * (r1 + r2)?
        #
        # (l1 * l2) * (r1 + r2) =
        #   expand (expand (expand l1 * expand l2) * expand r1) +
        #   expand (expand (expand l1 * expand l2) * expand r2) +
        #
        if lhs_mul and rhs_add:
            # el gets deconstructed and reconstructed in mcomb
            el = mcomb(self.lhs.lhs, self.lhs.rhs)
            return self.rhs.create(mcomb(el, self.rhs.lhs),
                                   mcomb(el, self.rhs.rhs))

        # Just expand the children. We don't use mcomb() as that
        # calls expand on the result, and we don't want that here,
        # to avoid infinite loops.
        #
        a = self.lhs.expand()
        b = self.rhs.expand()
        return self.create(a, b)

    def separate(self):
        """Separate out the 'additive' terms in the model expression.

        Returns
        -------
        models : list of Model instances
            The additive terms of the model. Summing the elements will
            return the equivalent model expression to that stored in the
            object.

        Notes
        -----
        The expansion is only at the top-most level, so an expression
        like 'a * (b - c)' will return only a list with a single
        component, whereas 'a * b - a * c' will return two components
        ('a * b' and '-(a * c)').
        """

        if self.op not in [numpy.add, numpy.subtract]:
            return [self]

        out = []
        out.extend(separate(self.lhs))

        # Do we have to convert 'a - b' to 'a + (-b)'?
        #
        rhs = self.rhs
        if self.op == numpy.subtract:
            rhs = -rhs

        out.extend(separate(rhs))
        return out


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
            s += ':%s' % filter.step

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

    def __init__(self, func):
        if isinstance(func, Model):
            raise ModelErr('badinstance', type(self).__name__)
        if not callable(func):
            raise ModelErr('noncall', type(self).__name__, type(func).__name__)
        self.func = func
        Model.__init__(self, func.__name__)

    def calc(self, p, *args, **kwargs):
        return self.func(*args, **kwargs)

    def startup(self, cache=False):
        pass

    def teardown(self):
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
    def wrapobj(obj):
        return _wrapobj(obj, ArithmeticFunctionModel)

    def __init__(self, outer, inner, *otherargs, **otherkwargs):
        self.outer = self.wrapobj(outer)
        self.inner = self.wrapobj(inner)
        self.otherargs = otherargs
        self.otherkwargs = otherkwargs
        CompositeModel.__init__(self,
                                f'{self.outer.name}({self.inner.name})',
                                (self.outer, self.inner))

    def startup(self, cache=False):
        self.inner.startup(cache)
        self.outer.startup(cache)
        CompositeModel.startup(self, cache)

    def teardown(self):
        self.inner.teardown()
        self.outer.teardown()
        CompositeModel.teardown(self)

    def calc(self, p, *args, **kwargs):
        nouter = len(self.outer.pars)
        return self.outer.calc(p[:nouter],
                               self.inner.calc(p[nouter:], *args, **kwargs),
                               *self.otherargs, **self.otherkwargs)


# TODO: can we remove this as it is now unused?
#
class MultigridSumModel(CompositeModel, ArithmeticModel):

    def __init__(self, models):
        self.models = tuple(models)
        name = '%s(%s)' % (type(self).__name__,
                           ','.join([m.name for m in models]))
        CompositeModel.__init__(self, name, self.models)

    def calc(self, p, arglist):
        vals = []
        for model, args in zip(self.models, arglist):
            # FIXME: we're not using p here (and therefore assuming that the
            # parameter values have already been updated to match the contents
            # of p)
            vals.append(model(*args))
        return sum(vals)


class RegridWrappedModel(CompositeModel, ArithmeticModel):

    def __init__(self, model, wrapper):
        self.model = self.wrapobj(model)
        self.wrapper = wrapper

        if hasattr(model, 'integrate'):
            self.wrapper.integrate = model.integrate

        CompositeModel.__init__(self,
                                f"{self.wrapper.name}({self.model.name})",
                                (self.model, ))

    def calc(self, p, *args, **kwargs):
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
    def wrapobj(obj):
        # TODO: choice of ArithmeticConstandModel or
        #       ArithmeticFunctionModel?
        return _wrapobj(obj, ArithmeticFunctionModel)


def _wrapobj(obj, wrapper):
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
    if isinstance(obj, (ArithmeticModel, ArithmeticConstantModel, ArithmeticFunctionModel)):
        return obj

    return wrapper(obj)


# Notebook representation
#
def modelcomponents_to_list(model):
    if hasattr(model, 'parts'):
        modellist = []
        for p in model.parts:
            modellist.extend(modelcomponents_to_list(p))
        return modellist
    else:
        return [model]


def html_model(mdl):
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
            else:
                this_comp_nrows +=1
        if this_comp_nrows > 0:
            complist.append(comp)
            nrows.append(this_comp_nrows)

    out = '<table class="model">'
    out += '<caption>Expression: {}</caption>'.format(formatting.clean_bracket(mdl.name))
    out += '<thead><tr>'
    cols = ['Component', 'Parameter', 'Thawed', 'Value',
            'Min', 'Max', 'Units']
    for col in cols:
        out += '<th>{}</th>'.format(col)

    out += '</tr></thead><tbody>'

    for mcount, (n, comp) in enumerate(zip(nrows, complist)):
        for i, par in enumerate(comp.pars):
            style = '' if ((i > 0) or (mcount == 0)) else ' class="block"'

            if par.hidden:
                continue

            def addtd(val):
                "Use the parameter to convert to HTML"
                return '<td>{}</td>'.format(par._val_to_html(val))

            out += '<tr{}>'.format(style)
            if i ==0 :
                cls = "model-" + ("even" if mcount % 2 == 1 else "odd")
                out += '<th class="{}" scope="rowgroup" '.format(cls)
                out += 'rowspan={}>{}</th>'.format(n, comp.name)

            out += '<td>{}</td>'.format(par.name)

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
                out += '<td colspan=2>&#8656; {}</td>'.format(formatting.clean_bracket(par.link.fullname))
            else:
                out += addtd(par.min)
                out += addtd(par.max)

            out += '<td>{}</td>'.format(par._units_to_html())
            out += '</tr>'

    out += '</tbody></table>'

    ls = ['<details open><summary>Model</summary>' + out + '</details>']
    return formatting.html_from_sections(mdl, ls)


# Model expansion
#
# Operators include
#    numpy.add, numpy.subtract
#    numpy.multiply, numpy.divide, numpy.floor_divide, numpy.true_divide
#
def is_add(mdl):
    """Does the model represent addition or subtraction?

    Paramerters
    -----------
    mdl : BinaryOpModel instance

    Returns
    -------
    flag : bool
    """
    return mdl.op in [numpy.add, numpy.subtract]


def is_mul(mdl):
    """Does the model represent multiplication or division?

    Paramerters
    -----------
    mdl : BinaryOpModel instance

    Returns
    -------
    flag : bool
    """
    return mdl.op in [numpy.multiply, numpy.divide,
                      numpy.floor_divide, numpy.true_divide]


# Experiment in "expanding out" model expressions which is useful
# when we want to try and combine all the "emission" components
# (that is, components combined with * or /).
#
def expand(model):
    """Expand the model expression.

    This is not guaranteed to support all CompositeModel sub-classes.

    Parameters
    ----------
    model : sherpa.models.model.Model instance
        The model expression to expand.

    Returns
    -------
    expanded : sherpa.models.model.Model instance
        The expanded expression. It can be the same as the input, and
        will reference the same model components.

    See Also
    --------
    separate

    Examples
    --------

    >>> box = sherpa.models.basic.Box1D()
    >>> line = sherpa.models.basic.Gauss1D()
    >>> bgnd = sherpa.models.basic.Polynom1D()
    >>> orig1 = box * (line + bgnd)
    >>> orig2 = box * (bgnd - line)
    >>> expand(box * line).name
    '(box1d * gauss1d)'
    >>> expand(line + bgnd).name
    '(gauss1d + polynom1d)'
    >>> expand(orig1).name
    '((box1d * gauss1d) + (box1d * polynom1d))'
    >>> expand(orig2).name
    '((box1d * polynom1d) - (box1d * gauss1d))'

    """

    try:
        return model.expand()
    except AttributeError:
        return model


def separate(mdl):
    """Separate out the additive terms of the model.

    Parameters
    ----------
    model : sherpa.models.model.Model instance
        The model expression. It is expected to be the output of
        `expand` but it does not have to be.

    Returns
    -------
    terms : list of sherpa.models.model.Model instances
        The individual additive terms. They are added together
        to create the full model expression.

    See Also
    --------
    expand

    Examples
    --------

    >>> box = sherpa.models.basic.Box1D()
    >>> line = sherpa.models.basic.Gauss1D()
    >>> bgnd = sherpa.models.basic.Polynom1D()
    >>> orig1 = box * (line + bgnd)
    >>> orig2 = box * (bgnd - line)
    >>> separate(box * line)
    [<BinaryOpModel model instance '(box1d * gauss1d)'>]
    >>> separate(line + bgnd)
    [<BinaryOpModel model instance '(gauss1d + polynom1d)'>]
    >>> separate(orig1)
    [<BinaryOpModel model instance '(box1d * (gauss1d + polynom1d))'>]
    >>> separate(expand(orig1))
    [<BinaryOpModel model instance '(box1d * gauss1d)'>,
     <BinaryOpModel model instance '(box1d * polynom1d)'>]
    >>> separate(expand(orig2))
    [<BinaryOpModel model instance '(box1d * polynom1d)'>,
     <UnaryOpModel model instance '-((box1d * gauss1d))'>]
    >>> [c.name for c in separate(expand(orig1))]
    ['(box1d * gauss1d)', '(box1d * polynom1d)']
    >>> [c.name for c in separate(expand(orig2))]
    ['(box1d * polynom1d)', '-((box1d * gauss1d))']

    """

    # Since we can't guarantee all model components will contain
    # a separate method. In fact only BinaryOpModels do at present,
    # which is why this is left as a top-level function rather
    # than an interface-like call.
    #
    try:
        return mdl.separate()
    except AttributeError:
        return [mdl]
