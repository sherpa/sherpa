from __future__ import absolute_import
#
#  Copyright (C) 2010, 2016, 2017, 2018, 2019, 2020
#      Smithsonian Astrophysical Observatory
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


import logging
import numpy
import hashlib
import warnings

from sherpa.models.regrid import EvaluationSpace1D, ModelDomainRegridder1D, EvaluationSpace2D, ModelDomainRegridder2D
from sherpa.utils import SherpaFloat, NoNewAttributesAfterInit
from sherpa.utils.err import ModelErr

from .parameter import Parameter

warning = logging.getLogger(__name__).warning


__all__ = ('Model', 'CompositeModel', 'SimulFitModel',
           'ArithmeticConstantModel', 'ArithmeticModel', 'RegriddableModel1D', 'RegriddableModel2D',
           'UnaryOpModel', 'BinaryOpModel', 'FilterModel', 'modelCacher1d',
           'ArithmeticFunctionModel', 'NestedModel', 'MultigridSumModel')


def boolean_to_byte(boolean_value):
    bmap = {True: b'1', False: b'0'}
    return bmap.get(boolean_value, b'0')


def modelCacher1d(func):

    def cache_model(cls, pars, xlo, *args, **kwargs):
        use_caching = cls._use_caching
        cache = cls._cache
        queue = cls._queue

        digest = ''
        if use_caching:

            data = [numpy.array(pars).tostring(), boolean_to_byte(kwargs.get('integrate', False)),
                    numpy.asarray(xlo).tostring()]
            if args:
                data.append(numpy.asarray(args[0]).tostring())

            token = b''.join(data)
            digest = hashlib.sha256(token).digest()
            if digest in cache:
                return cache[digest].copy()

        vals = func(cls, pars, xlo, *args, **kwargs)

        if use_caching:
            # remove first item in queue and remove from cache
            key = queue.pop(0)
            cache.pop(key, None)

            # append newest model values to queue
            queue.append(digest)
            cache[digest] = vals.copy()

        return vals

    cache_model.__name__ = func.__name__
    cache_model.__doc__ = func.__doc__
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

    Attributes
    ----------
    name : str
        The name given to the instance.
    pars : tuple of sherpa.parameter.Parameter objects
        The parameters of the model instance.

    Notes
    -----
    Parameters can be accessed via the ``pars`` attribute, but it is
    expected that they will generally be accessed directly, as the
    class provides case-insensitive access to the parameter names
    as object attributes. That is, if the model contains parameters
    called ``breakFreq`` and ``norm``, and the instance is stored in
    the variable ``mdl``, then the following can be used to access
    the parameters::

        print("Break frequency = {}".format(mdl.breakfreq))

        mdl.norm = 1.2e-3

    """

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
        par = getattr(self, name.lower(), None)
        if (par is not None) and isinstance(par, Parameter):
            # When setting an attribute that is a Parameter, set the parameter's
            # value instead.
            par.val = val
        else:
            NoNewAttributesAfterInit.__setattr__(self, name, val)
            if isinstance(val, Parameter):
                # Update parameter index
                self._par_index[val.name.lower()] = val
                if val.aliases:
                    # Update index of aliases, if necessary
                    for alias in val.aliases:
                        self._par_index[alias] = val

    def startup(self):
        """Called before a model may be evaluated multiple times.

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
        setup
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

    thawedpars = property(_get_thawed_pars, _set_thawed_pars)

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

    thawedparmins = property(_get_thawed_par_mins, _set_thawed_pars_mins)

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

    thawedparmaxes = property(_get_thawed_par_maxes, _set_thawed_pars_maxes)

    def _get_thawed_par_hardmins(self):
        return [p.hard_min for p in self.pars if not p.frozen]

    thawedparhardmins = property(_get_thawed_par_hardmins)

    def _get_thawed_par_hardmaxes(self):
        return [p.hard_max for p in self.pars if not p.frozen]

    thawedparhardmaxes = property(_get_thawed_par_hardmaxes)

    def reset(self):
        for p in self.pars:
            p.reset()


class CompositeModel(Model):

    def __init__(self, name, parts):
        self.parts = tuple(parts)
        allpars = []
        for part in self.parts:
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
            except:
                warning("Could not determine whether the model is discrete.\n"+
                        "This probably means that you have restored a session saved with a previous version of Sherpa.\n"+
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

            parts.append(p)
            if isinstance(p, CompositeModel):
                parts.extend(p._get_parts())

        # FIXME: do we want to remove duplicate components from parts?

        return parts

    def startup(self, cache):
        pass

    def teardown(self):
        pass


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

    def startup(self, cache):
        for part in self:
            part.startup(cache)
        CompositeModel.startup(self, cache)

    def teardown(self):
        for part in self:
            part.teardown()
        CompositeModel.teardown(self)


class ArithmeticConstantModel(Model):

    def __init__(self, val, name=None):
        if name is None:
            name = str(val)
        self.name = name
        self.val = SherpaFloat(val)
        Model.__init__(self, self.name)

    def startup(self, cache):
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

    def __init__(self, name, pars=()):
        self.integrate = True

        # Model caching ability
        # queue memory of maximum size
        self.cache = 5
        self._use_caching = True  # FIXME: reduce number of variables?
        self._queue = ['']
        self._cache = {}
        Model.__init__(self, name, pars)

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

        if 'cache' not in state:
            self.__dict__['cache'] = 5

    def __getitem__(self, filter):
        return FilterModel(self, filter)

    def startup(self, cache):
        self._queue = ['']
        self._cache = {}
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


class RegriddableModel1D(ArithmeticModel):
    def regrid(self, *arrays, **kwargs):
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
        valid_keys = ('interp')
        for key in kwargs.keys():
            if key not in valid_keys:
                raise TypeError("unknown keyword argument: '%s'" % key)
        eval_space = EvaluationSpace1D(*arrays)
        regridder = ModelDomainRegridder1D(eval_space, **kwargs)
        regridder._make_and_validate_grid(arrays)
        return regridder.apply_to(self)


class RegriddableModel2D(ArithmeticModel):
    def regrid(self, *arrays):
        eval_space = EvaluationSpace2D(*arrays)
        regridder = ModelDomainRegridder2D(eval_space)
        return regridder.apply_to(self)


class UnaryOpModel(CompositeModel, ArithmeticModel):

    def __init__(self, arg, op, opstr):
        self.arg = arg
        self.op = op
        CompositeModel.__init__(self, ('%s(%s)' % (opstr, self.arg.name)),
                                (self.arg,))

    def calc(self, p, *args, **kwargs):
        return self.op(self.arg.calc(p, *args, **kwargs))


class BinaryOpModel(CompositeModel, ArithmeticModel):

    @staticmethod
    def wrapobj(obj):
        if isinstance(obj, ArithmeticModel):
            return obj
        return ArithmeticConstantModel(obj)

    def __init__(self, lhs, rhs, op, opstr):
        self.lhs = self.wrapobj(lhs)
        self.rhs = self.wrapobj(rhs)
        self.op = op
        CompositeModel.__init__(self,
                                ('(%s %s %s)' %
                                 (self.lhs.name, opstr, self.rhs.name)),
                                (self.lhs, self.rhs))

    def startup(self, cache):
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


class FilterModel(CompositeModel, ArithmeticModel):

    def __init__(self, model, filter):
        self.model = model
        self.filter = filter

        if isinstance(filter, tuple):
            filter_str = ','.join([self._make_filter_str(f) for f in filter])
        else:
            filter_str = self._make_filter_str(filter)

        CompositeModel.__init__(self,
                                ('(%s)[%s]' % (self.model.name, filter_str)),
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

    def __init__(self, func):
        if isinstance(func, Model):
            raise ModelErr('badinstance', type(self).__name__)
        if not callable(func):
            raise ModelErr('noncall', type(self).__name__, type(func).__name__)
        self.func = func
        Model.__init__(self, func.__name__)

    def calc(self, p, *args, **kwargs):
        return self.func(*args, **kwargs)

    def startup(self):
        pass

    def teardown(self):
        pass


class NestedModel(CompositeModel, ArithmeticModel):

    @staticmethod
    def wrapobj(obj):
        if isinstance(obj, ArithmeticModel):
            return obj
        return ArithmeticFunctionModel(obj)

    def __init__(self, outer, inner, *otherargs, **otherkwargs):
        self.outer = self.wrapobj(outer)
        self.inner = self.wrapobj(inner)
        self.otherargs = otherargs
        self.otherkwargs = otherkwargs
        CompositeModel.__init__(self,
                                ('%s(%s)' %
                                 (self.outer.name, self.inner.name)),
                                (self.outer, self.inner))

    def startup(self, cache):
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
                                "{}({})".format(self.wrapper.name,
                                                self.model.name),
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
        if isinstance(obj, ArithmeticModel):
            return obj
        else:
            return ArithmeticFunctionModel(obj)
