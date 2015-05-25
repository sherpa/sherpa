# 
#  Copyright (C) 2010  Smithsonian Astrophysical Observatory
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


from itertools import izip
import logging
import numpy
import hashlib
from sherpa.utils import SherpaFloat, NoNewAttributesAfterInit
from sherpa.utils.err import ModelErr

from parameter import Parameter

warning = logging.getLogger(__name__).warning


__all__ = ('Model', 'CompositeModel', 'SimulFitModel',
           'ArithmeticConstantModel', 'ArithmeticModel',
           'UnaryOpModel', 'BinaryOpModel', 'FilterModel', 'modelCacher1d',
           'ArithmeticFunctionModel', 'NestedModel', 'MultigridSumModel')

def modelCacher1d(func):

    def cache_model(cls, pars, xlo, *args, **kwargs):
        use_caching = cls._use_caching
        cache = cls._cache
        queue = cls._queue

        digest = ''
        if use_caching:

            data = [ numpy.array(pars).tostring(), str(kwargs.get('integrate',0)), 
                     numpy.asarray(xlo).tostring() ]
            if args:
                data.append( numpy.asarray(args[0]).tostring() )

            token = ''.join(data)
            digest = hashlib.sha256(token).digest()
            if digest in cache:
                return cache[digest]

        vals = func(cls, pars, xlo, *args, **kwargs)

        if use_caching:
            # remove first item in queue and remove from cache
            key = queue.pop(0)
            cache.pop(key, None)

            # append newest model values to queue
            queue.append(digest)
            cache[digest] = vals

        return vals

    cache_model.__name__ = func.__name__
    cache_model.__doc__ = func.__doc__
    return cache_model


class Model(NoNewAttributesAfterInit):

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
        hfmt = '\n   %-12s %-6s %12s %12s %12s %10s'
        s += hfmt % ('Param', 'Type', 'Value', 'Min', 'Max', 'Units')
        s += hfmt % ('-'*5, '-'*4, '-'*5, '-'*3, '-'*3, '-'*5)
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

    # Make parameter access case insensitive
    def __getattr__(self, name):
        par = None
        for key in self.__dict__.keys():
            if (type(key) == str):
                if (name.lower() == key.lower()):
                    par = self.__dict__.get(key)
                    break
        if (par is not None) and isinstance(par, Parameter):
            return par
        # this must be AttributeError for 'getattr' to work
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, name))

    def __setattr__(self, name, val):
        par = getattr(self, name.lower(), None)
        if (par is not None) and isinstance(par, Parameter):
            par.val = val
        else:
            NoNewAttributesAfterInit.__setattr__(self, name, val)

    def startup(self):
        raise NotImplementedError

    def calc(self, p, *args, **kwargs):
        raise NotImplementedError

    def teardown(self):
        raise NotImplementedError

    def guess(self, dep, *args, **kwargs):
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

        for p, v in izip(tpars, vals):
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

        for p, v in izip(tpars, vals):
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

        for p, v in izip(tpars, vals):
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

    def startup(self):
        #print 'Starting up %s...' % type(self).__name__
        pass

    def teardown(self):
        #print 'Tearing down %s...' % type(self).__name__
        pass


class SimulFitModel(CompositeModel):

    def __iter__(self):
        return iter(self.parts)


    def startup(self):
        #print 'Starting up %s...' % type(self).__name__
        for part in self:
            part.startup()
        CompositeModel.startup(self)


    def teardown(self):
        #print 'Tearing down %s...' % type(self).__name__
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

    def startup(self):
        #print 'Starting up %s...' % type(self).__name__
        pass

    def calc(self, p, *args, **kwargs):
        return self.val

    def teardown(self):
        #print 'Tearing down %s...' % type(self).__name__
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
        self.integrate=True

        # Model caching ability
        # queue memory of maximum size
        self.cache = 5
        self._use_caching = False  # FIXME: reduce number of variables?
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

        if not state.has_key('_use_caching'):
            self.__dict__['_use_caching'] = False

        if not state.has_key('_queue'):
            self.__dict__['_queue'] = ['']

        if not state.has_key('_cache'):
            self.__dict__['_cache'] = {}

        if not state.has_key('cache'):
            self.__dict__['cache'] = 5


    def __getitem__(self, filter):
        return FilterModel(self, filter)

    def startup(self):
        self._queue = ['']
        self._cache = {}
        if int(self.cache) > 0:
            self._queue = ['']*int(self.cache)
            frozen = numpy.array([par.frozen for par in self.pars], dtype=bool)
            if len(frozen) > 0 and frozen.all():
                self._use_caching = True

    def teardown(self):
        self._use_caching = False

    def apply(self, outer, *otherargs, **otherkwargs):
        return NestedModel(outer, self, *otherargs, **otherkwargs)

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


    def startup(self):
        self.lhs.startup()
        self.rhs.startup()
        CompositeModel.startup(self)


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
        #print 'Starting up %s...' % type(self).__name__
        pass

    def teardown(self):
        #print 'Tearing down %s...' % type(self).__name__
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


    def startup(self):
        self.inner.startup()
        self.outer.startup()
        CompositeModel.startup(self)


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
        for model, args in izip(self.models, arglist):
            # FIXME: we're not using p here (and therefore assuming that the
            # parameter values have already been updated to match the contents
            # of p)
            vals.append(model(*args))
        return sum(vals)
