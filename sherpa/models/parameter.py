# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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
from sherpa.utils import SherpaFloat, NoNewAttributesAfterInit
from sherpa.utils.err import ParameterErr

warning = logging.getLogger(__name__).warning


__all__ = ('Parameter', 'CompositeParameter', 'ConstantParameter',
           'UnaryOpParameter', 'BinaryOpParameter')


# Default minimum and maximum magnitude for parameters
#tinyval = 1.0e-120
#hugeval = 1.0e+120
tinyval = numpy.float(numpy.finfo(numpy.float32).tiny) # FLT_TINY
hugeval = numpy.float(numpy.finfo(numpy.float32).max)  # FLT_MAX
#tinyval = 1.0e-38
#hugeval = 1.0e+38

def _make_set_limit(name):
    def _set_limit(self, val):
        val = SherpaFloat(val)
        # Ensure that we don't try to set any value that is outside
        # the hard parameter limits.
        if val < self._hard_min:
            raise ParameterErr('edge', self.fullname, 'hard minimum', self._hard_min)
        if val > self._hard_max:
            raise ParameterErr('edge', self.fullname, 'hard maximum', self._hard_max)
        # Ensure that we don't try to set a parameter range, such that
        # the minimum will be greater than the current parameter value,
        # or that the maximum will be less than the current parameter value.

        # But we only want to do this check *after* parameter has been
        # created and fully initialized; we are doing this check some time
        # *later*, when the user is trying to reset a parameter range
        # such that the new range will leave the current value
        # *outside* the new range.  We want to warn against and disallow that.

        # Due to complaints about having to rewrite existing user scripts,
        # downgrade the ParameterErr issued here to mere warnings.  Also,
        # set the value to the appropriate soft limit.
        if (hasattr(self, "_NoNewAttributesAfterInit__initialized") == True and
            self._NoNewAttributesAfterInit__initialized == True):
            if (name == "_min"):
                if (val > self.val):
                    self.val = val
                    warning(('parameter %s less than new minimum; %s reset to %g') % (self.fullname, self.fullname, self.val))
            if (name == "_max"):
                if (val < self.val):
                    self.val = val
                    warning(('parameter %s greater than new maximum; %s reset to %g') % (self.fullname, self.fullname, self.val))

        setattr(self, name, val)
    return _set_limit

def _make_unop(op, opstr):
    def func(self):
        return UnaryOpParameter(self, op, opstr)
    return func

def _make_binop(op, opstr):
    def func(self, rhs):
        return BinaryOpParameter(self, rhs, op, opstr)
    def rfunc(self, lhs):
        return BinaryOpParameter(lhs, self, op, opstr)
    return (func, rfunc)


class Parameter(NoNewAttributesAfterInit):

    #
    # Read-only properties
    #

    def _get_alwaysfrozen(self):
        return self._alwaysfrozen
    alwaysfrozen = property(_get_alwaysfrozen)

    def _get_hard_min(self):
        return self._hard_min
    hard_min = property(_get_hard_min)

    def _get_hard_max(self):
        return self._hard_max
    hard_max = property(_get_hard_max)

    #
    # 'val' property
    #
    
    def _get_val(self):
        if hasattr(self, 'eval'):
            return self.eval()
        if self.link is not None:
            return self.link.val
        return self._val    

    def _set_val(self, val):
        if isinstance(val, Parameter):
            self.link = val
        else:
            # Reset link
            self.link = None

            # Validate new value
            val = SherpaFloat(val)
            if val < self.min:
                raise ParameterErr('edge', self.fullname, 'minimum', self.min)
            if val > self.max:
                raise ParameterErr('edge', self.fullname, 'maximum', self.max)
            
            self._val = val
            self._default_val = val

    val = property(_get_val, _set_val)

    #
    # '_default_val' property
    #
    
    def _get_default_val(self):
        if hasattr(self, 'eval'):
            return self.eval()
        if self.link is not None:
            return self.link.default_val
        return self._default_val    

    def _set_default_val(self, default_val):
        if isinstance(default_val, Parameter):
            self.link = default_val
        else:
            # Reset link
            self.link = None

            # Validate new value
            default_val = SherpaFloat(default_val)
            if default_val < self.min:
                raise ParameterErr('edge', self.fullname, 'minimum', self.min)
            if default_val > self.max:
                raise ParameterErr('edge', self.fullname, 'maximum', self.max)
            
            self._default_val = default_val

    default_val = property(_get_default_val, _set_default_val)

    #
    # 'min' and 'max' properties
    #

    def _get_min(self):
        return self._min
    min = property(_get_min, _make_set_limit('_min'))

    def _get_max(self):
        return self._max
    max = property(_get_max, _make_set_limit('_max'))

    #
    # 'default_min' and 'default_max' properties
    #

    def _get_default_min(self):
        return self._default_min
    default_min = property(_get_default_min, _make_set_limit('_default_min'))

    def _get_default_max(self):
        return self._default_max
    default_max = property(_get_default_max, _make_set_limit('_default_max'))


    #
    # 'frozen' property
    #
    
    def _get_frozen(self):
        if self.link is not None:
            return True
        return self._frozen
    def _set_frozen(self, val):
        val = bool(val)
	if self._alwaysfrozen and (not val):
	    raise ParameterErr('alwaysfrozen', self.fullname)
        self._frozen = val
    frozen = property(_get_frozen, _set_frozen)

    #
    # 'link' property'
    #

    def _get_link(self):
        return self._link
    def _set_link(self, link):
        if link is not None:
            if self._alwaysfrozen:
                raise ParameterErr('frozennolink', self.fullname)
            if not isinstance(link, Parameter):
                raise ParameterErr('notlink')

            # Short cycles produce error
            # e.g. par = 2*par+3
            if self in link:
               raise ParameterErr('linkcycle')

            # Correctly test for link cycles in long trees.
            cycle = False
            ll = link
            while isinstance(ll, Parameter):
                if ll == self or self in ll:
                    cycle = True
                ll = ll.link

            # Long cycles are overwritten BUG #12287
            if cycle and isinstance(link, Parameter):
                link.link = None

        self._link = link
    link = property(_get_link, _set_link)

    #
    # Methods
    #

    def __init__(self, modelname, name, val, min=-hugeval, max=hugeval,
		 hard_min=-hugeval, hard_max=hugeval, units='',
                 frozen=False, alwaysfrozen=False, hidden=False):
	self.modelname = modelname
	self.name = name
	self.fullname = '%s.%s' % (modelname, name)

	self._hard_min = SherpaFloat(hard_min)
	self._hard_max = SherpaFloat(hard_max)
	self.units = units

	self._alwaysfrozen = bool(alwaysfrozen)
	if alwaysfrozen:
	    self._frozen = True
	else:
	    self._frozen = frozen

        self.hidden = hidden

	# Set validated attributes.  Access them via their properties so that
	# validation takes place.
	self.min = min
	self.max = max
	self.val = val
        self.default_min = min
        self.default_max = max
        self.default_val = val
        self.link = None
        self._guessed = False

        NoNewAttributesAfterInit.__init__(self)

    def __iter__(self):
        return iter([self])

    def __repr__(self):
        r = "<%s '%s'" % (type(self).__name__, self.name)
        if self.modelname:
            r += " of model '%s'" % self.modelname
        r += '>'
        return r

    def __str__(self):
        if self.link is not None:
            linkstr = self.link.fullname
        else:
            linkstr = str(None)

	return (('val         = %s\n' +
		 'min         = %s\n' +
		 'max         = %s\n' +
		 'units       = %s\n' +
		 'frozen      = %s\n' +
                 'link        = %s\n'
                 'default_val = %s\n' +
                 'default_min = %s\n' +
                 'default_max = %s') %
		(str(self.val), str(self.min), str(self.max), self.units,
                 self.frozen, linkstr, str(self.default_val), str(self.default_min),
                 str(self.default_max)))

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

    def freeze(self):
        self.frozen = True

    def thaw(self):
        self.frozen = False

    def unlink(self):
        self.link = None

    def reset(self):
        # circumvent the attr checks for simplicity, as the defaults have
        # already passed (defaults either set by user or through self.set).
        if self._guessed:
            self._min = self.default_min
            self._max = self.default_max
            self._guessed = False
        self._val = self.default_val


    def set(self, val=None, min=None, max=None, frozen=None,
            default_val=None, default_min=None, default_max=None):

        if max is not None and max > self.max:
            self.max = max
        if default_max is not None and default_max > self.default_max:
            self.default_max = default_max


        if min is not None and min < self.min:
            self.min = min
        if default_min is not None and default_min < self.default_min:
            self.default_min = default_min


        if val is not None:
            self.val = val
        if default_val is not None:
            self.default_val = default_val

        if min is not None:
            self.min = min
        if max is not None:
            self.max = max

        if default_min is not None:
            self.default_min = default_min
        if default_max is not None:
            self.default_max = default_max

        if frozen is not None:
            self.frozen = frozen


class CompositeParameter(Parameter):

    def __init__(self, name, parts):
        self.parts = tuple(parts)
        Parameter.__init__(self, '', name, 0.0)
        self.fullname = name

    def __iter__(self):
        return iter(self._get_parts())

    def _get_parts(self):
        parts = []

        for p in self.parts:
            # A CompositeParameter should not hold a reference to itself
            assert (p is not self), (("'%s' object holds a reference to " +
                                      "itself") % type(self).__name__)

            parts.append(p)
            if isinstance(p, CompositeParameter):
                parts.extend(p._get_parts())

        # FIXME: do we want to remove duplicate components from parts?

        return parts

    def eval(self):
        raise NotImplementedError


class ConstantParameter(CompositeParameter):

    def __init__(self, value):
        self.value = SherpaFloat(value)
        CompositeParameter.__init__(self, str(value), ())

    def eval(self):
        return self.value


class UnaryOpParameter(CompositeParameter):

    def __init__(self, arg, op, opstr):
        self.arg = arg
	self.op = op
	CompositeParameter.__init__(self,
                                    '%s(%s)' % (opstr, self.arg.fullname),
                                    (self.arg,))

    def eval(self):
        return self.op(self.arg.val)


class BinaryOpParameter(CompositeParameter):

    @staticmethod
    def wrapobj(obj):
	if isinstance(obj, Parameter):
	    return obj
	return ConstantParameter(obj) 

    def __init__(self, lhs, rhs, op, opstr):
        self.lhs = self.wrapobj(lhs)
        self.rhs = self.wrapobj(rhs)
	self.op = op
	CompositeParameter.__init__(self, '(%s %s %s)' %
                                    (self.lhs.fullname, opstr,
                                     self.rhs.fullname), (self.lhs, self.rhs))

    def eval(self):
        return self.op(self.lhs.val, self.rhs.val)
