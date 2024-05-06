#
#  Copyright (C) 2007, 2017, 2020 - 2024
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

"""Support for model parameter values.

Parameter creation, evaluation, and combination are normally done as
part of the model interface provided by
sherpa.models.model.ArithmeticModel.

In the following the p variable corresponds to

    >>> from sherpa.models.parameter import Parameter
    >>> p = Parameter('model', 'eta', 2)
    >>> print(p)
    val         = 2.0
    min         = -3.4028234663852886e+38
    max         = 3.4028234663852886e+38
    units       =
    frozen      = False
    link        = None
    default_val = 2.0
    default_min = -3.4028234663852886e+38
    default_max = 3.4028234663852886e+38

Naming
======

The first two arguments to a parameter are the model name and then the
parameter name, so `"model"` and `"eta"` here. They are used to
display the `name` and `fullname` attributes:

    >>> p.name
    'eta'
    >>> p.fullname
    'model.eta'

A units field can be attached to the parameter, but this is used purely
for screen output by the sherpa.models.model.ArithmeticModel class and
is not used when changing or evaluating the parameter.

Changing parameter values
=========================

The `val` attribute is used to retrieve or change the parameter value:

    >>> p.val
    2.0
    >>> p.val = 3

Parameter limits
================

The parameter is forced to lie within the `min` and `max` attributes
of the parameter (known as the "soft" limits). The default values
for these are the 32-bit floating point maximum value, and it's
negative:

    >>> p.max
    3.4028234663852886e+38
    >>> p.min
    -3.4028234663852886e+38

Setting a value outside this range will raise a
`sherpa.utils.err.ParameterErr` exception:

    >>> p.val = 1e40
    Traceback (most recent call last):
      ...
    sherpa.utils.err.ParameterErr: parameter model.eta has a maximum of 3.40282e+38

These limits can be changed, as shown below, but they must lie
within the `hard_min` to `hard_max` range of the parameter (the
"hard" limits), which can not be changed:

    >>> p.min = 0
    >>> p.max = 10

Freezing and thawing
====================

When fitting a model expression it is useful to be able to restrict
the fit to a subset of parameters. This is done by only selecting
those parameters which are not "frozen". This can be indicated by
calling the `freeze` and `thaw` methods, or changing the `frozen`
attribute directly:

    >>> p.frozen
    False
    >>> p.freeze()
    >>> p.frozen
    True
    >>> p.thaw()
    >>> p.frozen
    False

Note that the `frozen` flag is used to indicate what parameters to
vary in a fit, but it is still possible to directly change a parameter
value when it is frozen:

    >>> p.freeze()
    >>> p.val = 6
    >>> print(p)
    val         = 6.0
    min         = 0.0
    max         = 10.0
    units       =
    frozen      = True
    link        = None
    default_val = 6.0
    default_min = -3.4028234663852886e+38
    default_max = 3.4028234663852886e+38

Changing multiple settings at once
==================================

The `set` method should be used when multiple settings need to be
changed at once, as it allows for changes to both the value and
limits, such as changing the value to be 20 and the limits to 8 to 30:

    >>> p.val = 20
    Traceback (most recent call last):
      ...
    sherpa.utils.err.ParameterErr: parameter model.eta has a maximum of 10
    >>> p.set(val=20, min=8, max=30)
    >>> print(p)
    val         = 20.0
    min         = 8.0
    max         = 30.0
    units       =
    frozen      = True
    link        = None
    default_val = 20.0
    default_min = -3.4028234663852886e+38
    default_max = 3.4028234663852886e+38

Linking parameters
==================

A parameter can be "linked" to another parameter, in which case the
value of the parameter is calculated based on the link expression,
such as being twice the other parameter:

    >>> q = Parameter('other', 'beta', 4)
    >>> p.val = 2 * q
    >>> print(p)
    val         = 8.0
    min         = 8.0
    max         = 30.0
    units       =
    frozen      = True
    link        = 2 * other.beta
    default_val = 8.0
    default_min = -3.4028234663852886e+38
    default_max = 3.4028234663852886e+38

The `link` attribute stores the expression:

    >>> p.link
    <BinaryOpParameter '2 * other.beta'>

A ParameterErr exception will be raised whenever the linked expression
is evaluated and the result lies outside the parameters soft limits
(the `min` to `max` range). For this example, p must lie between 8 and
30 and so changing parameter q to a value of 3 will cause an error,
but only when parameter p is checked, not when the related parameter
(here q) is changed:

    >>> q.val = 3
    >>> print(p)
    Traceback (most recent call last):
      ...
    sherpa.utils.err.ParameterErr: parameter model.eta has a minimum of 8
    >>> p.val
    Traceback (most recent call last):
      ...
    sherpa.utils.err.ParameterErr: parameter model.eta has a minimum of 8

Resetting a parameter
=====================

The `reset` method is used to restore the parameter value and soft
limits to a known state. The idea is that if the values were changed
in a fit then `reset` will change them back to the values before a
fit.

"""

from __future__ import annotations

import logging
from typing import Any, Callable, Iterator, Optional, \
    Sequence, SupportsFloat, Union

import numpy as np

from sherpa.utils import NoNewAttributesAfterInit, formatting
from sherpa.utils.err import ParameterErr
from sherpa.utils.numeric_types import SherpaFloat

from .op import get_precedences_op, get_precedence_expr, \
    get_precedence_lhs, get_precedence_rhs

warning = logging.getLogger(__name__).warning


__all__ = ('Parameter', 'CompositeParameter', 'ConstantParameter',
           'UnaryOpParameter', 'BinaryOpParameter')


# Default minimum and maximum magnitude for parameters
# tinyval = 1.0e-120
# hugeval = 1.0e+120
# tinyval = 1.0e-38
# hugeval = 1.0e+38
#
# Use FLT_TINY and FLT_MAX
tinyval = float(np.finfo(np.float32).tiny)
hugeval = float(np.finfo(np.float32).max)


def _make_set_limit(name: str) -> Callable[[Any, SupportsFloat], None]:
    def _set_limit(self: Any,
                   val: SupportsFloat) -> None:
        val = SherpaFloat(val)
        # Ensure that we don't try to set any value that is outside
        # the hard parameter limits.
        if val < self._hard_min:
            raise ParameterErr('edge', self.fullname,
                               'hard minimum', self._hard_min)
        if val > self._hard_max:
            raise ParameterErr('edge', self.fullname,
                               'hard maximum', self._hard_max)

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
        if hasattr(self, "_NoNewAttributesAfterInit__initialized") and \
           self._NoNewAttributesAfterInit__initialized:
            if name == "_min" and (val > self.val):
                self.val = val
                warning('parameter %s less than new minimum; '
                        '%s reset to %g', self.fullname, self.fullname,
                        self.val)

            if name == "_max" and (val < self.val):
                self.val = val
                warning('parameter %s greater than new maximum; '
                        '%s reset to %g', self.fullname, self.fullname,
                        self.val)

        setattr(self, name, val)

    return _set_limit


# It's hard to come up with a sensible typing rule for the operator
# in _make_unop/binop so we just use Callable.
#
def _make_unop(op: Callable,
               opstr: str,
               strformat: Optional[str] = None) -> Callable:

    if strformat is None:
        def func(self):
            return UnaryOpParameter(self, op, opstr)
    else:
        def func(self):
            return UnaryOpParameter(self, op, opstr, strformat=strformat)

    return func


def _make_binop(op: Callable,
                opstr: str) -> tuple[Callable, Callable]:
    def func(self, rhs):
        return BinaryOpParameter(self, rhs, op, opstr)

    def rfunc(self, lhs):
        return BinaryOpParameter(lhs, self, op, opstr)

    return (func, rfunc)


class Parameter(NoNewAttributesAfterInit):
    """Represent a model parameter.

    Parameters
    ----------
    modelname : str
        The name of the model component containing the parameter.
    name : str
        The name of the parameter. It should be considered to be
        matched in a case-insensitive manner.
    val : number
        The default value for the parameter.
    min, max, hard_min, hard_max: number, optional
        The soft and hard limits for the parameter value.
    units : str, optional
        The units for the parameter value.
    frozen : bool, optional
        Does the parameter default to being frozen?
    alwaysfrozen : bool, optional
        If set then the parameter can never be thawed.
    hidden : bool, optional
        Should the parameter be included when displaying the model
        contents?
    aliases : None or list of str
        If not None then alternative names for the parameter (these
        are expected to be matched in a case-insensitive manner).

    """

    #
    # Read-only properties
    #

    def _get_alwaysfrozen(self) -> bool:
        return self._alwaysfrozen
    alwaysfrozen = property(_get_alwaysfrozen,
                            doc='Is the parameter always frozen?')

    def _get_hard_min(self) -> SherpaFloat:
        return self._hard_min
    hard_min = property(_get_hard_min,
                        doc="""The hard minimum of the parameter.

See Also
--------
hard_max
""")

    def _get_hard_max(self) -> SherpaFloat:
        return self._hard_max
    hard_max = property(_get_hard_max,
                        doc="""The hard maximum of the parameter.

See Also
--------
hard_min
""")

    # 'val' property
    #
    # Note that _get_val has to check the parameter value when it
    # is a link, to ensure that it isn't outside the parameter's
    # min/max range. See issue #742.
    #
    _val: SherpaFloat  # needed for typing

    def _get_val(self) -> SherpaFloat:
        if hasattr(self, 'eval'):
            return self.eval()
        if self.link is None:
            return self._val

        val = self.link.val
        if val < self.min:
            raise ParameterErr('edge', self.fullname, 'minimum', self.min)
        if val > self.max:
            raise ParameterErr('edge', self.fullname, 'maximum', self.max)

        return val

    def _set_val(self, val: Union[Parameter, SupportsFloat]) -> None:
        if isinstance(val, Parameter):
            self.link = val
            return

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

    val = property(_get_val, _set_val,
                   doc="""The current value of the parameter.

If the parameter is a link then it is possible that accessing
the value will raise a ParameterErr in cases where the link
expression falls outside the soft limits of the parameter.

See Also
--------
default_val, link, max, min
""")

    #
    # '_default_val' property
    #

    def _get_default_val(self) -> SherpaFloat:
        if hasattr(self, 'eval'):
            return self.eval()
        if self.link is not None:
            return self.link.default_val
        return self._default_val

    def _set_default_val(self,
                         default_val: Union[Parameter, SupportsFloat]) -> None:
        if isinstance(default_val, Parameter):
            self.link = default_val
            return

        # Reset link
        self.link = None

        # Validate new value
        default_val = SherpaFloat(default_val)
        if default_val < self.min:
            raise ParameterErr('edge', self.fullname, 'minimum', self.min)
        if default_val > self.max:
            raise ParameterErr('edge', self.fullname, 'maximum', self.max)

        self._default_val = default_val

    default_val = property(_get_default_val, _set_default_val,
                           doc="""The default value of the parameter.

See Also
--------
val
""")

    #
    # 'min' and 'max' properties
    #

    def _get_min(self) -> SupportsFloat:
        return self._min
    min = property(_get_min, _make_set_limit('_min'),
                   doc="""The minimum value of the parameter.

The minimum must lie between the hard_min and hard_max limits.

See Also
--------
max, val
""")

    def _get_max(self) -> SupportsFloat:
        return self._max
    max = property(_get_max, _make_set_limit('_max'),
                   doc="""The maximum value of the parameter.

The maximum must lie between the hard_min and hard_max limits.

See Also
--------
min, val
""")

    #
    # 'default_min' and 'default_max' properties
    #
    _default_min: SupportsFloat  # needed for typnig
    _default_max: SupportsFloat

    def _get_default_min(self) -> SupportsFloat:
        return self._default_min
    default_min = property(_get_default_min, _make_set_limit('_default_min'))

    def _get_default_max(self) -> SupportsFloat:
        return self._default_max

    default_max = property(_get_default_max, _make_set_limit('_default_max'))

    #
    # 'frozen' property
    #
    _frozen: bool  # needed for typing

    def _get_frozen(self) -> bool:
        if self.link is not None:
            return True
        return self._frozen

    def _set_frozen(self, val: bool) -> None:
        val = bool(val)
        if self._alwaysfrozen and (not val):
            raise ParameterErr('alwaysfrozen', self.fullname)
        self._frozen = val

    frozen = property(_get_frozen, _set_frozen,
                      doc="""Is the parameter currently frozen?

Those parameters created with `alwaysfrozen` set can not
be changed.

See Also
--------
alwaysfrozen
""")

    #
    # 'link' property'
    #

    def _get_link(self) -> Optional[Parameter]:
        return self._link

    def _set_link(self, link: Optional[Parameter]) -> None:
        if link is None:
            self._link = None
            return

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
    link = property(_get_link, _set_link,
                    doc="""The link expression to other parameters, if set.

The link expression defines if the parameter is not
a free parameter but is actually defined in terms of
other parameters.

See Also
--------
val

Examples
--------

>>> a = Parameter("mdl", "a", 2)
>>> b = Parameter("mdl", "b", 1)
>>> b.link = 10 - a
>>> a.val
2.0
>>> b.val
8.0
""")

    #
    # Methods
    #

    def __init__(self,
                 modelname: str,
                 name: str,
                 val: SupportsFloat,
                 min: SupportsFloat = -hugeval,
                 max: SupportsFloat = hugeval,
                 hard_min: SupportsFloat = -hugeval,
                 hard_max: SupportsFloat = hugeval,
                 units: str = '',
                 frozen: bool = False,
                 alwaysfrozen: bool = False,
                 hidden: bool = False,
                 aliases: Optional[Sequence[str]] = None) -> None:
        self.modelname = modelname
        self.name = name
        self.fullname = f"{modelname}.{name}"

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

        self.aliases = [a.lower() for a in aliases] if aliases is not None else []

        NoNewAttributesAfterInit.__init__(self)

    def __iter__(self) -> Iterator[Parameter]:
        return iter([self])

    def __repr__(self) -> str:
        r = f"<{type(self).__name__} '{self.name}'"
        if self.modelname:
            r += f" of model '{self.modelname}'"
        r += '>'
        return r

    def __str__(self) -> str:
        if self.link is not None:
            linkstr = self.link.fullname
        else:
            linkstr = str(None)

        out = [f'val         = {self.val}',
               f'min         = {self.min}',
               f'max         = {self.max}',
               f'units       = {self.units}',
               f'frozen      = {self.frozen}',
               f'link        = {linkstr}',
               f'default_val = {self.default_val}',
               f'default_min = {self.default_min}',
               f'default_max = {self.default_max}']
        return "\n".join(out)

    # Support 'rich display' representations
    #
    def _val_to_html(self, v: SupportsFloat) -> str:
        """Convert a value to a string for use by the HTML output.

        The conversion to a string uses the Python defaults for most
        cases. The units field is currently used to determine whether
        to convert angles to factors of pi (this probably should be
        done by a subclass or mixture).
        """

        # The use of equality rather than some form of tolerance
        # should be okay here.
        #
        if v == hugeval:
            return 'MAX'
        if v == -hugeval:
            return '-MAX'

        if v == tinyval:
            return 'TINY'
        if v == -tinyval:
            return '-TINY'

        if self.units in ['radian', 'radians']:
            tau = 2 * np.pi

            if v == tau:
                return '2&#960;'
            if v == -tau:
                return '-2&#960;'

            if v == np.pi:
                return '&#960;'
            if v == -np.pi:
                return '-&#960;'

        return str(v)

    def _units_to_html(self) -> str:
        """Convert the unit to HTML.

        This is provided for future expansion/experimentation,
        and is not guaranteed to remain.
        """

        return self.units

    def _repr_html_(self) -> str:
        """Return a HTML (string) representation of the parameter
        """
        return html_parameter(self)

    # Unary operations
    # It is safest to always say -(..) even if arg is a single field
    __neg__ = _make_unop(np.negative, '-', strformat='-({arg})')
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

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if not method == '__call__':
            return NotImplemented
        if hasattr(np, ufunc.__name__):
            name = f"numpy.{ufunc.__name__}"
        else:
            # Unfortunately, there is no ufunc.__module__ we could use
            # This could be a ufunc from e.g. scipy or generated from an
            # arbitrary Python function with numpy.frompyfunc.
            # In the latter case, the name will be something like
            # "func (vectorized)", which looks confusing in out string
            # representation, so we simplify it.
            name = ufunc.__name__.replace(' (vectorized)', '')
        if ufunc.nin == 1:
            return UnaryOpParameter(inputs[0], ufunc, name)
        if ufunc.nin == 2:
            return BinaryOpParameter(inputs[0], inputs[1], ufunc, name,
                                     strformat='{opstr}({lhs}, {rhs})')
        return NotImplemented

    def freeze(self) -> None:
        """Set the `frozen` attribute for the parameter.

        See Also
        --------
        thaw
        """
        self.frozen = True

    def thaw(self) -> None:
        """Unset the `frozen` attribute for the parameter.

        Raises
        ------
        ParameterErr
           The parameter is marked as always frozen.

        See Also
        --------
        frozen
        """
        self.frozen = False

    def unlink(self) -> None:
        """Remove any link to other parameters."""
        self.link = None

    def reset(self) -> None:
        """Reset the parameter value and limits to their default values."""
        # circumvent the attr checks for simplicity, as the defaults have
        # already passed (defaults either set by user or through self.set).
        if self._guessed:
            # TODO: It is not clear the logic for when _guessed gets set
            # (see sherpa.utils.param_apply_limits) so we do not
            # describe the logic in the docstring yet.
            self._min = self.default_min
            self._max = self.default_max
            self._guessed = False

        self._val = self.default_val

    def set(self,
            val: Optional[SupportsFloat] = None,
            min: Optional[SupportsFloat] = None,
            max: Optional[SupportsFloat] = None,
            frozen: Optional[bool] = None,
            default_val: Optional[SupportsFloat] = None,
            default_min: Optional[SupportsFloat] = None,
            default_max: Optional[SupportsFloat] = None) -> None:
        """Change a parameter setting.

        Parameters
        ----------
        val : number or None, optional
            The new parameter value.
        min, max : number or None, optional
            The new parameter range.
        frozen : bool or None, optional
            Should the frozen flag be set?
        default_val : number or None, optional
            The new default parameter value.
        default_min, default_max : number or None, optional
            The new default parameter limits.
        """

        # The validation checks are left to the individual properties.
        # However, it means that the logic here has to handle cases
        # of 'set(val=1, min=0, max=2)' but a value of 1 lies
        # outside the min/max of the object before the call, and
        # we don't want the call to fail because of this.
        #
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
    """Represent a parameter with composite parts.

    This is the base class for representing expressions that combine
    multiple parameters and values.

    Parameters
    ----------
    name : str
        The name for the collection.
    parts : sequence of Parameter objects
        The parameters.

    Notes
    -----
    Composite parameters can be iterated through to find their
    components:

       >>> p = Parameter('m', 'p', 2)
       >>> q = Parameter('m', 'q', 4)
       >>> c = (p + q) / 2
       >>> c
       <BinaryOpParameter '(m.p + m.q) / 2'>
       >>> for cpt in c:
       ...     print(type(cpt))
       ...
       <class 'sherpa.models.parameter.BinaryOpParameter'>
       <class 'sherpa.models.parameter.Parameter'>
       <class 'sherpa.models.parameter.Parameter'>
       <class 'sherpa.models.parameter.ConstantParameter'>

    """

    def __init__(self,
                 name: str,
                 parts: Sequence[Parameter]) -> None:
        self.parts = tuple(parts)
        Parameter.__init__(self, '', name, 0.0)
        self.fullname = name

    def __iter__(self) -> Iterator[Parameter]:
        return iter(self._get_parts())

    def _get_parts(self) -> list[Parameter]:
        parts = []

        for p in self.parts:
            # A CompositeParameter should not hold a reference to itself
            assert (p is not self), (f"'{type(self).__name__}' object "
                                     "holds a reference to itself")

            parts.append(p)
            if isinstance(p, CompositeParameter):
                parts.extend(p._get_parts())

        # FIXME: do we want to remove duplicate components from parts?

        return parts

    def eval(self) -> SupportsFloat:
        """Evaluate the composite expression."""
        raise NotImplementedError


class ConstantParameter(CompositeParameter):
    """Represent an expression containing 1 or more parameters."""

    def __init__(self, value: SupportsFloat) -> None:
        self.value = SherpaFloat(value)
        CompositeParameter.__init__(self, str(value), ())

        # Ensure the constant is considered frozen
        self._alwaysfrozen = True
        self.frozen = True

    def eval(self) -> SherpaFloat:
        return self.value


class UnaryOpParameter(CompositeParameter):
    """Apply an operator to a parameter expression.

    Parameters
    ----------
    arg : Parameter instance
    op : function reference
        The ufunc to apply to the parameter value.
    opstr : str
        The symbol used to represent the operator.
    strformat : str
        Format string for printing this operation. Elements
        that can be used are `arg` and `opstr`, e.g.
        `strformat='{opstr}({arg})'`.

    See Also
    --------
    BinaryOpParameter
    """

    def __init__(self,
                 arg: Parameter,
                 op: Callable[[SupportsFloat], SupportsFloat],
                 opstr: str,
                 strformat: str = '{opstr}({arg})') -> None:
        self.arg = arg
        self.op = op
        self.opprec = get_precedences_op(op)[0]

        fullname = self.arg.fullname
        name = strformat.format(opstr=opstr, arg=fullname)
        CompositeParameter.__init__(self, name, (self.arg,))

    def eval(self) -> SupportsFloat:
        return self.op(self.arg.val)


class BinaryOpParameter(CompositeParameter):
    """Combine two parameter expressions.

    Parameters
    ----------
    lhs : Parameter instance
        The left-hand side of the expression.
    rhs : Parameter instance
        The right-hand side of the expression.
    op : function reference
        The ufunc to apply to the two parameter values.
    opstr : str
        The symbol used to represent the operator.
    strformat : str
        Format string for printing this operation. Elements
        that can be used are `lhs`, `rhs`, and `opstr`, e.g.
        `strformat='{opstr}({lhs}, {rhs})'`.

    See Also
    --------
    UnaryOpParameter
    """

    @staticmethod
    def wrapobj(obj: Any) -> Parameter:
        if isinstance(obj, Parameter):
            return obj
        return ConstantParameter(obj)

    def __init__(self,
                 lhs: Any,
                 rhs: Any,
                 op: Callable[[SupportsFloat, SupportsFloat], SupportsFloat],
                 opstr: str,
                 strformat: str = '{lhs} {opstr} {rhs}') -> None:
        self.lhs = self.wrapobj(lhs)
        self.rhs = self.wrapobj(rhs)
        self.op = op

        p, a = get_precedences_op(op)
        self.opprec = p

        # Is this an infix or prefix operator? This could be specified
        # explicitly (and, in fact, could replace the use of the
        # strformat argument), but for now the behvaiour is inferred
        # by assuming that a prefix operator has strformat beginning
        # with '{opstr}(. This heuristic is not perfect, but should be
        # sufficient for our needs.
        #
        if strformat.startswith("{opstr}("):
            # This is a prefix form so we do not need to worry about
            # adding brackets around the lhs and rhs terms.
            #
            lstr = self.lhs.fullname
            rstr = self.rhs.fullname

        else:
            # Simplify the expression if possible.
            #
            lp = get_precedence_expr(self.lhs)
            rp = get_precedence_expr(self.rhs)

            lstr = get_precedence_lhs(self.lhs.fullname, lp, p, a)
            rstr = get_precedence_rhs(self.rhs.fullname, opstr, rp, p)

        name = strformat.format(lhs=lstr, rhs=rstr, opstr=opstr)
        CompositeParameter.__init__(self, name, (self.lhs, self.rhs))

    def eval(self) -> SupportsFloat:
        return self.op(self.lhs.val, self.rhs.val)


def expand_par(expr: Parameter) -> list[Parameter]:
    """Return the individual parameter components referenced in expt

    The aim is that changing any of the return values should
    change the supplied parameter expression.

    Parameters
    ----------
    expr : Parameter

    Returns
    -------
    pars : list of Parameter
       The unique parameters in the expression (thawed and frozen).

    Examples
    --------

    >>> p1 = Parameter("m", "a", 5)
    >>> expand_par(p1)
    [<Parameter 'a' of model 'm'>]

    >>> expand_par(p1 + p1)
    [<Parameter 'a' of model 'm'>]

    >>> p2 = Parameter("n", "b", 2)
    >>> expand_par(p1 + p2 + 2)
    [<Parameter 'a' of model 'm'>, <Parameter 'b' of model 'n'>]

    >>> p2.freeze()
    >>> expand_par(p1 + p2 + 2)
    [<Parameter 'a' of model 'm'>, <Parameter 'b' of model 'n'>]

    >>> p3 = Parameter("m", "x", 2)
    >>> p2.link = 3 * p3
    >>> expand_par(p1 + p2 + 2)
    [<Parameter 'a' of model 'm'>, <Parameter 'x' of model 'm'>]

    """

    # The link expression may just be a parameter, or a composite
    # parameter. However iter should handle this.
    #
    out = []
    for par in iter(expr):
        # For some reason we include the composite terms too in the
        # iteration, so just skip then. This also skips
        # ConstantParameter terms, as they happen to extend
        # CompositeParameter.
        #
        if isinstance(par, CompositeParameter):
            continue

        # If par is a link we need to expand the link.
        #
        if par.link:
            try:
                out.extend(expand_par(par.link))
            except RecursionError:
                # Not sure what to do here so include the
                # "problematic" parameter.
                #
                if par not in out:
                    out.append(par)

            continue

        if par in out:
            continue

        out.append(par)

    return out


# Notebook representation
#
def html_parameter(par: Parameter) -> str:
    """Construct the HTML to display the parameter."""

    # Note that as this is a specialized table we do not use
    # formatting.html_table but create everything directly.
    #
    def addtd(val):
        "Use the parameter to convert to HTML"
        return f'<td>{par._val_to_html(val)}</td>'

    out = '<table class="model">'
    out += '<thead><tr>'
    cols = ['Component', 'Parameter', 'Thawed', 'Value',
            'Min', 'Max', 'Units']
    for col in cols:
        out += f'<th>{col}</th>'

    out += '</tr></thead><tbody><tr>'

    out += f'<th class="model-odd">{par.modelname}</th>'
    out += f'<td>{par.name}</td>'

    linked = par.link is not None
    if linked:
        out += "<td>linked</td>"
    else:
        out += '<td><input disabled type="checkbox"'
        if not par.frozen:
            out += ' checked'
        out += '></input></td>'

    out += addtd(par.val)
    if linked:
        # 8592 is single left arrow
        # 8656 is double left arrow
        #
        val = formatting.clean_bracket(par.link.fullname)
        out += f'<td colspan="2">&#8656; {val}</td>'

    else:
        out += addtd(par.min)
        out += addtd(par.max)

    out += f'<td>{par._units_to_html()}</td>'
    out += '</tr>'

    out += '</tbody></table>'

    ls = ['<details open><summary>Parameter</summary>' + out + '</details>']
    return formatting.html_from_sections(par, ls)
