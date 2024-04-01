#
#  Copyright (C) 2024
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

"""
Operator-related routines for the model and parameter classes.
"""

from typing import Any, Callable

import numpy as np


__all__ = ('get_precedences_op', 'get_precedence_expr',
           'get_precedence_lhs', 'get_precedence_rhs',
           )


# We make remainder and floor_divide have the same precedence
# as power (aka **) just to make things clear.
#
lprec : dict[Callable, int]
lprec = {np.power: 4,
         np.remainder: 4,
         np.floor_divide: 4,
         np.multiply: 3,
         np.divide: 3,
         np.true_divide: 3,
         np.add: 2,
         np.subtract: 2
         }
"""Represent the precedence for the left-hand side of a binary operator.

All values must be below the `sherpa.models.op.NO_PRECEDENCE` value.
"""


rprec : dict[Callable, bool]
rprec = {np.power: True,
         np.remainder: True,
         np.floor_divide: True
         }
"""Should the left-hand side of a binary operator be wrapped by brackets?"""


NO_PRECEDENCE : int
NO_PRECEDENCE = 9
"""Represent "no precedence"."""


# Expression terms (used to remove excess brackets from model and
# parameter expressions).
#
def get_precedences_op(op: Callable) -> tuple[int, bool]:
    """Return precedences for the operation.

    Unrecognized parameters are mapped to
    (`sherpa.models.op.NO_PRECEDENCE`, False).

    Parameters
    ----------
    op : callable
       The operator (e.g. np.multiply or np.power).

    Returns
    -------
    lprec, rprec : (int, bool)
       The left and right "precedences" (the right case only cares
       about being set or not hence we use a boolean).

    See Also
    --------
    get_precedence_expr, get_precedence_lhs, get_precedence_rhs

    Notes
    -----
    This uses the module-level symbols: `sherpa.models.op.lprec`
    and `sherpa.models.op.rprec`.

    """

    return lprec.get(op, NO_PRECEDENCE), rprec.get(op, False)


# Typing is hard to get right here given that we want to allow this
# for both models and parameters.
#
def get_precedence_expr(expr: Any) -> int:
    """Return precedence for the expression.

    This is the precedence of the operator in the expression,
    if one exists.

    Parameters
    ----------
    expr : Model or Parameter
       The expression. It may have a .opprec field.

    Returns
    -------
    prec : int
       The "precedence" value `sherpa.models.op.NO_PRECEDENCE` is
       returned if expr has an unknown operator or it does not
       contain an operator.

    See Also
    --------
    get_precedences_op, get_precedence_lhs, get_precedence_rhs

    """

    try:
        return expr.opprec
    except AttributeError:
        return NO_PRECEDENCE


def get_precedence_lhs(lstr: str, lp: int, p: int, a: bool) -> str:
    """Return the string to use for the left side of a binary operator.

    This is only needed when the operator (represented by the p value)
    is written in infix form - such as 'a + b'- rather than prefix
    form like 'np.add(a, b)'.

    Parameters
    ----------
    lstr : str
       The term to the left of the operator.
    lp : int
       Precedences of any operator in lstr.
    p : int
       Precedences of the current operator.
    a : bool
       Do we care about power-like terms.

    Returns
    -------
    term : str
       Either lstr or (lstr).

    See Also
    --------
    get_precedences_op, get_precedence_expr, get_precedence_rhs

    Examples
    --------

    >>> get_precedence_lhs("a", NO_PRECEDENCE, lprec[np.divide], False)
    'a'

    >>> get_precedence_lhs("a + b", lprec[np.add], lprec[np.divide], False)
    '(a + b)'

    >>> get_precedence_lhs("a * b", lprec[np.multiply], lprec[np.add], False)
    'a * b'

    """

    if lp < p:
        return f"({lstr})"

    if not a:
        return lstr

    # We could combine all these into one, but for now keep
    # them separate.
    #
    if lp == p:
        # For now ensure the left term is bracketed just to be
        # clear.
        #
        return f"({lstr})"

    if lstr[0] == "-":
        # If the term is negative then ensure it is included
        # in a bracket before passed through to the power
        # term.  Python has "-a ** 2" actually mapping to "-(a
        # ** 2)" so we need to say "(-a) ** 2" if we really
        # want the unary operator before the power term.
        #
        return f"({lstr})"

    return lstr


def get_precedence_rhs(rstr: str, opstr: str, rp: int, p: int) -> str:
    """Return the string to use for the right side of a binary operator.

    This is only needed when the operator (represented by the p value)
    is written in infix form - such as 'a + b'- rather than prefix
    form like 'np.add(a, b)'.

    Parameters
    ----------
    rstr : str
       The term to the right of the operator.
    opstr : str
       The string representing the operator.
    rp : int
       Precedences of any operator in rstr.
    p : int
       Precedences of the current operator.

    Returns
    -------
    term : str
       Either rstr or (rstr).

    See Also
    --------
    get_precedences_op, get_precedence_expr, get_precedence_lhs

    Examples
    --------

    >>> get_precedence_rhs("a", "/", NO_PRECEDENCE, lprec[np.divide])
    'a'

    >>> get_precedence_rhs("a + b", "/", lprec[np.add], lprec[np.divide])
    '(a + b)'

    """

    if opstr in ["+", "*"]:
        condition = rp < p
    else:
        condition = rp <= p

    if condition:
        return f"({rstr})"

    return rstr
