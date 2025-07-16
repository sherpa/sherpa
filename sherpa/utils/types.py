#
#  Copyright (C) 2024-2025
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

"""Useful types for Sherpa.

This module should be considered an internal module as its contents are
likely to change as types get added to Sherpa and the typing ecosystem
in Python matures.

"""

from collections.abc import Callable, Sequence
from typing import Any, Concatenate, ParamSpec, Protocol

import numpy as np


# Some of these may change to TypeAlias (once Python 3.10 is the
# minimum) and then type/TypeAliasType once Python 3.12 is the
# minimum.
#

P = ParamSpec('P')

# Represent identifiers; mainly used in the UI code.
#
IdType = int | str

# Try to be generic when using arrays as input or output. There is no
# attempt to encode the data type or shape for ndarrays at this time.
#
ArrayType = Sequence[float] | np.ndarray

# Return the statistic and a per-bin value.
#
StatResults = tuple[float, np.ndarray]

# Represent statistic evaluation.
#
StatFunc = Callable[[ArrayType], StatResults]

# What is the best typing rule here?
#
StatErrFunc = Callable[[ArrayType], ArrayType]

# What do the optimization functions return?
#
OptReturn = tuple[bool,            # did the optimiser succeed
                  np.ndarray,      # final parameter values
                  float,           # final statistic value
                  str,             # message
                  dict[str, Any]]  # information to pass back

# The optimizers are sent the statistic, data, and other
# arguments, and return OptReturn.
#
OptFunc = Callable[Concatenate[StatFunc,
                               ArrayType,  # starting position
                               ArrayType,  # minimum values
                               ArrayType,  # maximum values
                               P],
                   OptReturn]


# Represent model evaluation. Using a Protocol may be better, but
# for now keep with a Callable. Ideally the model would just
# return an ndarray but
#
# - an ArithmeticConstantModel can return a scalar
# - models could return a sequence rather than a ndarray
#
ModelFunc = Callable[..., ArrayType]


# Fit-like function used by both the optimisation and fit modules.
#
class FitFunc(Protocol):
    def __call__(self,
                 statfunc: StatFunc,
                 pars: ArrayType,
                 parmins: ArrayType,
                 parmaxes: ArrayType,
                 statargs: Any = None,
                 statkwargs: Any = None
                 ) -> OptReturn:
        ...
