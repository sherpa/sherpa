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

"""Useful types for Sherpa.

This module should be considered an internal module as its contents are
likely to change as types get added to Sherpa and the typing ecosystem
in Python matures.

"""

from typing import Callable, Sequence, Union

import numpy as np


# Some of these may change to TypeAlias (once Python 3.10 is the
# minimum) and then type/TypeAliasType once Python 3.12 is the
# minimum.
#

# Represent identifiers; mainly used in the UI code.
#
IdType = Union[int, str]

# Try to be generic when using arrays as input or output. There is no
# attempt to encode the data type or shape for ndarrays at this time.
#
ArrayType = Union[Sequence[float], np.ndarray]

# Represent statistic evaluation.
#
# Ideally the callable routines would be labelled as accepting
# [ArrayType] and return np.ndarray, but we do not enforce this yet
# (and it is not entirely clear whether the stat functions may accept
# multiple arguments).
#
StatErrFunc = Callable[..., ArrayType]
StatResults = tuple[float, np.ndarray]
StatFunc = Callable[..., StatResults]


# Represent model evaluation. Using a Protocol may be better, but
# for now keep with a Callable. Ideally the model would just
# return an ndarray but
#
# - an ArithmeticConstantModel can return a scalar
# - models could return a sequence rather than a ndarray
#
ModelFunc = Callable[..., ArrayType]
