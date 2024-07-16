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

"""Useful types for Sherpa Astronomy I/O.

This module should be considered an internal module as its contents is
likely to change as types get added to Sherpa and the typing ecosystem
in Python matures.

"""

from typing import Any, Mapping, Sequence, Union

import numpy as np


# Some variants are named xxx and xxxArg, where the former is the
# return value (an invariant type, like dict) and the latter is
# covariant (such as Mapping) as it's used as an argument to a
# function.
#
KeyType = Union[str, bool, int, float]
NamesType = Sequence[str]
HdrTypeArg = Mapping[str, KeyType]
HdrType = dict[str, KeyType]

# Note that ColumnsTypeArg allows more for the values than does
# ColumnsType.
#
ColumnsTypeArg = Mapping[str, Union[np.ndarray, list, tuple]]
ColumnsType = dict[str, np.ndarray]

# It's hard to type the arguments to the Data constructors
DataTypeArg = Mapping[str, Any]
DataType = dict[str, Any]
