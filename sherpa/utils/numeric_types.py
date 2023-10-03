#
#  Copyright (C) 2023
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

"""Basic numeric types.

The idea is to match the definitions used for the C/C++ extension
code.

"""

import numpy

__all__ = ("SherpaInt", "SherpaUInt", "SherpaFloat")

# Default numeric types (these match the typedefs in extension.hh)
#
SherpaInt = numpy.intp
"""The default signed integer type."""

SherpaUInt = numpy.uintp
"""The default unsigner-integer type."""

SherpaFloat = numpy.float64
"""The default floating-point type."""
