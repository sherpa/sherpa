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
"""Represent FITS-like concepts for the backends.

This is an experimental module, and is likely to change.

"""

from dataclasses import dataclass
from typing import Any, Optional, Union


__all__ = ("HeaderItem", "Column", "TableHDU")


# Perhaps the HeaderItem and Column types can be used in
# sherpa.astro.io to pass information to and from the backends, rather
# than the current system. However, that is for later work.
#
@dataclass
class HeaderItem:
    """Represent a FITS header card.

    This does not support all FITS features.
    """
    name: str
    """The keyword name (case insensitive)"""
    value: Union[str, bool, int, float]  # will we need more?
    """The keyword value"""
    desc: Optional[str] = None
    """The description for the keyword"""
    unit: Optional[str] = None
    """The units of the value"""


@dataclass
class Column:
    """Represent a FITS column.

    This does not support all FITS features.
    """
    name: str
    """The name of the column (will be converted to upper case and leading/trailing spaces removed)."""
    values: Any  # should be typed
    """The values for the column, as a ndarray.

    Variable-field arrays are represented as ndarrays with an object
    type.
    """
    desc: Optional[str] = None
    """The column description"""
    unit: Optional[str] = None
    """The units of the column"""

    def __post_init__(self):
        self.name = self.name.strip().upper()
        if self.name == "":
            raise ValueError("name can not be empty")


# To be generic this should probably be split into Primary, Table, and
# Image HDUs.
#
@dataclass
class TableHDU:
    """Represent a HDU: header and optional columns"""
    name: str
    """The name of the HDU (will be converted to upper case and leading/trailing spaces removed)."""
    header: list[HeaderItem]
    """The header values"""
    data: Optional[list[Column]] = None
    """The column data.

    This should be empty for a primary header, and have at least one
    entry for a table HDU.
    """

    def __post_init__(self):
        self.name = self.name.strip().upper()
        if self.name == "":
            raise ValueError("name can not be empty")
