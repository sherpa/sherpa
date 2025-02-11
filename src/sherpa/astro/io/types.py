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
in Python matures. These types are primarily for communicating
between the front and back ends, as well as providing some structure
to FITS-like data structures.

"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Mapping, Optional, Sequence, Union

import numpy as np

if TYPE_CHECKING:
    from sherpa.astro.io.wcs import WCS

# Useful names to export
#
__all__ = ("KeyType", "HeaderItem", "Header", "Column", "Block",
           "TableBlock", "ImageBlock", "SpectrumBlock",
           "SpecrespBlock", "MatrixBlock", "EboundsBlock",
           "BlockList", "BlockType")


# Some variants are named xxx and xxxArg, where the former is the
# return value (an invariant type, like dict) and the latter is
# covariant (such as Mapping) as it's used as an argument to a
# function.
#
# Note that there are issues with whether we also want to include
# NumPy types - e.g. np.bool_, np.integer, np.floating - in some of
# these types.
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
#
DataTypeArg = Mapping[str, Any]
DataType = dict[str, Any]

# More-specific types designed to handle FITS-like data. The aim is to
# provide value/type checking at object construction time.
#
@dataclass
class HeaderItem:
    """Represent a FITS header card.

    This does not support all FITS features.
    """

    name: str
    """The keyword name (case insensitive, no spaces)"""

    value: KeyType
    """The keyword value"""

    desc: Optional[str] = None
    """The description for the keyword"""

    unit: Optional[str] = None
    """The units of the value"""

    def __post_init__(self) -> None:

        if not isinstance(self.name, str):
            raise ValueError(f"Invalid HDU name: '{self.name}' ({type(self.name)})")

        if " " in self.name:
            raise ValueError(f"Invalid key name: '{self.name}'")

        if not isinstance(self.value, (bool, int, str, float, np.bool_, np.integer, np.floating)):
            raise ValueError(f"Invalid key: {self.name} = {self.value} ({type(self.value)})")

        self.name = self.name.upper()


@dataclass
class Header:
    """Represent multiple FITS header cards.

    This does not support all FITS features.
    """

    values: list[HeaderItem]
    """The stored header items."""

    def add(self, item: HeaderItem) -> None:
        """Add the item (not checking for duplicates)."""
        self.values.append(item)

    def get(self, name: str) -> Optional[HeaderItem]:
        """Return the first occurrence of the key (case insensitive)."""

        uname = name.upper()
        for val in self.values:
            if uname == val.name.upper():
                return val

        return None

    def delete(self, name: str) -> Optional[HeaderItem]:
        """Remove the first occurrence of the key, if it exists (case insensitive)."""

        store = []
        out = None
        uname = name.upper()
        for val in self.values:
            if out is None and uname == val.name.upper():
                out = val
            else:
                store.append(val)

        if out is not None:
            self.values = store

        return out


@dataclass
class Column:
    """Represent a FITS column.

    This does not support all FITS features.
    """

    name: str
    """The column name (case insensitive)"""

    values: np.ndarray
    """The values for the column, as a ndarray.

    Variable-field arrays are represented as ndarrays with an object
    type.
    """

    desc: Optional[str] = None
    """The column description"""

    unit: Optional[str] = None
    """The units of the column"""

    minval: Optional[Union[int, float]] = None
    """The minimum value (corresponds to FITS TLMIN setting)."""

    maxval: Optional[Union[int, float]] = None
    """The maximum value (corresponds to FITS TLMAX setting)."""

    def __post_init__(self) -> None:

        if not isinstance(self.name, str):
            raise ValueError(f"Invalid column name: '{self.name}' ({type(self.name)})")

        if " " in self.name:
            raise ValueError(f"Invalid column name: '{self.name}'")

        self.name = self.name.upper()

        if not isinstance(self.values, np.ndarray):
            raise ValueError(f"Invalid column '{self.name}': values not ndarray")


# Represent a block of data: we have
#    - header only
#    - header + columns
#    - header + image
#
@dataclass
class Block:
    """Represent a block (header only)"""

    name: str
    """The name of the HDU (case insensitive)"""

    header: Header
    """The header values"""

    def __post_init__(self) -> None:

        if not isinstance(self.name, str):
            raise ValueError(f"Invalid HDU name: '{self.name}' ({type(self.name)})")

        if not isinstance(self.header, Header):
            raise ValueError(f"header is not set correctly for {self.name}")

        # HDU names can contain spaces, such as "SPECRESP MATRIX"
        self.name = self.name.strip().upper()


@dataclass
class TableBlock(Block):
    """Represent header and columns"""

    columns: list[Column]
    """The column data. This list must not be empty."""

    def __post_init__(self) -> None:

        super().__post_init__()
        if self.columns is None or len(self.columns) == 0:
            raise ValueError(f"Columns are missing or empty for {self.name}")

        for idx, col in enumerate(self.columns, 1):
            if not isinstance(col, Column):
                raise ValueError(f"Column {idx} is not a Column object in {self.name}")

    def get(self, colname: str) -> Optional[Column]:
        """Return the column (case insensitive) if it exists."""

        uname = colname.upper()
        for col in self.columns:
            if uname == col.name.upper():
                return col

        return None

    def rget(self, colname: str) -> Column:
        """Return a required column (case insensitive).

        Raise a ValueError if colname does not exist.
        """

        col = self.get(colname)
        if col is None:
            raise ValueError(f"column {colname} does not exist in {self.name}")

        return col


@dataclass
class SpectrumBlock(TableBlock):
    """Represent a PHA dataset.

    This ensures that the column CHANNEL exists and is 1 or 2D, and
    that either the COUNTS column exists, or RATE column and EXPOSURE
    keyword exists.

    """

    def __post_init__(self) -> None:

        super().__post_init__()
        chan = self.get("CHANNEL")
        if chan is None:
            raise ValueError(f"The PHA SPECTRUM block {self.name} "
                             "is missing the column: 'CHANNEL'")
        if chan.values.ndim not in [1, 2]:
            raise ValueError("Unable to handle CHANNEL shape: "
                             f"{chan.values.shape}")

        if self.get("COUNTS") is not None:
            return

        if self.get("RATE") is None:
            raise ValueError(f"The PHA SPECTRUM block {self.name} "
                             "is missing one of: 'COUNTS' or 'RATE'")

        if self.header.get("EXPOSURE") is None:
            raise ValueError(f"The PHA SPECTRUM block {self.name} "
                             "is missing the EXPOSURE keyword")


@dataclass
class SpecrespBlock(TableBlock):
    """Represent an ARF.

    This ensures that the columns SPECRESP, ENERG_LO, and ENERG_HI
    exist. It currently does not enforce any header settings.

    """

    def __post_init__(self) -> None:

        super().__post_init__()
        for name in ["SPECRESP", "ENERG_LO", "ENERG_HI"]:
            col = self.get(name)
            if col is None:
                raise ValueError(f"The ARF SPECRESP block {self.name} "
                                 "is missing the column: '{name}'")


@dataclass
class MatrixBlock(TableBlock):
    """Represent the MATRIX block of a RMF.

    This ensures that the columns ENERG_LO, ENERG_HI, N_GRP, F_CHAN,
    N_CHAN, and MATRIX exist. It currently does not enforce any header
    settings.

    """

    def __post_init__(self) -> None:

        super().__post_init__()
        for name in ["ENERG_LO", "ENERG_HI", "N_GRP", "F_CHAN",
                     "N_CHAN", "MATRIX"]:
            col = self.get(name)
            if col is None:
                raise ValueError(f"The RMF MATRIX block {self.name} is "
                                 f"missing the column: '{name}'")


@dataclass
class EboundsBlock(TableBlock):
    """Represent the EBOUNDS block of a RMF.

    This ensures that the columns CHANNEL, E_MIN, and E_MAX exist.  It
    currently does not enforce any header settings.

    """

    def __post_init__(self) -> None:

        super().__post_init__()
        for name in ["CHANNEL", "E_MIN", "E_MAX"]:
            col = self.get(name)
            if col is None:
                raise ValueError(f"The RMF EBOUNDS block {self.name} is "
                                 f"missing the column: '{name}'")


# The sky and eqpos field depends on whether the WCS code is
# available.
#
@dataclass
class ImageBlock(Block):
    """Represent header and image data."""

    image: np.ndarray
    """The image data."""

    sky: Optional[WCS] = None
    """The WCS for the physical/sky coordinate system."""

    eqpos: Optional[WCS] = None
    """The WCS for the WCS coordinate system."""


BlockType = Union[TableBlock, ImageBlock]


@dataclass
class BlockList:
    """Represent a set of blocks.

    This follows the FITS approach where the first block may be header
    only. It however does not require this.

    """

    blocks: list[BlockType]
    """The data for the blocks."""

    header: Optional[Header] = None
    """An optional header.

    If set this is used to create the first block with no data.
    """

    def __post_init__(self) -> None:

        # If there's no header then the first  block must be an image.
        #
        if self.header is not None:
            return

        if len(self.blocks) == 0 or \
           not isinstance(self.blocks[0], ImageBlock):
            raise ValueError("If header is empty the first block must be an image.")
