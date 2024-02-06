#
#  Copyright (C) 2021 - 2024
#  MIT
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
'''A dummy backend for I/O.

This backend provides no functionality and raises an error if any of its
functions are used. It is provided as a model for what is needed in a
backend, even if it does nothing, and to allow `sherpa.astro.io` to be
imported even if no usable backend is available.

'''

import logging
from typing import Any, Mapping, Optional, Sequence, Union

import numpy as np

from ..data import Data1D

from .backends import BaseBackend
from .xstable import TableHDU


__all__ = ('Backend', )


warning = logging.getLogger(__name__).warning
warning("""Cannot import usable I/O backend.
    If you are using CIAO, this is most likely an error and you should contact the CIAO helpdesk.
    If you are using Standalone Sherpa, please install astropy.""")


# Some variants are named xxx and xxxArg, where the former is the
# return value (an invariant type, like dict) and the latter is
# covariant (such as Mapping) as it's used as an argument to a
# function.
#
KeyType = Union[str, bool, int, float]
NamesType = Sequence[str]
HdrTypeArg = Mapping[str, KeyType]
HdrType = dict[str, KeyType]
ColumnsType = Mapping[str, Union[np.ndarray, list, tuple]]

# It's hard to type the arguments to the Data constructors
DataTypeArg = Mapping[str, Any]
DataType = dict[str, Any]


class Backend(BaseBackend):
    """A dummy backend for I/O.

    If used, any attempt to call an I/O routine (read or write) will
    raise a NotImplementedError exception.

    """

    name  = "dummy"

    def get_table_data(self, arg,
                       ncols: int = 1,
                       colkeys: Optional[NamesType] = None,
                       make_copy: bool = True,
                       fix_type: bool = True,
                       blockname: Optional[str] = None,
                       hdrkeys: Optional[NamesType] = None
                       ) -> tuple[list[str], list[np.ndarray], str, HdrType]:
        """Read columns from a file or object."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def get_header_data(self, arg,
                        blockname: Optional[str] = None,
                        hdrkeys: Optional[NamesType] = None
                        ) -> HdrType:
        """Read the metadata."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def get_image_data(self, arg,
                       make_copy: bool = True,
                       fix_type: bool = True
                       ) -> tuple[DataType, str]:
        """Read image data."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def get_column_data(self, *args) -> list[np.ndarray]:
        """Extract the column data."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def get_ascii_data(self,
                       filename: str,
                       ncols: int = 1,
                       colkeys: Optional[NamesType] = None,
                       sep: str = ' ',
                       dstype: type = Data1D,
                       comment: str = '#',
                       require_floats: bool = True
                       ) -> tuple[list[str], list[np.ndarray], str]:
        """Read columns from an ASCII file."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def get_arf_data(self, arg,
                     make_copy: bool = True
                     ) -> tuple[DataType, str]:
        """Read in the ARF."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def get_rmf_data(self, arg,
                     make_copy: bool = True
                     ) -> tuple[DataType, str]:
        """Read in the RMF."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def get_pha_data(self, arg,
                     make_copy: bool = True,
                     use_background: bool = False
                     ) -> tuple[list[DataType], str]:
        """Read in the PHA."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def pack_table_data(self,
                        data: ColumnsType,
                        col_names: NamesType,
                        header: Optional[HdrTypeArg] = None) -> Any:
        """Create the tabular data."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def pack_image_data(self,
                        data: DataTypeArg,
                        header: HdrTypeArg) -> Any:
        """Create the image data."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def pack_pha_data(self,
                      data: ColumnsType,
                      col_names: NamesType,
                      header: Optional[HdrTypeArg] = None) -> Any:
        """Create the PHA."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def pack_arf_data(self,
                      data: ColumnsType,
                      col_names: NamesType,
                      header: Optional[HdrTypeArg] = None) -> Any:
        """Create the ARF."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def pack_rmf_data(self, blocks) -> Any:
        """Create the RMF."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def pack_hdus(self,
                  blocks: Sequence[TableHDU]) -> Any:
        """Create a dataset."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def set_table_data(self,
                       filename: str,
                       data: ColumnsType,
                       col_names: NamesType,
                       header: Optional[HdrTypeArg] = None,
                       ascii: bool = False,
                       clobber: bool = False) -> None:
        """Write out the tabular data."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def set_image_data(self,
                       filename: str,
                       data: DataTypeArg,
                       header: HdrTypeArg,
                       ascii: bool = False,
                       clobber: bool = False) -> None:
        """Write out the image data."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def set_pha_data(self,
                     filename: str,
                     data: ColumnsType,
                     col_names: NamesType,
                     header: Optional[HdrTypeArg] = None,
                     ascii: bool = False,
                     clobber: bool = False) -> None:
        """Write out the PHA."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def set_arf_data(self,
                     filename: str,
                     data: ColumnsType,
                     col_names: NamesType,
                     header: Optional[HdrTypeArg] = None,
                     ascii: bool = False,
                     clobber: bool = False) -> None:
        """Write out the ARF."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def set_rmf_data(self,
                     filename: str,
                     blocks,
                     clobber: bool = False) -> None:
        """Write out the RMF."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def set_hdus(self,
                 filename: str,
                 blocks: Sequence[TableHDU],
                 clobber: bool = False) -> None:
        """Write out (possibly multiple) blocks."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def read_table_blocks(self, arg,
                          make_copy: bool = False
                          ) -> tuple[str,
                                     dict[int, dict[str, np.ndarray]],
                                     dict[int, HdrType]]:
        """Read in tabular data with no restrictions on the columns."""
        raise NotImplementedError('No usable I/O backend was imported.')


    def write_arrays(self,
                     filename: str,
                     args: Sequence[np.ndarray],
                     fields: Optional[NamesType] = None,
                     ascii: bool = True) -> None:
        """Write out columns."""
        raise NotImplementedError('No usable I/O backend was imported.')
