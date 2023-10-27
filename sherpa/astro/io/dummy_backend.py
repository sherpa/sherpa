#
#  Copyright (C) 2021, 2022, 2023
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
from typing import Any, Optional, Sequence, Union

import numpy as np

from ..data import Data1D


__all__ = ('get_table_data', 'get_header_data', 'get_image_data',
           'get_column_data', 'get_ascii_data', 'get_arf_data',
           'get_rmf_data', 'get_pha_data',
           #
           'pack_table_data', 'pack_image_data', 'pack_pha_data',
           'pack_arf_data', 'pack_rmf_data', 'pack_hdus',
           #
           'set_table_data', 'set_image_data', 'set_pha_data',
           'set_arf_data', 'set_rmf_data', 'set_hdus')


lgr = logging.getLogger(__name__)

warning = lgr.warning
warning("""Cannot import usable I/O backend.
    If you are using CIAO, this is most likely an error and you should contact the CIAO helpdesk.
    If you are using Standalone Sherpa, please install astropy.""")


KeyType = Union[str, bool, int, float]
NamesType = Sequence[str]
HdrType = dict[str, KeyType]


def get_table_data(arg,
                   ncols: int = 1,
                   colkeys: Optional[NamesType] = None,
                   make_copy: bool = True,
                   fix_type: bool = True,
                   blockname: Optional[str] = None,
                   hdrkeys: Optional[NamesType] = None
                   ) -> tuple[list[str], list[np.ndarray], str, HdrType]:
    """Read columns."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_header_data(arg,
                    blockname: Optional[str] = None,
                    hdrkeys: Optional[NamesType] = None
                    ) -> dict[str, KeyType]:
    """Read the metadata."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_image_data(arg,
                   make_copy: bool = True,
                   fix_type: bool = True
                   ) -> tuple[dict[str, Any], str]:
    """Read image data."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_column_data(*args) -> list[np.ndarray]:
    """Extract the column data."""
    raise NotImplementedError('No usable I/O backend was imported.')


# Follow sherpa.io.get_ascii_data API.
#
def get_ascii_data(filename: str,
                   ncols: int = 2,
                   colkeys: Optional[NamesType] = None,
                   sep: str = ' ',
                   dstype: type = Data1D,
                   comment: str = '#',
                   require_floats: bool = True
                   ) -> tuple[list[str], list[np.ndarray], str]:
    """Read columns from an ASCII file"""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_arf_data(arg,
                 make_copy: bool = False
                 ) -> tuple[dict[str, Any], str]:
    """Read in the ARF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_rmf_data(arg,
                 make_copy: bool = False
                 ) -> tuple[dict[str, Any], str]:
    """Read in the RMF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_pha_data(arg,
                 make_copy: bool = False,
                 use_background: bool = False
                 ) -> tuple[list[dict[str, Any]], str]:
    """Read in the PHA."""
    raise NotImplementedError('No usable I/O backend was imported.')


def pack_table_data(data, col_names, header=None) -> Any:
    """Create the tabular data."""
    raise NotImplementedError('No usable I/O backend was imported.')


def pack_image_data(data, header, ascii=False) -> Any:
    """Create the image data."""
    raise NotImplementedError('No usable I/O backend was imported.')


def pack_pha_data(data, col_names, header=None) -> Any:
    """Create the PHA."""
    raise NotImplementedError('No usable I/O backend was imported.')


def pack_arf_data(data, col_names, header=None) -> Any:
    """Create the ARF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def pack_rmf_data(blocks) -> Any:
    """Create the ARF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def pack_hdus(blocks) -> Any:
    """Create a dataset."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_table_data(filename: str,
                   data, col_names, header=None,
                   ascii: bool = False,
                   clobber: bool = False) -> None:
    """Write out the tabular data."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_image_data(filename: str,
                   data, header,
                   ascii: bool = False,
                   clobber: bool = False) -> None:
    """Write out the image data."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_pha_data(filename: str,
                 data, col_names, header=None,
                 ascii: bool = False,
                 clobber: bool = False) -> None:
    """Write out the PHA."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_arf_data(filename: str,
                 data, col_names, header=None,
                 ascii: bool = False,
                 clobber: bool = False) -> None:
    """Write out the ARF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_rmf_data(filename: str,
                 blocks,
                 clobber: bool = False) -> None:
    """Write out the RMF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_hdus(filename: str,
             blocks,
             clobber: bool = False) -> None:
    """Write out (possibly multiple) blocks."""
    raise NotImplementedError('No usable I/O backend was imported.')


def read_table_blocks(arg,
                      make_copy: bool = False
                      ) -> tuple[str,
                                 dict[int, dict[str, np.ndarray]],
                                 dict[int, HdrType]]:
    """Read in tabular data with no restrictions on the columns."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_arrays(filename: str,
               args: Sequence[np.ndarray],
               fields: Optional[NamesType] = None,
               ascii: bool = True,
               clobber: bool = False) -> None:
    """Write out columns."""
    raise NotImplementedError('No usable I/O backend was imported.')
