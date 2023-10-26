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

from ..data import Data1D


__all__ = ('get_table_data', 'get_header_data', 'get_image_data',
           'get_column_data', 'get_ascii_data',
           'get_arf_data', 'get_rmf_data', 'get_pha_data',
           'set_table_data', 'set_image_data', 'set_pha_data',
           'set_arf_data', 'set_rmf_data', 'set_hdus')


lgr = logging.getLogger(__name__)

warning = lgr.warning
warning("""Cannot import usable I/O backend.
    If you are using CIAO, this is most likely an error and you should contact the CIAO helpdesk.
    If you are using Standalone Sherpa, please install astropy.""")


def get_table_data(arg, ncols=1, colkeys=None, make_copy=True, fix_type=True,
                   blockname=None, hdrkeys=None):
    """Read colums."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_header_data(arg, blockname=None, hdrkeys=None):
    """Read the metadata."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_image_data(arg, make_copy=True, fix_type=True):
    """Read image data."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_column_data(*args):
    """Extract the column data."""
    raise NotImplementedError('No usable I/O backend was imported.')


# Follow sherpa.io.get_ascii_data API.
#
def get_ascii_data(filename, ncols=2, colkeys=None, sep=' ',
                   dstype=Data1D, comment='#', require_floats=True):
    """Read columns from an ASCII file"""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_arf_data(arg, make_copy=False):
    """Read in the ARF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_rmf_data(arg, make_copy=False):
    """Read in the RMF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def get_pha_data(arg, make_copy=False, use_background=False):
    """Read in the PHA."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_table_data(filename, data, col_names, header=None,
                   ascii=False, clobber=False, packup=False):
    """Create the tabular data."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_image_data(filename, data, header, ascii=False, clobber=False,
                   packup=False):
    """Create the image data."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_pha_data(filename, data, col_names, header=None,
                 ascii=False, clobber=False, packup=False):
    """Create the PHA."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_arf_data(filename, data, col_names, header=None,
                 ascii=False, clobber=False, packup=False):
    """Create the ARF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_rmf_data(filename, blocks, clobber=False):
    """Create the ARF."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_hdus(filename, hdulist, clobber=False):
    """Write out (possibly multiple) blocks."""
    raise NotImplementedError('No usable I/O backend was imported.')


def read_table_blocks(arg, make_copy=False):
    """Read in tabular data with no restrictions on the columns."""
    raise NotImplementedError('No usable I/O backend was imported.')


def set_arrays(filename, args, fields=None, ascii=True, clobber=False):
    """Write out columns."""
    raise NotImplementedError('No usable I/O backend was imported.')
