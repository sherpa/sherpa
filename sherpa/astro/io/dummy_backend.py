#
#  Copyright (C) 2021, 2022
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
functions are used. It is just here to ensure that `sherpa.astro.io` can be
imported, even if no FITS reader is installed.
'''
import logging

__all__ = ('get_table_data', 'get_header_data', 'get_image_data',
           'get_column_data', 'get_ascii_data',
           'get_arf_data', 'get_rmf_data', 'get_pha_data',
           'set_table_data', 'set_image_data', 'set_pha_data')


lgr = logging.getLogger(__name__)

warning = lgr.warning
warning("""Cannot import usable I/O backend.
    If you are using CIAO, this is most likely an error and you should contact the CIAO helpdesk.
    If you are using Standalone Sherpa, please install astropy.""")


def get_table_data(*args, **kwargs):
    """A do-nothing operation"""
    raise NotImplementedError('No usable I/O backend was imported.')


get_header_data = get_table_data
get_image_data = get_table_data
get_column_data = get_table_data
get_ascii_data = get_table_data
get_arf_data = get_table_data
get_rmf_data = get_table_data
get_pha_data = get_table_data
set_table_data = get_table_data
set_image_data = get_table_data
set_pha_data = get_table_data
