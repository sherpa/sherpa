#
# Copyright (C) 2015, 2016, 2019, 2020, 2021
#               Smithsonian Astrophysical Observatory
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

from importlib import import_module

from sherpa.plot import plotter
from sherpa.utils.logging import config_logger

logger = config_logger(__name__)

name = plotter.name

backend_map = {
    'pylab': 'plot_matplotlib',
    'dummy': 'plot_dummy'
}


def _update_globals(module):
    globals().update((k, v)
                     for k, v in module.__dict__.items() if k not in globals())


try:
    backend_module = import_module("." + backend_map[name], __name__)
except ImportError:
    backend_module = import_module(".plot_dummy", __name__)

_update_globals(backend_module)
