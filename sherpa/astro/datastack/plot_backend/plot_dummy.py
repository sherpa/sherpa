#
#  Copyright (C) 2010, 2014, 2015  Smithsonian Astrophysical Observatory
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

"""
Plotting routines for the data stack module. These are dummy
routines that do nothing, and are for when no supported
plotting module is available.
"""
from sherpa.utils.logging import config_logger

logger = config_logger(__name__)

name = "dummy_backend"


def dummy(*args, **kwargs):
    logger.warning("using dummy plotting backend")

initialize_backend = dummy
initialize_plot = dummy
select_plot = dummy
save_plot = dummy
plot_savefig = dummy
plot_xlabel = dummy
plot_ylabel = dummy
plot_title = dummy
plot_xlim = dummy
plot_ylim = dummy
plot_set_xscale = dummy
plot_set_yscale = dummy
