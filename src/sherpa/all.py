# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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

Import all standard Sherpa subpackages

This is a convenience module that allows one to import all of the
standard Sherpa subpackages with a single 'import sherpa.all'
statement.  Note that discipline-specific subpackages
(e.g. sherpa.astro) are not imported by this module.

"""


import logging
warning = logging.getLogger(__name__).warning

from sherpa.data import *
from sherpa.estmethods import *
from sherpa.fit import *

try:
    from sherpa.image import *
except ImportError:
    warning('failed to import sherpa.image; imaging routines will not be ' +
            'available')

from sherpa.io import *
from sherpa.models import *
from sherpa.optmethods import *

try:
    from sherpa.plot import *
except ImportError:
    warning('failed to import sherpa.plot; plotting routines will not be ' +
            'available')

from sherpa.stats import *
from sherpa.instrument import *
from sherpa.utils import *

from sherpa.sim import *

# We don't want these listed as part of the module
del logging, warning
