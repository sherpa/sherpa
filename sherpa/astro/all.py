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

import logging
warning = logging.getLogger(__name__).warning

from sherpa.astro.background import *
from sherpa.astro.data import *
from sherpa.astro.instrument import *
from sherpa.astro.flux import *

try:
    from sherpa.astro.io import *
except ImportError:
    warning('failed to import sherpa.astro.io; FITS I/O routines will not ' +
            'be available')

from sherpa.astro.models import *
from sherpa.astro.optical import *

try:
    from sherpa.astro.plot import *
except ImportError:
    warning('failed to import sherpa.astro.plot; astronomical plotting ' +
            'routines will not be available')

from sherpa.astro.utils import *
from sherpa.astro.sim import *

try:
    from sherpa.astro.xspec import *
except ImportError:
    warning('failed to import sherpa.astro.xspec; XSPEC models will not be ' +
            'available')

# We don't want these listed as part of the module
del logging, warning
