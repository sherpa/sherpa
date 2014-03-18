#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
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
