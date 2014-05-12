#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
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
