#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2011)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.sim import MCMC as _MCMC
from sherpa.sim import _samplers as __samplers
from sherpa.sim import _walkers as __walkers
from sherpa.astro.sim.pragbayes import *
from sherpa.astro.sim.fullbayes import *

_samplers = __samplers.copy()
_samplers.update(dict(pragbayes=PragBayes, fullbayes=FullBayes))

_walkers = __walkers.copy()
_walkers.update(dict(pragbayes=WalkWithSubIters, fullbayes=WalkWithSubIters))

class MCMC(_MCMC):
    """

    High-level UI to pyBLoCXS that joins the loop in 'Walk' with the jumping
    rule in 'Sampler'.  Implements a user interface for configuration.  This
    class implements a calc_stat() function using the Sherpa interface to 'Fit'.

    Overrides the base class with an astronomy specific class for incorporating
    uncertainties in the calibration information in Ancillary Response files
    (ARF).

    """
    __samplers = _samplers.copy()

    __walkers = _walkers.copy()
