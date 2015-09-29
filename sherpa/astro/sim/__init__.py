#
#  Copyright (C) 2011, 2015  Smithsonian Astrophysical Observatory
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

"""Include calibration uncertainties in Bayesian Low-Count X-ray Spectral (BLoCXS) analysis.

This module adds the "pragbayes" and "fullbayes" jumping rules for use with
the pyBLoCXS MCMC code. These rules support calibration uncertainties in the
ARF of PHA files, as discussed in [1]_:

- ``PragBayes`` is used when effective area calibration uncertainty is to be
  included in the calculation. In this case, at eacn nominal MCMC iteration, a
  new calibration product is generated, and a series of ``N`` steps (which can
  be changed using ``set_sampler_opt``) are carried out with this calibration
  product, using the ``MetropolisMH`` rule. When the steps are finished, the
  *last* iteration is stored in the chain and the process repeated. This
  can *only* be used with PHA data sets, and requires a special formulation
  of the ARF to include calibration uncertainties.

- ``FullBayes``, is an experimental rule that extends the behavior of
  the ``PragBayes`` rule.

References
----------

.. [1] "Accounting for Calibration Uncertainties in X-ray Analysis:
       Effective Areas in Spectral Fitting", Lee et al., 2011, ApJ, 731, 126
       http://adsabs.harvard.edu/abs/2011ApJ...731..126L

"""

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
