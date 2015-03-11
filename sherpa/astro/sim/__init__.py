# 
#  Copyright (C) 2011  Smithsonian Astrophysical Observatory
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
