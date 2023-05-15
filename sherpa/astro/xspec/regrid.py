#
#  Copyright (C) 2020, 2023
#  Smithsonian Astrophysical Observatory
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

"""Regrid support for XSPEC models.

We need to handle multiplicative models specially, since
they have integrate=False but are sent the low and high
edges. One solution would be to just treat the models
as evaluating the data at the left-edge of each bin and
then regrid that, but tests did not seem to work well
(but should be reviewed).

"""

import numpy as np

from sherpa.models.regrid import ModelDomainRegridder1D
from sherpa.utils.err import ModelErr


__all__ = ('XSMultiplicativeRegridder', )


class XSMultiplicativeRegridder(ModelDomainRegridder1D):
    """XSPEC Multiplicative regridder

    XSPEC multiplicative models are always treated as non-integrated,
    in that we use interpolation, not re-binning, when regridding.
    """

    def _evaluate(self, data_space, pars, modelfunc, **kwargs):
        """Evaluate the grid, treating as non-integrated for interpolation."""

        eval_space = self.evaluation_space
        if not eval_space.is_integrated:
            raise ModelErr('needsint')

        # Note that unlike the standard evaluation for interpolation by regrid
        # we do not join the grids (since it's hard to do cleanly here as we
        # have bins not points). We do have to handle the case where the
        # evaluation falls outside the evaluation_space.
        #
        # We also treat each bin being at its central energy, which is not
        # guaranteed for any particular model.
        #
        elo, ehi = eval_space.grid[0], eval_space.grid[1]
        y = modelfunc(pars, elo, ehi, **kwargs)

        emid = eval_space.midpoint_grid
        dmid = data_space.midpoint_grid

        if eval_space.start <= data_space.start and eval_space.end >= data_space.end:
            return self.method(dmid, emid, y)

        # Should this check the lo/hi bins instead of the mid-points?
        #
        yout = np.zeros(dmid.size, dtype=y.dtype)
        idx = np.where((dmid >= eval_space.start) &
                       (dmid <= eval_space.end))
        yout[idx] = self.method(dmid[idx], emid, y)
        return yout
