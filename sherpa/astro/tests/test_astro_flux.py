#
#  Copyright (C) 2025
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

"""Low-level tests of sherpa.astro.flux"""

import numpy as np
import pytest

from sherpa.astro.data import DataPHA
from sherpa.astro.flux import sample_flux
from sherpa.astro.instrument import Response1D, create_arf, create_delta_rmf
from sherpa.data import Data1D
from sherpa.fit import Fit
from sherpa.models.basic import Polynom1D, Scale1D
from sherpa.models.parameter import Parameter
from sherpa.models.model import ArithmeticModel
from sherpa.stats import Chi2DataVar


def validate_fit_linked_par(fit: Fit) -> None:
    """Fit the data and check it is at the expected location."""

    fres = fit.fit()
    assert fres.succeeded
    assert fit.model.thawedpars == pytest.approx([17.076190476182614,
                                                  6.571428571477132])


def setup_data1d_linked_par(link: bool
                            ) -> tuple[Fit, ArithmeticModel, ArithmeticModel]:
    """Create simple setup for sample_flux testing

    Fit a straight line but using a "complex" expression so that
    we can have multiple components with free paramaters. The
    linked parameter, if selected, is the last free parameter, to
    make it possible to do direct comparisons between both versions.

    """

    x = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
    y = [20, 20, 25, 24, 24, 27]
    staterr = [2] * 6
    d = Data1D("test-link" if link else "test-unlinked", x, y,
               staterror=staterr)

    cpt1 = Scale1D("cpt1")
    cpt2 = Polynom1D("cpt2")
    cpt2.c0.freeze()
    cpt2.c1.thaw()

    if link:
        lmdl = Scale1D("linkedpar")
        # copy over the original starting value
        lmdl.c0 = 0
        cpt2.c1 = lmdl.c0

    mdl = cpt1 + cpt2

    assert len(mdl.get_thawed_pars()) == 2
    assert len(mdl.thawedpars) == 2

    tvals = [p.val for p in mdl.pars if not p.frozen]
    if link:
        assert len(tvals) == 1
    else:
        assert len(tvals) == 2

    f = Fit(d, mdl, stat=Chi2DataVar())
    validate_fit_linked_par(f)

    return f, mdl, cpt2


def setup_datapha_linked_par(link: bool
                             ) -> tuple[Fit, ArithmeticModel, ArithmeticModel]:
    """Create simple setup for sample_flux testing

    setup_data1d_linked_pars but with a "perfect" response
    included.

    """

    x = [1, 2, 3, 4, 5, 6]
    y = [20, 20, 25, 24, 24, 27]
    staterr = [2] * 6
    pha = DataPHA("test-link" if link else "test-unlinked", x, y,
                  staterror=staterr)
    pha.exposure = 1.0

    cpt1 = Scale1D("cpt1")
    cpt2 = Polynom1D("cpt2")
    cpt2.c0.freeze()
    cpt2.c1.thaw()

    if link:
        lmdl = Scale1D("linkedpar")
        # copy over the original starting value
        lmdl.c0 = 0
        cpt2.c1 = lmdl.c0

    mdl = cpt1 + cpt2

    # Ensure we evaluate this "un integrated". Note that
    # polynom1d uses the left edge if given an integrated
    # dataset and integrate is False.
    #
    # cpt1.integrate = False  scale1d is not integrated by default
    cpt2.integrate = False

    egrid = np.asarray([0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi)

    arf = create_arf(elo, ehi)

    pha.set_rmf(rmf)
    pha.set_arf(arf)

    rsp = Response1D(pha)
    full_mdl = rsp(mdl)

    assert len(full_mdl.get_thawed_pars()) == 2
    assert len(full_mdl.thawedpars) == 2

    tvals = [p.val for p in full_mdl.pars if not p.frozen]
    if link:
        assert len(tvals) == 1
    else:
        assert len(tvals) == 2

    f = Fit(pha, full_mdl, stat=Chi2DataVar())
    validate_fit_linked_par(f)

    # NOTE: do not apply the response to the returned models
    return f, mdl, cpt2


def setup_seed():
    """We want a consistent set of random numbers.

    It is not important if it is not the most-random set.
    """
    return np.random.RandomState(962355427)


# The same errors AND parameters are used for the _full and _subset
# runs, so the results should be the same.
#
SAMPLE_C0 = np.asarray([16.36343231,
                        19.61766136,
                        18.19143038,
                        18.4230697,
                        13.46291245])
SAMPLE_C1 = np.asarray([2.10661396,
                        5.30479342,
                        5.49025147,
                        8.85666698,
                        7.01016254])

# Ensure we can't accidentally change them
SAMPLE_C0.setflags(write=False)
SAMPLE_C1.setflags(write=False)


@pytest.mark.parametrize("lo,hi",
                         [(0, 40),
                          (None, None)
                          ])
@pytest.mark.parametrize("link", [False, True])
@pytest.mark.parametrize("setup,fluxes1,fluxes2",
                         [(setup_data1d_linked_par,
                           [0] * 5,
                           [7.44279111e-08, 9.85727388e-08, 9.37630615e-08, 1.06897220e-07, 8.11086789e-08]),
                          (setup_datapha_linked_par,
                           [1.67167473e-07, 2.21044981e-07, 2.10197208e-07, 2.39277125e-07, 1.81511770e-07],
                           [1.67167473e-07, 2.21044981e-07, 2.10197208e-07, 2.39277125e-07, 1.81511770e-07])
                          ])
def test_sf_full_linked(setup, fluxes1, fluxes2, link, lo, hi):
    """Check sample_flux with a basic dataset/model: Data1D/DataPHA

    Three things that are checked
      - behaviour with an upper-limit
      - checks some, but not all, cases of lo/hi set or unset
        (the behaviour when only one is given is "interesting" for
         historical reasons and is not tested here)
      - using a linked parameter gives the same results
        (this requires the free par order to match for the linked and
        unlinked cases)

    This uses the full model expression for the sample_flux
    call.

    """

    # This is known to fail with
    #   AssertionError: No mask has been found with limits None None
    # under these conditions. As they come from different parametrize
    # calls apply the xfail here.
    #
    if lo is None and setup == setup_data1d_linked_par:
        pytest.xfail("test conditions known to fail")

    f, mdl, _ = setup(link)
    rng = setup_seed()

    res = sample_flux(f, f.data, mdl, lo=lo, hi=hi, num=5, rng=rng,
                      clip="hard")

    # Check the sampled values.
    #
    assert res[:, 1] == pytest.approx(SAMPLE_C0)
    assert res[:, 2] == pytest.approx(SAMPLE_C1)

    # There's no clipping.
    #
    assert res[:, 3] == pytest.approx([0, 0, 0, 0, 0])

    # These fluxes are not meaningful, but can be checked. We
    # could work out the expected values based on the model and
    # parameters, but just treat as a regression test.
    #
    # The behaviour when lo/hi is None is different for the
    # Data1D and DataPHA cases, hence sending in fluxes1 and
    # fluxes2.
    #
    if lo is None:
        assert hi is None  # if this changes need to change the test
        assert res[:, 0] == pytest.approx(fluxes1)
    else:
        assert res[:, 0] == pytest.approx(fluxes2)


@pytest.mark.parametrize("lo,hi",
                         [(0, 40),
                          (None, None)
                          ])
@pytest.mark.parametrize("link", [False, True])
@pytest.mark.parametrize("setup,fluxes1,fluxes2",
                         [(setup_data1d_linked_par,
                           [0] * 5,
                           [1.15068538e-08, 2.31384429e-08, 2.38129428e-08, 3.60563965e-08, 2.93407693e-08]),
                          (setup_datapha_linked_par,
                           [2.55950939e-08, 5.13178152e-08, 5.28094405e-08, 7.98852720e-08, 6.50339731e-08],
                           [2.55950939e-08, 5.13178152e-08, 5.28094405e-08, 7.98852720e-08, 6.50339731e-08])
                          ])
def test_sf_subset_linked(setup, fluxes1, fluxes2, link, lo, hi):
    """Check sample_flux with a basic dataset/model: Data1D/DataPHA

    Three things that are checked
      - behaviour with an upper-limit
      - checks some, but not all, cases of lo/hi set or unset
        (the behaviour wheh only one is given is "interesting" for
         historical reasons and is not tested here)
      - using a linked parameter gives the same results
        (this requires the free par order to match for the linked and
        unlinked cases)

    This uses a subset of the model expression for the sample_flux
    call.

    """

    # This is known to fail with
    #   AssertionError: No mask has been found with limits None None
    # under these conditions. As they come from different parametrize
    # calls apply the xfail here.
    #
    if lo is None and setup == setup_data1d_linked_par:
        pytest.xfail("test conditions known to fail")

    f, _, mdl = setup(link)
    rng = setup_seed()

    res = sample_flux(f, f.data, mdl, lo=lo, hi=hi, num=5, rng=rng,
                      clip="hard")

    # Check the sampled values.
    #
    assert res[:, 1] == pytest.approx(SAMPLE_C0)
    assert res[:, 2] == pytest.approx(SAMPLE_C1)

    # There's no clipping.
    #
    assert res[:, 3] == pytest.approx([0, 0, 0, 0, 0])

    # Flux checks.
    #
    if lo is None:
        assert hi is None  # if this changes need to change the test
        assert res[:, 0] == pytest.approx(fluxes1)
    else:
        assert res[:, 0] == pytest.approx(fluxes2)
