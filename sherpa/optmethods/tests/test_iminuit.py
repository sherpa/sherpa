#
#  Copyright (C) 2022  Smithsonian Astrophysical Observatory
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

import pytest
import numpy

from sherpa.utils.testing import requires_iminuit
from sherpa.optmethods import _tstoptfct

try:
    from sherpa.optmethods.myminuit import myminuit
except ImportError:
    pass

def init(name, npar):
    x0, _, _, fmin = _tstoptfct.init(name, npar)
    # make sure minuit does not transform
    xmin = npar * [-numpy.inf]
    xmax = npar * [numpy.inf]
    return x0, xmin, xmax, fmin

def tst_opt(fct, migrad, npar, reltol=1.0e-3, abstol=1.0e-3):
    """The central function for all the optimization test
    1) Out of the 35 tests from:
    J. MORE', B. GARBOW & K. HILLSTROM,
    "Algorithm 566: Fortran Subroutines for Testing Unconstrained
    Optimization Software.", ACM TOMS, VOL. 7, PAGES 14-41 AND 136-140, 1981
    xfail: lmdif(6), minim(5), neldermead(3), moncar(2)
    2) The remaining random 32 'global' func tests:
    xfail: minim(14), montecarlo(0), neldermead(10)
    """
    x0, xmin, xmax, fmin = init(fct.__name__, npar)
    status, x, fval, msg, xtra = myminuit(fct, x0, xmin, xmax, migrad=migrad)
    assert fval == pytest.approx(fmin, rel=reltol, abs=abstol)


###############################################################################
@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(True, marks=pytest.mark.xfail)])
def test_rosenbrock(migrad, npar=4):
    tst_opt(_tstoptfct.rosenbrock, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_freudensteinroth(migrad, npar=4):
    tst_opt(_tstoptfct.freudenstein_roth, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_powell_badly_scaled(migrad, npar=2):
    tst_opt(_tstoptfct.powell_badly_scaled, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_brown_badly_scaled(migrad, npar=2):
    tst_opt(_tstoptfct.brown_badly_scaled, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_beale(migrad, npar=2):
    tst_opt(_tstoptfct.beale, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_jennrich_sampson(migrad, npar=2):
    tst_opt(_tstoptfct.jennrich_sampson, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_helical_valley(migrad, npar=3):
    tst_opt(_tstoptfct.helical_valley, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_bard(migrad, npar=3):
    tst_opt(_tstoptfct.bard, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, False])
def test_gaussian(migrad, npar=3):
    tst_opt(_tstoptfct.gaussian, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_meyer(migrad, npar=3):
    tst_opt(_tstoptfct.meyer, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_gulf_research_development(migrad, npar=3):
    tst_opt(_tstoptfct.gulf_research_development, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_box3d(migrad, npar=3):
    tst_opt(_tstoptfct.box3d, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_powell_singular(migrad, npar=4):
    tst_opt(_tstoptfct.powell_singular, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_wood(migrad, npar=4):
    tst_opt(_tstoptfct.wood, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_kowalik_osborne(migrad, npar=4):
    tst_opt(_tstoptfct.kowalik_osborne, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_brown_dennis(migrad, npar=4):
    tst_opt(_tstoptfct.brown_dennis, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_osborne1(migrad, npar=5):
    tst_opt(_tstoptfct.osborne1, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_biggs(migrad, npar=6):
    tst_opt(_tstoptfct.biggs, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_osborne2(migrad, npar=11):
    tst_opt(_tstoptfct.osborne2, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_watson(migrad, npar=6):
    tst_opt(_tstoptfct.watson, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_extended_rosenbrock(migrad, npar=12):
    tst_opt(_tstoptfct.rosenbrock, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_extended_powell_singular(migrad, npar=8):
    tst_opt(_tstoptfct.powell_singular, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_penaltyI(migrad, npar=4):
    tst_opt(_tstoptfct.penaltyI, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_penaltyII(migrad, npar=4):
    tst_opt(_tstoptfct.penaltyII, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_variably_dimensioned(migrad, npar=6):
    tst_opt(_tstoptfct.variably_dimensioned, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_trigonometric(migrad, npar=4):
    tst_opt(_tstoptfct.trigonometric, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_brown_almost_linear(migrad, npar=3):
    tst_opt(_tstoptfct.brown_almost_linear, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_discrete_boundary(migrad, npar=4):
    tst_opt(_tstoptfct.discrete_boundary, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_discrete_integral(migrad, npar=4):
    tst_opt(_tstoptfct.discrete_integral, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_broyden_tridiagonal(migrad, npar=4):
    tst_opt(_tstoptfct.broyden_tridiagonal, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_broyden_banded(migrad, npar=4):
    tst_opt(_tstoptfct.broyden_banded, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_linear_fullrank(migrad, npar=4):
    tst_opt(_tstoptfct.linear_fullrank, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_linear_fullrank1(migrad, npar=4):
    tst_opt(_tstoptfct.linear_fullrank1, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_linear_fullrank0cols0rows0(migrad, npar=4):
    tst_opt(_tstoptfct.linear_fullrank0cols0rows, migrad, npar)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [True, pytest.param(False, marks=pytest.mark.xfail)])
def test_linear_chebyquad(migrad, npar=9):
    tst_opt(_tstoptfct.chebyquad, migrad, npar)

# ###############################################################################


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Ackley(migrad, npar=4):
    tst_opt(_tstoptfct.Ackley, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Booth(migrad, npar=6):
    tst_opt(_tstoptfct.Booth, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Bohachevsky1(migrad, npar=2):
    tst_opt(_tstoptfct.Bohachevsky1, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Bohachevsky2(migrad, npar=2):
    tst_opt(_tstoptfct.Bohachevsky2, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Bohachevsky3(migrad, npar=2):
    tst_opt(_tstoptfct.Bohachevsky3, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Branin(migrad, npar=2):
    tst_opt(_tstoptfct.Branin, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Branin2(migrad, npar=2):
    tst_opt(_tstoptfct.Branin2, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Chichinadze(migrad, npar=2):
    tst_opt(_tstoptfct.Chichinadze, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Cola(migrad, npar=17):
    tst_opt(_tstoptfct.Cola, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Colville(migrad, npar=4):
    tst_opt(_tstoptfct.Colville, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_dcs(migrad, npar=5):
    tst_opt(_tstoptfct.dcs, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_decanom(migrad, npar=2):
    tst_opt(_tstoptfct.decanom, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_dodecal(migrad, npar=3):
    tst_opt(_tstoptfct.dodecal, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_DixonPrice(migrad, npar=5):
    tst_opt(_tstoptfct.DixonPrice, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Easom(migrad, npar=5):
    tst_opt(_tstoptfct.Easom, migrad, npar, False)

@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_factor(migrad, npar=5):
    tst_opt(_tstoptfct.factor, migrad, npar, False)

@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Func1(migrad, npar=2):
    tst_opt(_tstoptfct.Func1, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Griewank(migrad, npar=2):
    tst_opt(_tstoptfct.Griewank, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Hansen(migrad, npar=2):
    tst_opt(_tstoptfct.Hansen, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Hartman6(migrad, npar=6):
    tst_opt(_tstoptfct.Hartman6, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Himmelblau(migrad, npar=2):
    tst_opt(_tstoptfct.Himmelblau, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Holzman1(migrad, npar=3):
    tst_opt(_tstoptfct.Holzman1, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Holzman2(migrad, npar=3):
    tst_opt(_tstoptfct.Holzman2, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Judge(migrad, npar=3):
    tst_opt(_tstoptfct.Judge, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Levy(migrad, npar=4):
    tst_opt(_tstoptfct.Levy, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_McCormick(migrad, npar=2):
    tst_opt(_tstoptfct.McCormick, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_McKinnon(migrad, npar=2):
    tst_opt(_tstoptfct.McKinnon, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Michalewicz(migrad, npar=2):
    tst_opt(_tstoptfct.Michalewicz, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Paviani(migrad, npar=10):
    tst_opt(_tstoptfct.Paviani, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Rastrigin(migrad, npar=4):
    tst_opt(_tstoptfct.Rastrigin, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_seqp(migrad, npar=2):
    tst_opt(_tstoptfct.seqp, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Shekel5(migrad, npar=4):
    tst_opt(_tstoptfct.Shekel5, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Shekel7(migrad, npar=4):
    tst_opt(_tstoptfct.Shekel7, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Shekel10(migrad, npar=4):
    tst_opt(_tstoptfct.Shekel10, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_ShekelModified(migrad, npar=2):
    tst_opt(_tstoptfct.ShekelModified, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Shubert(migrad, npar=2):
    tst_opt(_tstoptfct.Shubert, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_SixHumpCamel(migrad, npar=2):
    tst_opt(_tstoptfct.SixHumpCamel, migrad, npar, False)

@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Trecanni(migrad, npar=2):
    tst_opt(_tstoptfct.Trecanni, migrad, npar, False)


@requires_iminuit
@pytest.mark.parametrize("migrad",
                         [pytest.param(True, marks=pytest.mark.xfail),
                          pytest.param(False, marks=pytest.mark.xfail)])
def test_Trefethen4(migrad, npar=2):
    tst_opt(_tstoptfct.Trefethen4, migrad, npar, False)
