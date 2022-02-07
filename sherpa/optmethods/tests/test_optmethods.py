#
#  Copyright (C) 2007, 2015, 2016, 2018, 2019, 2020, 2021, 2022
#     Smithsonian Astrophysical Observatory
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

from sherpa.optmethods import _tstoptfct
from sherpa.optmethods.optfcts import lmdif, minim, montecarlo, neldermead, \
    grid_search
from sherpa.utils import _ncpus


def init(name, npar):
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    return x0, xmin, xmax, fmin


def tst_opt(opt, fct, npar, reltol=1.0e-3, abstol=1.0e-3):
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
    status, x, fval, msg, xtra = opt(fct, x0, xmin, xmax)
    assert fmin == pytest.approx(fval, rel=reltol, abs=abstol)
    if opt == lmdif and _ncpus > 1:
        status, x, fval, msg, xtra = opt(fct, x0, xmin, xmax, numcores=_ncpus)
        assert fmin == pytest.approx(fval, rel=reltol, abs=abstol)
        assert xtra.get('num_parallel_map') != 0


###############################################################################
def test_gridsearch_maxfev():
    fct = _tstoptfct.rosenbrock
    x0, xmin, xmax, fmin = init(fct.__name__, 2)
    result = grid_search(fct, x0, xmin, xmax, maxfev=4)
    assert result[4]['nfev'] == 4

def test_gridsearch_maxfev():
    fct = _tstoptfct.rosenbrock
    x0, xmin, xmax, fmin = init(fct.__name__, 2)
    result = grid_search(fct, x0, xmin, xmax, num=2, method='nEldErmEad')
    assert result[4]['nfev'] == 313

@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_rosenbrock(opt, npar=4):
    tst_opt(opt, _tstoptfct.rosenbrock, npar)

@pytest.mark.parametrize("opt", [pytest.param(lmdif, marks=pytest.mark.xfail),
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 pytest.param(montecarlo,
                                              marks=pytest.mark.xfail),
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_freudensteinroth(opt, npar=4):
    tst_opt(opt, _tstoptfct.freudenstein_roth, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_powell_badly_scaled(opt, npar=2):
    tst_opt(opt, _tstoptfct.powell_badly_scaled, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_brown_badly_scaled(opt, npar=2):
    tst_opt(opt, _tstoptfct.brown_badly_scaled, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_beale(opt, npar=2):
    tst_opt(opt, _tstoptfct.beale, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_jennrich_sampson(opt, npar=2):
    tst_opt(opt, _tstoptfct.jennrich_sampson, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_helical_valley(opt, npar=3):
    tst_opt(opt, _tstoptfct.helical_valley, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_bard(opt, npar=3):
    tst_opt(opt, _tstoptfct.bard, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_gaussian(opt, npar=3):
    tst_opt(opt, _tstoptfct.gaussian, npar)


@pytest.mark.parametrize("opt", [pytest.param(lmdif, marks=pytest.mark.xfail),
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_meyer(opt, npar=3):
    tst_opt(opt, _tstoptfct.meyer, npar)


@pytest.mark.slow
@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_gulf_research_development(opt, npar=3):
    tst_opt(opt, _tstoptfct.gulf_research_development, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_box3d(opt, npar=3):
    tst_opt(opt, _tstoptfct.box3d, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_powell_singular(opt, npar=4):
    tst_opt(opt, _tstoptfct.powell_singular, npar)


@pytest.mark.parametrize("opt", [pytest.param(lmdif, marks=pytest.mark.xfail),
                                 minim, montecarlo, neldermead])
def test_wood(opt, npar=4):
    tst_opt(opt, _tstoptfct.wood, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_kowalik_osborne(opt, npar=4):
    tst_opt(opt, _tstoptfct.kowalik_osborne, npar)


@pytest.mark.slow
@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_brown_dennis(opt, npar=4):
    tst_opt(opt, _tstoptfct.brown_dennis, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_osborne1(opt, npar=5):
    tst_opt(opt, _tstoptfct.osborne1, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_biggs(opt, npar=6):
    tst_opt(opt, _tstoptfct.biggs, npar)


@pytest.mark.parametrize("opt", [lmdif,
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_osborne2(opt, npar=11):
    tst_opt(opt, _tstoptfct.osborne2, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_watson(opt, npar=6):
    tst_opt(opt, _tstoptfct.watson, npar)


@pytest.mark.parametrize("opt", [lmdif,
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_extended_rosenbrock(opt, npar=12):
    tst_opt(opt, _tstoptfct.rosenbrock, npar)


@pytest.mark.slow
@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_extended_powell_singular(opt, npar=8):
    tst_opt(opt, _tstoptfct.powell_singular, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_penaltyI(opt, npar=4):
    tst_opt(opt, _tstoptfct.penaltyI, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_penaltyII(opt, npar=4):
    tst_opt(opt, _tstoptfct.penaltyII, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_variably_dimensioned(opt, npar=6):
    tst_opt(opt, _tstoptfct.variably_dimensioned, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_trigonometric(opt, npar=4):
    tst_opt(opt, _tstoptfct.trigonometric, npar)


@pytest.mark.parametrize("opt", [pytest.param(lmdif, marks=pytest.mark.xfail),
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 pytest.param(montecarlo,
                                              marks=pytest.mark.xfail),
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_brown_almost_linear(opt, npar=3):
    tst_opt(opt, _tstoptfct.brown_almost_linear, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_discrete_boundary(opt, npar=4):
    tst_opt(opt, _tstoptfct.discrete_boundary, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_discrete_integral(opt, npar=4):
    tst_opt(opt, _tstoptfct.discrete_integral, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_broyden_tridiagonal(opt, npar=4):
    tst_opt(opt, _tstoptfct.broyden_tridiagonal, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_broyden_banded(opt, npar=4):
    tst_opt(opt, _tstoptfct.broyden_banded, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_linear_fullrank(opt, npar=4):
    tst_opt(opt, _tstoptfct.linear_fullrank, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_linear_fullrank1(opt, npar=4):
    tst_opt(opt, _tstoptfct.linear_fullrank1, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_linear_fullrank0cols0rows0(opt, npar=4):
    tst_opt(opt, _tstoptfct.linear_fullrank0cols0rows, npar)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_linear_chebyquad(opt, npar=9):
    tst_opt(opt, _tstoptfct.chebyquad, npar)

###############################################################################


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Ackley(opt, npar=4):
    tst_opt(opt, _tstoptfct.Ackley, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Booth(opt, npar=6):
    tst_opt(opt, _tstoptfct.Booth, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Bohachevsky1(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky1, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Bohachevsky2(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky2, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Bohachevsky3(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky3, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Branin(opt, npar=2):
    tst_opt(opt, _tstoptfct.Branin, npar)


# why the following test fails for AMD64?
# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  montecarlo, neldermead])
# def test_Branin2(opt, npar=2):
#     tst_opt(opt, _tstoptfct.Branin2, npar)

# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_Chichinadze(opt, npar=2):
#     tst_opt(opt, _tstoptfct.Chichinadze, npar)

# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_Cola(opt, npar=17):
#     tst_opt(opt, _tstoptfct.Cola, npar)

@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Colville(opt, npar=4):
    tst_opt(opt, _tstoptfct.Colville, npar)

def test_minim_no_reflect(reltol=1.0e-3, abstol=1.0e-3):
     fct = _tstoptfct.Colville
     wrong_fval = 174.28569111739617
     x0, xmin, xmax, fmin = init( fct.__name__, 4)
     status, x, fval, msg, xtra = minim(fct, x0, xmin, xmax, reflect=False)
     assert fval == pytest.approx(wrong_fval, rel=reltol, abs=abstol)
     assert fmin != wrong_fval

# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_dcs(opt, npar=5):
#     tst_opt(opt, _tstoptfct.dcs, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_decanom(opt, npar=2):
    tst_opt(opt, _tstoptfct.decanom, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_dodecal(opt, npar=3):
    tst_opt(opt, _tstoptfct.dodecal, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_DixonPrice(opt, npar=5):
    tst_opt(opt, _tstoptfct.DixonPrice, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_Easom(opt, npar=5):
    tst_opt(opt, _tstoptfct.Easom, npar)


# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_factor(opt, npar=5):
#     tst_opt(opt, _tstoptfct.factor, npar)

@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Func1(opt, npar=2):
    tst_opt(opt, _tstoptfct.Func1, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Griewank(opt, npar=2):
    tst_opt(opt, _tstoptfct.Griewank, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_Hansen(opt, npar=2):
    tst_opt(opt, _tstoptfct.Hansen, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Hartman6(opt, npar=6):
    tst_opt(opt, _tstoptfct.Hartman6, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Himmelblau(opt, npar=2):
    tst_opt(opt, _tstoptfct.Himmelblau, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Holzman1(opt, npar=3):
    tst_opt(opt, _tstoptfct.Holzman1, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Holzman2(opt, npar=3):
    tst_opt(opt, _tstoptfct.Holzman2, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_Judge(opt, npar=3):
    tst_opt(opt, _tstoptfct.Judge, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_Levy(opt, npar=4):
    tst_opt(opt, _tstoptfct.Levy, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_McCormick(opt, npar=2):
    tst_opt(opt, _tstoptfct.McCormick, npar)


# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_McKinnon(opt, npar=2):
#     tst_opt(opt, _tstoptfct.McKinnon, npar)

@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Michalewicz(opt, npar=2):
    tst_opt(opt, _tstoptfct.Michalewicz, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_Paviani(opt, npar=10):
    tst_opt(opt, _tstoptfct.Paviani, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Rastrigin(opt, npar=4):
    tst_opt(opt, _tstoptfct.Rastrigin, npar)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_seqp(opt, npar=2):
    tst_opt(opt, _tstoptfct.seqp, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Shekel5(opt, npar=4):
    tst_opt(opt, _tstoptfct.Shekel5, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Shekel7(opt, npar=4):
    tst_opt(opt, _tstoptfct.Shekel7, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_Shekel10(opt, npar=4):
    tst_opt(opt, _tstoptfct.Shekel10, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 pytest.param(montecarlo,
                                              marks=pytest.mark.xfail),
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_ShekelModified(opt, npar=2):
    tst_opt(opt, _tstoptfct.ShekelModified, npar)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_Shubert(opt, npar=2):
    tst_opt(opt, _tstoptfct.Shubert, npar)


# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_SixHumpCamel(opt, npar=2):
#     tst_opt(opt, _tstoptfct.SixHumpCamel, npar)

@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Trecanni(opt, npar=2):
    tst_opt(opt, _tstoptfct.Trecanni, npar)


# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_Trefethen4(opt, npar=2):
#     tst_opt(opt, _tstoptfct.Trefethen4, npar)
