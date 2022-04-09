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
import operator
import numpy

from sherpa.data import Data1D
from sherpa.fit import Fit
from sherpa.models import Polynom1D
from sherpa.optmethods import NelderMead, LevMar
from sherpa.optmethods.transformation import Transformation
from sherpa.stats import Chi2Gehrels
from sherpa.optmethods import _tstoptfct
from sherpa.optmethods.optfcts import lmdif, minim, montecarlo, neldermead
from sherpa.utils import _ncpus

EPSILON64 = numpy.float_(numpy.finfo(numpy.float64).eps)


class TestTransformation(Transformation):
    def __init__(self):
        Transformation.__init__(self)
        return
    def init(self, max_range=8192):
        super().init(max_range)
        return


def expand(arg, factor, val, op1, op2 ):
    if arg > 0.0:
        return arg * op1(1.0, factor)
    elif arg < 0.0:
        return arg * op2(1.0, factor)
    else:
        return factor * val

def init(name, npar, transform, factor):
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    if transform:
        xmin = numpy.empty_like(x0)
        xmax = numpy.empty_like(x0)
        for ii, arg in enumerate(x0):
            xmin[ii] = expand(arg, factor, -1, operator.sub, operator.add)
            xmax[ii] = expand(arg, factor, 1, operator.add, operator.sub)
        assert sum(xmax - xmin) / len(x0) <= 8192
    return x0, xmin, xmax, fmin

def tst_opt(opt, fct, npar, transform, factor=4, reltol=1.0e-3, abstol=1.0e-3):
    """The central function for all the optimization test
    1) Out of the 35 tests from:
    J. MORE', B. GARBOW & K. HILLSTROM,
    "Algorithm 566: Fortran Subroutines for Testing Unconstrained
    Optimization Software.", ACM TOMS, VOL. 7, PAGES 14-41 AND 136-140, 1981
    xfail: lmdif(6), minim(5), neldermead(3), moncar(2)
    2) The remaining random 32 'global' func tests:
    xfail: minim(14), montecarlo(0), neldermead(10)
    """
    x0, xmin, xmax, fmin = init(fct.__name__, npar, transform, factor)
    if opt == montecarlo:
        status, x, fval, msg, xtra = opt(fct, x0, xmin, xmax)
    else:
        status, x, fval, msg, xtra = opt(fct, x0, xmin, xmax,
                                         transformation=transform)
    assert fmin == pytest.approx(fval, rel=reltol, abs=abstol)
    if opt == lmdif and _ncpus > 1:
        status, x, fval, msg, xtra = opt(fct, x0, xmin, xmax, numcores=_ncpus,
                                         transformation=transform)
        assert fmin == pytest.approx(fval, rel=reltol, abs=abstol)
        assert xtra.get('num_parallel_map') != 0


###############################################################################
@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_rosenbrock(opt, npar=4):
    tst_opt(opt, _tstoptfct.rosenbrock, npar, False)
    tst_opt(opt, _tstoptfct.rosenbrock, npar, True)


@pytest.mark.parametrize("opt", [pytest.param(lmdif, marks=pytest.mark.xfail),
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 pytest.param(montecarlo,
                                              marks=pytest.mark.xfail),
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_freudensteinroth(opt, npar=4):
    tst_opt(opt, _tstoptfct.freudenstein_roth, npar, False)
    tst_opt(opt, _tstoptfct.freudenstein_roth, npar, Transformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_powell_badly_scaled(opt, npar=2):
    tst_opt(opt, _tstoptfct.powell_badly_scaled, npar, False)
    tst_opt(opt, _tstoptfct.powell_badly_scaled, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_brown_badly_scaled(opt, npar=2):
    tst_opt(opt, _tstoptfct.brown_badly_scaled, npar, False)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_beale(opt, npar=2):
    tst_opt(opt, _tstoptfct.beale, npar, False)
    tst_opt(opt, _tstoptfct.beale, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_jennrich_sampson(opt, npar=2):
    tst_opt(opt, _tstoptfct.jennrich_sampson, npar, False)
    tst_opt(opt, _tstoptfct.jennrich_sampson, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_helical_valley(opt, npar=3):
    tst_opt(opt, _tstoptfct.helical_valley, npar, False)
    tst_opt(opt, _tstoptfct.helical_valley, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_bard(opt, npar=3):
    tst_opt(opt, _tstoptfct.bard, npar, False)
    tst_opt(opt, _tstoptfct.bard, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_gaussian(opt, npar=3):
    tst_opt(opt, _tstoptfct.gaussian, npar, False)
    tst_opt(opt, _tstoptfct.gaussian, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(lmdif, marks=pytest.mark.xfail),
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_meyer(opt, npar=3):
    tst_opt(opt, _tstoptfct.meyer, npar, False)


@pytest.mark.slow
@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_gulf_research_development(opt, npar=3):
    tst_opt(opt, _tstoptfct.gulf_research_development, npar, False)
    tst_opt(opt, _tstoptfct.gulf_research_development, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_box3d(opt, npar=3):
    tst_opt(opt, _tstoptfct.box3d, npar, False)
    tst_opt(opt, _tstoptfct.box3d, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_powell_singular(opt, npar=4):
    tst_opt(opt, _tstoptfct.powell_singular, npar, False)
    tst_opt(opt, _tstoptfct.powell_singular, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(lmdif, marks=pytest.mark.xfail),
                                 minim, montecarlo, neldermead])
def test_wood(opt, npar=4):
    tst_opt(opt, _tstoptfct.wood, npar, False)
    tst_opt(opt, _tstoptfct.wood, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_kowalik_osborne(opt, npar=4):
    tst_opt(opt, _tstoptfct.kowalik_osborne, npar, False)
    tst_opt(opt, _tstoptfct.kowalik_osborne, npar, TestTransformation(), 2)


@pytest.mark.slow
@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_brown_dennis(opt, npar=4):
    tst_opt(opt, _tstoptfct.brown_dennis, npar, False)
    tst_opt(opt, _tstoptfct.brown_dennis, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_osborne1(opt, npar=5):
    tst_opt(opt, _tstoptfct.osborne1, npar, False)
    tst_opt(opt, _tstoptfct.osborne1, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_biggs(opt, npar=6):
    tst_opt(opt, _tstoptfct.biggs, npar, False)
    tst_opt(opt, _tstoptfct.biggs, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif,
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_osborne2(opt, npar=11):
    tst_opt(opt, _tstoptfct.osborne2, npar, False)
    tst_opt(opt, _tstoptfct.osborne2, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_watson(opt, npar=6):
    tst_opt(opt, _tstoptfct.watson, npar, False)
    tst_opt(opt, _tstoptfct.watson, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif,
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_extended_rosenbrock(opt, npar=12):
    tst_opt(opt, _tstoptfct.rosenbrock, npar, False)
    tst_opt(opt, _tstoptfct.rosenbrock, npar, TestTransformation())


@pytest.mark.slow
@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_extended_powell_singular(opt, npar=8):
    tst_opt(opt, _tstoptfct.powell_singular, npar, False)
    tst_opt(opt, _tstoptfct.powell_singular, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_penaltyI(opt, npar=4):
    tst_opt(opt, _tstoptfct.penaltyI, npar, False)
    tst_opt(opt, _tstoptfct.penaltyI, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_penaltyII(opt, npar=4):
    tst_opt(opt, _tstoptfct.penaltyII, npar, False)
    tst_opt(opt, _tstoptfct.penaltyII, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_variably_dimensioned(opt, npar=6):
    tst_opt(opt, _tstoptfct.variably_dimensioned, npar, False)
    tst_opt(opt, _tstoptfct.variably_dimensioned, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_trigonometric(opt, npar=4):
    tst_opt(opt, _tstoptfct.trigonometric, npar, False)
    tst_opt(opt, _tstoptfct.trigonometric, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(lmdif, marks=pytest.mark.xfail),
                                 pytest.param(minim, marks=pytest.mark.xfail),
                                 pytest.param(montecarlo,
                                              marks=pytest.mark.xfail),
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_brown_almost_linear(opt, npar=3):
    tst_opt(opt, _tstoptfct.brown_almost_linear, npar, False)
    tst_opt(opt, _tstoptfct.brown_almost_linear, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_discrete_boundary(opt, npar=4):
    tst_opt(opt, _tstoptfct.discrete_boundary, npar, False)
    tst_opt(opt, _tstoptfct.discrete_boundary, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_discrete_integral(opt, npar=4):
    tst_opt(opt, _tstoptfct.discrete_integral, npar, False)
    tst_opt(opt, _tstoptfct.discrete_integral, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_broyden_tridiagonal(opt, npar=4):
    tst_opt(opt, _tstoptfct.broyden_tridiagonal, npar, False)
    tst_opt(opt, _tstoptfct.broyden_tridiagonal, npar, TestTransformation(), 2)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_broyden_banded(opt, npar=4):
    tst_opt(opt, _tstoptfct.broyden_banded, npar, False)
    tst_opt(opt, _tstoptfct.broyden_banded, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_linear_fullrank(opt, npar=4):
    tst_opt(opt, _tstoptfct.linear_fullrank, npar, False)
    tst_opt(opt, _tstoptfct.linear_fullrank, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_linear_fullrank1(opt, npar=4):
    tst_opt(opt, _tstoptfct.linear_fullrank1, npar, False)
    tst_opt(opt, _tstoptfct.linear_fullrank1, npar, TestTransformation())


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_linear_fullrank0cols0rows0(opt, npar=4):
    tst_opt(opt, _tstoptfct.linear_fullrank0cols0rows, npar, False)


@pytest.mark.parametrize("opt", [lmdif, minim, montecarlo, neldermead])
def test_linear_chebyquad(opt, npar=9):
    tst_opt(opt, _tstoptfct.chebyquad, npar, False)
    tst_opt(opt, _tstoptfct.chebyquad, npar, TestTransformation())

###############################################################################


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Ackley(opt, npar=4):
    tst_opt(opt, _tstoptfct.Ackley, npar, False)
    tst_opt(opt, _tstoptfct.Ackley, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Booth(opt, npar=6):
    tst_opt(opt, _tstoptfct.Booth, npar, False)
    tst_opt(opt, _tstoptfct.Booth, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Bohachevsky1(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky1, npar, False)
    tst_opt(opt, _tstoptfct.Bohachevsky1, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Bohachevsky2(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky2, npar, False)
    tst_opt(opt, _tstoptfct.Bohachevsky2, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Bohachevsky3(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky3, npar, False)
    tst_opt(opt, _tstoptfct.Bohachevsky3, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Branin(opt, npar=2):
    tst_opt(opt, _tstoptfct.Branin, npar, False)


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
    tst_opt(opt, _tstoptfct.Colville, npar, False)


def test_minim_no_reflect(reltol=1.0e-3, abstol=1.0e-3):
     fct = _tstoptfct.Colville
     wrong_fval = 174.28569111739617
     x0, xmin, xmax, fmin = init( fct.__name__, 4, None, False)
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
    tst_opt(opt, _tstoptfct.decanom, npar, False)
    tst_opt(opt, _tstoptfct.decanom, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_dodecal(opt, npar=3):
    tst_opt(opt, _tstoptfct.dodecal, npar, False)
    tst_opt(opt, _tstoptfct.dodecal, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_DixonPrice(opt, npar=5):
    tst_opt(opt, _tstoptfct.DixonPrice, npar, False)
    tst_opt(opt, _tstoptfct.DixonPrice, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_Easom(opt, npar=5):
    tst_opt(opt, _tstoptfct.Easom, npar, False)
    tst_opt(opt, _tstoptfct.Easom, npar, TestTransformation())


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
    tst_opt(opt, _tstoptfct.Func1, npar, False)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Griewank(opt, npar=2):
    tst_opt(opt, _tstoptfct.Griewank, npar, False)
    tst_opt(opt, _tstoptfct.Griewank, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_Hansen(opt, npar=2):
    tst_opt(opt, _tstoptfct.Hansen, npar, False)
    tst_opt(opt, _tstoptfct.Hansen, npar, TestTransformation(), 6)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Hartman6(opt, npar=6):
    tst_opt(opt, _tstoptfct.Hartman6, npar, False)
    tst_opt(opt, _tstoptfct.Hartman6, npar, TestTransformation(), 6)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Himmelblau(opt, npar=2):
    tst_opt(opt, _tstoptfct.Himmelblau, npar, False)
    tst_opt(opt, _tstoptfct.Himmelblau, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Holzman1(opt, npar=3):
    tst_opt(opt, _tstoptfct.Holzman1, npar, False)


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Holzman2(opt, npar=3):
    tst_opt(opt, _tstoptfct.Holzman2, npar, False)
    tst_opt(opt, _tstoptfct.Holzman2, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_Judge(opt, npar=3):
    tst_opt(opt, _tstoptfct.Judge, npar, False)
    tst_opt(opt, _tstoptfct.Judge, npar, TestTransformation(), 2)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_Levy(opt, npar=4):
    tst_opt(opt, _tstoptfct.Levy, npar, False)
    # tst_opt(opt, _tstoptfct.Levy, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_McCormick(opt, npar=2):
    tst_opt(opt, _tstoptfct.McCormick, npar, False)
    tst_opt(opt, _tstoptfct.McCormick, npar, TestTransformation(), 2)


# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_McKinnon(opt, npar=2):
#     tst_opt(opt, _tstoptfct.McKinnon, npar)

@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Michalewicz(opt, npar=2):
    tst_opt(opt, _tstoptfct.Michalewicz, npar, False)
    tst_opt(opt, _tstoptfct.Michalewicz, npar, TestTransformation(), 2)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_Paviani(opt, npar=10):
    tst_opt(opt, _tstoptfct.Paviani, npar, False)
    tst_opt(opt, _tstoptfct.Paviani, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Rastrigin(opt, npar=4):
    tst_opt(opt, _tstoptfct.Rastrigin, npar, False)
    tst_opt(opt, _tstoptfct.Rastrigin, npar, TestTransformation())


@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_seqp(opt, npar=2):
    tst_opt(opt, _tstoptfct.seqp, npar, False)
    tst_opt(opt, _tstoptfct.seqp, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Shekel5(opt, npar=4):
    tst_opt(opt, _tstoptfct.Shekel5, npar, False)


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_Shekel7(opt, npar=4):
    tst_opt(opt, _tstoptfct.Shekel7, npar, False)
    tst_opt(opt, _tstoptfct.Shekel7, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo,
                                 pytest.param(neldermead, marks=pytest.mark.xfail)])
def test_Shekel10(opt, npar=4):
    tst_opt(opt, _tstoptfct.Shekel10, npar, False)
    tst_opt(opt, _tstoptfct.Shekel10, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 pytest.param(montecarlo,
                                              marks=pytest.mark.xfail),
                                 pytest.param(neldermead,
                                              marks=pytest.mark.xfail)])
def test_ShekelModified(opt, npar=2):
    tst_opt(opt, _tstoptfct.ShekelModified, npar, False)
    tst_opt(opt, _tstoptfct.ShekelModified, npar, TestTransformation())


@pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
                                 montecarlo, neldermead])
def test_Shubert(opt, npar=2):
    tst_opt(opt, _tstoptfct.Shubert, npar, False)
    tst_opt(opt, _tstoptfct.Shubert, npar, TestTransformation())


# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_SixHumpCamel(opt, npar=2):
#     tst_opt(opt, _tstoptfct.SixHumpCamel, npar)

@pytest.mark.parametrize("opt", [minim, montecarlo, neldermead])
def test_Trecanni(opt, npar=2):
    tst_opt(opt, _tstoptfct.Trecanni, npar, False)
    tst_opt(opt, _tstoptfct.Trecanni, npar, TestTransformation())


# @pytest.mark.parametrize("opt", [pytest.param(minim, marks=pytest.mark.xfail),
#                                  pytest.param(montecarlo,
#                                               marks=pytest.mark.xfail),
#                                  pytest.param(neldermead,
#                                               marks=pytest.mark.xfail)])
# def test_Trefethen4(opt, npar=2):
#     tst_opt(opt, _tstoptfct.Trefethen4, npar)

################################Transform###################################
def setup_4_transform(opt):
    x = [ 0.5,  1.5,  2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5]
    y = [1.6454, 1.7236, 1.9472, 2.2348, 2.6187, 2.8642, 3.1263, 3.2073, \
         3.2852, 3.3092, 3.4496]
    z = 11 * [0.04114]
    data = Data1D('1', x, y, z)
    poly = Polynom1D()
    poly.pars[1].thaw()
    poly.pars[2].thaw()
    poly.pars[3].thaw()
    method=opt()
    method.config['transformation'] = TestTransformation()
    if opt == LevMar:
        method.config['epsfcn'] = EPSILON64
    fit = Fit(data, poly, method=method, stat=Chi2Gehrels())
    parvals = (1.498430848959475, 0.14469980089323808, 0.032293618881299224,
               -0.0027772921523028238)
    return parvals, poly, fit


@pytest.mark.parametrize("opt", [LevMar, NelderMead])
def test_transform_hi(opt, factor=4):
    parvals, poly, fit = setup_4_transform(opt)
    poly.c0.max = expand(poly.c0.val, factor, 1, operator.add, operator.sub)
    poly.c1.max = expand(poly.c1.val, factor, 1, operator.add, operator.sub)
    poly.c2.max = expand(poly.c2.val, factor, 1, operator.add, operator.sub)
    poly.c3.max = expand(poly.c3.val, factor, 1, operator.add, operator.sub)
    result = fit.fit()
    assert parvals == pytest.approx(result.parvals, abs=1.0e-6, rel=1.0e-6)
    if opt == LevMar:
        expect = [0.052084507576, 0.041377269137, 0.008749967892,
                  0.000523424920]
        covarerr = numpy.sqrt(result.extra_output['covar'].diagonal())
        assert expect == pytest.approx(covarerr, abs=1.0e-6, rel=1.0e-6)


@pytest.mark.parametrize("opt", [LevMar, NelderMead])
def test_transform_lo(opt, factor=4):
    parvals, poly, fit = setup_4_transform(opt)
    poly.c0.min = expand(poly.c0.val, factor, -1, operator.sub, operator.add)
    poly.c1.min = expand(poly.c1.val, factor, -1, operator.sub, operator.add)
    poly.c2.min = expand(poly.c2.val, factor, -1, operator.sub, operator.add)
    poly.c3.min = expand(poly.c3.val, factor, -1, operator.sub, operator.add)
    result = fit.fit()
    assert parvals == pytest.approx(result.parvals, abs=1.0e-6, rel=1.0e-6)
    if opt == LevMar:
        expect = [0.052084557616, 0.041377269241, 0.008749967662,
                  0.000523424913]
        covarerr = numpy.sqrt(result.extra_output['covar'].diagonal())
        assert expect == pytest.approx(covarerr, abs=1.0e-6, rel=1.0e-6)


@pytest.mark.parametrize("opt", [LevMar, NelderMead])
def test_transform_lohi(opt, factor=4):
    parvals, poly, fit = setup_4_transform(opt)
    poly.c0.min = expand(poly.c0.val, factor, -1, operator.sub, operator.add)
    poly.c0.max = expand(poly.c0.val, factor, 1, operator.add, operator.sub)
    poly.c1.min = expand(poly.c1.val, factor, -1, operator.sub, operator.add)
    poly.c1.max = expand(poly.c1.val, factor, 1, operator.add, operator.sub)
    poly.c2.min = expand(poly.c2.val, factor, -1, operator.sub, operator.add)
    poly.c2.max = expand(poly.c2.val, factor, 1, operator.add, operator.sub)
    poly.c3.min = expand(poly.c3.val, factor, -1, operator.sub, operator.add)
    poly.c3.max = expand(poly.c3.val, factor, 1, operator.add, operator.sub)
    result = fit.fit()
    assert parvals == pytest.approx(result.parvals, abs=1.0e-6, rel=1.0e-6)
    if opt == LevMar:
        expect = [0.052084634545, 0.041377280240, 0.008749970667,
                  0.000523421758]
        covarerr = numpy.sqrt(result.extra_output['covar'].diagonal())
        assert expect == pytest.approx(covarerr, abs=1.0e-6, rel=1.0e-6)


@pytest.mark.parametrize("opt", [LevMar, NelderMead])
def test_transform_partial(opt, factor=4):
    parvals, poly, fit = setup_4_transform(opt)
    poly.c0.min = expand(poly.c0.val, factor, -1, operator.sub, operator.add)
    poly.c0.max = expand(poly.c0.val, factor, 1, operator.add, operator.sub)
    poly.c1.min = expand(poly.c1.val, factor, -1, operator.sub, operator.add)
    poly.c2.max = expand(poly.c2.val, factor, 1, operator.add, operator.sub)
    result = fit.fit()
    assert parvals == pytest.approx(result.parvals, abs=1.0e-6, rel=1.0e-6)
    if opt == LevMar:
        expect = [0.052084632179, 0.041377269482, 0.008749967999,
                  0.000523424925]
        covarerr = numpy.sqrt(result.extra_output['covar'].diagonal())
        assert expect == pytest.approx(covarerr, abs=1.0e-6, rel=1.0e-6)


@pytest.mark.parametrize("opt", [lmdif, minim, neldermead])
def test_will_not_work(opt, npar=2):


    class WillNotWork:
        def __init__(self):
            pass

    with pytest.raises(TypeError) as excinfo:
        tst_opt(opt, _tstoptfct.rosenbrock, npar, WillNotWork())

    errmsg = "transformation must inherit from class Transformation"
    assert errmsg in str(excinfo.value)
################################Transform###################################
