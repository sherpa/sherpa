#
#  Copyright (C) 2019, 2020  Smithsonian Astrophysical Observatory
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

import numpy as np
import pytest

from sherpa.optmethods import _tstoptfct
from sherpa.optmethods.ncoresnm import ncoresNelderMead
from sherpa.optmethods.ncoresde import ncoresDifEvo
from sherpa.utils import _ncpus

def init(name, npar):
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    return x0, xmin, xmax, fmin

def tst_opt(opt, fcn, npar, reltol=1.0e-3, abstol=1.0e-3):
    def func(arg):
        return fcn(arg)[0]
    x0, xmin, xmax, fmin = init(fcn.__name__, npar)
    nfev, fval, par = opt(func, x0, xmin, xmax)
    assert fmin == pytest.approx(fval, rel=reltol, abs=abstol)

NCORES_NM = ncoresNelderMead()
NCORES_DE = ncoresDifEvo()

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Ackley(opt, npar=4):
    tst_opt(opt, _tstoptfct.Ackley, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Bohachevsky1(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky1, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Bohachevsky2(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky2, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Bohachevsky3(opt, npar=2):
    tst_opt(opt, _tstoptfct.Bohachevsky3, npar)

@pytest.mark.slow
@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Booth(opt, npar=6):
    tst_opt(opt, _tstoptfct.Booth, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Branin(opt, npar=2):
    tst_opt(opt, _tstoptfct.Branin, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Branin2(opt, npar=2):
    tst_opt(opt, _tstoptfct.Branin2, npar)

@pytest.mark.slow
@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Colville(opt, npar=4):
    tst_opt(opt, _tstoptfct.Colville, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_decanom(opt, npar=2):
    tst_opt(opt, _tstoptfct.decanom, npar)

@pytest.mark.slow
@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_dodecal(opt, npar=3):
    tst_opt(opt, _tstoptfct.dodecal, npar)

@pytest.mark.slow    
@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_DixonPrice(opt, npar=5):
    tst_opt(opt, _tstoptfct.DixonPrice, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Hansen(opt, npar=2):
    tst_opt(opt, _tstoptfct.Hansen, npar)

@pytest.mark.slow    
@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Holzman2(opt, npar=3):
    tst_opt(opt, _tstoptfct.Holzman2, npar)

@pytest.mark.slow
@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Judge(opt, npar=3):
    tst_opt(opt, _tstoptfct.Judge, npar)

@pytest.mark.parametrize("opt", [NCORES_NM])
def test_McCormick(opt, npar=2):
    tst_opt(opt, _tstoptfct.McCormick, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Michalewicz(opt, npar=2):
    tst_opt(opt, _tstoptfct.Michalewicz, npar)

@pytest.mark.parametrize("opt", [NCORES_NM])
def test_Paviani(opt, npar=10):
    tst_opt(opt, _tstoptfct.Paviani, npar)

@pytest.mark.slow
def test_Paviani(npar=10):
    tst_opt(NCORES_DE, _tstoptfct.Paviani, npar)

@pytest.mark.parametrize("opt", [NCORES_NM])
def test_Rastrigin(opt, npar=4):
    tst_opt(opt, _tstoptfct.Rastrigin, npar)

@pytest.mark.slow
def test_Rastrigin(npar=4):
    tst_opt(NCORES_DE, _tstoptfct.Rastrigin, npar)
    
@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_seqp(opt, npar=2):
    tst_opt(opt, _tstoptfct.seqp, npar)

@pytest.mark.slow
def test_seqp(npar=2):
    tst_opt(NCORES_DE, _tstoptfct.seqp, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Shubert(opt, npar=2):
    tst_opt(opt, _tstoptfct.Shubert, npar)

@pytest.mark.parametrize("opt", [NCORES_NM, NCORES_DE])
def test_Trecanni(opt, npar=2):
    tst_opt(opt, _tstoptfct.Trecanni, npar)

    
