#
#  Copyright (C) 2010, 2016, 2018, 2019, 2020, 2021, 2022, 2023
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

import time

import numpy as np

import pytest

from sherpa.data import Data1D
from sherpa.models.basic import Gauss1D
from sherpa.optmethods import LevMar
from sherpa.stats import LeastSq
from sherpa.fit import Fit, DataSimulFit, SimulFitModel
from sherpa.utils.parallel import parallel_map, parallel_map_funcs


@pytest.mark.parametrize("ntasks", [1, 2, 8])
def test_parallel_map_on_error(ntasks, caplog):
    """What happens if one of the processes raises an error?"""

    def func(x):
        if x == 2:
            # Assume the other tasks will have completed within 0.1
            # second, so this failure is at, or near the end, of the
            # tasks being processed.
            #
            time.sleep(0.1)
            raise ValueError("x can not be 2")
        return x

    args = [-3, 0, 1, 2, 3, 4]
    with pytest.raises(ValueError):
        parallel_map(func, args, numcores=ntasks)


@pytest.mark.parametrize("num_tasks, num_segments",
                         [
                            (1, 1),
                            (8, 1),
                            (1, 8),
                            (10, 5),
                            (5, 10),
                            (5, 5)
                         ])
def test_parallel_map(num_tasks, num_segments):
    """Simple tests of test_parallel_map"""

    iterable = [np.arange(1, 2 + 2 * i) for i in range(num_segments)]

    result = list(map(np.sum, iterable))
    result = np.asarray(result)

    pararesult = parallel_map(np.sum, iterable, num_tasks)

    assert np.asarray(pararesult) == pytest.approx(result)


def test_parallel_map_funcs1():
    for arg in range(1, 5):
        func = [np.sum]
        funcs = arg * func
        datas = [np.arange(1, 2+2*i) for i in range(arg)]
        result = []
        for func, data in zip(funcs, datas):
            result.extend(np.asarray(list(map(func, data))))

        assert parallel_map_funcs(funcs, datas, arg) == pytest.approx(result)


def test_parallel_map_funcs2():
    # TODO: should this pass through ncores or not?
    def tst(ncores, sg, stat, opt):
        sd = DataSimulFit('sd', [d, d], numcores=2)
        f = Fit(sd, sg, stat, opt)
        result = f.fit()
        return result

    def cmp_results(result, tol=1.0e-3):
        assert result.succeeded
        parvals = (1.7555670572301785, 1.5092728216164186, 4.893136872267538)
        assert result.numpoints == 200

        # use tol in approx?
        assert result.parvals == pytest.approx(parvals)

    np.random.seed(0)
    x = np.linspace(-5., 5., 100)
    ampl = 5
    pos = 1.5
    sigma = 0.75
    err = 0.25
    y = ampl * np.exp(-0.5 * (x - pos)**2 / sigma**2)
    y += np.random.normal(0., err, x.shape)
    d = Data1D('junk', x, y)
    g = Gauss1D()
    opt = LevMar()
    stat = LeastSq()
    sg = SimulFitModel('sg', [g, g])

    result = tst(1, sg, stat, opt)
    cmp_results(result)

    result = tst(2, sg, stat, opt)
    cmp_results(result)
