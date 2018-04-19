from __future__ import print_function

#
#  Copyright (C) 2018  Smithsonian Astrophysical Observatory
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

import numpy

from sherpa.utils import requires_data, requires_xspec, requires_fits
from sherpa.utils import SherpaTestCase

from sherpa.models import Polynom1D, SimulFitModel
from sherpa.astro.models import Lorentz1D
from sherpa.data import Data1D, DataSimulFit
from sherpa.optmethods import NelderMead
from sherpa.stats import Cash
from sherpa.fit import Fit

import logging
logger = logging.getLogger("sherpa")

class test_cache(SherpaTestCase):

    _fit_same_poly_bench = {
        'numpoints': 153,
        'dof': 151,
        'istatval': 306.0,
        'statval': -19417.55547762454,
        'parnames': ('polynom1d.c0', 'polynom1d.c1'),
        'parvals': numpy.array([1.6874869548492573, 0.4796635263358297])
        }

    _fit_diff_poly_bench = {
        'numpoints': 153,
        'dof': 147,
        'istatval': 306.0,
        'statval': -19418.614329430166,
        'parnames': ('polynom1d.c0', 'polynom1d.c1', 'polynom1d.c0',
                     'polynom1d.c1', 'polynom1d.c0', 'polynom1d.c1'),
        'parvals': numpy.array([1.755035876061136, 0.46847103779088173,
                                1.9361738167035432, 0.4791446119040638,
                                1.342017385707699, 0.49194814027685246])
        }
    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        x1 = [  1.,  3.,  5.,  7.,  9., 11., 13., 15., 17., 19., 21., 23.,
                25., 27., 29., 31., 33., 35., 37., 39., 41., 43., 45.,
                47., 49., 51., 53., 55., 57., 59., 61., 63., 65., 67.,
                69., 71., 73., 75., 77., 79., 81., 83., 85., 87., 89.,
                91., 93., 95., 97., 99.,101.]
        y1 = [ 1., 5., 2., 4., 7.,11., 9., 8.,12.,18.,12.,11.,13.,12.,13.,
               13.,20.,23.,16.,20.,24.,17.,21.,26.,22.,24.,24.,21.,28.,
               28.,26.,25.,34.,26.,34.,33.,25.,38.,31.,43.,35.,42.,50.,
               41.,43.,47.,57.,53.,60.,46.,54.]

        x2 = [ 1.,  3.,  5.,  7.,  9., 11., 13., 15., 17., 19., 21., 23.,
               25., 27., 29., 31., 33., 35., 37., 39., 41., 43., 45., 47.,
               49., 51., 53., 55., 57., 59., 61., 63., 65., 67., 69., 71.,
               73., 75., 77., 79., 81., 83., 85., 87., 89., 91., 93., 95.,
               97., 99.,101.]
        y2 = [ 0., 7., 6., 3., 5., 5., 9.,11.,13., 8.,14.,13.,14.,18.,11.,
               15.,17.,26., 15.,19.,25.,30.,15.,29.,16.,25.,27.,29.,36.,
               41.,22.,27.,33.,32.,45.,37.,38.,38.,34.,52.,40.,41.,31.,
               47.,38.,52.,57.,33.,48.,53.,45.]

        x3 = [ 1.,  3.,  5.,  7.,  9., 11., 13., 15., 17., 19., 21., 23.,
               25., 27., 29., 31., 33., 35., 37., 39., 41., 43., 45., 47.,
               49., 51., 53., 55., 57., 59., 61., 63., 65., 67., 69., 71.,
               73., 75., 77., 79., 81., 83., 85., 87., 89., 91., 93., 95.,
               97., 99.,101.]
        y3 = [ 1., 2., 4., 2., 5., 8.,15.,10.,13.,10.,16.,10.,13.,12.,16.,
               17.,17.,20., 23.,16.,25.,22.,19.,31.,26.,24.,21.,29.,36.,
               30.,33.,30.,37.,27.,36.,32., 42.,44.,39.,30.,40.,33.,39.,
               49.,56.,47.,46.,35.,63.,40.,57.]

        self.d1 = Data1D('1', x1, y1)
        self.d2 = Data1D('2', x2, y2)
        self.d3 = Data1D('3', x3, y3)

    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    def compare_results(self, arg1, arg2, tol=1e-6):

        for key in ["numpoints", "dof"]:
            assert arg1[key] == int(getattr(arg2, key))

        for key in ["istatval", "statval"]:
            numpy.testing.assert_allclose(float(arg1[key]),
                                          float(getattr(arg2, key)), tol)

        for key in ["parvals"]:
            try:
                numpy.testing.assert_allclose(arg1[key],
                                              getattr(arg2, key), tol)
            except AssertionError:
                print('parvals bench: ', arg1[key])
                print('parvals fit:   ', getattr(arg2, key))
                print('results', arg2)
                raise

    def test_diff_cache(self):
        poly1 = Polynom1D()
        poly2 = Polynom1D()
        poly3 = Polynom1D()
        poly1.pars[1].thaw()
        poly2.pars[1].thaw()
        poly3.pars[1].thaw()
        sdata = DataSimulFit('all', (self.d1, self.d2, self.d3))
        smodel = SimulFitModel('diff', (poly1, poly2, poly3))
        sfit = Fit(sdata, smodel, method=NelderMead(), stat=Cash())
        result = sfit.fit()
        self.compare_results(self._fit_diff_poly_bench, result)

    def test_same_cache(self):
        poly = Polynom1D()
        poly.pars[1].thaw()
        sdata = DataSimulFit('all', (self.d1, self.d2, self.d3))
        smodel = SimulFitModel('same', (poly, poly, poly))
        sfit = Fit(sdata, smodel, method=NelderMead(), stat=Cash())
        result = sfit.fit()
        self.compare_results(self._fit_same_poly_bench, result)
