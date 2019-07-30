from __future__ import print_function
#
#  Copyright (C) 2019  Smithsonian Astrophysical Observatory
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


from sherpa.utils import _ncpus
from sherpa import get_config
from sherpa.utils.testing import SherpaTestCase
from sherpa.optmethods import optfcts
from sherpa.optmethods import _tstoptfct
from sherpa.optmethods.ncoresnm import ncoresNelderMead
from sherpa.optmethods.ncoresde import ncoresDifEvo
import sherpa.optmethods.opt as tstopt

class test_optmethods(SherpaTestCase):

    def setUp(self):
        self.tolerance = 1.0e-2
        self.ncores_nm = ncoresNelderMead()
        self.ncores_de = ncoresDifEvo()
        self.numpar = 10
        
    def print_result( self, name, f, x, nfev ):
        print('%s(%s) = %g in %d nfev' % (name, x, f, nfev))

    def tst_opt(self, opt, fcn, x0, xmin, xmax, fmin):
        def func(arg):
            return fcn(arg)[0]
        nfev, fval, par = opt(func, x0, xmin, xmax)
        self.assertEqualWithinTol(fval, fmin, self.tolerance)        

    def tst_func(self, opt, func, x0, xmin, xmax, fmin):
        nfev, fval, par = opt(func, x0, xmin, xmax)
        self.assertEqualWithinTol(fval, fmin, self.tolerance)        

    def test_Ackley(self):
        xmin = self.numpar * [-32.768]
        xmax = self.numpar * [32.768]
        x0 = self.numpar * [12.3]
        if _ncpus != 1:
            self.tst_func(self.ncores_nm, tstopt.Ackley, x0, xmin, xmax, 0.0)
        self.tst_func(self.ncores_de, tstopt.Ackley, x0, xmin, xmax, 0.0)

    def test_Bohachevsky1(self):
        xmin = [-100, -100]
        xmax = [100, 100]
        x0 = [-12, 10]
        self.tst_func(self.ncores_nm, tstopt.Bohachevsky1, x0, xmin, xmax, 0.0)
        self.tst_func(self.ncores_de, tstopt.Bohachevsky1, x0, xmin, xmax, 0.0)

    def test_Bohachevsky2(self):
        xmin = [-100, -100]
        xmax = [100, 100]
        x0 = [12, 10]
        self.tst_func(self.ncores_nm, tstopt.Bohachevsky2, x0, xmin, xmax, 0.0)
        self.tst_func(self.ncores_de, tstopt.Bohachevsky2, x0, xmin, xmax, 0.0)

    def test_Bohachevsky3(self):
        xmin = [-100, -100]
        xmax = [100, 100]
        x0 = [-61.2, 51.0]
        self.tst_func(self.ncores_nm, tstopt.Bohachevsky3, x0, xmin, xmax, 0.0)
        self.tst_func(self.ncores_de, tstopt.Bohachevsky3, x0, xmin, xmax, 0.0)

    def test_Booth(self):
        xmin = [-10, -10]
        xmax = [10, 10]
        x0 = [-6.2, 5.0]
        self.tst_func(self.ncores_nm, tstopt.Booth, x0, xmin, xmax, 0.0)
        self.tst_func(self.ncores_de, tstopt.Booth, x0, xmin, xmax, 0.0)

    def test_BoxBetts(self):
        xmin = [0.9, 9.0, 0.9]
        xmax = [1.2, 11.2, 1.2]
        x0 = [(xmin[0] + xmax[0]) * 0.5, (xmin[1] + xmax[1]) * 0.5,
              (xmin[2] + xmax[2]) * 0.5]
        self.tst_func(self.ncores_nm, tstopt.BoxBetts, x0, xmin, xmax, 0.0)
        self.tst_func(self.ncores_de, tstopt.BoxBetts, x0, xmin, xmax, 0.0)

    def test_Branin(self):
        xmin = [-5, 0]
        xmax = [10, 15]
        x0 = [-3.2, 5.0]
        self.tst_func(self.ncores_nm, tstopt.Branin, x0, xmin, xmax, 0.397887)
        self.tst_func(self.ncores_de, tstopt.Branin, x0, xmin, xmax, 0.397887)

    def test_DixonPrixe(self):
        xmin = self.numpar * [-10]
        xmax = self.numpar * [10]
        x0 = self.numpar * [-1.]
        self.tst_func(self.ncores_de, tstopt.DixonPrice, x0, xmin, xmax, 0.0)

    def test_Easom(self):
        xmin = [-100, -100]
        xmax = [100, 100]
        x0 = [25.0, 25.0]
        self.tst_func(self.ncores_de, tstopt.Easom, x0, xmin, xmax, -1.0)

    def test_FreudensteinRoth(self):
        xmin = self.numpar * [-1000, -1000]
        xmax = self.numpar * [1000, 1000]
        x0 = self.numpar * [0.5, -2]
        self.tst_func(self.ncores_nm, tstopt.FreudensteinRoth, x0, xmin, xmax,
                      0.0)
        self.tst_func(self.ncores_de, tstopt.FreudensteinRoth, x0, xmin, xmax,
                      0.0)

    def test_GoldsteinPrice(self):
        xmin = [-2, -2]
        xmax = [2, 2]
        x0 = [-1, 1]
        self.tst_func(self.ncores_nm, tstopt.GoldsteinPrice, x0, xmin, xmax,
                      3.0)
        self.tst_func(self.ncores_de, tstopt.GoldsteinPrice, x0, xmin, xmax,
                      3.0)

    def test_Hump(self):
        xmin = [-5, -5]
        xmax = [5, 5]
        x0 = [-3.2, 5.0]
        self.tst_func(self.ncores_nm, tstopt.Hump, x0, xmin, xmax, 0.0)
        self.tst_func(self.ncores_de, tstopt.Hump, x0, xmin, xmax, 0.0)        
        
    def test_Matyas(self):
        xmin = [-10, -10]
        xmax = [10, 10]
        x0 = [-3.2, 5.0]
        self.tst_func(self.ncores_nm, tstopt.Matyas, x0, xmin, xmax, 0.0)
        self.tst_func(self.ncores_de, tstopt.Matyas, x0, xmin, xmax, 0.0)

    def test_McCormick(self):
        xmin = [-1.5, -3.0]
        xmax = [4.0, 4.0]
        x0 = [0.0, 0.0]
        if _ncpus != 1:
            self.tst_func(self.ncores_nm, tstopt.McCormick, x0, xmin, xmax,
                          -1.91)
        # self.tst_func(self.ncores_de, tstopt.McCormick, x0, xmin, xmax, -1.91)
        
    def test_Paviani(self):
        xmin = self.numpar * [2.001]
        xmax = self.numpar * [9.999]
        x0 = self.numpar * [5.0]
        self.tst_func(self.ncores_nm, tstopt.Paviani, x0, xmin, xmax, -45.7)
        self.tst_func(self.ncores_de, tstopt.Paviani, x0, xmin, xmax, -45.7)

    def test_Rastrigin(self):
        xmin = self.numpar * [-5.12]
        xmax = self.numpar * [5.12]
        x0	 = self.numpar * [-2.0]
        if _ncpus != 1:
            self.tst_func(self.ncores_nm, tstopt.Rastrigin, x0, xmin, xmax,
                          0.0)
        self.tst_func(self.ncores_de, tstopt.Rastrigin, x0, xmin, xmax, 0.0)

    def test_Shubert(self):
        xmin = [-10, -10]
        xmax = [10, 10]
        x0	 = [-2.0, 5.0]
        if _ncpus != 1:
            self.tst_func(self.ncores_nm, tstopt.Shubert, x0, xmin, xmax,
                          -186.7309)
        self.tst_func(self.ncores_de, tstopt.Shubert, x0, xmin, xmax,
                      -186.7309)

    #
    # comment out a few time consuming tests, anyone
    # modifying the ncores opt should uncomment the tests
    #
    # def test_Levy(self):
    #     xmin = self.numpar * [-10]
    #     xmax = self.numpar * [10]
    #     x0 = self.numpar * [-5.]
    #     self.tst_func(self.ncores_nm, tstopt.Levy, x0, xmin, xmax, 0.0)
    #     self.tst_func(self.ncores_de, tstopt.Levy, x0, xmin, xmax, 0.0)
    # def test_Sphere(self):
    #     xmin = self.numpar * [-5.12]
    #     xmax = self.numpar * [5.12]
    #     x0 = self.numpar * [-2.0]
    #     self.tst_func(self.ncores_nm, tstopt.Sphere, x0, xmin, xmax, 0.0)
    #     self.tst_func(self.ncores_de, tstopt.Sphere, x0, xmin, xmax, 0.0)

    # def test_SumSquares(self):
    #     xmin = self.numpar * [-10]
    #     xmax = self.numpar * [10]
    #     x0 = self.numpar * [-2.0]
    #     self.tst_func(self.ncores_nm, tstopt.SumSquares, x0, xmin, xmax, 0.0)
    #     self.tst_func(self.ncores_de, tstopt.SumSquares, x0, xmin, xmax, 0.0)

    # def test_Zakharov(self):
    #     xmin = self.numpar * [-5, -5]
    #     xmax = self.numpar * [10, 10]
    #     x0   = self.numpar * [0.5, -2]
    #     self.tst_func(self.ncores_nm, tstopt.Zakharov, x0, xmin, xmax, 0.0)
    #     self.tst_func(self.ncores_de, tstopt.Zakharov, x0, xmin, xmax, 0.0)
