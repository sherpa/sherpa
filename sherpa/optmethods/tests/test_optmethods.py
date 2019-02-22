from __future__ import print_function
#
#  Copyright (C) 2007, 2015, 2016, 2018, 2019  Smithsonian Astrophysical Observatory
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
from sherpa.utils.testing import SherpaTestCase
from sherpa.optmethods import optfcts
from sherpa.optmethods import _tstoptfct


class test_optmethods(SherpaTestCase):

    def setUp(self):
        self.tolerance = 1.0e-2
        self.mc = '_montecarlo'
        self.nm = '_neldermead'
        self.lm = '_lmdif'
        self.verbose = False

    def print_result( self, name, f, x, nfev ):
        print('%s(%s) = %g in %d nfev' % (name, x, f, nfev))

    def tst_all( self, name, fct, fmin, x0, xmin, xmax,
                 iprint=False ):
        self.tst( optfcts.neldermead, name + self.nm, fct, fmin,
                  x0, xmin, xmax, iprint=iprint )
        self.tst( optfcts.montecarlo,  name + self.mc , fct, fmin,
                  x0, xmin, xmax, iprint=iprint, maxfev=8192*len(x0) )
        self.tst( optfcts.montecarlo,  name + self.mc , fct, fmin,
                  x0, xmin, xmax, iprint=iprint, maxfev=8192*len(x0),
                  numcores=_ncpus )
        self.tst( optfcts.lmdif, name + self.lm, fct, fmin,
                  x0, xmin, xmax, iprint=iprint )

    def tst( self, optmethod, name, fct, fmin, x0, xmin, xmax,
             maxfev=4096, iprint=False, numcores=1 ):
        if 1 == numcores:
            status, x, fval, msg, stuff = optmethod( fct, x0, xmin, xmax,
                                                     maxfev=maxfev*len(x0))
        else:
            status, x, fval, msg, stuff = optmethod( fct, x0, xmin, xmax,
                                                     maxfev=maxfev*len(x0),
                                                     numcores=numcores)
        nfev = stuff.get('nfev')
        if iprint:
            print('fmin = %g vs fval = %g' % ( fmin, fval ))
        if self.verbose or iprint:

            self.print_result( name, fval, x, nfev )
        self.assertEqualWithinTol( fval, fmin, self.tolerance )

    def test_rosenbrock(self):
        name = 'rosenbrock'
        npar = 4
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.rosenbrock, fmin, x0, xmin, xmax )

    def test_powell_badly_scaled(self):
        name = 'powell_badly_scaled'
        npar = 2
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst( optfcts.neldermead, name + self.nm,
                  _tstoptfct.powell_badly_scaled, fmin, x0, xmin, xmax )
        self.tst( optfcts.montecarlo, name + self.mc,
                  _tstoptfct.powell_badly_scaled, fmin, x0, xmin, xmax )

    def test_brown_badly_scaled(self):
        name = 'brown_badly_scaled'
        npar = 2
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.brown_badly_scaled, fmin, x0, xmin,
                      xmax )

    def test_beale(self):
        name = 'beale'
        npar = 2
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name,_tstoptfct.beale, fmin, x0, xmin, xmax )

    def test_jennrich_sampson(self):
        name = 'jennrich_sampson'
        npar = 2
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.jennrich_sampson, fmin, x0, xmin, xmax )

    def test_helical_valley(self):
        name = 'helical_valley'
        npar = 3
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.helical_valley, fmin, x0, xmin, xmax )

    def test_bard(self):
        name = 'bard'
        npar = 3
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.bard, fmin, x0, xmin, xmax )

    def test_gaussian(self):
        name = 'gaussian'
        npar = 3
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst( optfcts.lmdif, name + self.lm, _tstoptfct.gaussian,
                  fmin, x0, xmin, xmax )

##     # This test actually passed, there is a bug with assertEqualWithinTol
##     def test_meyer(self):
##         name = 'meyer'
##         npar = 3
##         x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
##         self.tst_all( name, _tstoptfct.meyer_fct, _tstoptfct.meyer,
##                       fmin, x0, xmin, xmax )

    def test_gulf_research_development(self):
        name = 'gulf_research_development'
        npar = 3
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.gulf_research_development,
                      fmin, x0, xmin, xmax )

    def test_box3d(self):
        name = 'box3d'
        npar = 3
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.box3d, fmin, x0, xmin, xmax )

    def test_wood(self):
        name = 'wood'
        npar = 4
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst( optfcts.neldermead, name + self.nm, _tstoptfct.wood,
                  fmin, x0, xmin, xmax )
        self.tst( optfcts.montecarlo, name + self.mc, _tstoptfct.wood,
                  fmin, x0, xmin, xmax )

    def test_kowalik_osborne(self):
        name = 'kowalik_osborne'
        npar = 4
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.kowalik_osborne, fmin, x0, xmin, xmax )

    def test_brown_dennis(self):
        name = 'brown_dennis'
        npar = 4
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.brown_dennis, fmin, x0, xmin, xmax )

##     # look at why it fails for monte carlo
    def test_biggs(self):
        name = 'biggs'
        npar = 6
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        #self.tst( optfcts.neldermead, name + self.nm, _tstoptfct.biggs_fct,
        #          fmin, x0, xmin, xmax )
        self.tst( optfcts.lmdif, name + self.lm, _tstoptfct.biggs,
                  fmin, x0, xmin, xmax )

    def test_osborne2(self):
        name = 'osborne2'
        npar = 11
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.osborne2, fmin, x0, xmin, xmax )

    def test_watson(self):
        name = 'watson'
        npar = 6
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst( optfcts.neldermead, name + self.nm, _tstoptfct.watson,
                  fmin, x0, xmin, xmax )
        self.tst( optfcts.montecarlo, name + self.mc, _tstoptfct.watson,
                  fmin, x0, xmin, xmax )

    def test_variably_dimensioned(self):
        name = 'variably_dimensioned'
        npar = 5
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.variably_dimensioned, fmin, x0, xmin,
                      xmax )

    def test_trigonometric(self):
        name = 'trigonometric'
        npar = 9
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst( optfcts.neldermead, name + self.nm,
                  _tstoptfct.trigonometric, fmin, x0, xmin, xmax )
        self.tst( optfcts.montecarlo, name + self.mc,
                  _tstoptfct.trigonometric, fmin, x0, xmin, xmax )

    def test_discrete_boundary(self):
        name = 'discrete_boundary'
        npar = 5
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.discrete_boundary, fmin, x0, xmin,
                      xmax )

    def test_discrete_integral(self):
        name = 'discrete_integral'
        npar = 5
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.discrete_integral, fmin, x0, xmin,
                      xmax )

    def test_broyden_tridiagonal(self):
        name = 'broyden_tridiagonal'
        npar = 16
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.broyden_tridiagonal, fmin, x0,
                      xmin, xmax )

    def test_broyden_banded(self):
        name = 'broyden_banded'
        npar = 18
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.broyden_banded, fmin, x0, xmin, xmax )

    def test_linear_fullrank(self):
        name = 'linear_fullrank'
        npar = 18
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.linear_fullrank, fmin, x0, xmin, xmax )

    def test_linear_fullrank1(self):
        name = 'linear_fullrank1'
        npar = 15
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.linear_fullrank1, fmin, x0, xmin, xmax )

    def test_linear_fullrank0cols0rows(self):
        name = 'linear_fullrank0cols0rows'
        npar = 13
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.linear_fullrank0cols0rows, fmin, x0,
                      xmin, xmax )

    def test_chebyquad(self):
        name = 'chebyquad'
        npar = 11
        x0, xmin, xmax, fmin = _tstoptfct.init( name, npar )
        self.tst_all( name, _tstoptfct.chebyquad, fmin, x0, xmin, xmax )


def tstme():
    from sherpa.utils import SherpaTest
    import sherpa.optmethods
    SherpaTest(sherpa.optmethods).test()

#
# pushd ../../../ ; makeme ; popd ; python test_optmethods.py
#
if __name__ == '__main__':
    tstme()
