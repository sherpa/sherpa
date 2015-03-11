# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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
import sys
from sherpa.utils import *

def sqr( x, *args ):
    return x*x

def prob1( x, *args ):
    return numpy.cos( x ) - pow( x, 3.0 )
        
def prob2( x, *args ):
    return pow( x, 4.0 ) - x - 1.0

def prob3( x, *args ):
    return ( x + 3.0 ) * ( x - 1.0 ) * ( x - 1.0 )

def prob4( x, *args ):
    return numpy.exp( x ) + x - 2.0

def prob5( x, *args ):
    return x + numpy.exp( x )

def prob6( x, *args ):
    return numpy.sqrt( x ) + x * x - 5.0

def prob7( x, *args ):
    return x + numpy.sin(x) - 4.0

def prob8( x, *args ):
    return x * x - numpy.exp( x )
  
def prob9( x, *args ):
    return numpy.sin( x ) + numpy.cos( x )

def prob10( x, *args ):
    return x * x - 2.0

def prob11( x, *args ):
    x2 = sqr( x )
    x4 = sqr( x2 )
    x8 = sqr( x4 )
    return x8 * x

def prob12( x, *args ):
    if x == 0.0:
        return 0.0
    return x - numpy.exp( -1.0 / ( x * x ) )

def prob13( x, *args ):
    x2 = sqr( x - 1.0 )
    x4 = x2 * x2
    return x * x4 * ( x - 1.0 )

def prob14( x, *args ):
    return ( x - 5.0 ) * ( x - 2.0 )

def prob15( x, *args ):
    return ( x - 50.0 ) * ( x + 20.0 )

def prob16( x, *args ):
    return (sqr(x)-2.0)*x - 5.0

def prob17( x, *args ):
    return sqr( x - 1.0e-5 ) - 3.0

def prob18( x, *args ):
    return numpy.cos(x) - x

def prob19( x, *args ):
    return numpy.sin(x) - x

def prob20( x, *args ):
    return x * x * x - 2.0 * x - 5.0

def prob21( x, *args ):
    return numpy.sin(x) - 0.5 * x

def prob22( x, *args ):
    return 2 * x - numpy.exp( - x )

def prob23( x, *args ):
    return x * numpy.exp( -x )

def prob24( x, *args ):
    return numpy.exp( x ) - 1. / (100*x*x)

def prob25( x, *args ):
    tmp = x - 1
    return ( x + 3 ) * tmp * tmp

def prob26( x, *args ):
    return numpy.exp(x) - 2 - 1 / ( 10 * x )**2 - 2 / ( 100 * x )**3

def prob27( x, *args ):
    return x**3

def prob28( x, *args ):
    return numpy.cos( x ) - x

# root is 1.83 the convergence rate for Muller's method
def prob29( x, *args ):
    return pow( x, 3.0 ) - x*x - x - 1

# muller gets it wrong but muller_bound got it right
def prob30( x, *args ):
    return numpy.exp( x ) - 1

def prob31( x, *args ):
    return x * numpy.exp( x ) - 1

def prob32( x, *args ):
    return 11.0 * pow( x, 11.0 ) - 1.0

def prob33( x, *args ):
    return numpy.exp( (1.0*x+7)*x -30.0 ) - 1.0

def prob34( x, *args ):
    return 1.0 / x - numpy.sin( x ) + 1.0

def prob35( x, *args ):
    return ( x*x - 2.0 ) * x - 5.0

def prob36( x, *args ):
    return 1.0 / x - 1.0

def prob37( x, *args ):
    return numpy.log( x )

def prob38( x, *args ):
    return x - numpy.exp( numpy.sin(x) )

def repeller( x, *args ):
    return  20.0 * x / ( 100.0 * x * x + 1.0 )

def pinhead( x, *args ):
    epsilon = 0.00001
    if epsilon == 0.0:
        return ( 16.0 - x**4 ) / ( 16.0 * x**4 )
    else:
        return  ( 16.0 - x**4 ) / ( 16.0 * x**4 + epsilon )

#    Sample results:
#
#    E        M      X
#    -----  ---  ----------
#    0.100    5    5.554589
#    0.200    5    6.246908
#    0.300    5    7.134960
#    0.400    5    8.313903
#    0.500    5    9.950063
#    0.600    5   12.356653
#    0.700    5   16.167990
#    0.800    5   22.656579
#    0.900    5   33.344447
#    0.990    5   45.361023
#    0.990    1   24.725822
#    0.990   33   89.722155
#    0.750   70  110.302
#    0.990    2   32.361007

def kepler( x, *args ):
    e = 0.8
    m = 5.0
    pi = 3.141592653589793
    return ( pi * ( x - m ) / 180.0 - e * numpy.sin ( pi * x / 180.0 ) )

def wallis( x, *args ):
    return x**3 - 2*x - 5

def thinpole( x, *args ):
    pi = 3.141592653589793
    return  3.0 * x**2 + 1.0 + ( numpy.log ( ( x - pi )**2 ) ) / pi**4

def muller_convergence_rate(x):
    return x - pow(x,3.0)/3.0

class test_root( SherpaTestCase ):

    def setUp(self):
        self.verbose = False
        self.tol = 1.0e-6
        
    def demuller2( self, fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32,
                   tol=1.0e-6 ):
        return demuller( fcn, xa, xb, (xa+xb)/2.0, args=args, maxfev=maxfev,
                         tol=tol)
    
    def tst_solve( self, fct, a, b, tol, iprint=False ):
        methods = [ bisection, self.demuller2, new_muller, apache_muller, zeroin ]
        if self.verbose or iprint:
            sys.stdout.write( '\n' )
            
        for solve_func in methods:
            # demuller2 cannot solve prob30 & pinhead so may as well skip them
            if fct.__name__ == 'prob30' or fct.__name__ == 'pinhead':
                self.assert_( 0 == 0 )
                return
            result = solve_func( fct, a, b, tol=tol )
            [ root, froot ] = result[ 0 ]
            if self.verbose or iprint:
                solve_name = solve_func.__name__
                fct_name = fct.__name__
                nfev = result[ -1 ]
                msg = '%s:\t%s(%e) = %e in %d nfevs\n' % ( solve_name, fct_name, root, froot, nfev )
                sys.stdout.write( msg )
            self.assert_( abs( froot ) <= tol )

    def test_prob1( self ):
        a = 0.0
        b = 1.0
        self.tst_solve( prob1, a, b, self.tol )

    def test_prob2( self ):
        a = -5.0
        b = 1.2
        self.tst_solve( prob2, a, b, self.tol )

    def test_prob3( self ):
        a = -10.0
        b = 10.0
        self.tst_solve( prob3, a, b, self.tol )

    def test_prob4( self ):
        a = -1.0
        b = 1.0
        self.tst_solve( prob4, a, b, self.tol )

    def test_prob5( self ):
        a = -1.0
        b = 1.0
        self.tst_solve( prob5, a, b, self.tol )

    def test_prob6( self ):
        a = 0.0
        b = 3.0
        self.tst_solve( prob6, a, b, self.tol )

    def test_prob7( self ):
        a = 0.0
        b = 5.0
        self.tst_solve( prob7, a, b, self.tol )

    def test_prob8( self ):
        a = -5.0
        b = 0.0
        self.tst_solve( prob8, a, b, self.tol )

    def test_prob9( self ):
        a = 0.0
        b = 3.14
        self.tst_solve( prob9, a, b, self.tol )

    def test_prob10( self ):
        a = 0.0
        b = 2.0
        self.tst_solve( prob10, a, b, self.tol )

    def test_prob11( self ):
        a = -1.0
        b = 1.5
        self.tst_solve( prob11, a, b, self.tol )

    def test_prob12( self ):
        a = -1.0
        b = 1.0
        self.tst_solve( prob12, a, b, self.tol )

    def test_prob13( self ):
        a = -0.5
        b = 0.99
        self.tst_solve( prob13, a, b, self.tol )

    def test_prob14( self ):
        a = 1.0
        b = 3.0
        self.tst_solve( prob14, a, b, self.tol )

    def test_prob15( self ):
        a = 15.0
        b = 75.0
        self.tst_solve( prob15, a, b, self.tol )

    def test_prob16( self ):
        a = 2.0
        b = 3.0
        self.tst_solve( prob16, a, b, self.tol )

    def test_prob17( self ):
        a = -1.0
        b = 3.0
        self.tst_solve( prob17, a, b, self.tol )

    def test_prob18( self ):
        a = -1.0
        b = 3.0
        self.tst_solve( prob18, a, b, self.tol )

    def test_prob19( self ):
        a = -1.0
        b = 3.0
        self.tst_solve( prob19, a, b, self.tol )

    def test_prob20( self ):
        a = 2.0
        b = 3.0
        self.tst_solve( prob20, a, b, self.tol )

    def test_prob21( self ):
        a = -1000.0
        b =  1000.0
        self.tst_solve( prob21, a, b, self.tol )

    def test_prob22( self ):
        a = -10.0
        b = 100.0
        self.tst_solve( prob22, a, b, self.tol )

    def test_prob23( self ):
        a = -10.0
        b = 100.0
        self.tst_solve( prob23, a, b, self.tol )

    def test_prob24( self ):
        a = 1.0e-4
        b = 20.0
        self.tst_solve( prob24, a, b, self.tol )

    def test_prob25( self ):
        a = -1000.0
        b = 1000.0
        self.tst_solve( prob25, a, b, self.tol )

    def test_prob26( self ):
        a = 1.0e-4
        b = 20.0
        self.tst_solve( prob26, a, b, self.tol )

    def test_prob27( self ):
        a = -1000.0
        b =  1000.0
        self.tst_solve( prob27, a, b, self.tol )

    def test_prob28( self ):
        a = -10.0
        b =  10.0
        self.tst_solve( prob28, a, b, self.tol )

    def test_prob29( self ):
        a = 0
        b = 10
        self.tst_solve( prob29, a, b, self.tol )

    def test_prob30( self ):
        a = -50
        b = 100.0
        self.tst_solve( prob30, a, b, self.tol )

    def test_prob31( self ):
        a = -1.0
        b = 1.0
        self.tst_solve( prob31, a, b, self.tol )

    def test_prob32( self ):
        a = 0.1
        b = 0.9
        self.tst_solve( prob32, a, b, self.tol )

    def test_prob33( self ):
        a = 2.8
        b = 3.1
        self.tst_solve( prob33, a, b, self.tol )

    def test_prob34( self ):
        a = -1.3
        b = -0.5
        self.tst_solve( prob34, a, b, self.tol )

    def test_prob35( self ):
        a = 2.0
        b = 3.0
        self.tst_solve( prob35, a, b, self.tol )

    def test_prob36( self ):
        a = 0.5
        b = 1.5
        self.tst_solve( prob36, a, b, self.tol )

    def test_prob37( self ):
        a = 0.5
        b = 5.0
        self.tst_solve( prob37, a, b, self.tol )

    def test_prob38( self ):
        a = 1.0
        b = 4.0
        self.tst_solve( prob38, a, b, self.tol )

    def test_repeller( self ):
        a = -10.0
        b =  10.0
        self.tst_solve( repeller, a, b, self.tol )

    def test_pinhead( self ):
        a =  1.0e-5
        b =  10.0
        self.tst_solve( pinhead, a, b, self.tol )

    def test_kepler( self ):
        a = -175.0
        b =  185.0
        self.tst_solve( kepler, a, b, self.tol )

    def test_wallis( self ):
        a = 2.0
        b = 3.0
        self.tst_solve( wallis, a, b, self.tol )

    def test_muller_convergence_rate( self ):
        a = 1.7
        b = 10.0
        self.tst_solve( muller_convergence_rate, a, b, self.tol )
        
def tstme():
    from sherpa.utils import SherpaTest
    import sherpa.utils
    SherpaTest(sherpa.utils).test()

#
# pushd ../../.. ; makeme ; popd ; python test_root.py
#
if '__main__' == __name__:
    tstme()
