# 
#  Copyright (C) 2010  Smithsonian Astrophysical Observatory
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
from sherpa.utils import *
from sherpa.utils import SherpaFloat

class test_utils(SherpaTestCase):

    def setUp(self):
        def f1(a, b, c):
            pass
        self.f1 = f1

        def f2(a, b=1, c=2, d=3, e=4):
            pass
        self.f2 = f2

    def test_NoNewAttributesAfterInit(self):
        class C(NoNewAttributesAfterInit):
            def __init__(self):
                self.x = 1
                self.y = 2
                self.z = 3
                del self.z
                NoNewAttributesAfterInit.__init__(self)

        c = C()
        self.assertEqual(c.x, 1)
        self.assertEqual(c.y, 2)
        self.assert_(not hasattr(c, 'z'))
        self.assertRaises(AttributeError, delattr, c, 'x')
        self.assertRaises(AttributeError, delattr, c, 'z')
        self.assertRaises(AttributeError, setattr, c, 'z', 5)
        c.x = 4
        self.assertEqual(c.x, 4)

    def test_export_method(self):
        class C(object):
            def m(self, x, y=2):
                'Instance method m()'
                return x * y

            @classmethod
            def cm(klass, x, y=2):
                'Class method cm()'
                return x * y

            @staticmethod
            def sm(x, y=2):
                'Static method sm()'
                return x * y

        c = C()

        # Basic usage
        for meth in (c.m, c.cm, c.sm):
            m = export_method(meth)
            self.assertEqual(m.__name__, meth.__name__)
            self.assert_(m.__doc__ is not None)
            self.assertEqual(m.__doc__, meth.__doc__)
            self.assertEqual(m(3), 6)
            self.assertEqual(m(3, 7), 21)
            e = None
            try:
                m()
            except TypeError, e:
                pass
            self.assertEqual(str(e),
                             '%s() takes at least 1 argument (0 given)' %
                             meth.__name__)

        # Non-method argument
        def f(x):  return 2*x
        self.assert_(export_method(f) is f)

        # Name and module name
        m = export_method(c.m, 'foo', 'bar')
        self.assertEqual(m.__name__, 'foo')
        self.assertEqual(m.__module__, 'bar')
        self.assertEqual(m(3), 6)
        self.assertEqual(m(3, 7), 21)

    def test_get_keyword_names(self):
        self.assertEqual(get_keyword_names(self.f1), [])
        l = ['b', 'c', 'd', 'e']
        self.assertEqual(get_keyword_names(self.f2), l)
        self.assertEqual(get_keyword_names(self.f2, 2), l[2:])

    def test_get_keyword_defaults(self):
        self.assertEqual(get_keyword_defaults(self.f1), {})
        d = {'b':1, 'c':2, 'd':3, 'e':4}
        self.assertEqual(get_keyword_defaults(self.f2), d)
        del d['b']
        del d['c']
        self.assertEqual(get_keyword_defaults(self.f2, 2), d)

    def test_print_fields(self):
        names = ['a', 'bb', 'ccc']
        vals = {'a': 3, 'bb': 'Ham', 'ccc': numpy.array([1.0, 2.0, 3.0])}
        self.assertEqual(print_fields(names, vals),
                         'a   = 3\nbb  = Ham\nccc = Float64[3]')

    def test_calc_total_error(self):
        stat = numpy.array([1,2])
        sys = numpy.array([3,4])
        self.assert_(calc_total_error(None, None) is None)
        self.assert_(calc_total_error(stat, None) is stat)
        self.assert_(calc_total_error(None, sys) is sys)
        self.assert_(numpy.all(calc_total_error(stat, sys) ==
                               numpy.sqrt(stat*stat + sys*sys)))

    def test_poisson_noise(self):
        out = poisson_noise(1000)
        self.assert_(type(out) is SherpaFloat)
        self.assert_(out > 0.0)

        for x in (-1000, 0):
            out = poisson_noise(x)
            self.assert_(type(out) is SherpaFloat)
            self.assert_(out == 0.0)

        out = poisson_noise([1001, 1002, 0.0, 1003, -1004])
        self.assert_(type(out) is numpy.ndarray)
        self.assert_(out.dtype.type is SherpaFloat)
        self.assert_(numpy.all(numpy.flatnonzero(out > 0.0) ==
                               numpy.array([0, 1, 3])))

        self.assertRaises(ValueError, poisson_noise, 'ham')
        self.assertRaises(TypeError, poisson_noise, [1, 2, 'ham'])


    def test_neville( self ):
        func = numpy.exp
        tol = 1.0e-6
        num = 10
        x = []
        y = []
        for ii in xrange( num ):
            x.append( ii / float( num ) )
            y.append( func( x[ ii ] ) )
        xx = numpy.array( x )
        yy = numpy.array( y )
        for ii in xrange( num ):
            tmp = 1.01 * ( ii/ float( num ) )
            answer = func( tmp )
            val = neville( tmp, xx, yy )
            self.assert_( Knuth_close( answer, val, tol ) )

    def test_neville2d( self ):
        funcx = numpy.sin
        funcy = numpy.exp
        nrow = 10
        ncol = 10
        tol = 1.0e-4
        x = numpy.zeros( (nrow,) )
        y = numpy.zeros( (ncol,) )
        fval = numpy.empty( ( nrow, ncol ) )
        row_tmp = numpy.pi / nrow
        col_tmp = 1.0 / float( ncol )
        for row in xrange( nrow ):
            x[ row ] = ( row + 1.0 ) * row_tmp
            for col in xrange( ncol ):
                y[ col ] = ( col + 1.0 ) / float( ncol )
                fval[ row ][ col ] = funcx( x[ row ] ) * funcy( y[ col ] )

        for row in xrange( ncol ):
            xx = ( -0.1 + ( row + 1.0 ) / float( nrow ) ) * numpy.pi
            for col in xrange( 4 ):
                yy = -0.1 +( col + 1.0 )/ float( ncol )
                answer = funcx( xx ) * funcy( yy )
                val = neville2d( xx, yy, x, y, fval )
                self.assert_( Knuth_close( answer, val, tol ) )
    
    def test_parallel_map(self):

        ncpus = 1
        try:
            import multiprocessing
            ncpus = multiprocessing.cpu_count()
            Pool = multiprocessing.Pool
        except:
            return

        numtasks = 8
        size = (64,64)
        vals = numpy.random.rand(*size)
        f = numpy.linalg.eigvals
        iterable = [vals]*numtasks

        result = map(f, iterable)

        pararesult = parallel_map(f, iterable, ncpus)

        pool = Pool(ncpus)
        poolresult = pool.map(f, iterable)

        self.assert_((numpy.asarray(result) == numpy.asarray(pararesult)).all())
        self.assert_((numpy.asarray(result) == numpy.asarray(poolresult)).all())



if __name__ == '__main__':

    import sherpa.utils as utils
    SherpaTest(utils).test()
