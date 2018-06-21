#
<<<<<<< HEAD
#  Copyright (C) 2010, 2016, 2018, 2019, 2020  Smithsonian Astrophysical Observatory
=======
#  Copyright (C) 2010, 2016, 2018  Smithsonian Astrophysical Observatory
>>>>>>> an initial release of simultaneous fit on multicores (slower for most, ie a lot, of cases :)
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

import pytest
import numpy
from numpy.testing import assert_allclose, assert_equal

from sherpa import utils
from sherpa.utils import SherpaFloat, NoNewAttributesAfterInit
from sherpa.utils.testing import SherpaTestCase


class test_utils(SherpaTestCase):

    def setUp(self):
        def f1(a, b, c):
            pass
        self.f1 = f1

        def f2(a, b=1, c=2, d=3, e=4):
            pass
        self.f2 = f2

        def f3(a=None, b=1, c=2, d=3, e=4):
            pass
        self.f3 = f3

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
        self.assertFalse(hasattr(c, 'z'))
        self.assertRaises(AttributeError, delattr, c, 'x')
        self.assertRaises(AttributeError, delattr, c, 'z')
        self.assertRaises(AttributeError, setattr, c, 'z', 5)
        c.x = 4
        self.assertEqual(c.x, 4)

    def test_export_method(self):
        class C():
            def m(self, x, y=2):
                'Instance method m()'
                return x * y

            def margs(self, x, y=2, *args):
                'Instance method margs() with *args'
                return x * y + len(args)

            def kwargs(self, x, y=2, **kwargs):
                'Instance method kwargs() with **kwargs'
                return x * y + 2 * len(kwargs)

            def bargs(self, x, y=2, *args, **kwargs):
                'Instance method bargs() with *args and **kwargs'
                return x * y + len(args) + 2 * len(kwargs)

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
        for meth in (c.m, c.margs, c.kwargs, c.bargs, c.cm, c.sm):
            m = utils.export_method(meth)
            self.assertEqual(m.__name__, meth.__name__)
            self.assertTrue(m.__doc__ is not None)
            self.assertEqual(m.__doc__, meth.__doc__)
            self.assertEqual(m(3), 6)
            self.assertEqual(m(3, 7), 21)
            e = None
            try:
                m()
            except TypeError as e:
                emsg2 = "{}() ".format(meth.__name__) + \
                        "takes at least 1 argument (0 given)"
                emsg3 = "{}() ".format(meth.__name__) + \
                        "missing 1 required positional argument: 'x'"
                self.assertIn(str(e), [emsg2, emsg3])

        # Check that *args/**kwargs are handled correctly for methods;
        # should perhaps be included above to avoid repeated calls
        # to export_method?
        #
        meth = utils.export_method(c.margs)
        self.assertTrue(meth(3, 7, "a", "b"), 23)
        try:
            meth(12, dummy=None)
        except TypeError as e:
            emsg = "margs() got an unexpected keyword argument 'dummy'"
            self.assertEqual(str(e), emsg)

        meth = utils.export_method(c.kwargs)
        self.assertTrue(meth(3, 7, foo="a", bar="b"), 25)
        try:
            meth(12, 14, 15)
        except TypeError as e:
            emsg2 = "kwargs() takes at most 2 arguments (3 given)"
            emsg3 = "kwargs() takes from 1 to 2 positional arguments " + \
                    "but 3 were given"
            self.assertIn(str(e), [emsg2, emsg3])

        meth = utils.export_method(c.bargs)
        self.assertTrue(meth(3, 7, 14, 15, foo=None), 25)

        # Non-method argument
        def f(x):
            return 2 * x
        self.assertTrue(utils.export_method(f) is f)

        # Name and module name
        m = utils.export_method(c.m, 'foo', 'bar')
        self.assertEqual(m.__name__, 'foo')
        self.assertEqual(m.__module__, 'bar')
        self.assertEqual(m(3), 6)
        self.assertEqual(m(3, 7), 21)

    def test_get_num_args(self):
        self.assertEqual(utils.get_num_args(self.f1), (3, 3, 0))
        self.assertEqual(utils.get_num_args(self.f2), (5, 1, 4))
        self.assertEqual(utils.get_num_args(self.f3), (5, 0, 5))

    def test_get_keyword_names(self):
        self.assertEqual(utils.get_keyword_names(self.f1), [])
        l = ['b', 'c', 'd', 'e']
        self.assertEqual(utils.get_keyword_names(self.f2), l)
        self.assertEqual(utils.get_keyword_names(self.f2, 2), l[2:])
        self.assertEqual(utils.get_keyword_names(self.f2, 7), [])
        l = ['a', 'b', 'c', 'd', 'e']
        self.assertEqual(utils.get_keyword_names(self.f3), l)
        self.assertEqual(utils.get_keyword_names(self.f3, 1), l[1:])
        self.assertEqual(utils.get_keyword_names(self.f3, 7), [])

    def test_get_keyword_defaults(self):
        self.assertEqual(utils.get_keyword_defaults(self.f1), {})
        d = {'b': 1, 'c': 2, 'd': 3, 'e': 4}
        self.assertEqual(utils.get_keyword_defaults(self.f2), d)
        del d['b']
        del d['c']
        self.assertEqual(utils.get_keyword_defaults(self.f2, 2), d)
        self.assertEqual(utils.get_keyword_defaults(self.f2, 7), {})
        d = {'a': None, 'b': 1, 'c': 2, 'd': 3, 'e': 4}
        self.assertEqual(utils.get_keyword_defaults(self.f3), d)
        del d['a']
        self.assertEqual(utils.get_keyword_defaults(self.f3, 1), d)
        self.assertEqual(utils.get_keyword_defaults(self.f3, 7), {})

    def test_print_fields(self):
        names = ['a', 'bb', 'ccc']
        vals = {'a': 3, 'bb': 'Ham', 'ccc': numpy.array([1.0, 2.0, 3.0])}
        self.assertEqual(utils.print_fields(names, vals),
                         'a   = 3\nbb  = Ham\nccc = Float64[3]')

    def test_calc_total_error(self):
        stat = numpy.array([1, 2])
        sys = numpy.array([3, 4])

        self.assertEqual(utils.calc_total_error(None, None), None)
        assert_equal(utils.calc_total_error(stat, None), stat)
        assert_equal(utils.calc_total_error(None, sys), sys)

        # Unlike the above tests, only look for equivalence within
        # a tolerance, since the numbers are manipulated rather than
        # copied in this case (although the equation should be the same
        # so the old approach of using equality should actually be okay
        # here).
        ans = numpy.sqrt(stat * stat + sys * sys)
        assert_allclose(utils.calc_total_error(stat, sys), ans)

    def test_poisson_noise(self):
        out = utils.poisson_noise(1000)
        self.assertEqual(type(out), SherpaFloat)
        self.assertGreater(out, 0.0)

        for x in (-1000, 0):
            out = utils.poisson_noise(x)
            self.assertEqual(type(out), SherpaFloat)
            self.assertEqual(out, 0.0)

        out = utils.poisson_noise([1001, 1002, 0.0, 1003, -1004])
        self.assertEqual(type(out), numpy.ndarray)
        self.assertEqual(out.dtype.type, SherpaFloat)
        ans = numpy.flatnonzero(out > 0.0)
        assert_equal(ans, numpy.array([0, 1, 3]))

        self.assertRaises(ValueError, utils.poisson_noise, 'ham')
        self.assertRaises(TypeError, utils.poisson_noise, [1, 2, 'ham'])

    def test_neville(self):
        func = numpy.exp
        tol = 1.0e-6
        num = 10
        # TODO: can we not just use vectorized code here?
        x = []
        y = []
        for ii in range(num):
            x.append(ii / float(num))
            y.append(func(x[ii]))
        xx = numpy.array(x)
        yy = numpy.array(y)
        for ii in range(num):
            tmp = 1.01 * (ii / float(num))
            answer = func(tmp)
            val = utils.neville(tmp, xx, yy)
            self.assertTrue(utils.Knuth_close(answer, val, tol))

    def test_neville2d(self):
        funcx = numpy.sin
        funcy = numpy.exp
        nrow = 10
        ncol = 10
        tol = 1.0e-4
        # TODO: As with test_neville; can this not be simplified with
        # vectorized code
        x = numpy.zeros((nrow, ))
        y = numpy.zeros((ncol, ))
        fval = numpy.empty((nrow, ncol))
        row_tmp = numpy.pi / nrow
        # col_tmp = 1.0 / float(ncol)
        for row in range(nrow):
            x[row] = (row + 1.0) * row_tmp
            for col in range(ncol):
                y[col] = (col + 1.0) / float(ncol)
                fval[row][col] = funcx(x[row]) * funcy(y[col])

        for row in range(ncol):
            xx = (-0.1 + (row + 1.0) / float(nrow)) * numpy.pi
            for col in range(4):
                yy = -0.1 + (col + 1.0) / float(ncol)
                answer = funcx(xx) * funcy(yy)
                val = utils.neville2d(xx, yy, x, y, fval)
                self.assertTrue(utils.Knuth_close(answer, val, tol))


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
        f = numpy.sum
        iterable = [numpy.arange(1, 2 + 2 * i) for i in range(num_segments)]

        result = list(map(f, iterable))
        result = numpy.asarray(result)

        pararesult = utils.parallel_map(f, iterable, num_tasks)

        assert_equal(result, numpy.asarray(pararesult))

<<<<<<< HEAD

@pytest.mark.parametrize("los, his, axis", [([], [], [0,1,2,3,4]),
                                            ([], [1], [0,1,2,3,4]),
                                            ([1], [], [0,1,2,3,4]),
                                            ([], [], [])])
def test_filter_bins_empty(los, his, axis):
    """Ensure filter_bins returns None if one input is empty."""

    # just check the input parameters include an empty array
    assert len(los) == 0 or len(his) == 0 or len(axis) == 0

    assert utils.filter_bins(los, his, [axis]) is None


def test_filter_bins_scalar_array_empty():
    """Edge case: do we care about this result?"""

    f = utils.filter_bins([1], [2], [[]])
    assert f.dtype == numpy.bool
    assert len(f) == 0


@pytest.mark.parametrize("axval, flag", [(0, False),
                                         (0.99, False),
                                         (1 - 1e-7, True),
                                         (1, True),
                                         (2.34, True),
                                         (5, True),
                                         (5 + 1e-7, True),
                                         (5.01, False),
                                         (7, False)])
def test_filter_bins_scalar_array(axval, flag):
    """Edge case: do we care about this result?"""

    f = utils.filter_bins([1], [5], [axval])
    assert f == flag


@pytest.mark.parametrize("lo, hi, res",
                         [(0, 10, [True] * 5),
                          (1, 5, [True] * 5),
                          (2, 5, [False, True, True, True, True]),
                          (2, None, [False, True, True, True, True]),
                          (1, 4, [True, True, True, True, False]),
                          (None, 4, [True, True, True, True, False]),
                          (1.1, 4.9, [False, True, True, True, False]),
                          (2, 4, [False, True, True, True, False]),
                          # Have minimum > maximum, so nothing is selected
                          (4, 3, [False, False, False, False, False]),
                          # Have minimum = maximum = bin value
                          (4, 4, [False, False, False, True, False]),
                          # Have minimum = maximum, not equal to a bin value
                          (3.1, 3.1, [False, False, False, False, False])
                         ])
def test_filter_bins_one(lo, hi, res):
    """Can we filter the array between lo and hi?"""

    dvals = numpy.asarray([1, 2, 3, 4, 5])
    flags = utils.filter_bins([lo], [hi], [dvals])
    assert (flags == res).all()


def test_filter_bins_two_none():
    """Use two different arrays for filtering with no filter.
    """

    y1 = [1, 2, 3, 4, 5]
    y2 = [10, 20, 30, 40, 50]

    flags = utils.filter_bins((None, None), (None, None), (y1, y2))
    assert flags is None


@pytest.mark.parametrize("lo1, lo2, hi1, hi2, expected",
                         [ (1.5, 21, 3.6, 44, [False, False, True, False, False]),
                           (1.5, None, None, 44, [False, True, True, True, False]),
                           (None, None, None, 44, [True, True, True, True, False]),
                           (1.5, None, None, None, [False, True, True, True, True])
                         ])
def test_filter_bins_two(lo1, lo2, hi1, hi2, expected):
    """Use two different arrays for filtering.

    This version uses tuples rather than arrays as the
    input arguments to filter_bins.

    """

    y1 = [1, 2, 3, 4, 5]
    y2 = [10, 20, 30, 40, 50]

    flags = utils.filter_bins((lo1, lo2), (hi1, hi2), (y1, y2))
    assert len(flags) == 5
    for i in range(5):
        assert flags[i] == expected[i], i


def test_filter_bins_unordered():
    """What happens if the array is unordered?"""

    flags = utils.filter_bins((3, ), (8, ), [[1,4,3,7,8,10,5]])

    expected = [False, True, True, True, True, False, True]

    assert len(flags) == len(expected)
    for got, exp in zip(flags, expected):
        assert got == exp
=======
def test_paralle_map_funcs():
    for arg in range(1, 5):
        func = [numpy.sum]
        funcs = arg * func
        datas = [numpy.arange(1, 2+2*i) for i in range(arg)]
        result = []
        for func, data in zip(funcs, datas):
            result.extend(numpy.asarray(list(map(func, data))))
        assert_equal(result, utils.parallel_map_funcs(funcs, datas, arg))
>>>>>>> an initial release of simultaneous fit on multicores (slower for most, ie a lot, of cases :)
