#
#  Copyright (C) 2010, 2016, 2018 - 2024
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

from io import StringIO

import numpy

import pytest

from sherpa import utils
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.utils.err import IOErr


def f1(a, b, c):
    pass


def f2(a, b=1, c=2, d=3, e=4):
    pass


def f3(a=None, b=1, c=2, d=3, e=4):
    pass


def test_NoNewAttributesAfterInit():
    class C(NoNewAttributesAfterInit):
        def __init__(self):
            self.x = 1
            self.y = 2
            self.z = 3
            del self.z
            NoNewAttributesAfterInit.__init__(self)

    c = C()

    assert c.x == 1
    c.x = 4
    assert c.x == 4

    assert c.y == 2
    assert not hasattr(c, 'z')

    # attribute that exists
    with pytest.raises(AttributeError):
        del c.x

    # attribute that not exist
    with pytest.raises(AttributeError):
        del c.z

    with pytest.raises(AttributeError):
        c.z = 5


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


def test_export_method_basic():

    c = C()

    # Basic usage. The error message depends on
    #  a) Python version
    #  b) what method is being wrapped
    # (before Python 3.10 it didn't).
    #
    for meth in (c.m, c.margs, c.kwargs, c.bargs, c.cm, c.sm):
        m = utils.export_method(meth)

        assert m.__name__ == meth.__name__
        assert m.__doc__ is not None
        assert m.__doc__ == meth.__doc__

        assert m(3) == 6
        assert m(3, 7) == 21

        with pytest.raises(TypeError) as exc:
            m()

        emsg = f"{meth.__name__}() " + \
            "missing 1 required positional argument: 'x'"

        if meth.__name__ == 'sm':
            # In Python 3.10 we see C.sm rather than sm
            # so we search for both. We could be more explicit
            # and check on the Python version (e.g. for >= 3.10)
            # but it doesn't seem worth it. It's interesting it's
            # only for the static method.
            #
            assert str(exc.value) in [emsg, 'C.' + emsg]
        else:
            assert str(exc.value) == emsg


def test_export_method_args_call():

    # Check that *args/**kwargs are handled correctly for methods;
    # should perhaps be included above to avoid repeated calls
    # to export_method?
    #
    c = C()
    meth = utils.export_method(c.margs)
    assert meth(3, 7, "a", "b") == 23

    meth = utils.export_method(c.bargs)
    assert meth(3, 7, 14, 15, foo=None) == 25


def test_export_method_args_errors():

    c = C()
    meth = utils.export_method(c.margs)

    with pytest.raises(TypeError) as exc:
        meth(12, dummy=None)

    emsg = "margs() got an unexpected keyword argument 'dummy'"
    assert str(exc.value) == emsg

    meth = utils.export_method(c.kwargs)
    assert meth(3, 7, foo="a", bar="b" == 25)

    with pytest.raises(TypeError) as exc:
        meth(12, 14, 15)

    emsg = "kwargs() takes from 1 to 2 positional arguments " + \
        "but 3 were given"
    assert str(exc.value) in emsg


def test_export_method_non_method():

    # Non-method argument
    def f(x):
        return 2 * x

    assert utils.export_method(f) is f


def test_export_method_names():

    # Name and module name
    c = C()
    m = utils.export_method(c.m, 'foo', 'bar')
    assert m.__name__ == 'foo'
    assert m.__module__ == 'bar'
    assert m(3) == 6
    assert m(3, 7) == 21


@pytest.mark.parametrize("func,expected",
                         [(f1, (3, 3, 0)),
                          (f2, (5, 1, 4)),
                          (f3, (5, 0, 5))])
def test_get_num_args(func, expected):
    assert utils.get_num_args(func) == expected


def test_get_keyword_names_f1():
    assert utils.get_keyword_names(f1) == []


def test_get_keyword_names_f2():
    l = ['b', 'c', 'd', 'e']
    assert utils.get_keyword_names(f2) == l
    assert utils.get_keyword_names(f2, 2) == l[2:]
    assert utils.get_keyword_names(f2, 7) == []


def test_get_keyword_names_f3():
    l = ['a', 'b', 'c', 'd', 'e']
    assert utils.get_keyword_names(f3) == l
    assert utils.get_keyword_names(f3, 1) == l[1:]
    assert utils.get_keyword_names(f3, 7) == []


def test_get_keyword_defaults_f1():
    assert utils.get_keyword_defaults(f1) == {}


def test_get_keyword_defaults_f2():
    d = {'b': 1, 'c': 2, 'd': 3, 'e': 4}
    assert utils.get_keyword_defaults(f2) == d

    del d['b']
    del d['c']
    assert utils.get_keyword_defaults(f2, 2) == d
    assert utils.get_keyword_defaults(f2, 7) == {}


def test_get_keyword_defaults_f3():
    d = {'a': None, 'b': 1, 'c': 2, 'd': 3, 'e': 4}
    assert utils.get_keyword_defaults(f3) == d

    del d['a']
    assert utils.get_keyword_defaults(f3, 1) == d
    assert utils.get_keyword_defaults(f3, 7) == {}


def test_print_fields():
    names = ['a', 'bb', 'ccc']
    vals = {'a': 3, 'bb': 'Ham', 'ccc': numpy.array([1.0, 2.0, 3.0])}

    out = 'a   = 3\nbb  = Ham\nccc = Float64[3]'
    assert utils.print_fields(names, vals) == out


def test_calc_total_error():
    stat = numpy.array([1, 2])
    sys = numpy.array([3, 4])

    assert utils.calc_total_error(None, None) is None
    assert (utils.calc_total_error(stat, None) == stat).all()
    assert (utils.calc_total_error(None, sys) == sys).all()

    # Unlike the above tests, only look for equivalence within
    # a tolerance, since the numbers are manipulated rather than
    # copied in this case (although the equation should be the same
    # so the old approach of using equality should actually be okay
    # here).
    ans = numpy.sqrt(stat * stat + sys * sys)
    assert utils.calc_total_error(stat, sys) == pytest.approx(ans)


def test_neville():
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
        assert utils.Knuth_close(answer, val, tol)


def test_neville2d():
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
            assert utils.Knuth_close(answer, val, tol)


@pytest.mark.parametrize("los, his, axis", [([], [], [0, 1, 2, 3, 4]),
                                            ([], [1], [0, 1, 2, 3, 4]),
                                            ([1], [], [0, 1, 2, 3, 4]),
                                            ([], [], [])])
def test_filter_bins_empty(los, his, axis):
    """Ensure filter_bins returns None if one input is empty."""

    # just check the input parameters include an empty array
    assert len(los) == 0 or len(his) == 0 or len(axis) == 0

    assert utils.filter_bins(los, his, [axis]) is None


def test_filter_bins_scalar_array_empty():
    """Edge case: do we care about this result?"""

    f = utils.filter_bins([1], [2], [[]])
    assert f.dtype == bool
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


@pytest.mark.parametrize("lo,hi,res",
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
                          (3.1, 3.1, [False, False, False, False, False])])
def test_filter_bins_one(lo, hi, res):
    """Can we filter the array between [lo,hi] [unintegrated]?

    This test is a regression test rather than a from-first-principles
    test: this is mainly relevant for the edge cases like: lo=hi and
    lo>hi.
    """

    expected = numpy.asarray(res)

    dvals = numpy.asarray([1, 2, 3, 4, 5])
    flags = utils.filter_bins([lo], [hi], [dvals])
    assert flags == pytest.approx(expected)

    # We can also check an identity: that
    #    a <= x <= b
    # is the same as
    #    a <= x
    #    x <= b
    #
    flags = utils.filter_bins([lo, None], [None, hi], [dvals, dvals])
    assert flags == pytest.approx(expected)



@pytest.mark.parametrize("lo,hi,res",
                         [(0, 10, [False] * 0 + [True] * 5 + [False] * 0),
                          (1, 5, [False] * 0 + [True] * 4 + [False] * 1),
                          (2, 5, [False] * 1 + [True] * 3 + [False] * 1),
                          (2, None, [False] * 1 + [True] * 4 + [False] * 0),
                          (1, 4, [False] * 0 + [True] * 3 + [False] * 2),
                          (None, 4, [False] * 0 + [True] * 3 + [False] * 2),
                          (1.1, 4.9, [False] * 0 + [True] * 4 + [False] * 1),
                          (2, 4, [False] * 1 + [True] * 2 + [False] * 2),
                          # Have minimum > maximum, which is technically invalid
                          (4, 3, [False] * 5 + [True] * 0 + [False] * 0),
                          # Have minimum = maximum = bin value
                          (4, 4, [False] * 5 + [True] * 0 + [False] * 0),
                          # Have minimum = maximum, not equal to a bin value
                          (3.1, 3.1, [False] * 2 + [True] * 1 + [False] * 2)])
def test_filter_bins_one_int(lo, hi, res):
    """Can we filter the array between [lo,hi) [integrated]?

    This test is a regression test rather than a from-first-principles
    test: this is mainly relevant for the edge cases like: lo=hi and
    lo>hi.

    The test replicates the logic of Data1DInt.notice
    """

    lovals = numpy.asarray([1, 2, 3, 4, 5])
    hivals = lovals + 1
    flags = utils.filter_bins([None, lo], [hi, None], [lovals, hivals],
                              integrated=True)
    assert flags == pytest.approx(numpy.asarray(res))



def test_filter_bins_two_none():
    """Use two different arrays for filtering with no filter.
    """

    y1 = [1, 2, 3, 4, 5]
    y2 = [10, 20, 30, 40, 50]

    flags = utils.filter_bins((None, None), (None, None), (y1, y2))
    assert flags is None


@pytest.mark.parametrize("lo1,lo2,hi1,hi2,expected",
                         [(1.5, 21, 3.6, 44, [False, False, True, False, False]),
                          (1.5, None, None, 44, [False, True, True, True, False]),
                          (None, None, None, 44, [True, True, True, True, False]),
                          (1.5, None, None, None, [False, True, True, True, True])])
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

    flags = utils.filter_bins((3, ), (8, ),
                              [[1, 4, 3, 7, 8, 10, 5]])

    expected = [False, True, True, True, True, False, True]

    assert len(flags) == len(expected)
    for got, exp in zip(flags, expected):
        assert got == exp


@pytest.mark.parametrize("arg", [None, '', '  \t  '])
def test_parse_expr_empty(arg):
    assert utils.parse_expr(arg) == [(None, None)]


@pytest.mark.parametrize("arg,expected",
                         [("2:", [(2.0, None)]),
                          ("1:4, 5:", [(1.0, 4.0), (5.0, None)]),
                          ("1:4, 5:6, 9:", [(1.0, 4.0), (5.0, 6.0), (9.0, None)]),
                          (":2", [(None, 2.0)]),
                          (":2, 3:4", [(None, 2.0), (3.0, 4.0)]),
                          (":2, 3:4, 5:6", [(None, 2.0), (3.0, 4.0), (5.0, 6.0)]),
                          (" 1:2,4:5  ", [(1.0, 2.0), (4.0, 5.0)]),
                          (":2, 3, 5:6", [(None, 2.0), (3.0, 3.0), (5.0, 6.0)]),
                          (" :2 ,4:  ", [(None, 2.0), (4.0, None)]),
                          ("1:2 , 4:5, 6:8", [(1.0, 2.0), (4.0, 5.0), (6.0, 8.0)]),
                          (" :2 , 4:5, 6: ", [(None, 2.0), (4.0, 5.0), (6.0, None)]),
                          (":", [(None, None)]),
                          ("2:3,:,4", [(2.0, 3.0), (None, None), (4.0, 4.0)]),
                          (",2", [(None, None), (2.0, 2.0)]),
                          ("2,", [(2.0, 2.0), (None, None)]),
                          ("2,,3:4", [(2.0, 2.0), (None, None), (3.0, 4.0)]),
                          (" , :", [(None, None), (None, None)])])
def test_parse_expr(arg, expected):
    """Check parse_expr with various conditions

    Would be nice to use something like hypothesis to use
    property-based testing.
    """
    assert utils.parse_expr(arg) == expected


@pytest.mark.parametrize("instr,bound",
                         [("None:1", "lower"),
                          ("1:None", "upper"),
                          ("None:None", "lower"),
                          (":2,None:4", "lower"),
                          (":2,3:None", "upper"),
                          (":2,None:", "lower"),
                          (":2,None:4,5:6", "lower"),
                          (":2,3:None,8:", "upper"),
                          ("1:2,3:5,None:None", "lower")])
def test_parse_expr_not_num(instr, bound):

    with pytest.raises(TypeError) as exc:
        utils.parse_expr(instr)

    emsg = "Invalid {} bound 'None'".format(bound)
    assert str(exc.value) == emsg


@pytest.mark.parametrize("instr",
                         ["1:2:4", "1:3 , 4 : 5:0.2"])
def test_parse_expr_unexpected_parses(instr):
    """You used to be able to say a:b:c:d:e and still have it parsed"""

    with pytest.raises(TypeError) as exc:
        utils.parse_expr(instr)

    assert str(exc.value) == "interval syntax requires a tuple, 'lo:hi'"


@pytest.mark.parametrize("instr,expected",
                         [("2:4,1:3", [(2.0, 4.0), (1.0, 3.0)]),
                          (":4, 2:3, 7:", [(None, 4.0), (2.0, 3.0), (7.0, None)])])
def test_parse_expr_no_range_checking(instr, expected):
    """Each range is separate, and there is no constraint."""

    assert utils.parse_expr(instr) == expected


@pytest.mark.parametrize("instr,expected",
                         [("4:2", [(4.0, 2.0)]),
                          ("2:1, :3", [(2.0, 1.0), (None, 3.0)]),
                          (":10, :12", [(None, 10.0), (None, 12.0)])])
def test_parse_expr_no_ordering(instr, expected):
    """There's no requirement that lo <= hi"""

    assert utils.parse_expr(instr) == expected


def test_create_expr_empty():
    """Simple test of create_expr with no mask."""

    out = utils.create_expr([])
    assert out == ""


def test_create_expr_empty_mask():
    """What happens if the mask is all masked out"""
    mask = numpy.zeros(5, dtype=bool)

    out = utils.create_expr([], mask)
    assert out == ""


@pytest.mark.parametrize("mask",
                         [numpy.zeros(5, dtype=bool),
                          numpy.ones(2, dtype=bool)])
def test_create_expr_mask_size_error(mask):
    """What happens when the mask][True] and vals arrays have different sizes?

    This should be an error but currently isn't.
    """

    with pytest.raises(ValueError) as exc:
        utils.create_expr([1, 2, 3, 4], mask)

    assert str(exc.value) == 'mask array mismatch with vals'


@pytest.mark.parametrize("val,expected",
                         [(numpy.int16(3), "3"),
                          (numpy.float32(3), "3.0")
                          ])
def test_create_expr_singleton(val, expected):
    """Simple test of create_expr with no mask."""

    out = utils.create_expr(numpy.asarray([val]))
    assert out == expected


def test_create_expr_nomask():
    """Simple test of create_expr with no mask."""

    chans = numpy.arange(1, 10, dtype=numpy.int16)
    out = utils.create_expr(chans)
    assert out == "1-9"


def test_create_expr_mask():
    """Simple test of create_expr with an identity mask."""

    chans = numpy.arange(1, 10, dtype=numpy.int16)
    mask = numpy.ones(9, dtype=bool)
    out = utils.create_expr(chans, mask)
    assert out == "1-9"


def test_create_expr_nomask_delim():
    """Simple test of create_expr with no mask."""

    chans = numpy.arange(1, 10, dtype=numpy.int16)
    out = utils.create_expr(chans, delim='::')
    assert out == "1::9"


def test_create_expr_mask_delim():
    """Simple test of create_expr with an identity mask."""

    chans = numpy.arange(1, 10, dtype=numpy.int16)
    mask = numpy.ones(9, dtype=bool)
    out = utils.create_expr(chans, mask, delim='::')
    assert out == "1::9"


def test_create_expr_missing_nomask():
    """Simple test of create_expr with no mask with missing data"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    chans = chans[(chans < 12) | (chans > 13)]
    out = utils.create_expr(chans)
    assert out == "1-11,14-19"


def test_create_expr_missing_nomask_nonchannels():
    """What happens if the input is not in channel units?"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = (chans < 12) | (chans > 13)

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], format='%4.2f')

    # I am not sure this is the correct output, but it's what
    # we create.
    #
    vs = ['{:4.2f}'.format(v) for v in vals[filt]]
    expected = ','.join(vs)
    assert out == expected


def test_create_expr_missing_mask():
    """Simple test of create_expr with an identity with missing data"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = (chans < 12) | (chans > 13)
    out = utils.create_expr(chans[filt], filt)
    assert out == "1-11,14-19"


def test_create_expr_missing_mask_nonchannels():
    """What happens if the input is not in channel units?"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = (chans < 12) | (chans > 13)

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')
    assert out == "0.40-0.50,0.53-0.58"


def test_create_expr_missing_start_nomask():
    """Simple test of create_expr with no mask with missing data

    Have a single bin, then missing data, then the rest
    """

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    chans = chans[chans != 2]
    out = utils.create_expr(chans)
    assert out == "1,3-19"


def test_create_expr_missing_start2_nomask():
    """Simple test of create_expr with no mask with missing data

    Have a single bin, then missing data, then the rest
    """

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    chans = chans[chans != 3]
    out = utils.create_expr(chans)

    assert out == "1-2,4-19"


def test_create_expr_missing_start2_mask():
    """How does this compare to nomask version?"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = chans != 3

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')

    assert out == "0.40-0.41,0.43-0.58"


def test_create_expr_missing_end_nomask():
    """Simple test of create_expr with no mask with missing data

    Ends with a single bin.
    """

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    chans = chans[chans != 18]
    out = utils.create_expr(chans)
    assert out == "1-17,19"


def test_create_expr_missing_multiple1_nomask():
    """Simple test of create_expr with no mask with missing data"""

    chans = numpy.asarray([1, 2, 4, 6, 7, 10, 18], dtype=numpy.int16)
    out = utils.create_expr(chans)

    assert out == "1-2,4,6-7,10,18"


def test_create_expr_missing_multiple2_nomask():
    """Simple test of create_expr with no mask with missing data"""

    chans = numpy.asarray([1, 2, 3, 6, 7, 10, 12, 13, 14, 17, 18], dtype=numpy.int16)
    out = utils.create_expr(chans)
    assert out == "1-3,6-7,10,12-14,17-18"


def test_create_expr_mask_dropstart():
    """Start is masked out"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = chans > 3

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')
    assert out == "0.43-0.58"


def test_create_expr_mask_dropstart2():
    """Start is masked out"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = (chans > 3) & ((chans < 8) | (chans > 10))

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')
    assert out == "0.43-0.46,0.50-0.58"


def test_create_expr_mask_dropend():
    """End is masked out"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = chans < 5

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')
    assert out == "0.40-0.43"


def test_create_expr_mask_dropend2():
    """End is masked out"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = (chans < 15) & ((chans < 8) | (chans > 11))

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')
    assert out == "0.40-0.46,0.51-0.53"


def test_create_expr_mask_drop():
    """mask all the things (but leave some data!)"""

    chans = numpy.arange(1, 20, dtype=numpy.int16)
    filt = (chans > 5) & (chans < 15) & ((chans < 10) | (chans > 12))

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')
    assert out == "0.45-0.48,0.52-0.53"


@pytest.mark.parametrize("idx", [0, 1, 5, 17, 18])
def test_create_expr_mask_singlebin(idx):
    """mask all the things (but leave one bin)"""

    filt = numpy.zeros(19, dtype=bool)
    filt[idx] = True

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')
    expected = "{:4.2f}".format(vals[filt][0])
    assert out == expected


def test_create_expr_mask_singlebins():
    """several single bin"""

    filt = numpy.zeros(19, dtype=bool)
    filt[0] = True
    filt[2] = True
    filt[16] = True
    filt[18] = True

    vals = numpy.arange(19) * 0.01 + 0.4
    out = utils.create_expr(vals[filt], filt, format='%4.2f')
    assert out == "0.40,0.42,0.56,0.58"


def test_create_expr_integrated_empty():
    """"No data"""

    out = utils.create_expr_integrated([], [])
    assert out == ''


def test_create_expr_integrated_size_error():
    """What happens when lovals and hivals have different sizes?"""

    with pytest.raises(ValueError) as exc:
        utils.create_expr_integrated([1, 2, 3, 4], [2, 3])

    assert str(exc.value) == 'hivals array mismatch with lovals'


def test_create_expr_integrated_mask_size_error():
    """What happens when lovals and mask have different sizes?"""

    with pytest.raises(ValueError) as exc:
        utils.create_expr_integrated([1, 2, 3, 4], [4, 5, 6, 7], [True, True])

    assert str(exc.value) == 'mask array mismatch with lovals'


@pytest.mark.parametrize("expected,lovals,hivals",
                         [('1-5', [1, 2, 3, 4], [2, 3, 4, 5]),
                          ('1-8', [1, 2, 4, 5, 7], [2, 3, 5, 6, 8]),
                          ('0.1-1.0', [0.1, 0.2, 0.4, 0.8], [0.2, 0.4, 0.8, 1.0]),
                          ('0.1-1.0', [0.1, 0.2, 0.4, 0.8], [0.2, 0.4, 0.6, 1.0])
                         ])
def test_create_expr_integrated_docstring_no_mask(expected, lovals, hivals):
    """"Check the docstring examples"""

    out = utils.create_expr_integrated(lovals, hivals)
    assert out == expected


@pytest.mark.parametrize("expected,lovals,hivals,mask",
                         [('0.1-1.0', [0.1, 0.2, 0.4, 0.8], [0.2, 0.4, 0.6, 1.0], [True, True, True, True]),
                          ('1-3,4-5', [1, 2, 4], [2, 3, 5], [True, True, False, True]),
                          ('1-4,5-6', [1, 3, 5], [2, 4, 6], [True, True, False, True]),
                          ('0.1-0.4,0.6-1.0', [0.1, 0.2, 0.6, 0.8], [0.2, 0.4, 0.8, 1.0],[True, True, False, True, True]),
                          ('0.1-0.3,0.4-0.5,0.8-1.0', [0.1, 0.2, 0.4, 0.8], [0.2, 0.3, 0.5, 1.0], [True, True, False, True, False, True]),
                          ('0.1-0.3,0.4-0.5,0.8-1.0', [0.1, 0.2, 0.4, 0.8], [0.2, 0.3, 0.5, 1.0], [False, True, True, False, True, False, True, False]),
                          ('1-2,2-5', [1, 2, 3, 4], [2, 3, 4, 5], [True, False, True, True, True]),
                         ])
def test_create_expr_integrated_docstring_mask(expected, lovals, hivals ,mask):
    """"Check the docstring examples"""

    out = utils.create_expr_integrated(lovals, hivals, mask)
    assert out == expected


def test_send_to_pager_file(tmpdir):
    """Check we can create a file"""

    # would like to use tmp_path but this needs a recent pytest

    fileexists = tmpdir.join('tmp.txt')
    fileexists.write_text('x', 'ascii')
    assert fileexists.read() == 'x'

    testdata = "test\ntest"
    utils.send_to_pager(testdata, str(fileexists), clobber=True)
    assert fileexists.read_text('ascii') == testdata + "\n"


@pytest.mark.parametrize('clobber', [False, True])
def test_send_to_pager_stringio(clobber):
    """Check we can create a StringIO object

    We ignore the clobber parameter in this approach.
    """

    out = StringIO()
    testdata = "test\ntest"
    utils.send_to_pager(testdata, out, clobber=clobber)
    assert out.getvalue() == testdata + "\n"


def test_send_to_pager_file_clobber(tmpdir):
    """Ensure clobber=False doesn't overwrite things"""

    # would like to use tmp_path but this needs a recent pytest

    fileexists = tmpdir.join('tmp.txt')
    fileexists.write_text('x', 'ascii')

    with pytest.raises(IOErr) as exc:
        utils.send_to_pager("test", str(fileexists))

    assert str(exc.value).startswith("file '")
    assert str(exc.value).endswith("/tmp.txt' exists and clobber is not set")
