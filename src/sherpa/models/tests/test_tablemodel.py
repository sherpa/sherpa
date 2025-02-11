#
#  Copyright (C) 2022
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

"""Test the tablemodel class.

The behavior of the tablemodel is rather specialized, which means that
some contorted tests need to be written.
"""

import pickle

import numpy as np

import pytest

from sherpa.data import Data1D
from sherpa.models.basic import TableModel
from sherpa.models.model import ArithmeticModel
from sherpa.utils import linear_interp, nearest_interp
from sherpa.utils.err import ModelErr


def test_create():
    """Basic creation test"""

    mdl = TableModel()
    assert isinstance(mdl, ArithmeticModel)


def test_name():
    """Basic creation test: name"""

    mdl = TableModel()
    assert mdl.name == "tablemodel"


def test_named():
    """Basic creation test: name"""

    mdl = TableModel("bobby")
    assert mdl.name == "bobby"


def test_filename():
    """Basic creation test: filename

    How is the filename attribute meant to be used?
    """

    mdl = TableModel()
    assert mdl.filename is None


def test_ndim():
    """What is the ndim?

    At the moment we can use this with nD data but only really if we do
    not set the "x" argument.
    """

    mdl = TableModel()
    assert mdl.ndim is None


def test_pars():
    """We have a single parameter."""

    mdl = TableModel()
    assert len(mdl.pars) == 1
    assert mdl.pars[0].name == "ampl"
    assert mdl.pars[0].val == pytest.approx(1)
    assert not mdl.pars[0].frozen


def test_show():
    """Basic show test"""

    mdl = TableModel("test")
    out = str(mdl).split("\n")
    assert out[0] == "test"
    assert out[1].startswith("   Param        Type ")
    assert out[2].startswith("   -----        ---- ")
    assert out[3].startswith("   test.ampl    thawed ")
    assert len(out) == 4


def test_empty():
    """Basic creation test: no data"""

    mdl = TableModel()
    assert mdl.get_x() is None
    assert mdl.get_y() is None


def test_method_change():
    """We can change the method"""

    mdl = TableModel()
    assert mdl.method == linear_interp

    mdl.method = nearest_interp
    assert mdl.method == nearest_interp


def test_method_must_be_callable():
    """method must be callable"""

    mdl = TableModel()

    # This comes from NoNewAttributesAfterInit via AttrithmeticModel
    #
    with pytest.raises(AttributeError,
                       match="^'TableModel' object attribute 'method' cannot be replaced with a non-callable attribute$"):
        mdl.method = 23


def test_load_no_x():
    """load: y only"""

    mdl = TableModel()
    mdl.load(None, [3, 2, 7])
    assert mdl.get_x() is None
    assert mdl.get_y() == pytest.approx([3, 2, 7])


def test_load_x_y_sorted():
    """load: x and y (sorted)"""

    mdl = TableModel()
    mdl.load([12, 14, 17], [3, 2, 7])
    assert mdl.get_x() == pytest.approx([12, 14, 17])
    assert mdl.get_y() == pytest.approx([3, 2, 7])


def test_load_x_y_unsorted():
    """load: x and y (unsorted)"""

    mdl = TableModel()
    mdl.load([17, 12, 14], [7, 3, 2])
    assert mdl.get_x() == pytest.approx([12, 14, 17])
    assert mdl.get_y() == pytest.approx([3, 2, 7])


def test_load_no_x_check_storage():
    """How are the x/y values stored?

    This is a regression test (that is, the behavior could be
    changed).

    """

    mdl = TableModel()
    mdl.load(None, [1, 3, 7])
    y = mdl.get_y()
    assert isinstance(y, np.ndarray)


def test_load_x_y_sorted_check_storage():
    """How are the x/y values stored?

    This is a regression test.
    """

    mdl = TableModel()
    mdl.load([1, 2, 3], [1, 3, 7])
    x = mdl.get_x()
    y = mdl.get_y()
    assert isinstance(x, np.ndarray)
    assert isinstance(y, np.ndarray)


def test_load_x_y_unsorted_check_storage():
    """How are the x/y values stored?

    This is a regression test.
    """

    mdl = TableModel()
    mdl.load([3, 2, 1], [1, 3, 7])
    x = mdl.get_x()
    y = mdl.get_y()
    assert isinstance(x, np.ndarray)
    assert isinstance(y, np.ndarray)


def test_load_no_y():
    """do we error out or not?

    This is a regression test.
    """

    mdl = TableModel()
    with pytest.raises(ModelErr,
                       match="^y must be set if x is set$"):
        mdl.load([12, 14, 17], None)


def test_load_x_longer_than_y():
    """do we error out or not?

    This is a regression test.
    """

    mdl = TableModel()
    with pytest.raises(ModelErr,
                       match="^size mismatch between x and y: 3 vs 2$"):
        mdl.load([12, 14, 17], [3, 2])


def test_load_x_shorter_than_y():
    """do we error out or not?

    This is a regression test.
    """

    mdl = TableModel()
    with pytest.raises(ModelErr,
                       match="^size mismatch between x and y: 3 vs 4$"):
        mdl.load([12, 17, 14], [3, 2, 4, 5])


def test_load_y_not_1d():
    """do we error out or not?

    This is a regression test.
    """

    x = np.arange(6) + 1
    y = np.arange(6).reshape(2, 3)
    mdl = TableModel()
    with pytest.raises(ModelErr,
                       match="Array must be 1D or None"):
        mdl.load(x, y)


# sherpa.utils.interpolate assumes inputs are ndarray, so work around
# this by forcing the input to be one.
#
@pytest.mark.parametrize("x2,iflag",
                         [(None, False),
                          ([14, 16, 20, 25], False),
                          ([14, 16, 20, 25], True)])
def test_eval_x_y(x2, iflag):
    """Check we can interpolate onto a grid after load(x,y)

    The application ignores the integrate setting/upper-edge of the bin.
    """

    tm = TableModel()
    tm.load([10, 15, 20], [4, 6, 2])
    tm.ampl = 10
    tm.integrate = iflag

    yout = tm(np.asarray([12, 14, 16, 20]), x2)
    exp = 10 * np.asarray([4.8, 5.6, 5.2, 2])
    assert yout == pytest.approx(exp)


@pytest.mark.parametrize("y", [[4, 6, 2], np.asarray([4, 6, 2])])
@pytest.mark.parametrize("x2,iflag",
                         [(None, False),
                          ([14, 16, 20], False),
                          ([14, 16, 20], True)])
def test_eval_no_x(y, x2, iflag):
    """Check we can return something after load(None,y)

    The application ignores the integrate setting/upper-edge of the
    bin.  We also ignore the x array values, presumably to allow users
    to load data with no dependent axis (or that we know we don't care
    about the x axis).

    """

    tm = TableModel()
    tm.load(None, y)
    tm.ampl = 10
    tm.integrate = iflag

    # fails when y is not a ndarray
    yout = tm([12, 14, 16], x2)
    assert yout == pytest.approx([40, 60, 20])


@pytest.mark.parametrize("y", [[4, 6, 2, 8, 7], np.asarray([4, 6, 2, 8, 7])])
@pytest.mark.parametrize("iflag", [False, True])
def test_eval_filtered(y, iflag):
    """Check we can return something after load(None, y) and filtering the dataset.
    """

    tm = TableModel()
    tm.load(None, y)
    tm.ampl = 10
    tm.integrate = iflag

    d = Data1D("ex", [100, 200, 300, 400, 500], [4, 5, 6, 7, 8])
    d.ignore(xhi=100)
    d.ignore(xlo=350, xhi=450)

    tm.fold(d)

    yout = d.eval_model_to_fit(tm)
    assert yout == pytest.approx([60, 20, 70])


def test_eval_pick_up_changed_load():
    """Corner case: does __filtered_y get cleaned up?

    There's no "nice" way to check if __filtered_y gets cleared after
    a second load call, so we test it by evaluating the model.

    """

    tm = TableModel()
    tm.load(None, [5, 7, 2, 8, 4])
    tm.ampl = 10

    d = Data1D('simple', [20, 30, 40, 50, 70], [1, 2, 3, 4, 5])
    d.notice(25)

    # Sets up tm.__filtered_y to be [7, 2, 8, 4]
    tm.fold(d)

    # check
    assert tm([1, 2, 3, 4]) == pytest.approx([70, 20, 80, 40])

    # Call load with a different y array
    #
    tm.load(None, [2, 3, 1, 4, 2])

    # The filter has been removed, so we can not call it with a
    # filtered array.
    #
    assert tm([1, 2, 3, 4, 5]) == pytest.approx([20, 30, 10, 40, 20])

    with pytest.raises(ModelErr,
                       match=r"^Mismatch between table model and data, \(5 vs 4\)$"):
        tm([1, 2, 3, 4])


@pytest.mark.parametrize("flag", [True, False])
def test_eval_filter_all_or_none(flag):
    """Check that .fold does nothing when the mask is True or False

    This is a regression test.
    """

    mdl = TableModel()
    mdl.load(None, np.asarray([5, 3, 4]))

    d = Data1D("ex", [1, 2, 3], [1, 1, 2])

    # Either ensure all the data is noticed or ignored.
    #
    d.mask = flag

    # If all the data is ignored, should we error out?
    #
    mdl.fold(d)

    y = mdl([1, 2, 3])
    assert y == pytest.approx([5, 3, 4])


def test_eval_pass_kwargs():
    """Check we can pass kwargs, even though we ignore them.

    This is a regression test.
    """

    mdl = TableModel()
    mdl.load(None, [10, 12 , 2])
    mdl.ampl = 5

    y = mdl([1, 2, 3], arg1=True, not_a_kwarg={"answer": 23})
    assert y == pytest.approx([50, 60, 10])


def test_eval_checks_length():
    """Check we error out when the lengths are wrong: calc via model evaluation

    This requires
      - the table data not being loaded

    """

    tm = TableModel('tbl')

    with pytest.raises(ModelErr,
                       match="^The tablemodel's load method must be called first$"):
        tm([1, 2, 3])


def test_fold_no_load_checks_length():
    """Check we error out when the lengths are wrong: fold.

    This requires
      - the table data having no X array
      - the data being partially filtered

    """

    ytbl = [5, 7, 3, 2]
    tm = TableModel('tbl')

    xdata = np.asarray([0, 10, 20, 25, 30, 40])
    ydata = np.ones_like(xdata)
    d = Data1D('ex', xdata, ydata)
    d.ignore(xhi=7)

    with pytest.raises(ModelErr,
                       match="^The tablemodel's load method must be called first$"):
        tm.fold(d)


def test_fold_load_checks_length():
    """Check we error out when the lengths are wrong: fold.

    This requires
      - the table data having no X array
      - the data being partially filtered

    """

    ytbl = [5, 7, 3, 2]
    tm = TableModel('tbl')
    tm.load(None, ytbl)

    xdata = np.asarray([0, 10, 20, 25, 30, 40])
    ydata = np.ones_like(xdata)
    d = Data1D('ex', xdata, ydata)
    d.ignore(xhi=7)

    with pytest.raises(ModelErr,
                       match=r"^Mismatch between table model and data, \(4 vs 6\)$"):
        tm.fold(d)


def test_cached_basic():
    """Basic check the cache works.

    This is mainly here so we can have confidence in
    test_cache_is_revoked.
    """

    mdl = TableModel()
    mdl.load([1, 2, 3], [5, 2, 12])

    # We could check the _cache_ctr all in one go but we do want false
    # positives if we ever decide to change how the cache works, so
    # test individual fields.
    #
    assert mdl._cache_ctr["check"] == 0

    y = mdl([1, 2, 3])
    assert y == pytest.approx([5, 2, 12])
    assert mdl._cache_ctr["check"] == 1
    assert mdl._cache_ctr["hits"] == 0
    assert mdl._cache_ctr["misses"] == 1

    # Is cached
    y = mdl([1, 2, 3])
    assert y == pytest.approx([5, 2, 12])
    assert mdl._cache_ctr["check"] == 2
    assert mdl._cache_ctr["hits"] == 1
    assert mdl._cache_ctr["misses"] == 1

    # Not cached
    y = mdl([1, 3])
    assert y == pytest.approx([5, 12])
    assert mdl._cache_ctr["check"] == 3
    assert mdl._cache_ctr["hits"] == 1
    assert mdl._cache_ctr["misses"] == 2


def test_cache_is_revoked():
    """Loading new data should clear the cache. Does it?
    """

    mdl = TableModel()
    mdl.load([1, 2, 3], [5, 2, 12])

    y = mdl([1, 2, 3])
    assert y == pytest.approx([5, 2, 12])

    mdl.load([1, 2, 3], [4, 9, -3])

    # We could check the "check" field here but I want an explicit
    # check that the model cache works (that is, the model evaluation
    # works), rather than an implicit check which can come later.
    #
    cache = mdl._cache_ctr.copy()

    y = mdl([1, 2, 3])
    assert y == pytest.approx([4, 9, -3])

    assert cache["check"] == 0

    assert mdl._cache_ctr["check"] == 1
    assert mdl._cache_ctr["hits"] == 0
    assert mdl._cache_ctr["misses"] == 1


def test_cache_method():
    """Does changing the method cause the cache to change?

    In this case we just test the model evaluation,
    and not the actual cache values.
    """

    mdl = TableModel()
    mdl.load([1, 2, 3], [5, 2, 12])
    assert mdl.method == linear_interp

    y = mdl([1.4, 2.6])
    assert y == pytest.approx([3.8, 8])

    mdl.method = nearest_interp
    y = mdl([1.4, 2.6])
    assert y == pytest.approx([5, 12])


def test_pickle_none(tmp_path):
    """Can we save/restore a TableModel with no data?"""

    x = [1, 2, 3]
    y = [4, -2, 13]
    xtest = [1.5, 2, 2.5, 3]
    ytest = [1, -2, 5.5, 13]

    mdl = TableModel("fred")

    out = tmp_path / "model.state"
    out = str(out)
    with open(out, "wb") as ofh:
        pickle.dump(mdl, ofh)

    with open(out, "rb") as ifh:
        nmdl = pickle.load(ifh)

    assert isinstance(nmdl, TableModel)
    assert nmdl.name == "fred"
    assert nmdl.filename is None
    assert nmdl.get_x() is None
    assert nmdl.get_y() is None
    assert nmdl.ampl.val == pytest.approx(1)
    assert not nmdl.ampl.frozen

    assert nmdl.method is linear_interp

    # Check we can use the restored object
    #
    nmdl.load(x, y)
    assert nmdl(xtest) == pytest.approx(ytest)


def test_pickle_x_y(tmp_path):
    """Can we save/restore a TableModel with x/y?"""

    x = [1, 2, 3]
    y = [4, -2, 13]
    xtest = [1.5, 2, 2.5, 3]
    ytest = [-5 * 1, -5 * -2, -5 * 5.5, -5 * 13]

    mdl = TableModel()
    mdl.load(x, y)
    mdl.ampl = -5
    mdl.ampl.freeze()
    assert mdl(xtest) == pytest.approx(ytest)

    out = tmp_path / "model.state"
    out = str(out)
    with open(out, "wb") as ofh:
        pickle.dump(mdl, ofh)

    with open(out, "rb") as ifh:
        nmdl = pickle.load(ifh)

    assert isinstance(nmdl, TableModel)
    assert nmdl.name == "tablemodel"
    assert nmdl.filename is None
    assert nmdl.get_x() == pytest.approx(x)
    assert nmdl.get_y() == pytest.approx(y)
    assert nmdl.ampl.val == pytest.approx(-5)
    assert nmdl.ampl.frozen

    assert nmdl(xtest) == pytest.approx(ytest)


def test_pickle_y(tmp_path):
    """Can we save/restore a TableModel with just y and filtered?"""

    x = [0, 1, 2, 3, 5]
    y = [0, 4, -2, 13, 2]
    xtest = [1, 2, 3]
    ytest = [-5 * 4, -5 * -2, -5 * 13]

    mdl = TableModel("test")
    mdl.load(None, y)
    mdl.ampl = -5
    mdl.ampl.freeze()

    # let's change the filename to test
    filename = "dummy.dat"
    mdl.filename = filename

    data = Data1D("dummy", x, [1] * 5)
    data.notice(1, 3)

    mdl.fold(data)

    assert mdl(xtest) == pytest.approx(ytest)

    out = tmp_path / "model.state"
    out = str(out)
    with open(out, "wb") as ofh:
        pickle.dump(mdl, ofh)

    with open(out, "rb") as ifh:
        nmdl = pickle.load(ifh)

    assert isinstance(nmdl, TableModel)
    assert nmdl.name == "test"
    assert nmdl.filename == filename
    assert nmdl.get_x() is None
    assert nmdl.get_y() == pytest.approx(y)
    assert nmdl.ampl.val == pytest.approx(-5)
    assert nmdl.ampl.frozen

    assert nmdl(xtest) == pytest.approx(ytest)


def test_pickle_interpolation(tmp_path):
    """Can we save/restore a TableModel with different interpolation

    This is only relevant for load(x, y) cases.
    """

    x = [0, 1, 2, 3, 5]
    y = [0, 4, -2, 13, 2]
    xtest = [0.6, 3.9]
    ytest = [-5 * 4, -5 * 13]

    mdl = TableModel("test")
    mdl.load(x, y)
    mdl.ampl = -5
    mdl.ampl.freeze()

    mdl.method = nearest_interp

    assert mdl(xtest) == pytest.approx(ytest)

    out = tmp_path / "model.state"
    out = str(out)
    with open(out, "wb") as ofh:
        pickle.dump(mdl, ofh)

    with open(out, "rb") as ifh:
        nmdl = pickle.load(ifh)

    assert isinstance(nmdl, TableModel)
    assert nmdl.name == "test"
    assert nmdl.filename is None
    assert nmdl.get_x() == pytest.approx(x)
    assert nmdl.get_y() == pytest.approx(y)
    assert nmdl.ampl.val == pytest.approx(-5)
    assert nmdl.ampl.frozen

    assert nmdl.method is nearest_interp

    assert nmdl(xtest) == pytest.approx(ytest)
