#
#  Copyright (C) 2019, 2020, 2021, 2022
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

import numpy

import pytest

from sherpa.data import Data, Data1D, DataSimulFit, Data1DInt, Data2D, \
    Data2DInt, BaseData
from sherpa.models import Polynom1D, Polynom2D
from sherpa.utils.err import NotImplementedErr, DataErr
from sherpa.ui.utils import Session
from sherpa.astro.ui.utils import Session as AstroSession


NAME = "data_test"
X_ARRAY = numpy.arange(0, 10, 1)
Y_ARRAY = numpy.arange(100, 110, 1)
X0_2D_RAW, X1_2D_RAW = numpy.meshgrid(X_ARRAY, X_ARRAY)
Y_2D_RAW = X0_2D_RAW + X1_2D_RAW
Y_2D = Y_2D_RAW.ravel()
X0_2D, X1_2D = X0_2D_RAW.ravel(), X1_2D_RAW.ravel()
SHAPE_2D = X_ARRAY.size, X_ARRAY.size
SYSTEMATIC_ERROR_ARRAY = numpy.arange(0, 0.10, 0.01)
STATISTICAL_ERROR_ARRAY = numpy.arange(0, 1, 0.1)
SYS_ERROR_2D = Y_2D / 10
STAT_ERROR_2D = Y_2D / 5
X_THRESHOLD = 3
MULTIPLIER = 2

DATA_1D_CLASSES = (Data1D, Data1DInt)
DATA_2D_CLASSES = (Data2D, Data2DInt)
ALL_DATA_CLASSES = DATA_1D_CLASSES + DATA_2D_CLASSES
REALLY_ALL_DATA_CLASSES = (Data, ) + ALL_DATA_CLASSES

DATA_ARGS = NAME, (X_ARRAY,), Y_ARRAY, STATISTICAL_ERROR_ARRAY, SYSTEMATIC_ERROR_ARRAY
DATA1D_ARGS = NAME, X_ARRAY, Y_ARRAY, STATISTICAL_ERROR_ARRAY, SYSTEMATIC_ERROR_ARRAY
DATA1DINT_ARGS = NAME, X_ARRAY - 0.5, X_ARRAY + 0.5, Y_ARRAY, STATISTICAL_ERROR_ARRAY, SYSTEMATIC_ERROR_ARRAY
DATA2D_ARGS = NAME, X0_2D, X1_2D, Y_2D, SHAPE_2D, STAT_ERROR_2D, SYS_ERROR_2D
DATA2DINT_ARGS = NAME, X0_2D - 0.5, X1_2D - 0.5, X0_2D + 0.5, X1_2D + 0.5, Y_2D, SHAPE_2D, STAT_ERROR_2D, SYS_ERROR_2D


INSTANCE_ARGS = {
    Data1D: DATA1D_ARGS,
    Data: DATA_ARGS,
    Data1DInt: DATA1DINT_ARGS,
    Data2D: DATA2D_ARGS,
    Data2DInt: DATA2DINT_ARGS
}


POS_X_ARRAY = {
    Data1D: 1,
    Data: 1,
    Data1DInt: 1,
    Data2D: 1,
    Data2DInt: 1,
    }

POS_Y_ARRAY = {
    Data1D: 2,
    Data: 2,
    Data1DInt: 3,
    Data2D: 3,
    Data2DInt: 5,
    }

POS_STATERR_ARRAY = {
    Data1D: 3,
    Data: 3,
    Data1DInt: 4,
    Data2D: 5,
    Data2DInt: 7,
    }


POS_SYSERR_ARRAY = {
    Data1D: 4,
    Data: 4,
    Data1DInt: 5,
    Data2D: 6,
    Data2DInt: 8,
    }


@pytest.fixture
def data(request):
    data_class = request.param

    return data_class(*INSTANCE_ARGS[data_class])


@pytest.fixture
def data_no_errors(request):
    data_class = request.param

    # Use the normal arguments but remove the error values
    all_args = INSTANCE_ARGS[data_class]
    no_errors = all_args[:POS_STATERR_ARRAY[data_class]]
    out = data_class(*no_errors)
    assert out.staterror is None
    assert out.syserror is None
    return out


@pytest.fixture
def data_simul_fit():
    data_one = Data1D("data_one", X_ARRAY, Y_ARRAY, STATISTICAL_ERROR_ARRAY, SYSTEMATIC_ERROR_ARRAY)
    data_two = Data1D("data_two", MULTIPLIER * X_ARRAY, MULTIPLIER * Y_ARRAY,
                      MULTIPLIER * STATISTICAL_ERROR_ARRAY, MULTIPLIER * SYSTEMATIC_ERROR_ARRAY)
    return DataSimulFit(NAME, (data_one, data_two))


@pytest.fixture
def data_simul_fit_no_errors():
    data_one = Data1D("data_one", X_ARRAY, Y_ARRAY)
    data_two = Data1D("data_two", MULTIPLIER * X_ARRAY, MULTIPLIER * Y_ARRAY)
    return DataSimulFit(NAME, (data_one, data_two))


@pytest.fixture
def data_simul_fit_some_errors():
    data_one = Data1D("data_one", X_ARRAY, Y_ARRAY, STATISTICAL_ERROR_ARRAY, SYSTEMATIC_ERROR_ARRAY)
    data_two = Data1D("data_two", MULTIPLIER * X_ARRAY, MULTIPLIER * Y_ARRAY)
    return DataSimulFit(NAME, (data_one, data_two))


@pytest.mark.xfail
def test_base_data_instantiation():
    with pytest.raises(NotImplementedErr):
        BaseData()


@pytest.mark.parametrize("data", (Data, Data2D, Data2DInt), indirect=True)
def test_data_get_x(data):
    with pytest.raises(AttributeError):
        data.get_x()


@pytest.mark.xfail
@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_get_x_special(data):
    # XFAIL: These classes still provide get_x
    with pytest.raises(AttributeError):
        data.get_x()


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_x0(data):
    with pytest.raises(AttributeError):
        data.get_x0()


@pytest.fixture
def data_for_load_arrays(request):
    data_class = request.param
    from sherpa.astro.ui.utils import Session
    session = Session()
    data_args = INSTANCE_ARGS[data_class]
    args = data_args + (data_class,)
    data = data_class(*data_args)
    return session, args, data


@pytest.mark.parametrize("data_for_load_arrays", ALL_DATA_CLASSES, indirect=True)
def test_load_arrays(data_for_load_arrays):
    session, args, data = data_for_load_arrays
    session.load_arrays(*args)
    new_data = session.get_data(data.name)
    assert new_data is not data
    # DATA-NOTE: Do we need an equality operator for data classes? These tests are very partial
    numpy.testing.assert_array_equal(new_data.get_indep(), data.get_indep())
    numpy.testing.assert_array_equal(new_data.get_dep(), data.get_dep())


# DATA-NOTE: In the current Sherpa Data cannot be correctly loaded using load_arrays
@pytest.mark.xfail
@pytest.mark.parametrize("data_for_load_arrays", (Data, ), indirect=True)
def test_load_arrays_data(data_for_load_arrays):
    session, args, _ = data_for_load_arrays
    session.load_arrays(*args)


@pytest.mark.parametrize("data_no_errors", ALL_DATA_CLASSES, indirect=True)
def test_load_arrays_no_errors(data_no_errors):
    from sherpa.astro.ui.utils import Session
    session = Session()
    data = data_no_errors
    data_class = data.__class__
    data_args = INSTANCE_ARGS[data_class]
    args = data_args + (data_class,)
    session.load_arrays(*args)
    new_data = session.get_data(data.name)
    assert new_data is not data
    # DATA-NOTE: Do we need an equality operator for data classes? These tests are very partial
    # Note that when they are created with load_arrays they seem to lose the name, which becomes the ID
    numpy.testing.assert_array_equal(new_data.get_indep(), data.get_indep())
    numpy.testing.assert_array_equal(new_data.get_dep(), data.get_dep())


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_x1(data):
    with pytest.raises(AttributeError):
        data.get_x1()


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_xlabel(data):
    assert data.get_xlabel() == "x"


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data_get_x0label(data):
    assert data.get_x0label() == "x0"


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data_get_x1label(data):
    assert data.get_x1label() == "x1"


@pytest.mark.parametrize("data", ALL_DATA_CLASSES, indirect=True)
def test_data_get_ylabel(data):
    assert data.get_ylabel() == "y"


@pytest.mark.parametrize("data", (Data, ), indirect=True)
def test_data_get_dims(data):
    assert data.get_dims() == ((X_ARRAY.size, ), X_ARRAY.size)


@pytest.mark.parametrize("data", (Data, ), indirect=True)
def test_data_str_repr(data):
    assert repr(data) == "<Data data set instance 'data_test'>"
    assert str(data) == 'name      = data_test\nindep     = (array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),)\ndep       ' \
                        '= Int64[10]\nstaterror = Float64[10]\nsyserror  = Float64[10]'


@pytest.mark.parametrize("data", (Data1D, ), indirect=True)
def test_data1d_str_repr(data):
    assert repr(data) == "<Data1D data set instance 'data_test'>"
    assert str(data) == 'name      = data_test\nx         = Int64[10]\ny         = Int64[10]\nstaterror = ' \
                        'Float64[10]\nsyserror  = Float64[10]'


@pytest.mark.parametrize("data", (Data, Data1D), indirect=True)
def test_data_get_indep(data):
    numpy.testing.assert_array_equal(data.get_indep(), [X_ARRAY, ])


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_get_indep(data):
    numpy.testing.assert_array_equal(data.get_indep(), (X_ARRAY-0.5, X_ARRAY+0.5))


@pytest.mark.parametrize("data", (Data1D, Data), indirect=True)
def test_data_get_indep_filter(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    numpy.testing.assert_array_equal(data.get_indep(filter=True), [X_ARRAY[:X_THRESHOLD + 1], ])


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_get_indep_filter(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    expected = (X_ARRAY-0.5)[:X_THRESHOLD + 1], (X_ARRAY+0.5)[:X_THRESHOLD + 1]
    numpy.testing.assert_array_equal(data.get_indep(filter=True), expected)


@pytest.mark.parametrize("data", (Data1D, ), indirect=True)
def test_data_1d_get_indep_ignore(data):
    data.ignore(0, X_THRESHOLD)
    numpy.testing.assert_array_equal(data.get_indep(filter=True), [X_ARRAY[X_THRESHOLD + 1:], ])


@pytest.mark.parametrize("data", (Data, ), indirect=True)
def test_data_get_indep_ignore(data):
    data.ignore((0, ), (X_THRESHOLD, ))
    numpy.testing.assert_array_equal(data.get_indep(filter=True), [X_ARRAY[X_THRESHOLD + 1:], ])


@pytest.mark.parametrize("data", (Data1D, ), indirect=True)
def test_data_1d_get_indep_ignore_string_lower(data):
    with pytest.raises(DataErr):
        data.ignore("0", 1)


@pytest.mark.parametrize("data", (Data, ), indirect=True)
def test_data_get_indep_ignore_string_lower(data):
    with pytest.raises(DataErr):
        data.ignore(("0", ), (1, ))


@pytest.mark.parametrize("data", (Data1D, ), indirect=True)
def test_data_1d_get_indep_ignore_string_upper(data):
    with pytest.raises(DataErr):
        data.ignore(0, "1")


@pytest.mark.parametrize("data", (Data, ), indirect=True)
def test_data_get_indep_ignore_string_upper(data):
    with pytest.raises(DataErr):
        data.ignore((0, ), ("1", ))


@pytest.mark.parametrize("data", (Data1D, ), indirect=True)
def test_data_1d_get_indep_notice(data):
    data.notice(0, X_THRESHOLD)
    numpy.testing.assert_array_equal(data.get_indep(filter=True), [X_ARRAY[:X_THRESHOLD + 1], ])


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_get_indep_notice(data):
    data.notice(0, X_THRESHOLD)
    expected = [(X_ARRAY-0.5)[:X_THRESHOLD + 1], (X_ARRAY+0.5)[:X_THRESHOLD + 1]]
    actual = data.get_indep(filter=True)
    numpy.testing.assert_array_equal(actual[0], expected[0])
    numpy.testing.assert_array_equal(actual[1], expected[1])


@pytest.mark.parametrize("data", (Data, ), indirect=True)
def test_data_get_indep_notice(data):
    data.notice((0, ), (X_THRESHOLD, ))
    numpy.testing.assert_array_equal(data.get_indep(filter=True), [X_ARRAY[:X_THRESHOLD + 1], ])


@pytest.mark.parametrize("data", (Data1D, Data), indirect=True)
def test_data_get_indep_mask(data):
    data.mask = X_ARRAY == 0
    numpy.testing.assert_array_equal(data.get_indep(filter=True)[0], X_ARRAY[0])


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_get_indep_mask(data):
    data.mask = X_ARRAY == 0
    numpy.testing.assert_array_equal(data.get_indep(filter=True), ([(X_ARRAY-0.5)[0]], [(X_ARRAY+0.5)[0]]))


@pytest.mark.parametrize("data", (Data1D, Data), indirect=True)
def test_data_get_indep_filter_mask(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    data.mask = X_ARRAY == 0
    numpy.testing.assert_array_equal(data.get_indep(filter=True)[0], [X_ARRAY[0]])


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_get_indep_filter_mask(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    data.mask = X_ARRAY == 0
    numpy.testing.assert_array_equal(data.get_indep(filter=True), ([(X_ARRAY-0.5)[0]], [(X_ARRAY+0.5)[0]]))


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_dep_filter(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    numpy.testing.assert_array_equal(data.get_dep(filter=True), Y_ARRAY[:X_THRESHOLD + 1])


@pytest.mark.parametrize("data", (Data1D, Data1DInt), indirect=True)
def test_data_set_dep_filter(data):
    # This used to be [0, 1] but why would we want this, so check we
    # can call set_dep with the expected argument size.
    #
    data.set_dep([0, 1] * 5)
    numpy.testing.assert_array_equal(data.get_dep(filter=True), [0, 1] * 5)

    # There's also support for scalar values.
    data.set_dep(0)
    numpy.testing.assert_array_equal(data.get_dep(filter=True), [0] * Y_ARRAY.size)


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_staterror(data):
    numpy.testing.assert_array_equal(data.get_staterror(), STATISTICAL_ERROR_ARRAY)


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_staterror_filter(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    numpy.testing.assert_array_equal(data.get_staterror(filter=True), STATISTICAL_ERROR_ARRAY[:X_THRESHOLD + 1])


@pytest.mark.parametrize("data_no_errors", DATA_1D_CLASSES, indirect=True)
def test_data_get_staterror_func(data_no_errors):
    data_no_errors.mask = X_ARRAY <= X_THRESHOLD
    stat_error = data_no_errors.get_staterror(filter=False, staterrfunc=lambda x: MULTIPLIER * x)  # type: numpy.ndarray
    numpy.testing.assert_array_equal(stat_error, MULTIPLIER * Y_ARRAY)


@pytest.mark.parametrize("data_no_errors", DATA_1D_CLASSES, indirect=True)
def test_data_get_staterror_filter_func(data_no_errors):
    data_no_errors.mask = X_ARRAY <= X_THRESHOLD
    stat_error = data_no_errors.get_staterror(filter=True, staterrfunc=lambda x: MULTIPLIER * x)  # type: numpy.ndarray
    numpy.testing.assert_array_equal(stat_error, MULTIPLIER * Y_ARRAY[:X_THRESHOLD + 1])


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_syserror(data):
    numpy.testing.assert_array_equal(data.get_syserror(), SYSTEMATIC_ERROR_ARRAY)


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_syserror_filter(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    numpy.testing.assert_array_equal(data.get_syserror(filter=True), SYSTEMATIC_ERROR_ARRAY[:X_THRESHOLD + 1])


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_error(data):
    error = data.get_error()
    expected_error = numpy.sqrt(SYSTEMATIC_ERROR_ARRAY ** 2 + STATISTICAL_ERROR_ARRAY ** 2)
    numpy.testing.assert_array_equal(error, expected_error)


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_yerr(data):
    error = data.get_yerr()
    expected_error = numpy.sqrt(SYSTEMATIC_ERROR_ARRAY ** 2 + STATISTICAL_ERROR_ARRAY ** 2)
    numpy.testing.assert_array_equal(error, expected_error)


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_dep(data):
    numpy.testing.assert_array_equal(data.get_dep(), Y_ARRAY)


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_y(data):
    numpy.testing.assert_array_equal(data.get_y(), Y_ARRAY)


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_get_y_filter(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    numpy.testing.assert_array_equal(data.get_y(filter=True), Y_ARRAY[:X_THRESHOLD + 1])


@pytest.mark.parametrize("data", (Data1D, Data), indirect=True)
def test_data_get_y_filter_func(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    y = data.get_y(filter=True, yfunc=lambda x: MULTIPLIER*x)
    expected_y = (Y_ARRAY[:X_THRESHOLD + 1], MULTIPLIER*X_ARRAY[:X_THRESHOLD + 1])
    numpy.testing.assert_array_equal(y[0], expected_y[0])
    numpy.testing.assert_array_equal(y[1], expected_y[1])


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_get_y_filter_func(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    y = data.get_y(filter=True, yfunc=lambda x, y: (MULTIPLIER*x, MULTIPLIER*y))
    expected_y = (Y_ARRAY[:X_THRESHOLD + 1], (MULTIPLIER*(X_ARRAY-0.5)[:X_THRESHOLD + 1],
                  MULTIPLIER*(X_ARRAY+0.5)[:X_THRESHOLD + 1]))
    numpy.testing.assert_array_equal(y[0], expected_y[0])
    numpy.testing.assert_array_equal(y[1], expected_y[1])


@pytest.mark.parametrize("data", (Data1D, Data), indirect=True)
def test_data_get_y_func(data):
    y = data.get_y(filter=True, yfunc=lambda x: MULTIPLIER*x)
    expected_y = (Y_ARRAY, MULTIPLIER*X_ARRAY)
    numpy.testing.assert_array_equal(y[0], expected_y[0])
    numpy.testing.assert_array_equal(y[1], expected_y[1])


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_get_y_func(data):
    y = data.get_y(filter=True, yfunc=lambda x, y: (MULTIPLIER*x, MULTIPLIER*y))
    expected_y = (Y_ARRAY, (MULTIPLIER*(X_ARRAY-0.5), MULTIPLIER*(X_ARRAY+0.5)))
    numpy.testing.assert_array_equal(y[0], expected_y[0])
    numpy.testing.assert_array_equal(y[1], expected_y[1])


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_eval_model(data):
    model = Polynom1D()
    model.c0 = 0
    model.c1 = MULTIPLIER
    evaluated_data = data.eval_model(model)
    numpy.testing.assert_array_equal(evaluated_data, MULTIPLIER * X_ARRAY)


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_eval_model_to_fit_no_filter(data):
    model = Polynom1D()
    model.c0 = 0
    model.c1 = MULTIPLIER
    evaluated_data = data.eval_model_to_fit(model)
    numpy.testing.assert_array_equal(evaluated_data, MULTIPLIER * X_ARRAY)


@pytest.mark.parametrize("data", (Data1D, Data), indirect=True)
def test_data_eval_model_to_fit_filter(data):
    model = Polynom1D()
    model.c0 = 0
    model.c1 = MULTIPLIER
    data.mask = X_ARRAY <= X_THRESHOLD
    evaluated_data = data.eval_model_to_fit(model)
    numpy.testing.assert_array_equal(evaluated_data, MULTIPLIER * X_ARRAY[:X_THRESHOLD + 1])


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_eval_model_to_fit_filter(data):
    model = Polynom1D()
    model.c0 = 0
    model.c1 = MULTIPLIER
    data.mask = X_ARRAY <= X_THRESHOLD
    evaluated_data = data.eval_model_to_fit(model)
    numpy.testing.assert_array_equal(evaluated_data, MULTIPLIER * X_ARRAY[:X_THRESHOLD + 1])


@pytest.mark.parametrize("data", (Data1D, Data), indirect=True)
def test_data_to_guess(data):
    actual = data.to_guess()
    expected = [Y_ARRAY, X_ARRAY]
    numpy.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_to_guess(data):
    actual = data.to_guess()
    expected = [Y_ARRAY, X_ARRAY-0.5]
    numpy.testing.assert_array_equal(actual[0], expected[0])
    numpy.testing.assert_array_equal(actual[1], expected[1])


@pytest.mark.parametrize("data", DATA_1D_CLASSES, indirect=True)
def test_data_1d_to_fit(data):
    actual = data.to_fit()
    expected = [Y_ARRAY, STATISTICAL_ERROR_ARRAY, SYSTEMATIC_ERROR_ARRAY]
    numpy.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("data", (Data1D, ), indirect=True)
def test_data_1d_to_plot(data):
    actual = data.to_plot()
    yerr = numpy.sqrt(SYSTEMATIC_ERROR_ARRAY ** 2 + STATISTICAL_ERROR_ARRAY ** 2)
    expected = [X_ARRAY, Y_ARRAY, yerr, None, "x", "y"]
    numpy.testing.assert_array_equal(actual[0], expected[0])
    numpy.testing.assert_array_equal(actual[1], expected[1])
    numpy.testing.assert_array_equal(actual[2], expected[2])
    numpy.testing.assert_array_equal(actual[3], expected[3])
    numpy.testing.assert_array_equal(actual[4], expected[4])
    numpy.testing.assert_array_equal(actual[5], expected[5])


@pytest.mark.parametrize("data", (Data1D, ), indirect=True)
def test_data_1d_to_component_plot(data):
    actual = data.to_component_plot()
    yerr = numpy.sqrt(SYSTEMATIC_ERROR_ARRAY ** 2 + STATISTICAL_ERROR_ARRAY ** 2)
    expected = [X_ARRAY, Y_ARRAY, yerr, None, "x", "y"]
    numpy.testing.assert_array_equal(actual[0], expected[0])
    numpy.testing.assert_array_equal(actual[1], expected[1])
    numpy.testing.assert_array_equal(actual[2], expected[2])
    numpy.testing.assert_array_equal(actual[3], expected[3])
    numpy.testing.assert_array_equal(actual[4], expected[4])
    numpy.testing.assert_array_equal(actual[5], expected[5])


@pytest.mark.parametrize("data", (Data, Data1D, Data1DInt), indirect=True)
def test_data_to_contour(data):
    with pytest.raises(AttributeError):
        data.to_contour()


@pytest.mark.parametrize("data", (Data, Data2D, Data2DInt), indirect=True)
def test_data_to_plot(data):
    with pytest.raises(AttributeError):
        data.to_plot()


@pytest.mark.parametrize("data", (Data, Data2D, Data2DInt), indirect=True)
def test_data_to_component_plot(data):
    with pytest.raises(AttributeError):
        data.to_component_plot()


def test_data_simul_fit(data_simul_fit):
    y, stat_error, systematic_error = data_simul_fit.to_fit()
    expected_y = numpy.concatenate((Y_ARRAY, MULTIPLIER * Y_ARRAY))
    expected_stat_error = numpy.concatenate((STATISTICAL_ERROR_ARRAY, MULTIPLIER * STATISTICAL_ERROR_ARRAY))
    expected_sys_error = numpy.concatenate((SYSTEMATIC_ERROR_ARRAY, MULTIPLIER * SYSTEMATIC_ERROR_ARRAY))
    numpy.testing.assert_array_equal(y, expected_y)
    numpy.testing.assert_array_equal(stat_error, expected_stat_error)
    numpy.testing.assert_array_equal(systematic_error, expected_sys_error)


def test_data_simul_fit_to_plot(data_simul_fit):
    actual = data_simul_fit.to_fit()
    expected_y = numpy.concatenate((Y_ARRAY, MULTIPLIER * Y_ARRAY))
    expected_stat_error = numpy.concatenate((STATISTICAL_ERROR_ARRAY, MULTIPLIER * STATISTICAL_ERROR_ARRAY))
    expected_sys_error = numpy.concatenate((SYSTEMATIC_ERROR_ARRAY, MULTIPLIER * SYSTEMATIC_ERROR_ARRAY))
    numpy.testing.assert_array_equal(actual[0], expected_y)
    numpy.testing.assert_array_equal(actual[1], expected_stat_error)
    numpy.testing.assert_array_equal(actual[2], expected_sys_error)


def test_data_simul_fit_no_errors(data_simul_fit_no_errors):
    y, stat_error, systematic_error = data_simul_fit_no_errors.to_fit()
    expected_y = numpy.concatenate((Y_ARRAY, MULTIPLIER * Y_ARRAY))
    numpy.testing.assert_array_equal(y, expected_y)
    assert stat_error is None
    assert systematic_error is None


def test_data_simul_fit_some_errors(data_simul_fit_some_errors):
    with pytest.raises(DataErr):
        data_simul_fit_some_errors.to_fit()


def test_data_simul_fit_eval_model_to_fit(data_simul_fit):
    model = Polynom1D()
    model.c0 = 0
    model.c1 = MULTIPLIER
    data_simul_fit.datasets[0].mask = X_ARRAY <= X_THRESHOLD
    data_simul_fit.datasets[1].mask = X_ARRAY <= X_THRESHOLD
    evaluated_data = data_simul_fit.eval_model_to_fit((model, model))
    expected_data = numpy.concatenate((MULTIPLIER * X_ARRAY[:X_THRESHOLD+1],
                                       MULTIPLIER**2 * X_ARRAY[:X_THRESHOLD+1]))
    numpy.testing.assert_array_equal(evaluated_data, expected_data)


@pytest.mark.parametrize("data", (Data1D,), indirect=True)
def test_data1d_get_dims(data):
    assert data.get_dims() == (X_ARRAY.size, )


@pytest.mark.parametrize("data", (Data1D,), indirect=True)
def test_data1d_get_filter(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    assert data.get_filter() == '0.0000:3.0000'


@pytest.mark.parametrize("data", (Data1D,), indirect=True)
def test_data1d_get_filter_mask(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    assert data.get_filter() == '0.0000:3.0000'


@pytest.mark.parametrize("data", (Data1D,), indirect=True)
def test_data1d_get_filter_expr(data):
    data.mask = X_ARRAY <= X_THRESHOLD
    assert data.get_filter_expr() == '0.0000-3.0000 x'


@pytest.mark.parametrize("data", (Data1D,), indirect=True)
def test_data1d_get_bounding_mask(data):
    mask = X_ARRAY <= X_THRESHOLD
    data.mask = mask
    actual = data.get_bounding_mask()
    numpy.testing.assert_array_equal(actual[0], mask)
    numpy.testing.assert_array_equal(actual[1], X_ARRAY.size)


@pytest.mark.parametrize("data", (Data1D,), indirect=True)
def test_data1d_get_img(data):
    numpy.testing.assert_array_equal(data.get_img(), [Y_ARRAY, ])


@pytest.mark.parametrize("data", (Data1D,), indirect=True)
def test_data1d_get_img_yfunc(data):
    actual = data.get_img(yfunc=lambda x: MULTIPLIER * x)
    expected = ([Y_ARRAY, ], [MULTIPLIER * X_ARRAY, ], )
    numpy.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("data", (Data1D,), indirect=True)
def test_data1d_get_imgerr(data):
    expected_error = numpy.sqrt(SYSTEMATIC_ERROR_ARRAY ** 2 + STATISTICAL_ERROR_ARRAY ** 2)
    numpy.testing.assert_array_equal(data.get_imgerr(), [expected_error, ])


@pytest.mark.parametrize("data", (Data1D, Data1DInt), indirect=True)
def test_data1d_get_x(data):
    numpy.testing.assert_array_equal(data.get_x(), X_ARRAY)


@pytest.mark.parametrize("data", (Data1D, ), indirect=True)
def test_data1d_get_xerr(data):
    assert data.get_xerr() is None


@pytest.mark.parametrize("data", (Data1DInt, ), indirect=True)
def test_data_1d_int_get_xerr(data):
    numpy.testing.assert_array_equal(data.get_xerr(), [1] * X_ARRAY.size)


@pytest.mark.parametrize("data", (Data1D, Data1DInt), indirect=True)
def test_data1d_get_y(data):
    numpy.testing.assert_array_equal(data.get_y(), Y_ARRAY)


@pytest.mark.parametrize("data", (Data2D, ), indirect=True)
def test_data2_get_x0(data):
    numpy.testing.assert_array_equal(data.get_x0(), X0_2D)


@pytest.mark.parametrize("data", (Data2DInt, ), indirect=True)
def test_data2_int_get_x0(data):
    numpy.testing.assert_array_equal(data.get_x0(), X0_2D)


@pytest.mark.parametrize("data", (Data2D, ), indirect=True)
def test_data2_get_x1(data):
    numpy.testing.assert_array_equal(data.get_x1(), X1_2D)


@pytest.mark.parametrize("data", (Data2DInt, ), indirect=True)
def test_data2_int_get_x1(data):
    numpy.testing.assert_array_equal(data.get_x1(), X1_2D)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_dims(data):
    assert data.get_dims() == (X_ARRAY.size, X_ARRAY.size)


# DATA-NOTE: Not sure this should work, really, as the 1D implementation does not account for the difference in 2D
#  data, but in 2D it is hard with the current implementation to figure out the shape is self.shape is None
@pytest.mark.xfail()
@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_dims_no_shape(data):
    data.shape = None
    assert data.get_dims() == (X_ARRAY.size, X_ARRAY.size)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_axes(data):
    numpy.testing.assert_array_equal(data.get_axes(), (X_ARRAY+1, X_ARRAY+1))


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_img(data):
    numpy.testing.assert_array_equal(data.get_img(), Y_2D_RAW)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_imgerr(data):
    expected_error = numpy.sqrt(STAT_ERROR_2D ** 2 + SYS_ERROR_2D ** 2).reshape(SHAPE_2D)
    numpy.testing.assert_array_equal(data.get_imgerr(), expected_error)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_xerr(data):
    with pytest.raises(AttributeError):
        data.get_xerr()


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_dep_filter(data):
    test_filter = X0_2D <= X_THRESHOLD
    data.mask = test_filter
    numpy.testing.assert_array_equal(data.get_dep(filter=True), Y_2D[test_filter])


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_staterror(data):
    numpy.testing.assert_array_equal(data.get_staterror(), STAT_ERROR_2D)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_staterror_filter(data):
    test_filter = X0_2D <= X_THRESHOLD
    data.mask = test_filter
    numpy.testing.assert_array_equal(data.get_staterror(filter=True), STAT_ERROR_2D[test_filter])


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_syserror(data):
    numpy.testing.assert_array_equal(data.get_syserror(), SYS_ERROR_2D)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_syserror_filter(data):
    test_filter = X0_2D <= X_THRESHOLD
    data.mask = test_filter
    numpy.testing.assert_array_equal(data.get_syserror(filter=True), SYS_ERROR_2D[test_filter])


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_error(data):
    error = data.get_error()
    expected_error = numpy.sqrt(SYS_ERROR_2D ** 2 + STAT_ERROR_2D ** 2)
    numpy.testing.assert_array_equal(error, expected_error)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_yerr(data):
    error = data.get_yerr()
    expected_error = numpy.sqrt(SYS_ERROR_2D ** 2 + STAT_ERROR_2D ** 2)
    numpy.testing.assert_array_equal(error, expected_error)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_dep(data):
    numpy.testing.assert_array_equal(data.get_dep(), Y_2D)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_y(data):
    numpy.testing.assert_array_equal(data.get_y(), Y_2D)


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_y_filter(data):
    test_filter = X0_2D <= X_THRESHOLD
    data.mask = test_filter
    numpy.testing.assert_array_equal(data.get_y(filter=True), Y_2D[test_filter])


@pytest.mark.parametrize("data", (Data2D, ), indirect=True)
def test_data2_get_y_filter_func(data):
    test_filter = X0_2D <= X_THRESHOLD
    data.mask = test_filter
    y = data.get_y(filter=True, yfunc=lambda x0, x1: MULTIPLIER*(x0 + x1))
    expected_y = Y_2D[test_filter], (MULTIPLIER * (X0_2D + X1_2D))[test_filter]
    numpy.testing.assert_array_equal(y[0], expected_y[0])
    numpy.testing.assert_array_equal(y[1], expected_y[1])


@pytest.mark.parametrize("data", (Data2D, ), indirect=True)
def test_data2_get_img_func(data):
    y = data.get_img(yfunc=lambda x0, x1: MULTIPLIER*(x0 + x1))
    expected_y = Y_2D_RAW, MULTIPLIER * (X0_2D + X1_2D).reshape(data.shape)
    numpy.testing.assert_array_equal(y[0], expected_y[0])
    numpy.testing.assert_array_equal(y[1], expected_y[1])


@pytest.mark.parametrize("data", (Data2DInt, ), indirect=True)
def test_data2_int_get_y_filter_func(data):
    test_filter = X0_2D <= X_THRESHOLD
    data.mask = test_filter
    y = data.get_y(filter=True, yfunc=lambda x0lo, x0hi, x1lo, x1hi: MULTIPLIER*((x0lo+x0hi)/2 + (x1lo+x1hi)/2))
    expected_y = Y_2D[test_filter], (MULTIPLIER * (X0_2D + X1_2D))[test_filter]
    numpy.testing.assert_array_equal(y[0], expected_y[0])
    numpy.testing.assert_array_equal(y[1], expected_y[1])


@pytest.mark.parametrize("data", (Data2DInt, ), indirect=True)
def test_data2_int_get_img_func(data):
    y = data.get_img(yfunc=lambda x0lo, x0hi, x1lo, x1hi: MULTIPLIER*((x0lo+x0hi)/2 + (x1lo+x1hi)/2))
    expected_y = Y_2D_RAW, MULTIPLIER * (X0_2D + X1_2D).reshape(data.shape)
    numpy.testing.assert_array_equal(y[0], expected_y[0])
    numpy.testing.assert_array_equal(y[1], expected_y[1])


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_eval_model(data):
    model = Polynom2D()
    model.c = 0
    model.cy1 = MULTIPLIER
    model.cx1 = MULTIPLIER
    evaluated_data = data.eval_model(model)
    numpy.testing.assert_array_equal(evaluated_data, MULTIPLIER * (X0_2D + X1_2D))


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_eval_model_to_fit_no_filter(data):
    model = Polynom2D()
    model.c = 0
    model.cy1 = MULTIPLIER
    model.cx1 = MULTIPLIER
    evaluated_data = data.eval_model_to_fit(model)
    numpy.testing.assert_array_equal(evaluated_data, MULTIPLIER * (X0_2D + X1_2D))


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_eval_model_to_fit_filter(data):
    model = Polynom2D()
    model.c = 0
    model.cy1 = MULTIPLIER
    model.cx1 = MULTIPLIER
    test_filter = X0_2D <= X_THRESHOLD
    data.mask = test_filter
    evaluated_data = data.eval_model_to_fit(model)
    numpy.testing.assert_array_equal(evaluated_data, (MULTIPLIER * (X0_2D + X1_2D))[test_filter])


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_max_pos(data):
    numpy.testing.assert_array_equal(data.get_max_pos(), (X_ARRAY.size-1, X_ARRAY.size-1))


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_data2_get_max_pos_dep(data):
    dep = 1/(Y_2D+1)  # +1 to avoid dividing by zero
    numpy.testing.assert_array_equal(data.get_max_pos(dep=dep), (0, 0))


# DATA-NOTE: This is failing because Data2D.notice isn't implemented correctly and it just combines the
# Masks on the two axes into one, i.e. mask_x0 && mask_x1 is applied to both axes.
# We probably never noticed because DataIMG defines a notice2d method which we always use.
@pytest.mark.xfail
@pytest.mark.parametrize("data", (Data2D, ), indirect=True)
def test_data2_get_indep_notice(data):
    test_filter_0 = X0_2D <= X_THRESHOLD
    test_filter_1 = X1_2D <= X_THRESHOLD + 1
    data.notice(0, X_THRESHOLD, 0, X_THRESHOLD + 1)
    expected = [X0_2D[test_filter_0], X1_2D[test_filter_1]]
    actual = data.get_indep(filter=True)
    numpy.testing.assert_array_equal(actual[0], expected[0])
    numpy.testing.assert_array_equal(actual[1], expected[1])


# DATA-NOTE: This is failing for a different reason (can't get_indep(filter=True) in the first place).
# Not sure whether I am doing something wrong, but it's unlikely, since the Data2DInt.notice()
# signature seems consistent with what I am doing. In any case the problem is that at some point the
# Data2DInt.mask is a (10, 10) array, while the shape of the data is (100, )
@pytest.mark.xfail
@pytest.mark.parametrize("data", (Data2DInt, ), indirect=True)
def test_data2_int_get_indep_notice(data):
    test_filter_0 = X0_2D <= X_THRESHOLD
    test_filter_1 = X1_2D <= X_THRESHOLD + 1
    data.notice(0, X_THRESHOLD, 0, X_THRESHOLD + 1)
    expected = [(X0_2D - 0.5)[test_filter_0],
                (X0_2D + 0.5)[test_filter_0],
                (X1_2D - 0.5)[test_filter_1],
                (X1_2D + 0.5)[test_filter_1]]
    actual = data.get_indep(filter=True)
    numpy.testing.assert_array_equal(actual[0], expected[0])
    numpy.testing.assert_array_equal(actual[1], expected[1])


# DATA-NOTE: This is just a notice call in disguise, so it's failing like just above.
@pytest.mark.xfail
@pytest.mark.parametrize("data", (Data2D, ), indirect=True)
def test_data2_get_indep_ignore(data):
    test_filter_0 = X0_2D > X_THRESHOLD
    test_filter_1 = X1_2D > X_THRESHOLD + 1
    data.ignore(0, X_THRESHOLD, 0, X_THRESHOLD + 1)
    expected = [X0_2D[test_filter_0], X1_2D[test_filter_1]]
    actual = data.get_indep(filter=True)
    numpy.testing.assert_array_equal(actual[0], expected[0])
    numpy.testing.assert_array_equal(actual[1], expected[1])


@pytest.mark.parametrize("data", (Data2D, ), indirect=True)
def test_data2_get_indep_mask(data):
    test_filter = X0_2D == 0
    data.mask = test_filter
    expected = (X0_2D[test_filter], X1_2D[test_filter])
    numpy.testing.assert_array_equal(data.get_indep(filter=True), expected)


# DATA-NOTE: this fails because get_indep() does not work. Either I am missing something fundamental
# or the Data2DInt methods are bogus
@pytest.mark.xfail
@pytest.mark.parametrize("data", (Data2DInt, ), indirect=True)
def test_data2_int_get_indep_mask(data):
    test_filter = X0_2D == 0
    data.mask = test_filter
    expected = (X0_2D[test_filter], X1_2D[test_filter])
    numpy.testing.assert_array_equal(data.get_indep(filter=True), expected)


@pytest.fixture
def array_sizes_fixture():
    x0low, x0high = 3000, 4000
    x1low, x1high = 4000, 4800
    dx = 500
    x1, x0 = numpy.mgrid[x1low:x1high:dx, x0low:x0high:dx]
    y = (x0 - 3500) ** 2 + (x1 - 4500) ** 2
    return x0, x1, dx, y


# https://github.com/sherpa/sherpa/issues/627
def test_data2d_wrong_array_size(array_sizes_fixture):
    x0, x1, dx, y = array_sizes_fixture

    with pytest.raises(TypeError):
        Data2D('name', x0, x1, y.flatten(), staterror=numpy.sqrt(y).flatten())


def test_data2d_wrong_y_array_size(array_sizes_fixture):
    x0, x1, dx, y = array_sizes_fixture

    with pytest.raises(TypeError):
        Data2D('name', x0.flatten(), x1.flatten(), y, staterror=numpy.sqrt(y).flatten())


def test_data2d_int_wrong_array_size(array_sizes_fixture):
    x0, x1, dx, y = array_sizes_fixture

    with pytest.raises(TypeError):
        Data2DInt('name', x0, x0, x1, x1, y.flatten(), staterror=numpy.sqrt(y).flatten())


def test_data2d_int_wrong_y_array_size(array_sizes_fixture):
    x0, x1, dx, y = array_sizes_fixture

    with pytest.raises(TypeError):
        Data2DInt('name', x0.flatten(), x0.flatten(), x1.flatten(), x1.flatten(), y, staterror=numpy.sqrt(y).flatten())


# https://github.com/sherpa/sherpa/issues/628
def test_data2d_int_eval_model_to_fit(array_sizes_fixture):
    from sherpa.fit import Fit
    from sherpa.optmethods import LevMar
    from sherpa.stats import Chi2
    from sherpa.models import Gauss2D

    x0, x1, dx, y = array_sizes_fixture
    data2 = Data2DInt('name', x0.flatten(), x1.flatten(),
                      x0.flatten() + dx, x1.flatten() + dx,
                      y.flatten(),
                      staterror=numpy.sqrt(y).flatten())

    model2 = Gauss2D()
    fitter = Fit(data2, model2, Chi2(), LevMar())
    fitter.fit()  # Failed in Sherpa 4.11.0


# https://github.com/sherpa/sherpa/issues/695
@pytest.mark.parametrize('arrpos', [POS_X_ARRAY])
@pytest.mark.parametrize("Dataclass", ALL_DATA_CLASSES)
def test_data_indep_masked_numpyarray(arrpos, Dataclass):
    i = arrpos[Dataclass]
    args = list(INSTANCE_ARGS[Dataclass])
    mask = numpy.random.rand(*(args[i].shape)) > 0.5
    args[i] = numpy.ma.array(args[i], mask=mask)
    with pytest.warns(UserWarning, match="for dependent variables only"):
        data = Dataclass(*args)
    assert len(data.get_dep(filter=True)) == len(args[POS_Y_ARRAY[Dataclass]])


@pytest.mark.parametrize('arrpos', [POS_STATERR_ARRAY, POS_SYSERR_ARRAY])
@pytest.mark.parametrize("Dataclass", ALL_DATA_CLASSES)
def test_data_err_masked_numpyarray(arrpos, Dataclass):
    i = arrpos[Dataclass]
    args = list(INSTANCE_ARGS[Dataclass])
    mask = numpy.random.rand(*(args[i].shape)) > 0.5
    args[i] = numpy.ma.array(args[i], mask=mask)
    with pytest.warns(UserWarning, match="differs from the mask of the dependent array"):
        data = Dataclass(*args)
    assert len(data.get_dep(filter=True)) == len(args[POS_Y_ARRAY[Dataclass]])


@pytest.mark.parametrize('arrpos', [POS_STATERR_ARRAY, POS_SYSERR_ARRAY])
@pytest.mark.parametrize("Dataclass", ALL_DATA_CLASSES)
def test_data_deperr_masked_numpyarray(arrpos, Dataclass):
    '''Error arrays can be masked as long as that mask is the same as the dependent array'''
    i = arrpos[Dataclass]
    j = POS_Y_ARRAY[Dataclass]
    args = list(INSTANCE_ARGS[Dataclass])
    mask = numpy.random.rand(*(args[i].shape)) > 0.5
    args[i] = numpy.ma.array(args[i], mask=mask)
    args[j] = numpy.ma.array(args[j], mask=mask)
    data = Dataclass(*args)
    assert len(data.get_dep(filter=True)) == (~mask).sum()


@pytest.mark.parametrize("Dataclass", REALLY_ALL_DATA_CLASSES)
def test_data_dep_masked_numpyarray(Dataclass):
    args = list(INSTANCE_ARGS[Dataclass])
    posy = POS_Y_ARRAY[Dataclass]
    mask = numpy.random.rand(*(args[posy].shape)) > 0.5
    args[posy] = numpy.ma.array(args[posy], mask=mask)
    data = Dataclass(*args)
    assert data.mask.shape == mask.shape
    assert numpy.all(data.mask == ~mask)
    assert len(data.get_dep(filter=True)) == (~mask).sum()


@pytest.mark.parametrize("Dataclass", REALLY_ALL_DATA_CLASSES)
def test_data_dep_masked_numpyarray_nomask(Dataclass):
    args = list(INSTANCE_ARGS[Dataclass])
    posy = POS_Y_ARRAY[Dataclass]
    # By default, numpy creates a masked array with no mask set
    args[posy] = numpy.ma.array(args[posy])
    data = Dataclass(*args)
    # Sherpa's way of saying "mask is not set"
    assert data.mask is True
    assert len(data.get_dep(filter=True)) == len(args[posy].flatten())


@pytest.mark.parametrize("Dataclass", ALL_DATA_CLASSES)
def test_data_indep_anyobj_with_mask(Dataclass):
    args = list(INSTANCE_ARGS[Dataclass])

    class DummyMask(list):
        mask = 'whatisthis'
    args[1] = DummyMask(args[1])
    with pytest.warns(UserWarning, match="for dependent variables only"):
        data = Dataclass(*args)
    assert data.mask is True
    assert len(data.get_dep(filter=True)) == len(args[POS_Y_ARRAY[Dataclass]])


@pytest.mark.parametrize("Dataclass", REALLY_ALL_DATA_CLASSES)
def test_data_dep_any_obj_with_mask(Dataclass):
    args = list(INSTANCE_ARGS[Dataclass])
    posy = POS_Y_ARRAY[Dataclass]

    class DummyMask(list):
        mask = 'whatisthis'
    args[posy] = DummyMask(args[posy])
    with pytest.warns(UserWarning, match="Set .mask"):
        data = Dataclass(*args)
    assert data.mask is True
    assert len(data.get_dep(filter=True)) == len(data.get_dep(filter=False))


# repeat set of tests except now by using ui
# Results should be idendical, but tests are fast, so we just test again
# To make sure that there is no heuristic in load_arrays or similar that
# interferes with the logic
@pytest.mark.parametrize('arrpos', [POS_X_ARRAY, POS_STATERR_ARRAY, POS_SYSERR_ARRAY])
@pytest.mark.parametrize('Session', [Session, AstroSession])
@pytest.mark.parametrize("data_for_load_arrays", ALL_DATA_CLASSES, indirect=True)
def test_data_indeperr_masked_numpyarray_ui(arrpos, data_for_load_arrays, Session):
    session, args, data = data_for_load_arrays
    session = Session()
    i = arrpos[type(data)]
    mask = numpy.random.rand(*(args[i].shape)) > 0.5
    args = list(args)
    args[1] = numpy.ma.array(args[i], mask=mask)
    with pytest.warns(UserWarning, match="for dependent variables only"):
        session.load_arrays(*args)
    new_data = session.get_data(data.name)
    assert len(new_data.get_dep(filter=True)) == len(args[i])


@pytest.mark.parametrize('Session', [Session, AstroSession])
@pytest.mark.parametrize("data_for_load_arrays", ALL_DATA_CLASSES, indirect=True)
def test_data_dep_masked_numpyarray_ui(data_for_load_arrays, Session):
    session, args, data = data_for_load_arrays
    session = Session()
    posy = POS_Y_ARRAY[type(data)]
    mask = numpy.random.rand(*(args[posy].shape)) > 0.5
    args = list(args)
    args[posy] = numpy.ma.array(args[posy], mask=mask)
    session.load_arrays(*args)
    new_data = session.get_data(data.name)
    assert new_data.mask.shape == mask.shape
    assert numpy.all(new_data.mask == ~mask)
    assert len(new_data.get_dep(filter=True)) == (~mask).sum()


@pytest.mark.parametrize('Session', [Session, AstroSession])
@pytest.mark.parametrize("data_for_load_arrays", ALL_DATA_CLASSES, indirect=True)
def test_data_dep_masked_numpyarray_nomask_ui(data_for_load_arrays, Session):
    session, args, data = data_for_load_arrays
    session = Session()
    posy = POS_Y_ARRAY[type(data)]
    args = list(args)
    args[posy] = numpy.ma.array(args[posy])
    session.load_arrays(*args)
    new_data = session.get_data(data.name)
    # Sherpa's way of saying "mask is not set"
    assert new_data.mask is True
    assert len(new_data.get_dep(filter=True)) == len(args[posy].flatten())


# https://github.com/sherpa/sherpa/issues/346
@pytest.mark.parametrize('Session', [Session, AstroSession])
def test_regression_346(Session):
    session = Session()
    x = numpy.arange(-5, 5.1)
    old_y = x*x + 23.2
    y = numpy.ma.masked_array(old_y, mask=old_y < 35)
    e = numpy.ones(x.size)
    session.load_arrays("mydata", x, y, e)
    filtered_y = session. get_dep("mydata", filter=True)
    assert numpy.allclose(filtered_y, [48.2, 39.2, 39.2, 48.2])


def test_manual_setting_mask():
    d = Data1D(name='test', x=[1, 2, 3], y=[0, 0, 0])
    d.mask = True
    assert len(d.get_dep(filter=True)) == 3

    d.mask = False
    # This test looks like it does not do anything, but in fact "mask"
    # is a property with complext logic, so the fact that setting it to
    # False makes is False is non-trivial.
    # I don't want to test for
    # len(d.get_dep(filter=True)) == 0
    # because the get_dep raises and error when no data is noticed
    # and I don't want to test get_dep here, but the fact that setting
    # the mask itself works.
    assert d.mask is False

    d.mask = [True, False, True]
    assert len(d.get_dep(filter=True)) == 2
    arr = numpy.ma.array([3, 4, 5])
    # aka numpy.ma.nomask, but used in a more natural way
    d.mask = arr.mask
    assert len(d.get_dep(filter=True)) == 3

    with pytest.raises(DataErr) as e:
        d.mask = None
    assert 'True, False, or a mask array' in str(e.value)


def test_data_filter_no_data():
    """Check we get a excludes-all-data error"""

    x = numpy.asarray([1, 2, 5])
    d = Data1D('x', x, x)
    assert d.mask
    d.ignore()
    assert d.mask is False

    with pytest.raises(DataErr) as de:
        d.apply_filter([1, 2, 3])

    assert str(de.value) == 'mask excludes all data'


@pytest.mark.parametrize("vals", [4, [4], [[2, 3, 4]], [[2, 3], [3, 2]]])
def test_data_filter_invalid_size(vals):
    """Check we get a size-mismatch error"""

    x = numpy.asarray([1, 2, 5])
    d = Data1D('x', x, x)
    d.ignore(None, 2)
    assert d.mask == pytest.approx([False, False, True])

    with pytest.raises(DataErr) as de:
        d.apply_filter(vals)

    assert str(de.value) == 'size mismatch between mask and data array'


@pytest.mark.parametrize("lo,hi,emsg", [("1:20", None, 'lower'), (None, "2", 'upper'), ("0.5", "7", 'lower')])
@pytest.mark.parametrize("ignore", [False, True])
def test_data1d_notice_errors_out_on_string_range(lo, hi, emsg, ignore):
    """Check we get an error if lo or hi are strings."""

    xlo = numpy.asarray([1, 2, 5])
    xhi = numpy.asarray([2, 3, 8])
    y = numpy.zeros(3)
    d = Data1D('tmp', xlo, xhi, y)
    with pytest.raises(DataErr) as de:
        d.notice(lo, hi, ignore=ignore)

    err = f'strings not allowed in {emsg} bound list'
    assert str(de.value) == err


@pytest.mark.parametrize("expected,args",
                         [('2.0:20.0', []),
                          ('', [(False, 1, 30)]),
                          ('10.0:17.0', [(True, 7.1, 18)]),
                          ('10.0:12.0,17.0', [(True, 7.1, 18), (False, 13, 16)]),
                          ('2.0:12.0,17.0', [(True, 7.1, 18), (False, 13, 16), (True, 0, 12)]),
                          ('10.0:12.0,17.0:20.0', [(True, 7.1, 18), (False, 13, 16), (True, 15.5, 30)]),
                          ('', [(True, 7.1, 18), (False, 13, 16), (True, 6, 17), (False, 1, 40)]),
                          ('2.0:20.0', [(True, 7.1, 18), (False, 13, 16), (True, 6, 17), (True, 1, 40)]),
                         ])
def test_data1d_get_filter_calls(expected, args):
    """Basic check of get_filter

    expected is the expected response
    args is a list of 3-tuples of (flag, loval, hival) where
    flag is True for notice and False for ignore; they define
    the filter to apply
    """

    xs = numpy.asarray([2, 5, 10, 12, 15, 17, 20])
    ys = numpy.ones(xs.size)

    d = Data1D('data', xs, ys)

    for (flag, lo, hi) in args:
        if flag:
            d.notice(lo, hi)
        else:
            d.ignore(lo, hi)

    assert d.get_filter(format='%.1f') == expected


@pytest.mark.parametrize("expected,args",
                         [('2.0:25.0', []),
                          ('', [(False, 1, 30)]),
                          ('5.0:20.0', [(True, 7.1, 18)]),
                          ('5.0:12.0,17.0:20.0', [(True, 7.1, 18), (False, 13, 16)]),
                          ('2.0:12.0,17.0:20.0', [(True, 7.1, 18), (False, 13, 16), (True, 0, 12)]),
                          ('5.0:12.0,15.0:25.0', [(True, 7.1, 18), (False, 13, 16), (True, 15.5, 30)]),
                          ('', [(True, 7.1, 18), (False, 13, 16), (True, 6, 17), (False, 1, 40)]),
                          ('2.0:25.0', [(True, 7.1, 18), (False, 13, 16), (True, 6, 17), (True, 1, 40)]),
                         ])
def test_data1dint_get_filter_calls(expected, args):
    """Basic check of get_filter

    expected is the expected response
    args is a list of 3-tuples of (flag, loval, hival) where
    flag is True for notice and False for ignore; they define
    the filter to apply
    """

    # Note this is not a contiguous grid
    xlos = numpy.asarray([2, 5, 10, 12, 15, 17, 20])
    xhis = numpy.asarray([5, 8, 12, 15, 16, 20, 25])

    ys = numpy.ones(xlos.size)

    d = Data1DInt('data', xlos, xhis, ys)

    for (flag, lo, hi) in args:
        if flag:
            d.notice(lo, hi)
        else:
            d.ignore(lo, hi)

    assert d.get_filter(format='%.1f') == expected


def test_data1dint_get_x_xerr():
    """Check get_x/get_xerr when filtering

    This was added because there was a bug when all data had been
    filtered out. It is essentially the same as
    test_data1dint_get_filter_calls since get_filter calls get_x,
    but it does add explicit checks and a check of get_xerr.

    """

    # Note this is not a contiguous grid
    xlos = numpy.asarray([2, 5, 10, 12, 15, 17, 20])
    xhis = numpy.asarray([5, 8, 12, 15, 16, 20, 25])

    ys = numpy.ones(xlos.size)

    d = Data1DInt('data', xlos, xhis, ys)

    x = [3.5, 6.5, 11, 13.5, 15.5, 18.5, 22.5]
    xerr = xhis - xlos
    assert d.get_x() == pytest.approx(x)
    assert d.get_xerr() == pytest.approx(xerr)

    assert d.get_x(True) == pytest.approx(x)
    assert d.get_xerr(True) == pytest.approx(xerr)

    # Ignore a few points at the start and end
    d.notice(11, 18)

    # Just check that the default behavior doesn't change with the filter
    assert d.get_x() == pytest.approx(x)
    assert d.get_xerr() == pytest.approx(xerr)

    assert d.get_x(True) == pytest.approx(x[2:-1])
    assert d.get_xerr(True) == pytest.approx(xerr[2:-1])

    # Now ignore all points
    d.ignore(0, 1000)

    assert d.get_x() == pytest.approx(x)
    assert d.get_xerr() == pytest.approx(xerr)

    assert d.get_x(True) == pytest.approx([])
    assert d.get_xerr(True) == pytest.approx([])


@pytest.mark.parametrize('ignore', [False, True])
@pytest.mark.parametrize('lo,hi,evals',
                         [(0.5, 2.3, (0, 10, 0)),
                          (0.7, 2.1, (1, 8, 1)),
                          (0.5, 0.7, (0, 2, 8)),
                          (1.1, 1.3, (3, 2, 5)),
                          (2.1, 2.3, (8, 2, 0)),
                          # special case filters that are within a single bin
                          (0.45, 0.55, (0, 1, 9)),
                          (0.65, 0.75, (1, 1, 8)),
                          (1.05, 1.15, (3, 1, 6)),
                          (2.25, 2.35, (9, 1, 0)),
                          # outside the limits
                          (0.1, 0.4, (0, 0, 10)),
                          (0.1, 0.5, (0, 1, 9)),
                          (2.41, 2.8, (10, 0, 0)),
                          # Now queries on the edge of each bin; these would ideally
                          # only match 1 bin
                          (0.4, 0.6, (0, 1, 9)),
                          (0.6, 0.8, (1, 1, 8)),
                          (0.8, 1.0, (2, 1, 7)),
                          (1.0, 1.2, (3, 1, 6)),
                          (1.2, 1.4, (4, 1, 5)),
                          (1.4, 1.6, (5, 1, 4)),
                          (1.6, 1.8, (6, 1, 3)),
                          (1.8, 2.0, (7, 1, 2)),
                          (2.0, 2.2, (8, 1, 1)),
                          (2.2, 2.4, (9, 1, 0)),
                          # check last upper limit
                          (2.4, 2.6, (10, 0, 0))
                         ])
def test_data1dint_check_limit(ignore, lo, hi, evals):
    """Does Data1DInt handle limits (in particular upper limits).

    This is based on sherpa/astro/tests/test_astro_data.py::test_pha_check_limit
    but without the need for an ARF. It selects different bins than the
    PHA case!
    """

    egrids = 0.2 + 0.2 * numpy.arange(1, 12)
    d = Data1DInt('exammple', egrids[:-1], egrids[1:],
                  numpy.ones(10))

    assert d.mask is True

    func = d.ignore if ignore else d.notice
    func(lo, hi)
    if ignore:
        vout = True
        vin = False
    else:
        vout = False
        vin = True

    c1, c2, c3 = evals
    expected = [vout] * c1 + [vin] * c2 + [vout] * c3
    assert d.mask == pytest.approx(expected)


@pytest.fixture
def make_data2dint():
    """Create a simple 2D Int data set."""

    # A 1 by 2 grid, which we make sure does not start at 1,1 to check
    # that this is handled correctly.
    #
    x1, x0 = numpy.mgrid[10:12, -5:-4]
    shape = x0.shape
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = numpy.asarray([10, 5])

    return Data2DInt("ival", x0, x1, x0 + 1, x1 + 1,
                     y, shape=shape)


def test_data2dint_create(make_data2dint):
    """Check we can create a basic integrated 2D data set.

    See issue #1379
    """

    x0 = numpy.asarray([-5, -5])
    x1 = numpy.asarray([10, 11])

    img = make_data2dint

    assert (img.dep == [10, 5]).all()

    assert len(img.indep) == 4
    assert (img.indep[0] == x0).all()
    assert (img.indep[1] == x1).all()
    assert (img.indep[2] == (x0 + 1)).all()
    assert (img.indep[3] == (x1 + 1)).all()

    # I was initially surprised there was no header, so just make sure
    # we check for it not being here.
    #
    assert not(hasattr(img, "header"))


def test_data2dint_show(make_data2dint):
    """Check we can show a basic integrated 2D data set.

    See issue #1379
    """

    img = make_data2dint

    # This fails because there's problems getting x0 and x0lo
    # attributes.
    #
    out = str(img).split("\n")

    assert out[0] == "name      = ival"
    assert out[1] == "x0lo      = Int64[2]"
    assert out[2] == "x1lo      = Int64[2]"
    assert out[3] == "x0hi      = Int64[2]"
    assert out[4] == "x1hi      = Int64[2]"
    assert out[5] == "y         = Int64[2]"
    assert out[6] == "staterror = None"
    assert out[7] == "syserror  = None"
    assert out[8] == "shape     = (2, 1)"
    assert len(out) == 9


def test_data2dint_get_x0(make_data2dint):
    x0 = numpy.asarray([-5, -5])
    x = (x0 + x0 + 1) / 2

    assert (make_data2dint.get_x0() == x).all()


def test_data2dint_x0(make_data2dint):
    x0 = numpy.asarray([-5, -5])
    x = (x0 + x0 + 1) / 2

    assert (make_data2dint.x0 == x).all()


def test_data2dint_get_x1(make_data2dint):
    x1 = numpy.asarray([10, 11])
    x = (x1 + x1 + 1) / 2

    assert (make_data2dint.get_x1() == x).all()


def test_data2dint_x1(make_data2dint):
    x1 = numpy.asarray([10, 11])
    x = (x1 + x1 + 1) / 2

    assert (make_data2dint.x1 == x).all()


def test_data2dint_x0lo(make_data2dint):
    assert (make_data2dint.x0lo == [-5, -5]).all()


def test_data2dint_x0hi(make_data2dint):
    assert (make_data2dint.x0hi == [-4, -4]).all()


def test_data2dint_x1lo(make_data2dint):
    assert (make_data2dint.x1lo == [10, 11]).all()


def test_data2dint_x1hi(make_data2dint):
    assert (make_data2dint.x1hi == [11, 12]).all()


def test_data2dint_get_y(make_data2dint):
    assert (make_data2dint.get_y() == [10, 5]).all()


def test_data2dint_y(make_data2dint):
    assert (make_data2dint.y == [10, 5]).all()


def test_data2dint_get_dep(make_data2dint):
    assert (make_data2dint.get_dep() == [10, 5]).all()


def test_data2dint_get_x0label(make_data2dint):
    assert make_data2dint.get_x0label() == "x0"


def test_data2dint_get_x1label(make_data2dint):
    assert make_data2dint.get_x1label() == "x1"


def test_data2dint_get_ylabel(make_data2dint):
    assert make_data2dint.get_ylabel() == "y"


def test_data2dint_get_axes(make_data2dint):
    axes = make_data2dint.get_axes()
    assert len(axes) == 2
    assert (axes[0] == [1]).all()
    assert (axes[1] == [1, 2]).all()


def test_data2dint_notice(make_data2dint):
    """basic notice call"""
    img = make_data2dint

    # The mask attribute can be True, False, or a ndarray. Fortunately
    # using an ndarray as a truthy value throws a ValueError.
    #
    assert img.mask

    # Data is defined on x0=-5, x1=10,11
    # so this excludes the second point.
    #
    img.notice(x1lo=10, x1hi=11)
    assert (img.mask == numpy.asarray([True, False])).all()


def test_data2dint_ignore(make_data2dint):
    """basic ignore call"""
    img = make_data2dint

    assert img.mask
    img.notice(x1lo=10, x1hi=11, ignore=True)
    assert (img.mask == numpy.asarray([False, True])).all()


def test_data2dint_ignore_get_filter(make_data2dint):
    """What exactly does get_filter return here?

    The current behavior does not look sensible.
    """
    img = make_data2dint

    assert img.mask
    img.notice(x1lo=10, x1hi=11, ignore=True)
    assert img.get_filter() == ''


def test_data2dint_ignore_get_filter_expr(make_data2dint):
    """What exactly does get_filter_expr return here?

    The current behavior does not look sensible.
    """
    img = make_data2dint

    assert img.mask
    img.notice(x1lo=10, x1hi=11, ignore=True)
    assert img.get_filter_expr() == ''


def test_data2dint_notice_get_x0(make_data2dint):
    """basic notice call + get_x0"""
    img = make_data2dint
    img.notice(x1lo=10, x1hi=11)
    assert (img.get_x0() == numpy.asarray([-4.5, -4.5])).all()
    assert (img.get_x0(True) == numpy.asarray([-4.5])).all()


def test_data2dint_notice_get_x1(make_data2dint):
    """basic notice call + get_x1"""
    img = make_data2dint
    img.notice(x1lo=10, x1hi=11)
    assert (img.get_x1() == numpy.asarray([10.5, 11.5])).all()
    assert (img.get_x1(True) == numpy.asarray([10.5])).all()


def test_data2dint_notice_get_y(make_data2dint):
    """basic notice call + get_y"""
    img = make_data2dint
    img.notice(x1lo=10, x1hi=11)
    assert (img.get_y() == numpy.asarray([10, 5])).all()
    assert (img.get_y(True) == numpy.asarray([10])).all()


def test_data2dint_get_dims(make_data2dint):
    assert make_data2dint.get_dims() == (1, 2)


def test_data2dint_get_img(make_data2dint):
    ival = make_data2dint.get_img()
    assert ival.shape == (2, 1)
    assert (ival == numpy.asarray([[10], [5]])).all()


def test_data2dint_get_img_model(make_data2dint):
    """Check we can evaluate a model AND we ignore a filter"""
    img = make_data2dint

    # This model evaluates
    #   mdl.c + mdl.cx1 * x0 + mdl.cy1 * x1
    #
    # which becomes, because we use the middle of the bin
    #
    #   10 + 1 * (-4.5) + 10 * (10.5, 11.5)
    #   = (110.5, 120.5)
    #
    mdl = Polynom2D()
    mdl.c = 10
    mdl.cy1 = 10
    mdl.cx1 = 1

    # This selects only one point.
    #
    img.notice(x1lo=10, x1hi=11)

    # This should ignore the filter.
    #
    ivals = img.get_img(mdl)
    assert len(ivals) == 2
    assert ivals[0].shape == (2, 1)
    assert ivals[1].shape == (2, 1)
    assert (ivals[0] == numpy.asarray([[10], [5]])).all()
    assert (ivals[1] == numpy.asarray([[110.5], [120.5]])).all()


def test_data2dint_get_max_pos(make_data2dint):
    assert make_data2dint.get_max_pos() == (-4.5, 10.5)


def test_data2dint_get_bounding_mask(make_data2dint):
    """Data2D/Data2DInt do not have get_bounding_mask"""
    assert not hasattr(make_data2dint, "get_bounding_mask")


@pytest.mark.parametrize("method",
                         ["get_error",
                          "get_imgerr",
                          "get_staterror",
                          "get_syserror",
                          "get_yerr"
                         ])
def test_data2dint_method_is_none(method, make_data2dint):
    """Check those methods that return None"""
    func = getattr(make_data2dint, method)
    assert func() is None


@pytest.mark.parametrize("attribute",
                         ["staterror",
                          "syserror"
                         ])
def test_data2dint_attribute_is_none(attribute, make_data2dint):
    """Check those attributes that return None"""
    attr = getattr(make_data2dint, attribute)
    assert attr is None


@pytest.mark.parametrize("data", (Data, ) + DATA_1D_CLASSES, indirect=True)
def test_reduce_axis_size_1d(data):
    """What happens if we reduce the independent axis?"""

    nindep = len(data.indep)
    for indep in data.indep:
        assert len(indep) == 10

    for attr in ["dep", "staterror", "syserror"]:
        aval = getattr(data, attr)
        assert len(aval) == 10

    # Sanity checks.
    #
    for a, b in zip(data.indep, data.get_indep()):
        assert numpy.all(a == b)

    # Let's make the independent axis smaller.
    #
    smaller = []
    for indep in data.indep:
        smaller.append(indep[1:-1])

    data.indep = tuple(smaller)

    # Check what has changed.
    #
    assert len(data.indep) == nindep
    for indep in data.indep:
        assert len(indep) == 8

    # Note that the other values have not changed.
    #
    for attr in ["dep", "staterror", "syserror"]:
        aval = getattr(data, attr)
        assert len(aval) == 10


@pytest.mark.parametrize("data", DATA_2D_CLASSES, indirect=True)
def test_reduce_axis_size_2d(data):
    """What happens if we reduce the independent axis?

    There is a shape attribute which could be changed, or
    maybe should be changed.
    """

    nindep = len(data.indep)
    for indep in data.indep:
        assert len(indep) == 100

    for attr in ["dep", "staterror", "syserror"]:
        aval = getattr(data, attr)
        assert len(aval) == 100

    assert data.shape == (10, 10)

    # Sanity checks.
    #
    for a, b in zip(data.indep, data.get_indep()):
        assert numpy.all(a == b)

    # Let's make the independent axis smaller.
    #
    smaller = []
    for indep in data.indep:
        smaller.append(indep[1:-1])

    data.indep = tuple(smaller)

    # Check what has changed.
    #
    assert len(data.indep) == nindep
    for indep in data.indep:
        assert len(indep) == 98

    # Note that the other values have not changed.
    #
    assert data.shape == (10, 10)

    for attr in ["dep", "staterror", "syserror"]:
        aval = getattr(data, attr)
        assert len(aval) == 100


@pytest.mark.parametrize("data", ALL_DATA_CLASSES, indirect=True)
def test_invalid_independent_axis(data):
    """What happens if we use the wrong number of independent axes?

    We just duplicate the current axes.
    """

    indep = data.indep
    with pytest.raises(TypeError):
        data.indep = tuple(list(indep) * 2)


@pytest.mark.parametrize("data", (Data1DInt, Data2D, Data2DInt), indirect=True)
def test_invalid_independent_axis_component_size(data):
    """What happens if we use mis-matched sizes?

    It only makes sense to do this for data classes with
    multiple components. We remove one entry from the
    second component,
    """

    indep = list(data.indep)
    indep[1] = indep[1][:-1]
    # At the moment this does not error out
    data.indep = tuple(indep)


@pytest.mark.parametrize("data", (Data1DInt, Data2D, Data2DInt), indirect=True)
def test_invalid_independent_axis_component_none(data):
    """What happens if we use mis-matched sizes (by setting one to None).

    See test_invalid_independent_axis_component_size.
    """

    indep = list(data.indep)
    indep[1] = None
    # At the moment this does not error out
    data.indep = tuple(indep)


@pytest.mark.parametrize("data", ALL_DATA_CLASSES, indirect=True)
def test_invalid_dependent_axis(data):
    """What happens if the dependent axis does not match the independent axis?
    """

    # This currently does not fail
    data.y = data.y[:-2]


@pytest.mark.parametrize("data_class", ALL_DATA_CLASSES)
def test_make_invalid_dependent_axis(data_class):
    """What happens if call constructor with invalid independent axis?
    """

    # Take the correct arguments and reduce the independent axis by one.
    # Use a copy of everything just in case.
    args = []
    for arg in INSTANCE_ARGS[data_class]:
        if isinstance(arg, numpy.ndarray):
            arg = arg.copy()

        args.append(arg)

    ypos = POS_Y_ARRAY[data_class]
    args[ypos] = args[ypos][:-1]

    # This does not raise an error
    data_class(*args)


@pytest.mark.parametrize("data", ALL_DATA_CLASSES, indirect=True)
def test_set_independent_axis_to_none(data):
    """What happens if we clear the independent axis?

    It's not entirely clear what we mean by setting the independent
    axis to None: do we need to set all components or is it sufficient
    to use a single None. Similarly, what should the response be. We
    currently only test the current behavior.

    """

    assert all(d is not None for d in data.indep)

    indep = [None for d in data.indep]
    data.set_indep(tuple(indep))
    assert all(d is None for d in data.indep)
