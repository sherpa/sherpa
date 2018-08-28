# Copyright 2018 Smithsonian Astrophysical Observatory
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
# following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
# disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
# products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
import numpy as np
import pytest
from pytest import approx

from sherpa.astro.data import DataIMG, DataIMGInt
from sherpa.astro.ui.utils import Session
from sherpa.data import Data1DInt
from sherpa.models import Const1D, ArithmeticModel, Parameter, Const2D, ArithmeticModel2D


@pytest.fixture
def setup():
    const = Const1D("const")
    const.c0 = 0
    const.c0.freeze()

    my_model = MyModel("my_model")
    my_model.integrate = False

    return Session(), my_model, const


@pytest.fixture
def setup2d():
    const = Const2D("const")
    const.c0 = 0
    const.c0.freeze()

    x = [2, 3, 2, 3]
    y = [2, 2, 3, 3]
    xhi = [2.1, 3.5, 2.1, 3.5]
    yhi = [2.1, 2.1, 3, 3.5]

    # This is the result when rebinning [100, ] * 4
    z = [225, ] * 4

    my_model = MyModel2D("my_model")

    return Session(), my_model, const, (x, y, xhi, yhi, z)


def test_evaluate_model_on_arbitrary_grid_point_list(setup):
    """
    The idea of the test is that the model will evaluate differently depending on the grid it is evaluated on.
    This is really arbitrary, it just exercises the high level API for a common workflow while making sure the results
    are the expected ones.
    """
    ui, my_model, const = setup

    # Load data
    ui.load_arrays(1, [1, 2, 3], [100, 100, 100])

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid([1, 2, 2.5, 4, 5])

    # The model will usually be part of a complex model expression, so let's pretend we add another component,
    # although that component is muted.
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, 1)

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = [1, 2, 3, 4, 5]

    assert_fit(ui, my_model, 0)


def test_evaluate_model_on_arbitrary_grid_point_list_2d(setup2d):
    """
    The idea of the test is that the model will evaluate differently depending on the grid it is evaluated on.
    This is really arbitrary, it just exercises the high level API for a common workflow while making sure the results
    are the expected ones.
    """
    ui, my_model, const, data = setup2d
    x, y, _, _, z = data

    # Load data
    ui.load_arrays(1, x, y, z, DataIMG)

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid([2, 2.5, 3], [2, 2.5, 3])

    # The model will usually be part of a complex model expression, so let's pretend we add another component,
    # although that component is muted.
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, (1, 1))

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = [2, 3], [2, 3]

    assert_fit(ui, my_model, (0, 0))


def test_evaluate_model_on_arbitrary_grid_integrated_list(setup):
    """
    Same as above, but with integrated models.
    """
    ui, my_model, const = setup

    # Load data
    ui.load_arrays(1, [1.5, 2.5, 3.5], [2.5, 3.5, 4.5], [100, 100, 100], Data1DInt)

    # Get a model that evaluates on a different grid
    # This is the important part.
    # Note that our model has the integrate flag disabled, so the midpoint of these values
    # is taken into account when calculating the model.
    regrid_model = my_model.regrid([0, 1, 2], [1, 2, 3])

    # The model will be part of a complex model expression, so let's pretend we add another component
    ui.set_source(regrid_model + const)

    # Fit and check the result
    with pytest.warns(UserWarning):
        assert_fit(ui, my_model, 1)

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = [1.5, 2.5, 3.5], [2.5, 3.5, 4.5]

    assert_fit(ui, my_model, 0)


def test_evaluate_model_on_arbitrary_grid_integrated_list_2d(setup2d):
    """
    Same as above, but with integrated models
    """
    ui, my_model, const, data = setup2d
    x, y, xhi, yhi, z = data

    # Load data
    ui.load_arrays(1, x, y, xhi, yhi, z, DataIMGInt)

    regrid_lo = [2, 2.5, 3]
    regrid_hi = np.array([2, 2.5, 3.5])

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid(regrid_lo, regrid_lo, regrid_hi, regrid_hi)

    # The model will usually be part of a complex model expression, so let's pretend we add another component,
    # although that component is muted.
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, (1, 1))

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = x, y, xhi, yhi

    assert_fit(ui, my_model, (0, 0))


def test_evaluate_model_on_arbitrary_grid_point_ndarray(setup):
    """
    The idea of the test is that the model will evaluate differently depending on the grid it is evaluated on.
    This is really arbitrary, it just exercises the high level API for a common workflow while making sure the results
    are the expected ones.
    """
    ui, my_model, const = setup

    # Load data
    ui.load_arrays(1, [1, 2, 3], [100, 100, 100])

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid(np.array([1, 2, 2.5, 4, 5]))

    # The model will be part of a complex model expression, so let's pretend we add another component
    ui.set_source(regrid_model + const)

    # Fit and check the result
    assert_fit(ui, my_model, 1)

    # Now fit with a different regrid.
    # This is also the important part.
    regrid_model.grid = np.array([1, 2, 3, 4, 5])

    assert_fit(ui, my_model, 0)


def test_evaluate_model_on_arbitrary_grid_integrated_ndarray(setup):
    """
    Same as above, but with integrated models.
    """
    ui, my_model, const = setup

    # Load data
    ui.load_arrays(1, [1.5, 2.5, 3.5], [2.5, 3.5, 4.5], [100, 100, 100], Data1DInt)

    # Get a model that evaluates on a different grid
    # This is the important part.
    regrid_model = my_model.regrid(np.array([0, 1, 2]), [1, 2, 3])

    # The model will be part of a complex model expression, so let's pretend we add another component
    ui.set_source(regrid_model + const)

    # Fit and check the result
    with pytest.warns(UserWarning):
        assert_fit(ui, my_model, 1)

    # Now fit with a different grid.
    # This is also the important part.
    regrid_model.grid = [1.5, 2.5, 3.5], np.array([2.5, 3.5, 4.5])

    assert_fit(ui, my_model, 0)


def test_evaluate_model_on_arbitrary_grid_no_overlap(setup):
    """
    If grids do not overlap, issue a warning and return zeros
    """
    ui, my_model, _ = setup

    # Get a model that evaluates on a different grid
    # This is the important part. Note that there is overlap, but
    # the start and end p
    regrid_model = my_model.regrid([2, 2.5], [2, 2.5])

    with pytest.warns(UserWarning):
        np.testing.assert_array_equal(regrid_model([1, 2], [1, 2]), [0, 0])


def test_evaluate_model_on_arbitrary_grid_no_overlap_2d(setup2d):
    """
    In the 2D case, the overlap is way more stringent than in the 1D case, due to the complexity of rebinning
    """
    ui, my_model, _, data = setup2d
    x, y, _, _, _ = data

    my_model.x_has_25 = 1  # To force the model to evaluate to something other than 0.

    # Get a model that evaluates on a different grid
    # This is the important part. Note that there is overlap, but
    # the start and end points are different.
    regrid_model = my_model.regrid([2, 2.5], [2, 2.5])

    with pytest.warns(UserWarning):
        np.testing.assert_array_equal(regrid_model(x, y), [0, 0, 0, 0])



class MyModel(ArithmeticModel):
    """
    A model that returns [100, ] * len(x) if 2.5 is in the input array x
    """
    def __init__(self, name):
        self.has_25 = Parameter(name, "has_25", 0, min=0, max=1)
        ArithmeticModel.__init__(self, name, (self.has_25,))

    def guess(self, dep, *args, **kwargs):
        raise NotImplementedError()

    def get_center(self):
        raise NotImplementedError()

    def set_center(self, *args, **kwargs):
        raise NotImplementedError()

    def calc(self, p, *args, **kwargs):
        x = args[0]

        if 2.5 not in x:
            if p[0] == 0:
                return [100, ] * len(x)
            else:
                return [100-p[0]*100,] * len(x)
        if p[0] == 1:
            return [100, ] * len(x)
        else:
            return [p[0]*100, ] * len(x)


class MyModel2D(ArithmeticModel2D):
    """
    A 2D model that returns [100, ] * len(x) * len(y) if 2.5 is in the input arrays x and y
    """
    def __init__(self, name):
        self.x_has_25 = Parameter(name, "x_has_25", 0, min=0, max=1)
        self.y_has_25 = Parameter(name, "y_has_25", 0, min=0, max=1)
        ArithmeticModel2D.__init__(self, name, (self.x_has_25, self.y_has_25))

    def guess(self, dep, *args, **kwargs):
        raise NotImplementedError()

    def get_center(self):
        raise NotImplementedError()

    def set_center(self, *args, **kwargs):
        raise NotImplementedError()

    def calc(self, p, *args, **kwargs):
        x, y, x_has_25, y_has_25 = args[0], args[1], p[0], p[1]

        x_eval = np.array(self._eval(x, x_has_25))
        y_eval = np.array(self._eval(y, y_has_25))

        return (x_eval + y_eval)/2

    def _eval(self, array, has_25):
        if 2.5 not in array:
            if has_25 == 0:
                return [100, ] * len(array)
            else:
                return [100-has_25*100,] * len(array)
        if has_25 == 1:
            return [100, ] * len(array)
        else:
            return [has_25*100, ] * len(array)


def assert_fit(ui, model, value):
    ui.fit()

    try:  # 2D, two values
        len(value)
        assert model.x_has_25.val == approx(value[0])
        assert model.y_has_25.val == approx(value[1])
    except TypeError:  # 1D, one value
        assert model.has_25.val == approx(value)
