#
# Copyright (C) 2015  Smithsonian Astrophysical Observatory
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
from sys import version_info
from sherpa.astro.datastack import DataStack

if version_info.major == 2:
    import __builtin__ as builtins
    import mock
if version_info.major == 3:
    import builtins
    from unittest import mock


import pytest

from sherpa.astro.datastack.utils import model_wrapper, create_stack_model,\
    simple_wrapper, fit_wrapper, plot_wrapper

# The following imports are included to simplify the mocking/stubbing.
# Implementation should not be used if not to infer the mocking specs
from sherpa.astro.instrument import RSPModelPHA, RMFModelPHA, ARFModelPHA
from sherpa.models.model import ArithmeticConstantModel
from sherpa.models import UserModel, Polynom1D

mockobj = mock.MagicMock()
base = {
    'pha': mockobj,
    'model': mockobj
}

model_classes = {
    RSPModelPHA: mock.MagicMock(spec=RSPModelPHA, rmf=mockobj, arf=mockobj, **base),
    RMFModelPHA: mock.MagicMock(spec=RMFModelPHA, rmf=mock.MagicMock(), **base),
    ARFModelPHA: mock.MagicMock(spec=ARFModelPHA, arf=mock.MagicMock(), **base),
    ArithmeticConstantModel: mock.MagicMock(spec=ArithmeticConstantModel,
                                            val=3.0)
}


@pytest.fixture(params=model_classes.keys())
def model_fixture(request):
    model_class = request.param
    model_class_name = "sherpa.astro.instrument." + model_class.__name__
    model_comps_len = 1
    part = mock.MagicMock(spec=model_class)
    part.name = "model.part"
    part.pars = []
    retval = model_classes[request.param]
    if request.param is not ArithmeticConstantModel:
        retval.parts = (part,)
    else:
        model_class = float
        model_class_name = "sherpa.models.model.ArithmeticConstantModel"
        model_comps_len = 0
    retval.name = "model.name"

    return retval, model_class, model_class_name, model_comps_len


def test_model_factory(model_fixture):
    model = model_fixture[0]
    model_class = model_fixture[1]
    model_class_name = model_fixture[2]
    model_comps_len = model_fixture[3]
    with mock.patch(model_class_name, autospec=True):
        new_model, model_comps = create_stack_model(model, "model.name")
    assert isinstance(new_model, model_class)
    assert len(model_comps) == model_comps_len


def test_model_factory_regular_model():
    model = mock.MagicMock(spec_set=["name"])
    model.name = "model.name"

    model, model_comps = create_stack_model(model, "1")

    assert len(model_comps) == 1
    assert "name" in model_comps
    assert model_comps["name"]["model_name"] == "name"


def test_model_factory_template():
    model = mock.MagicMock(spec=Polynom1D)
    model.name = "polynom1d.p__ID"

    model, model_comps = create_stack_model(model, "1")

    assert len(model_comps) == 1
    assert "p" in model_comps
    assert model_comps["p"]["model_name"] == "p1"


def test_model_factory_usermodel():
    model = mock.MagicMock(spec=UserModel)
    model.name = "model.usermodel"

    model, model_comps = create_stack_model(model, "1")

    assert len(model_comps) == 1
    assert "usermodel" in model_comps

    assert model_comps["usermodel"]["model_name"] == "usermodel"


@mock.patch("sherpa.ui.create_model_component")
def test_model_factory_template_usermodel(_):
    model = mock.MagicMock(spec=UserModel)
    model.name = "somemodel.mymodel__ID"

    model, model_comps = create_stack_model(model, "1")

    assert len(model_comps) == 1
    assert "mymodel" in model_comps
    assert model_comps["mymodel"]["model_name"] == "mymodel1"


def test_model_factory_error_unexpected():
    model = mock.MagicMock(spec_set=["parts"])

    with pytest.raises(ValueError):
        create_stack_model(model, "model.name")


def test_model_factory_error_model_name():
    model = mock.MagicMock(spec_set=["name"])

    with pytest.raises(ValueError):
        create_stack_model(model, "model")


@mock.patch.object(builtins, 'eval')
def test_model_wrapper(eval, model_fixture):
    model = model_fixture[0]
    model_class_name = model_fixture[2]

    eval.return_value = model

    func = mock.MagicMock()
    func.__name__ = "sample function"
    func.__doc__ = "sample doc"

    dataset = mock.MagicMock()
    dataset.__getitem__.side_effect = lambda x: "TESTID" if x == "id" else dict()

    datastack = mock.MagicMock(spec=DataStack)
    datastack.filter_datasets.return_value = (dataset,)

    wrapper = model_wrapper(func)
    assert wrapper.__name__ == func.__name__
    assert wrapper.__doc__ == func.__doc__

    with mock.patch(model_class_name, autospec=True):
        wrapper(datastack, "model")

    expected = model.val if hasattr(model, 'val') else model
    func.assert_called_once_with("TESTID", mock.ANY)
    assert isinstance(func.call_args[0][1], expected.__class__)


def test_simple_wrapper():
    func = mock.MagicMock()
    func.side_effect = ["OK1", "OK2"]
    func.__name__ = "sample function"
    func.__doc__ = "sample doc"

    datasets = []

    dataset = mock.MagicMock()
    dataset.__getitem__.side_effect = lambda x: "TESTID1" if x == "id" else dict()
    datasets.append(dataset)

    dataset = mock.MagicMock()
    dataset.__getitem__.side_effect = lambda x: "TESTID2" if x == "id" else dict()
    datasets.append(dataset)

    datastack = mock.MagicMock(spec=DataStack)
    datastack.filter_datasets.return_value = datasets

    wrapper = simple_wrapper(func)
    assert wrapper.__name__ == func.__name__
    assert wrapper.__doc__ == func.__doc__
    args = ["a", 1, []]
    kwargs = {"one": 1, "two": (2,)}
    ans = wrapper(datastack, *args, **kwargs)
    assert ans == ["OK1", "OK2"]
    calls = [mock.call("TESTID1", *args, **kwargs),
             mock.call("TESTID2", *args, **kwargs)]
    func.assert_has_calls(calls)


def test_fit_wrapper():
    func = mock.MagicMock()
    func.side_effect = ["OK1", "OK2"]
    func.__name__ = "sample function"
    func.__doc__ = "sample doc"

    datasets = []

    dataset = mock.MagicMock()
    dataset.__getitem__.side_effect = lambda x: "TESTID1" if x == "id" else dict()
    datasets.append(dataset)

    dataset = mock.MagicMock()
    dataset.__getitem__.side_effect = lambda x: "TESTID2" if x == "id" else dict()
    datasets.append(dataset)

    datastack = mock.MagicMock(spec=DataStack)
    datastack.filter_datasets.return_value = datasets

    wrapper = fit_wrapper(func)
    assert wrapper.__name__ == func.__name__
    assert wrapper.__doc__ != func.__doc__
    args = ["a", 1, []]
    kwargs = {"one": 1, "two": (2,)}
    ans = wrapper(datastack, *args, **kwargs)
    assert ans is None
    func.assert_called_once_with("TESTID1", "TESTID2", *args, **kwargs)


def test_plot_wrapper():
    func = mock.MagicMock()
    func.side_effect = ["OK1", "OK2"]
    func.__name__ = "sample function"
    func.__doc__ = "sample doc"

    datasets = []

    dataset1 = mock.MagicMock()
    dataset1.__getitem__.side_effect = lambda x: "TESTID1" if x == "id" else dict()
    datasets.append(dataset1)

    dataset2 = mock.MagicMock()
    dataset2.__getitem__.side_effect = lambda x: "TESTID2" if x == "id" else dict()
    datasets.append(dataset2)

    datastack = mock.MagicMock(spec=DataStack)
    datastack.filter_datasets.return_value = datasets
    datastack.ids = ["TESTID1", "TESTID2"]

    with mock.patch("sherpa.astro.datastack.utils.backend") as backend:
        wrapper = plot_wrapper(func)
        assert wrapper.__name__ == func.__name__
        assert wrapper.__doc__ != func.__doc__
        args = ["a", 1, []]
        kwargs = {"one": 1, "two": (2,)}
        ans = wrapper(datastack, *args, **kwargs)
        backend.initialize_backend.assert_called_once_with()
        calls = [mock.call(dataset1, ["TESTID1", "TESTID2"]), mock.call(
            dataset2, ["TESTID1", "TESTID2"])]
        backend.initialize_plot.assert_has_calls(calls)
        assert ans is None
        calls = [mock.call("TESTID1", *args, **kwargs),
                 mock.call("TESTID2", *args, **kwargs)]
        func.assert_has_calls(calls)
