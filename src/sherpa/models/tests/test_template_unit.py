#
#  Copyright (C) 2016, 2018, 2022, 2023
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

import os

import numpy as np

import pytest

from sherpa.astro.ui.utils import Session as AstroSession
from sherpa.data import Data1D
from sherpa.models.basic import TableModel, Gauss1D
from sherpa.models.parameter import Parameter
from sherpa.models.template import Template, TemplateModel, \
    create_template_model
from sherpa.ui.utils import Session
from sherpa.utils.err import ModelErr
from sherpa.utils.testing import has_package_from_list, requires_data


@pytest.fixture
def skip_if_no_io(request):
    """For tests with a session fixture that need to be skipped if
    the I/O backend is not present ONLY for AstroSession.
    """

    if request.getfixturevalue("session") == Session:
        return

    # This requires knowledge of requires_fits, but we don't make it
    # possible to access this so we hard code the knowledge here. If
    # we ever get another I/O backend we need to revisit this.
    #
    if has_package_from_list("astropy.io.fits", "pycrates"):
        return

    pytest.skip(reason="FITS backend required")


@requires_data
@pytest.mark.parametrize("session", [Session, AstroSession])
def test_309(session, make_data_path, skip_if_no_io):

    idval = 'bug309'

    # have values near unity for the data
    ynorm = 1e9

    s = session()

    dname = make_data_path('load_template_with_interpolation-bb_data.dat')

    s.load_data(idval, dname)
    s.get_data(idval).y *= ynorm

    indexname = 'bb_index.dat'
    datadir = make_data_path('')

    # Need to load the data from the same directory as the index
    basedir = os.getcwd()
    os.chdir(datadir)
    try:
        s.load_template_model('bbtemp', indexname)
    finally:
        os.chdir(basedir)

    bbtemp = s.get_model_component('bbtemp')
    s.set_source(idval, bbtemp * ynorm)

    s.set_method('gridsearch')
    s.set_method_opt('sequence', None)
    s.fit(idval)


def test_create_template_model_parnames_does_not_match():
    """What happens if we send in too many parnames?"""

    templates = [TableModel("foo")]
    parvals = np.asarray([[2.3]])
    with pytest.raises(ValueError,
                       match="number of parvals and names do not match: 1 vs 2"):
        create_template_model("foo", ["a", "b"], parvals, templates)


def test_create_template_model_parvals_not_2d():
    """What happens if we send parvals not 2D"""

    # These arguments do not really make sense. They are enough to
    # check the code.
    #
    templates = [TableModel("foo")]
    parvals = np.asarray([2.3])
    with pytest.raises(ValueError,
                       match="parvals must be 2D, sent 1D"):
        create_template_model("foo", ["a"], parvals, templates)


def test_create_template_model_parvals_wrong_size():
    """What happens if we send parvals not the right size"""

    # These arguments do not really make sense. They are enough to
    # check the code.
    #
    templates = [TableModel("foo")]
    parvals = np.arange(12).reshape(3, 4)
    with pytest.raises(ValueError,
                       match="number of parvals and names do not match: 4 vs 1"):
        create_template_model("foo", ["a"], parvals, templates)


def test_create_template_model_unknown_interpolator():
    """What happens if we send an unknown interpolator?"""

    # These arguments do not really make sense. They are enough to
    # check the code.
    #
    templates = [TableModel("foo"), TableModel("bar")]
    parvals = np.asarray([[12], [13]])
    with pytest.raises(ModelErr,
                       match="Unknown template_interpolator_name=made_up"):
        create_template_model("foo", ["a"], parvals, templates,
                              template_interpolator_name="made_up")


def test_templatemodel_wrong_number_of_pars():
    """Check we error out"""

    templates = [TableModel("foo"), TableModel("bar")]
    with pytest.raises(ModelErr,
                       match="Number of parameter values and templates do not match"):
        TemplateModel(templates=templates)


def test_templatemodel_wrong_parvals_element():
    """Check we error out"""

    templates = [TableModel("foo"), TableModel("bar")]
    pars = [Parameter("fooy", "x", 1), Parameter("fooy", "y", 2)]
    parvals = [[2, 3], [1], [4, 5, 54]]
    with pytest.raises(ModelErr,
                       match="Number of parameter values and templates do not match"):
        TemplateModel(pars=pars, parvals=parvals,
                      templates=templates)


def test_templatemodel_no_arguments():
    """What happens if we have no parameters?"""

    with pytest.raises(ModelErr,
                       match="^parvals is empty or not set"):
        empty = TemplateModel("empty")


def test_templatemodel_pars_no_templates():
    """What happens if we have pars but no templates"""

    pars = [Parameter("empty", "bob", 12)]
    with pytest.raises(ModelErr,
                       match="^parvals is empty or not set"):
        TemplateModel("empty", pars=pars)


def test_templatemodel_templates_not_tablemodel():
    """What happens if we send in templates that are not TableModels?"""

    templates = [Gauss1D("m1"), Gauss1D("m2")]
    woop = create_template_model("woop", ["woop"], np.asarray([[5], [10]]),
                                 templates)

    with pytest.raises(AttributeError,
                       match="'Gauss1D' object has no attribute 'fold'"):
        woop.fold(Data1D("x", [1, 2, 3], [1, 2, 3]))

    # Not sure why this actually fails. Ideally the error would be clearer,
    # but really this should have been caught at object creation.
    #
    with pytest.raises(TypeError,
                       match="1D model evaluation input array sizes do not match, xlo: 3 vs xhi: 1"):
        woop([1, 2, 3])


def setup_basic(iname=None, no_x=False):
    """Can we get a TemplateModel set up?"""

    x = np.arange(1, 5)
    xdata = None if no_x else x

    templates = []
    for scale, name in enumerate(["m1", "m2", "m3"], 1):
        template = TableModel(name)
        template.load(xdata, x * scale)
        templates.append(template)

    return create_template_model("bob", ["pa"], np.asarray([[10], [20], [30]]),
                                 templates, template_interpolator_name=iname)


@pytest.mark.parametrize("iname,cls", [(None, TemplateModel), ("default", Template)])
def test_templatemodel_basic_returns_class(iname, cls):
    """Check we can get the correct class"""

    assert isinstance(setup_basic(iname), cls)


@pytest.mark.parametrize("iname", [None, "default"])
def test_templatemodel_basic_show(iname):
    """Check we can get a reasonable show value"""

    tmpl = setup_basic(iname)
    out = str(tmpl).split("\n")
    assert out[0] == "bob"
    assert out[1].startswith("   Param  ")
    assert out[2].startswith("   -----  ")
    assert out[3].startswith("   bob.pa ")
    assert len(out) == 4


@pytest.mark.parametrize("iname", [None, "default"])
def test_templatemodel_basic_has_pars(iname):
    """Check we can get expected parameters"""

    tmpl = setup_basic(iname)
    assert len(tmpl.pars) == 1

    par = tmpl.pars[0]
    assert isinstance(par, Parameter)
    assert par.name == "pa"
    assert par.fullname == "bob.pa"
    assert not par.frozen
    assert par.val == pytest.approx(10)
    assert par.min == pytest.approx(10)
    assert par.max == pytest.approx(30)


@pytest.mark.parametrize("iname", [None, pytest.param("default", marks=pytest.mark.xfail)])
@pytest.mark.parametrize("pval,pname", [(10, "m1"), (20, "m2"), (30, "m3")])
def test_templatemodel_basic_can_query(iname, pval, pname):
    """Check we can query

    This does not work if interpolation is selected, but perhaps it
    should, which is why the test is marked xfail.

    """

    tmpl = setup_basic(iname)
    tbl = tmpl.query([pval])
    assert isinstance(tbl, TableModel)
    assert tbl.name == pname


@pytest.mark.parametrize("iname", [None, "default"])
@pytest.mark.parametrize("pval", [10, 20, 30])
def test_templatemodel_basic_can_evaluate(iname, pval):
    """Check we can evaluate the model for the parameter values."""

    tmpl = setup_basic(iname)

    # We know that the scaling is pval / 10 * [1, 2, 3, 4]
    #
    x = np.arange(1, 5)
    expected = x * pval / 10

    tmpl.pa = pval
    assert tmpl(x) == pytest.approx(expected)


def test_templatemodel_basic_can_evaluate_not_at_par_none():
    """There is no interpolation."""

    tmpl = setup_basic(None)
    x = np.arange(1, 5)
    tmpl.pa = 17.5

    with pytest.raises(ModelErr,
                       match="^Interpolation of template parameters was disabled for this model, "):
        tmpl(x)


def test_templatemodel_basic_can_evaluate_not_at_par_default():
    """There is interpolation within a parameter"""

    tmpl = setup_basic("default")
    x = np.arange(1, 5)
    tmpl.pa = 17.5

    # This is linear scaling, so it "should" be easy;
    #   pa = 10  -> 1, 2, 3, 4
    #   pa = 20  -> 2, 4, 6, 8
    #
    expected = 1.75 * x
    assert tmpl(x) == pytest.approx(expected)


@pytest.mark.parametrize("iname", [None, "default"])
@pytest.mark.parametrize("pval", [10, 20, 30])
def test_templatemodel_basic_can_interpolate_x(iname, pval):
    """Check we can evaluate the model for the parameter values withinterpolating x"""

    tmpl = setup_basic(iname)

    # We know that the scaling is pval / 10 * x when x is in the range
    # 1-4, so let's see how the interpolation does.
    #
    x = np.asarray([1.4, 2.5, 2.8, 3.9])
    expected = x * pval / 10

    tmpl.pa = pval
    assert tmpl(x) == pytest.approx(expected)


def test_templatemodel_basic_can_evaluate_not_at_par_default_interpolate_x():
    """There is interpolation within a parameter and within x"""

    tmpl = setup_basic("default")
    x = np.arange(1, 5)
    tmpl.pa = 23.5

    # We need to assume the interpolation both between the parameter
    # values and then within x works.
    #
    x = np.asarray([1.4, 2.5, 2.8, 3.9])
    expected = x * 2.35
    assert tmpl(x) == pytest.approx(expected)


@pytest.mark.parametrize("iname", [None, "default"])
def test_templatemodel_basic_fold(iname):
    """We can fold the dataset and use this correctly"""

    tmpl = setup_basic(iname, no_x=True)

    data = Data1D("orig", [1, 2, 3, 4], [1] * 4)
    data.ignore(xlo=3, xhi=3)

    # check we have a subset
    assert data.mask == pytest.approx([1, 1, 0, 1])

    tmpl.fold(data)

    # Pick a known template so we can handle both iname=None and default
    tmpl.pa = 20

    assert data.eval_model_to_fit(tmpl) == pytest.approx([2, 4, 8])
    assert data.eval_model(tmpl) == pytest.approx([2, 4, 6, 8])


def weighted(xs):
    """Sum up the weight values."""
    w = sum(x[0] for x in xs)
    return sum(x[0] * x[1] for x in xs) / w


# k=1 means nearest match
# k=2 means linear interpolation
# k=3 means weighted nearest 3 matches
#
# It is unclear what the code will do when values are equidistant, so
# these tests may need to check for "a or b" rather than just a, but
# wait until that becomes a problem before worrying about it.
#
@pytest.mark.parametrize("order", [1, 2, 3])
@pytest.mark.parametrize("k,pval,expected",
                         [(1, 11, 5),
                          (1, 15, 5),
                          (1, 27, 20),
                          (1, 31, 20),
                          (1, 39, 40),

                          (2, 11, weighted([(0.9, 5), (0.1, 15)])),
                          (2, 15, weighted([(0.5, 5), (0.5, 15)])),
                          (2, 27, weighted([(0.3, 15), (0.7, 20)])),
                          (2, 31, weighted([(0.9, 20), (0.1, 40)])),
                          (2, 39, weighted([(0.1, 20), (0.9, 40)])),

                          # Weights are 1/<distance>
                          (3, 11, weighted([(1/1, 5), (1/9, 15), (1/19, 20)])),
                          (3, 15, weighted([(1/5, 5), (1/5, 15), (1/15, 20)])),
                          (3, 27, weighted([(1/7, 15), (1/3, 20), (1/13, 40)])),
                          (3, 31, weighted([(1/11, 15), (1/1, 20), (1/9, 40)])),
                          (3, 39, weighted([(1/9, 20), (1/1, 40), (1/11, 44)])),
                          ])
def test_interpolate_par1(order, k, pval, expected):
    """Check KNNInterpolator/Template distance calculation: 1 par

    The idea is that we use a "constant" model so we can know what
    the interpolated value should be.

    I was surprised that order=1, 2, 3 all gave the same results,
    but we only have one parameter so there's no real "power"
    for the order value to address.

    """

    x = np.asarray([2, 4, 6])
    templates = []
    for scale, name in zip([5, 15, 20, 40, 44], ["m1", "m2", "m3", "m4", "m5"]):
        template = TableModel(name)
        template.load(x, np.ones_like(x) * scale)
        templates.append(template)

    tmpl = create_template_model("bob", ["pa"], np.asarray([[10], [20], [30], [40], [50]]),
                                 templates)
    assert isinstance(tmpl, Template)
    assert tmpl.k == 2
    assert tmpl.order == 2

    # The x axis values don't really matter here, but keep them inside
    # the original x=2, 4, 6 setting. A single value could be used,
    # but include 2 just in case.
    #
    tmpl.order = order
    tmpl.k = k
    tmpl.pa = pval
    assert tmpl([3, 5]) == pytest.approx([expected, expected])


def check3(p1, p2, p3):
    """What do we expect the answer to be?"""
    x = np.asarray([2, 4, 6])
    return p1 * x * x + p2 * x + p3


@pytest.mark.parametrize("order,k,pval1,pval2,pval3,expected",
                         # First checks are linear interpolation along
                         # p3 axis at points where we have p1 and p2.
                         #
                         [(1, 2, 2, 3, 25, check3(2, 3, 25)),
                          (2, 2, 2, 3, 25, check3(2, 3, 25)),
                          (3, 2, 2, 3, 25, check3(2, 3, 25)),
                          (1, 2, 0.5, 1, 15, check3(0.5, 1, 15)),
                          (2, 2, 0.5, 1, 15, check3(0.5, 1, 15)),
                          (3, 2, 0.5, 1, 15, check3(0.5, 1, 15)),

                          # Actual values are [39, 57, 81]
                          # Why does the order setting not matter?
                          #
                          (1, 2, 0.75, 4.5, 27, [41, 58, 81]),
                          (2, 2, 0.75, 4.5, 27, [41, 58, 81]),
                          (3, 2, 0.75, 4.5, 27, [41, 58, 81])
                          ])
def test_interpolate_par3(order, k, pval1, pval2, pval3, expected):
    """Check KNNInterpolator/Template distance calculation: 3 pars

    The underlying model is

        y(p1, p2, p3) = p1 x^2 + p2 x + p3

    """

    x = np.asarray([2, 4, 6])
    templates = []

    # Here we use the same parameter values for the models
    # and the parameter values, unlike the test_interpolate_par1
    # case, where we had different scaling.
    #
    pvals1 = [0.5, 1, 1.5, 2, 2.5]
    pvals2 = [1, 2, 3, 4, 5]
    pvals3 = [10, 20, 30, 40, 50]

    pvals = []
    for p3 in pvals3:
        for p2 in pvals2:
            for p1 in pvals1:
                template = TableModel()
                template.load(x, p1 * x * x + p2 * x + p3)
                templates.append(template)
                pvals.append([p1, p2, p3])

    pvals = np.asarray(pvals)
    tmpl = create_template_model("bob", ["pa", "pb", "pc"], pvals,
                                 templates)
    assert isinstance(tmpl, Template)
    assert tmpl.k == 2
    assert tmpl.order == 2

    # To simplify comparison, we use the actual x array rather than
    # add on more interpolation.
    #
    tmpl.order = order
    tmpl.k = k
    tmpl.pa = pval1
    tmpl.pb = pval2
    tmpl.pc = pval3
    assert tmpl(x) == pytest.approx(expected)


@pytest.mark.parametrize("pa,expected",
                         [(1, [4.8, 6, 7.6, 8.2, 10]),
                          (1.5, [10.4, 11, 16.466667, 18.516667, 24.66667]),
                          (2, [16, 16, 25.333333, 28.83333, 39.33333])])
def test_template_with_different_independent_axes(pa, expected):
    """Check we can have different x values."""

    t1 = TableModel("m1")
    t1.load([10, 20, 30], [4, 6, 8])

    t2 = TableModel("m2")
    t2.load([12, 20, 32], [16, 16, 30])

    templates = [t1, t2]

    mdl = create_template_model("bob", ["pa"], np.asarray([[1], [2]]),
                                templates)

    # Evaluate in the grid 14, 20, 28, 31, 40
    #
    mdl.pa = pa
    got = mdl([14, 20, 28, 31, 40])
    assert got == pytest.approx(expected)
