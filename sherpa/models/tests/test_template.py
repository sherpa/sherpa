#
#  Copyright (C) 2011, 2015, 2016, 2018, 2021, 2022, 2023
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

from sherpa.models import TableModel, Gauss1D
from sherpa.models.template import create_template_model
from sherpa.utils.testing import requires_data
from sherpa.utils.err import ModelErr
from sherpa import ui


# Regression reported by Aneta (original report by Fabrizio Nicastro)
# When restoring a file saved with an older version of sherpa,
# we need to make sure models are assigned an is_discrete field.
# For model that do not have the is_discrete field, fallback to False.
#
@requires_data
def test_restore_is_discrete(run_thread, clean_ui):
    # TODO: test.py test is a dummy test. We need to implement a real
    # test. Take a look at the old testdata/ciao4.3/template_restore
    # directory for help.
    run_thread('template_restore', 'test.py')


# TestCase 1 load_template_model enables interpolation by default
@requires_data
def test_load_template_with_interpolation(run_thread, clean_ui):
    run_thread('load_template_with_interpolation')
    pvals = ui.get_fit_results().parvals
    pmin = pvals[0]
    pmax = pvals[1]
    if pmax < pmin:
        (pmin, pmax) = (pmax, pmin)

    tol = 0.001
    assert pmin == pytest.approx(2023.46, rel=tol)
    assert pmax == pytest.approx(2743.47, rel=tol)


@requires_data
def test_load_template_interpolator(run_thread, clean_ui):
    run_thread('load_template_interpolator')
    pval = ui.get_fit_results().parvals[0]
    assert pval == pytest.approx(2743.91, rel=0.001)


# TestCase 2 load_template_model with template_interpolator_name=None
# disables interpolation
#
# TestCase 3.1 discrete templates fail when probed for values they do
# not represent (gridsearch with wrong values)
@requires_data
def test_load_template_model_without_interpolation_2_31(run_thread, clean_ui):
    emsg = 'Interpolation of template parameters was disabled for this model, but parameter values not in the template library have been requested. Please use gridsearch method and make sure the sequence option is consistent with the template library'

    with pytest.raises(ModelErr, match=emsg):
        run_thread('load_template_without_interpolation',
                   scriptname='test_case_2_and_3.1.py')


# TestCase 3.2 discrete templates fail when probed for values they do
# not represent (continuous method with discrete template)
#
@requires_data
def test_load_template_model_without_interpolation_32(run_thread, clean_ui):
    emsg = r"You are trying to fit a model which has a discrete template model component with a continuous optimization method. Since CIAO4.6 this is not possible anymore. Please use gridsearch as the optimization method and make sure that the 'sequence' option is correctly set, or enable interpolation for the templates you are loading \(which is the default behavior\)."
    with pytest.raises(ModelErr, match=emsg):
        run_thread('load_template_without_interpolation',
                   scriptname='test_case_3.2.py')


# TestCase 4 gridsearch with right values succeeds
@requires_data
def test_grid_search_with_discrete_template(run_thread, clean_ui):
    run_thread('load_template_without_interpolation',
               scriptname='test_case_4.py')


# TestCase 5 user can access interpolators' parvals
@requires_data
def test_grid_search_with_discrete_template_parvals(run_thread, clean_ui):
    run_thread('load_template_with_interpolation',
               scriptname='test_case_5.py')


@pytest.fixture
def setUp():
    x = numpy.linspace(0.1, 5, 50)

    num = 4
    ntemplates = 2**num
    g1 = Gauss1D('g1')

    # create a 4-dimensional grid from 0 to 1 inclusive, shape = (16,4)
    grid = numpy.mgrid[[slice(0, 2, 1) for ii in range(num)]]
    grid = numpy.asarray(list(map(numpy.ravel, grid))).T
    coords = numpy.linspace(0.01, 6, 100)
    names = [f"p{i}" for i in range(num)]

    # Since we randomize the parameters use for the gaussian models,
    # use a fixed seed.
    #
    rng = numpy.random.default_rng(9427205)

    templates = []
    for ii in range(ntemplates):
        t = TableModel()
        g1.fwhm = rng.uniform(0.5, 2.0)
        g1.pos = rng.uniform(1.0, 4.5)
        g1.ampl = rng.uniform(1.0, 50.)
        t.load(coords, g1(coords))
        templates.append(t)

    return x, create_template_model("mdl", names, grid, templates)


def test_template_model_evaluation(setUp):
    x, model = setUp
    model.thawedpars = [0, 1, 0, 1]

    # Evaluate on a very-restricted grid and check the result (this is
    # a regression test and so the predicted values may change).
    #
    y = model([0.55, 1.25, 3.65, 4.95])
    assert y == pytest.approx([0.006695045904199339,
                               0.33558934132320956,
                               5.274350209830429,
                               0.024472491880120843])
