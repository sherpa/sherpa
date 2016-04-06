#
#  Copyright (C) 2011, 2015, 2016  Smithsonian Astrophysical Observatory
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


from sherpa.models import TableModel, Gauss1D
from sherpa.models.template import create_template_model
from sherpa.utils import SherpaTestCase, requires_data
from sherpa.utils.err import ModelErr
from sherpa import ui
import numpy
import logging
import os

logger = logging.getLogger("sherpa")


@requires_data
class test_new_templates_ui(SherpaTestCase):
    def assign_model(self, name, obj):
        self.locals[name] = obj

    def run_thread(self, name, scriptname='fit.py'):
        ui.clean()
        ui.set_model_autoassign_func(self.assign_model)
        super(test_new_templates_ui, self).run_thread(name, scriptname=scriptname)

    def setUp(self):
        self.loggingLevel = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

    def tearDown(self):
        logger.setLevel(self.loggingLevel)

    # Regression reported by Aneta (original report by Fabrizio Nicastro)
    # When restoring a file saved with an older version of sherpa,
    # we need to make sure models are assigned an is_discrete field.
    # For model that do not have the is_discrete field, fallback to False.
    def test_restore_is_discrete(self):
        # TODO: test.py test is a dummy test. We need to implement a real
        # test. Take a look at the old testdata/ciao4.3/template_restore
        # directory for help.
        self.run_thread('template_restore', 'test.py')

    # TestCase 1 load_template_model enables interpolation by default
    def test_load_template_with_interpolation(self):
        self.run_thread('load_template_with_interpolation')
        pvals = ui.get_fit_results().parvals
        pmin = pvals[0]
        pmax = pvals[1]
        if pmax < pmin:
            (pmin, pmax) = (pmax, pmin)

        tol = 0.001
        self.assertEqualWithinTol(2023.46, pmin, tol)
        self.assertEqualWithinTol(2743.47, pmax, tol)

    def test_load_template_interpolator(self):
        self.run_thread('load_template_interpolator')
        pval = ui.get_fit_results().parvals[0]
        self.assertEqualWithinTol(2743.91, pval, 0.001)

    # TestCase 2 load_template_model with template_interpolator_name=None
    # disables interpolation
    #
    # TestCase 3.1 discrete templates fail when probed for values they do
    # not represent (gridsearch with wrong values)
    def test_load_template_model_without_interpolation_2_31(self):
        try:
            self.run_thread('load_template_without_interpolation',
                            scriptname='test_case_2_and_3.1.py')
        except ModelErr:
            return
        self.fail('Fit should have failed: gridsearch with wrong parvals')

    # TestCase 3.2 discrete templates fail when probed for values they do
    # not represent (continuous method with discrete template)
    #
    def test_load_template_model_without_interpolation_32(self):
        try:
            self.run_thread('load_template_without_interpolation',
                            scriptname='test_case_3.2.py')
        except ModelErr:
            return
        self.fail('Fit should have failed: gridsearch with wrong parvals')

    # TestCase 4 gridsearch with right values succeeds
    def test_grid_search_with_discrete_template(self):
        self.run_thread('load_template_without_interpolation',
                        scriptname='test_case_4.py')

    # TestCase 5 user can access interpolators' parvals
    def test_grid_search_with_discrete_template_parvals(self):
        self.run_thread('load_template_with_interpolation',
                        scriptname='test_case_5.py')


class test_template(SherpaTestCase):

    def setUp(self):
        self.num = 4
        self.ncoords = 100
        self.ntemplates = 2**self.num
        self.x = numpy.linspace(0.1, 5, 50)
        g1 = Gauss1D('g1')

        # create a 4-dimensional grid from 0 to 1 inclusive, shape = (16,4)
        grid = numpy.mgrid[[slice(0, 2, 1) for ii in range(self.num)]]
        grid = numpy.asarray(map(numpy.ravel, grid)).T
        coords = numpy.linspace(0.01, 6, 100)
        names = ["p%i" % i for i in range(self.num)]
        templates = []
        for ii in range(self.ntemplates):
            t = TableModel()
            g1.fwhm = numpy.random.uniform(0.5, 2.0)
            g1.pos  = numpy.random.uniform(1.0, 4.5)
            g1.ampl = numpy.random.uniform(1.0, 50.)
            t.load(coords, g1(coords))
            templates.append(t)

        self.model = create_template_model("mdl", names, grid, templates)

    def tearDown(self):
        self.model = None

    def test_template_model_evaluation(self):
        self.model.thawedpars = [0, 1, 0, 1]
        # We want to evaluate the model, but do not check the result
        self.model(self.x)

#    def test_template_query_index(self):
#        expected = 5
#        result = self.model.query_index([0,1,0,1])
#        self.assertEqual(expected, result)
#
#    def test_template_query(self):
#        result = self.model.query([0,1,0,1])
