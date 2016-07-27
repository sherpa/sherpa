from __future__ import print_function
#
#  Copyright (C) 2012, 2015  Smithsonian Astrophysical Observatory
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

import six

from sherpa.utils import SherpaTestCase, requires_data
from sherpa.models import ArithmeticModel, Parameter
from sherpa.models.basic import PowLaw1D
from sherpa.utils.err import StatErr, SessionErr
from sherpa import ui
import numpy
import logging

logger = logging.getLogger("sherpa")


class UserModel(ArithmeticModel):
    def __init__(self, name='usermodel'):
        self.param1 = Parameter(name, 'param1', 1, min=0, max=100)
        self.param2 = Parameter(name, 'param2', 1, min=-100, max=100)

        ArithmeticModel.__init__(self, name, (self.param1,
                                              self.param2))

    def calc(self, p, x, *args, **kwargs):
        return p[0] * x + p[1]


@requires_data
class test_get_draws(SherpaTestCase):
    """
    Tests for PR #155

    TODO: Some test cases would be more readable if using pytest to parameterize a test case
    TODO: Some test cases cannot be implemented (easily) without mock
    In particular one cannot test that the correct matrix is used when one is not
    provided and that an error is thrown when, for any reason, the computed covariance matrix
    is None.
    """

    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        ui.clean()

        self.ascii = self.make_path('sim.poisson.1.dat')

        self.wrong_stat_msg = "Fit statistic must be cash, cstat or wstat, not {}"
        self.wstat_err_msg = "No background data has been supplied. Use cstat"
        self.no_covar_msg = "covariance has not been performed"
        self.fail_msg = "Call should not have succeeded"
        self.right_stats = {'cash', 'cstat', 'wstat'}
        self.model = PowLaw1D("p1")

        ui.load_data(self.ascii)
        ui.set_model(self.model)

    def tearDown(self):
        if hasattr(self, '_old_logger_level'):
            logger.setLevel(self._old_logger_level)
        ui.clean()

    # Test an exception is thrown is the proper stat is not set
    def test_covar_wrong_stat(self):
        ui.covar()
        fail = False
        wrong_stats = set(ui.list_stats()) - self.right_stats
        for stat in wrong_stats:
            ui.set_stat(stat)
            try:
                ui.get_draws()
            except ValueError as ve:
                self.assertEqual(self.wrong_stat_msg.format(stat), str(ve))
                continue
            fail = True
            break
        if fail:
            self.fail(self.fail_msg)

    # Test an exception is thrown when wstat is used without background
    def test_covar_wstat_no_background(self):
        ui.covar()
        ui.set_stat("wstat")
        try:
            ui.get_draws()
        except StatErr as ve:
            self.assertEqual(self.wstat_err_msg, str(ve))
            return
        self.fail(self.fail_msg)

    # Test an exception is thrown if covar is not run
    def test_no_covar(self):
        for stat in self.right_stats:
            ui.set_stat(stat)
            try:
                ui.get_draws()
            except SessionErr as ve:
                self.assertEqual(self.no_covar_msg, str(ve))
                return
        self.fail(self.fail_msg)

    # Test get_draws returns a valid response when the covariance matrix is provided
    # Note the accuracy of the returned values is not assessed here
    def test_covar_as_argument(self):
        for stat in self.right_stats - {'wstat'}:
            ui.set_stat(stat)
            ui.fit()
            matrix = [[0.00064075, 0.01122127], [0.01122127, 0.20153251]]
            niter = 10
            stat, accept, params = ui.get_draws(niter=niter, covar_matrix=matrix)
            self.assertEqual(niter + 1, stat.size)
            self.assertEqual(niter + 1, accept.size)
            self.assertEqual((2, niter + 1), params.shape)
            self.assertTrue(numpy.any(accept))

    # Test get_draws returns a valid response when the covariance matrix is not provided
    # Note the accuracy of the returned values is not assessed here
    def test_covar_as_none(self):
        for stat in self.right_stats - {'wstat'}:
            ui.set_stat(stat)
            ui.fit()
            ui.covar()
            niter = 10
            stat, accept, params = ui.get_draws(niter=niter)
            self.assertEqual(niter + 1, stat.size)
            self.assertEqual(niter + 1, accept.size)
            self.assertEqual((2, niter + 1), params.shape)
            self.assertTrue(numpy.any(accept))


@requires_data
class test_ui(SherpaTestCase):
    def setUp(self):
        ui.clean()
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        self.ascii = self.make_path('sim.poisson.1.dat')
        self.single = self.make_path('single.dat')
        self.double = self.make_path('double.dat')
        self.filter = self.make_path('filter_single_integer.dat')
        self.func = lambda x: x

        ui.dataspace1d(1, 1000, dstype=ui.Data1D)

    def tearDown(self):
        if hasattr(self, '_old_logger_level'):
            logger.setLevel(self._old_logger_level)
        ui.clean()

    def test_ascii(self):
        ui.load_data(1, self.ascii)
        ui.load_data(1, self.ascii, 2)
        ui.load_data(1, self.ascii, 2, ("col2", "col1"))

    # Test table model
    def test_table_model_ascii_table(self):
        ui.load_table_model('tbl', self.single)
        ui.load_table_model('tbl', self.double)

    # Test user model
    def test_user_model_ascii_table(self):
        ui.load_user_model(self.func, 'mdl', self.single)
        ui.load_user_model(self.func, 'mdl', self.double)

    def test_filter_ascii(self):
        ui.load_filter(self.filter)
        ui.load_filter(self.filter, ignore=True)

    def test_add_model(self):
        ui.add_model(UserModel)
        ui.set_model('usermodel.user1')

    def test_set_full_model(self):
        ui.load_psf('psf1', 'gauss2d.g1')
        ui.set_full_model('psf1(gauss2d.g2)+const2d.c1')
        ui.get_model()

    #        ui.get_source()

    # Bug 12644
    def test_source_methods_with_full_model(self):
        from sherpa.utils.err import IdentifierErr

        ui.load_data('full', self.ascii)
        ui.set_full_model('full', 'powlaw1d.p1')

        # Test Case 1
        try:
            ui.get_source('full')
        except IdentifierErr as e:
            six.assertRegex(self, str(e),
                                     "Convolved model\n.*\n is set for dataset full. You should use get_model instead.",
                                     str(e))
        try:
            ui.plot_source('full')
        except IdentifierErr as e:
            six.assertRegex(self, str(e),
                                     "Convolved model\n.*\n is set for dataset full. You should use plot_model instead.",
                                     str(e))

        # Test Case 2
        ui.set_source('full', 'powlaw1d.p2')
        ui.get_source('full')

        # Test Case 3
        ui.load_data('not_full', self.ascii)
        try:
            ui.get_source('not_full')
        except IdentifierErr as e:
            self.assertEqual('source not_full has not been set, consider using set_source() or set_model()', str(e))


class test_psf_ui(SherpaTestCase):
    models1d = ['gauss1d', 'delta1d', 'normgauss1d']
    models2d = ['gauss2d', 'delta2d', 'normgauss2d']

    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

    def tearDown(self):
        if hasattr(self, '_old_logger_level'):
            logger.setLevel(self._old_logger_level)

    def test_psf_model1d(self):
        ui.dataspace1d(1, 10)
        for model in self.models1d:
            try:
                ui.load_psf('psf1d', model + '.mdl')
                ui.set_psf('psf1d')
                mdl = ui.get_model_component('mdl')
                self.assertTrue((numpy.array(mdl.get_center()) ==
                                 numpy.array([4])).all())
            except:
                print(model)
                raise

    def test_psf_model2d(self):
        ui.dataspace2d([216, 261])
        for model in self.models2d:
            try:
                ui.load_psf('psf2d', model + '.mdl')
                ui.set_psf('psf2d')
                mdl = ui.get_model_component('mdl')
                self.assertTrue((numpy.array(mdl.get_center()) ==
                                 numpy.array([108, 130])).all())
            except:
                print(model)
                raise
