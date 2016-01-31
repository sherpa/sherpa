#
#  Copyright (C) 2012, 2015, 2016  Smithsonian Astrophysical Observatory
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
import sys

from sherpa.utils import SherpaTest, SherpaTestCase, requires_data
from sherpa.utils import linear_interp, nearest_interp, neville
from sherpa.utils.err import IdentifierErr, IOErr, ModelErr
from sherpa.models import ArithmeticModel, Parameter
from sherpa.models.basic import TableModel
from sherpa import ui

import numpy as np
from numpy.testing import assert_allclose


class UserModel(ArithmeticModel):

    def __init__(self, name='usermodel'):
        self.param1 = Parameter(name, 'param1', 1, min=0, max=100)
        self.param2 = Parameter(name, 'param2', 1, min=-100, max=100)

        ArithmeticModel.__init__(self, name, (self.param1,
                                              self.param2))

    def calc(self, p, x, *args, **kwargs):
        return p[0] * x + p[1]


@requires_data
class test_ui(SherpaTestCase):

    def setUp(self):
        self.ascii = self.make_path('threads/ascii_table/sim.poisson.1.dat')
        self.single = self.make_path('single.dat')
        self.double = self.make_path('double.dat')
        self.filter = self.make_path('filter_single_integer.dat')
        self.func = lambda x: x

        ui.dataspace1d(1, 1000, dstype=ui.Data1D)

    def test_ascii(self):
        ui.load_data(1, self.ascii)
        ui.load_data(1, self.ascii, 2)
        ui.load_data(1, self.ascii, 2, ("col2", "col1"))

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
        # ui.get_source()

    # Bug 12644
    def test_source_methods_with_full_model(self):

        ui.load_data('full', self.ascii)
        ui.set_full_model('full', 'powlaw1d.p1')

        # depending on how the test is run the model name can be
        # reported as 'p1' or 'powlaw1d.p1', so use a regexp.
        #
        def mk_regexp(func_head):
            return "Convolved model\n.*\n is set for dataset full. " + \
                "You should use {}_model instead.".format(func_head)

        # Test Case 1
        try:
            ui.get_source('full')
        except IdentifierErr as e:
            re = mk_regexp('get')
            self.assertRegexpMatches(str(e), re, msg=str(e))

        try:
            ui.plot_source('full')
        except IdentifierErr as e:
            re = mk_regexp('plot')
            self.assertRegexpMatches(str(e), re, msg=str(e))

        # Test Case 2
        ui.set_source('full', 'powlaw1d.p2')
        ui.get_source('full')

        # Test Case 3
        ui.load_data('not_full', self.ascii)
        try:
            ui.get_source('not_full')
        except IdentifierErr as e:
            emsg = 'source not_full has not been set, consider ' + \
                   'using set_source() or set_model()'
            self.assertEquals(emsg, str(e))


# For now have the data files as part of the Sherpa
# repository, rather than the sherpa-test-data submodule
#
# @requires_data
class BaseTableModelTestCase:   # important not to derive from (SherpaTestCase):

    # Where should the functions used to set/query/change the
    # Sherpa be accessed.
    #
    state = None

    # What interpolation functions should be used. Note that
    # the code is somewhat hard-coded to these values (e.g. for
    # selecting tolerances), so perhaps more, or less, abstraction
    # is needed.
    #
    interp1d = [linear_interp, nearest_interp, neville]

    # would like to over-ride make_path, but it is used
    # for the FITS tests, so need a separate one.
    #
    def make_local_path(self, fname):
        """Use local data directory"""
        raise NotImplementedError

    def setUp(self):
        self.ascii_onecol = self.make_local_path('gauss1d-onecol.dat')
        self.ascii_twocol = self.make_local_path('gauss1d.dat')
        self.ascii_threecol = self.make_local_path('gauss1d-error.dat')

        # the values in self.ascii_threecol
        dtype = np.float32
        self.x = np.asarray([80, 95, 110, 125, 140, 155, 170,
                             185, 200], dtype=dtype)
        self.y = np.asarray([0, 0, 9, 35, 93, 96, 49, 15, 0],
                            dtype=dtype)
        self.dy = np.asarray([1.86603, 1.86603, 4.1225, 6.97913,
                              10.6825, 10.8362, 8.05337, 4.96863,
                              1.86603], dtype=dtype)

    # should really be in SherpaTestCase
    def tearDown(self):
        self.state.clean()

    def test_basic(self):
        """Check that a table can be created and has expected properties"""

        flag = 'tbl' in self.state.list_model_components()
        self.assertFalse(flag)

        self.state.load_table_model('tbl', self.ascii_onecol)

        flag = 'tbl' in self.state.list_model_components()
        self.assertTrue(flag,
                        msg='Model component is added to system database')

        mdl = self.state.get_model_component('tbl')
        self.assertTrue(isinstance(mdl, TableModel), msg='Is a Table Model')

        # There is something strange with the name field, as it can be
        # either, depending on how the test is run. So skip the test for
        # now (although it suggests that there is something flaky with
        # the name logic for models).
        #
        # self.assertEqual('tablemodel.tbl', mdl.name)
        # self.assertEqual('tbl', mdl.name)

        # The x/y values are checked in later tests, so do not
        # include that here.

        pars = mdl.pars
        self.assertEqual(1, len(pars), msg='One parameter')

        par = pars[0]
        self.assertEqual("ampl", par.name)
        self.assertEqual(False, par.frozen, msg='Is ampl frozen')

        self.assertAlmostEqual(1.0, par.val)

    # def test_fail_on_missing_col(self):
    #     """Error out if a column is missing."""
    #
    #     self.assertRaises(IOErr, self.state.load_table_model, 'failed',
    #                       self.ascii_twocol, colkeys=['a', 'b'])

    def _checkcol(self, expvals, gotvals):
        """Check a column"""

        if expvals is None:
            self.assertEqual(None, gotvals)
            return

        # There is no promise about the data types of the columns,
        # so just check that it is a real, floating-point value.
        self.assertIn(gotvals.dtype, (np.float32, np.float64))
        assert_allclose(expvals, gotvals)

    def test_ascii_table1(self):
        """Read in a one-column ASCII file"""

        self.state.load_table_model('tbl1', self.ascii_onecol)
        mdl = self.state.get_model_component('tbl1')
        self._checkcol(None, mdl.get_x())
        self._checkcol(self.y, mdl.get_y())

    def test_ascii_table2(self):
        """Read in a two-column ASCII file"""
        self.state.load_table_model('tbl2', self.ascii_twocol)
        mdl = self.state.get_model_component('tbl2')
        self._checkcol(self.x, mdl.get_x())
        self._checkcol(self.y, mdl.get_y())

    def test_ascii_table3_col12(self):
        """Read in a three-column ASCII file: col1 col2"""
        self.state.load_table_model('tbl3', self.ascii_threecol)
        mdl = self.state.get_model_component('tbl3')
        self._checkcol(self.x, mdl.get_x())
        self._checkcol(self.y, mdl.get_y())

    def test_ascii_table3_col13(self):
        """Read in a three-column ASCII file: col1 col3"""
        self.state.load_table_model('tbl3', self.ascii_threecol,
                                    colkeys=['X', 'STAT_ERR'])
        mdl = self.state.get_model_component('tbl3')
        self._checkcol(self.x, mdl.get_x())
        self._checkcol(self.dy, mdl.get_y())

    def test_ascii_table3_col23(self):
        """Read in a three-column ASCII file: col2 col3"""
        self.state.load_table_model('tbl3', self.ascii_threecol,
                                    colkeys=['Y', 'STAT_ERR'])
        mdl = self.state.get_model_component('tbl3')
        idx = np.argsort(self.y)
        self._checkcol(self.y[idx], mdl.get_x())
        self._checkcol(self.dy[idx], mdl.get_y())

    def test_ascii_table3_ncols1(self):
        """Read in a three-column ASCII file: ncols=1"""
        self.state.load_table_model('tbl1', self.ascii_threecol, ncols=1)
        mdl = self.state.get_model_component('tbl1')
        self._checkcol(None, mdl.get_x())
        self._checkcol(self.x, mdl.get_y())

    def test_ascii_table3_ncols2(self):
        """Read in a three-column ASCII file: ncols=2"""
        self.state.load_table_model('tbl3', self.ascii_threecol, ncols=2)
        mdl = self.state.get_model_component('tbl3')
        self._checkcol(self.x, mdl.get_x())
        self._checkcol(self.y, mdl.get_y())

    # Since the previous tests have checked that the data has
    # been read in, the evaluation tests do not have to cover
    # all these cases, just some basic ones:
    #   - 2 column with and with-out filtering
    #   - integrated (probably not really needed)
    #
    def test_ascii_eval(self):
        """Basic test of evaluating the table model"""

        norm = 2.1e6
        self.state.set_stat('leastsq')

        for interp in self.interp1d:
            self.state.load_table_model('tbl', self.ascii_twocol,
                                        method=interp)
            mdl = self.state.get_model_component('tbl')
            mdl.ampl = norm

            self.state.load_arrays(1, self.x, self.y * norm,
                                   self.state.Data1D)
            self.state.set_source(mdl)

            res = self.state.get_stat_info()[0]
            self.assertAlmostEqual(0.0, res.statval)
            self.assertEqual(9, res.numpoints)

            self.state.ignore_id(1, None, self.x[1])
            self.state.ignore_id(1, self.x[3], self.x[4])

            res = self.state.get_stat_info()[0]
            self.assertAlmostEqual(0.0, res.statval)
            self.assertEqual(5, res.numpoints)

    def test_ascii_eval_interp(self):
        """Basic test of evaluating the table model: interpolation"""

        norm = 2.1e6
        self.state.set_stat('leastsq')

        for interp in self.interp1d:
            self.state.load_table_model('tbl', self.ascii_twocol,
                                        method=interp)
            mdl = self.state.get_model_component('tbl')
            mdl.ampl = norm

            self.state.load_arrays(1, self.x[:-1], self.x[1:],
                                   self.y[:-1] * norm,
                                   self.state.Data1DInt)
            self.state.set_source(mdl)

            res = self.state.get_stat_info()[0]
            self.assertAlmostEqual(0.0, res.statval)
            self.assertEqual(8, res.numpoints)

            self.state.ignore_id(1, None, self.x[1])
            self.state.ignore_id(1, self.x[3], self.x[4])

            res = self.state.get_stat_info()[0]
            self.assertAlmostEqual(0.0, res.statval)
            self.assertEqual(3, res.numpoints)

    def test_ascii_onecol_fail(self):
        """Model evaluation fails if #rows does not match."""

        self.state.load_table_model('fmodel', self.ascii_onecol)
        self.state.load_arrays(99, self.x[1:], self.y[1:])
        self.state.set_source(99, 'fmodel')
        self.assertRaises(ModelErr, self.state.calc_stat, 99)


class test_table_model(BaseTableModelTestCase, SherpaTestCase):

    state = ui

    def make_local_path(self, fname):
        """Use local data directory"""
        # Is there a better way than this?
        thisfile = sys.modules[self.__module__].__file__
        thisdir = os.path.dirname(thisfile)
        return os.path.join(thisdir, 'data', fname)

    # Is it worth checking that FITS files are not supported (it
    # does trigger an error condition so improves code coverage)?
    #
    @requires_data
    def test_fits_errors_out(self):
        """FITS files are unsupported"""
        fname = self.make_path('double.fits')
        self.assertRaises(IOErr,
                          self.state.load_table_model, 'ftbl', fname)

    # TODO: why does this test fail in sherpa.astro.ui - i.e. no
    #       error is thrown? Have moved out of Base...TestCase for now
    def test_fail_on_missing_col(self):
        """Error out if a column is missing."""

        self.assertRaises(IOErr, self.state.load_table_model, 'failed',
                          self.ascii_twocol, colkeys=['a', 'b'])


class test_psf_ui(SherpaTestCase):

    models1d = ['gauss1d', 'delta1d', 'normgauss1d']
    models2d = ['gauss2d', 'delta2d', 'normgauss2d']

    # Commented out setUp/tearDown as they do nothing
    #
    # def setUp(self):
    #     pass

    def tearDown(self):
        ui.clean()

    def test_psf_model1d(self):
        ui.dataspace1d(1, 10)
        for model in self.models1d:
            try:
                ui.load_psf('psf1d', model + '.mdl')
                ui.set_psf('psf1d')
                mdl = ui.get_model_component('mdl')
                self.assertTrue((np.array(mdl.get_center()) ==
                                 np.array([4])).all())
            except:
                print model
                raise

    def test_psf_model2d(self):
        ui.dataspace2d([216, 261])
        for model in self.models2d:
            try:
                ui.load_psf('psf2d', model + '.mdl')
                ui.set_psf('psf2d')
                mdl = ui.get_model_component('mdl')
                self.assertTrue((np.array(mdl.get_center()) ==
                                 np.array([108, 130])).all())
            except:
                print model
                raise


if __name__ == '__main__':

    if len(sys.argv) > 1:
        datadir = sys.argv[1]
    else:
        datadir = None

    SherpaTest(ui).test(datadir=datadir)
