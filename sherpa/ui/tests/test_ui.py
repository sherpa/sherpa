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
class test_table_model(SherpaTestCase):

    interp1d = [linear_interp, nearest_interp, neville]

    # would like to over-ride make_path, but this method is used
    # for the FITS tests
    def make_local_path(self, fname):
        """Use local data directory, rather than sherpa-test-data"""
        thisfile = sys.modules[self.__module__].__file__
        thisdir = os.path.dirname(thisfile)
        return os.path.join(thisdir, 'data', fname)

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
        ui.clean()

    def test_basic(self):
        """Check that a table can be created and has expected properties"""

        self.assertTrue('tbl' not in ui.list_model_components())
        ui.load_table_model('tbl', self.ascii_onecol)
        self.assertTrue('tbl' in ui.list_model_components(),
                        msg='Model component is added to system database')
        mdl = ui.get_model_component('tbl')
        self.assertTrue(isinstance(mdl, TableModel), msg='Is a Table Model')

        # There is something strange with the name field, as it can be
        # either, depending on how the test is run. So skip the test for
        # now (although it suggests that there is something flaky with
        # the name logic for models).
        #
        # self.assertEqual('tablemodel.tbl', mdl.name)
        # self.assertEqual('tbl', mdl.name)

        pars = mdl.pars
        self.assertEqual(1, len(pars), msg='One parameter')

        par = pars[0]
        self.assertEqual("ampl", par.name)
        self.assertEqual(False, par.frozen, msg='Is ampl frozen')

        self.assertAlmostEqual(1.0, par.val)

    def test_fail_on_missing_col(self):
        """Error out if a column is missing."""

        self.assertRaises(IOErr, ui.load_table_model, 'failed',
                          self.ascii_twocol, colkeys=['a', 'b'])

    def _test_table1(self, tname):
        """Tests for a one-column file"""

        mdl = ui.get_model_component(tname)
        x = mdl.get_x()
        y = mdl.get_y()

        self.assertEqual(None, x, msg='No X axis for table')

        y = mdl.get_y()
        self.assertEqual(9, y.size,
                         msg='Correct #rows for table')

        # I do not think we guarantee the type of the column, so
        # support both 32 and 64 bit data types
        self.assertIn(y.dtype, (np.float32, np.float64))

    def _test_table2(self, tname):
        """Tests for a two-column file"""
        mdl = ui.get_model_component(tname)

        x = mdl.get_x()
        y = mdl.get_y()

        self.assertEqual(9, x.size,
                         msg='Correct #rows for table (X)')
        self.assertEqual(9, y.size,
                         msg='Correct #rows for table (Y)')

        # I do not think we guarantee the type of the column, so
        # support both 32 and 64 bit data types. I would expect
        # both columns to have the same data type, but there
        # could be a reason that this is not true, so do not
        # enforce it here.
        self.assertIn(x.dtype, (np.float32, np.float64))
        self.assertIn(y.dtype, (np.float32, np.float64))

    def _test_table3(self, tname):
        """Tests for a three-column file: col1 col2"""
        mdl = ui.get_model_component(tname)

        x = mdl.get_x()
        y = mdl.get_y()

        self.assertEqual(9, x.size,
                         msg='Correct #rows for table (X)')
        self.assertEqual(9, y.size,
                         msg='Correct #rows for table (Y)')

        # I do not think we guarantee the type of the column, so
        # support both 32 and 64 bit data types. I would expect
        # both columns to have the same data type, but there
        # could be a reason that this is not true, so do not
        # enforce it here.
        self.assertIn(x.dtype, (np.float32, np.float64))
        self.assertIn(y.dtype, (np.float32, np.float64))

    def test_ascii_table1(self):
        """Read in a one-column ASCII file"""
        ui.load_table_model('tbl1', self.ascii_onecol)
        self._test_table1('tbl1')

        m = ui.get_model_component('tbl1')
        y = m.get_y()
        assert_allclose(self.y, y)

    def test_ascii_table2(self):
        """Read in a two-column ASCII file"""
        ui.load_table_model('tbl2', self.ascii_twocol)
        self._test_table2('tbl2')

        m = ui.get_model_component('tbl2')
        x = m.get_x()
        y = m.get_y()
        assert_allclose(self.x, x)
        assert_allclose(self.y, y)

    def test_ascii_table3_col12(self):
        """Read in a three-column ASCII file: col1 col2"""
        ui.load_table_model('tbl3', self.ascii_threecol)
        self._test_table2('tbl3')

        m = ui.get_model_component('tbl3')
        x = m.get_x()
        y = m.get_y()
        assert_allclose(self.x, x)
        assert_allclose(self.y, y)

    def test_ascii_table3_col13(self):
        """Read in a three-column ASCII file: col1 col3"""
        ui.load_table_model('tbl3', self.ascii_threecol,
                            colkeys=['X', 'STAT_ERR'])
        self._test_table2('tbl3')

        m = ui.get_model_component('tbl3')
        x = m.get_x()
        y = m.get_y()
        assert_allclose(self.x, x)
        assert_allclose(self.dy, y)

    def test_ascii_table3_col23(self):
        """Read in a three-column ASCII file: col2 col3"""
        ui.load_table_model('tbl3', self.ascii_threecol,
                            colkeys=['Y', 'STAT_ERR'])
        self._test_table2('tbl3')

        # Note: the values are sorted on read
        m = ui.get_model_component('tbl3')
        x = m.get_x()
        y = m.get_y()

        idx = np.argsort(self.y)
        assert_allclose(self.y[idx], x)
        assert_allclose(self.dy[idx], y)

    def test_ascii_table3_ncols1(self):
        """Read in a three-column ASCII file: ncols=1"""
        ui.load_table_model('tbl1', self.ascii_threecol, ncols=1)
        self._test_table1('tbl1')

        m = ui.get_model_component('tbl1')
        y = m.get_y()
        assert_allclose(self.x, y)

    def test_ascii_table3_ncols2(self):
        """Read in a three-column ASCII file: ncols=2"""
        ui.load_table_model('tbl3', self.ascii_threecol, ncols=2)
        self._test_table2('tbl3')

        m = ui.get_model_component('tbl3')
        x = m.get_x()
        y = m.get_y()
        assert_allclose(self.x, x)
        assert_allclose(self.y, y)

    # Is it worth checking that FITS files are not supported (it
    # does trigger an error condition so improves code coverage)?
    #
    @requires_data
    def test_fits_errors_out(self):
        """FITS files are unsupported"""
        fname = self.make_path('double.fits')
        self.assertRaises(IOErr, ui.load_table_model, 'ftbl', fname)

    def test_eval_basic1_fail(self):
        """Model evaluation fails if #rows does not match."""

        ui.load_table_model('fmodel', self.ascii_onecol)
        ui.load_arrays(99, self.x[1:], self.y[1:])
        ui.set_source(99, 'fmodel')
        self.assertRaises(ModelErr, ui.calc_stat, 99)

    def test_eval_basic1(self):
        """Check the one-column file evaluates sensibly."""

        ui.set_stat('leastsq')
        nval = 2.1e6
        for interp in self.interp1d:
            ui.load_table_model('tbl1', self.ascii_onecol,
                                method=interp)
            m = ui.get_model_component('tbl1')
            ui.load_arrays('tbl', self.x, self.y * nval, ui.Data1D)
            ui.set_source('tbl', 'tbl1')
            m.ampl = nval

            sval = ui.calc_stat('tbl')
            self.assertAlmostEqual(0.0, sval)

    def test_eval_basic2(self):
        """Check the two-column file evaluates sensibly.

        Here the X axis values match the model, so
        theoretically no interpolation is needed. Try out
        the interpolation models.
        """

        ui.set_stat('leastsq')
        nval = 2.1e6
        for interp in self.interp1d:
            ui.load_table_model('tbl2', self.ascii_twocol,
                                method=interp)
            m = ui.get_model_component('tbl2')
            ui.load_arrays('tbl', self.x, self.y * nval, ui.Data1D)
            ui.set_source('tbl', 'tbl2')
            m.ampl = nval

            sval = ui.calc_stat('tbl')
            self.assertAlmostEqual(0.0, sval)

    def test_eval_interp2(self):
        """Check the two-column file evaluates sensibly.

        Differences to test_eval_basic2 include:

        * it uses different X-axis values
        * model evaluation is explicit, rather than
          implicit
        * uses hard-coded values for the expected values
        """

        # Pick some values outside the data range; the expected values
        # are based on a visual check of the interpolation results,
        # and so are a regression test. The first/last points of
        # nevile are exclued from the final comparison, since the
        # tolerance on the comparison on these values is not really
        # well defined.
        #
        xvals = [75, 100, 130, 150, 190, 210]
        yexp = {}
        yexp[linear_interp] = [0, 3, 54.3, 95, 10, -10]
        yexp[nearest_interp] = [0, 0, 35, 96, 15, 0]
        yexp[neville] = [165.55, 9.5, 55.3, 103.6, 5.3, 199.6]

        for interp in self.interp1d:
            ui.load_table_model('tbl2', self.ascii_twocol,
                                method=interp)
            m = ui.get_model_component('tbl2')
            yint = m(xvals)

            ydiff = yint - yexp[interp]
            if interp == neville:
                ydiff = ydiff[1:-1]

            ymax = np.abs(ydiff).max()
            self.assertLess(ymax, 0.1)

    def test_eval_filter2(self):
        """Can we filter the data?

        Based on test_eval_interp2:

        * uses Sherpa data set

        * adds in a filter

        """

        # 120 and 140 will be filtered out
        xvals = [75, 100, 120, 130, 140, 150, 190, 210]
        yexp = {}
        yexp[linear_interp] = [0, 3, -99, 54.3, -99, 95, 10, -10]
        yexp[nearest_interp] = [0, 0, -99, 35, -99, 96, 15, 0]
        yexp[neville] = [165.55, 9.5, -99, 55.3, -99, 103.6, 5.3, 199.6]

        # These tolerances have been chosen to allow the tests to pass
        # on a single machine; they may need to be adjusted for other
        # machines.
        smax = {linear_interp: 0.1, nearest_interp: 0.01, neville: 1.0}
        ui.set_stat('leastsq')
        for interp in self.interp1d:
            ui.load_arrays('filt', xvals, yexp[interp])

            ui.load_table_model('filt2', self.ascii_twocol,
                                method=interp)
            ui.set_source('filt', 'filt2')

            ui.ignore_id('filt', 110, 125)
            ui.ignore_id('filt', 135, 149)

            # Use calc_stat_info in order to check dof
            svals = ui.get_stat_info()
            sval = svals[0]
            self.assertEqual(5, sval.dof)
            self.assertLess(sval.statval, smax[interp])

    def test_eval_integrated(self):
        """Test handling with integrated data sets.

        This is just to check that the current documentation,
        which says that the table x-axis is the left edge of
        the bin.
        """

        # make sure the grid is not evenly spaced, since this
        # should show up any issues with trying to integrate
        # the model over the grid.
        #
        xgrid = np.asarray([10, 17, 20, 22, 31, 35, 40])
        xlo = xgrid[:-1]
        xhi = xgrid[1:]

        ytbl = np.asarray([-10, 10, 15, 20, 20, 10])
        ymdl = ytbl / 2.0

        # Use low-level routines to create the table model,
        # to avoid having a new data set. This side-steps any
        # of the code to handle filtering, but that is not
        # relevant here.
        #
        def make_tm(method):
            tmdl = TableModel('tmp')
            tmdl.filename = 'tmp'
            tmdl.method = method
            tmdl.load(xlo, ytbl)
            return tmdl

        for interp in self.interp1d:
            # Just to be sure that there is no cache issues,
            # re-create the table each iteration
            tmdl = make_tm(interp)

            tmdl.ampl = 0.5
            ycalc = tmdl(xlo, xhi)
            assert_allclose(ymdl, ycalc)


class test_psf_ui(SherpaTestCase):

    models1d = ['gauss1d', 'delta1d', 'normgauss1d']
    models2d = ['gauss2d', 'delta2d', 'normgauss2d']

    # Commented out setUp/tearDown as they do nothing
    #
    # def setUp(self):
    #     pass

    # def tearDown(self):
    #     pass

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
