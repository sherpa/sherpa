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

from sherpa.utils import SherpaTest, SherpaTestCase, requires_data
from sherpa.utils.err import IdentifierErr, IOErr
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


@requires_data
class test_table_model(SherpaTestCase):

    def setUp(self):
        self.ascii_onecol = self.make_path('single.dat')
        self.ascii_twocol = self.make_path('double.dat')
        # note: There is no FITS support in sherpa.ui.load_table_model,
        #       it is only available in sherpa.astro.ui.load_table_model
        self.fits_onecol = self.make_path('single.fits')
        self.fits_twocol = self.make_path('double.fits')

        # TODO: add tests to check that the model evaluation/interpolation
        #       works
        # ui.dataspace1d(1, 1000, dstype=ui.Data1D)

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

    def _test_table1(self, tname):
        """Tests for a one-column file"""

        mdl = ui.get_model_component(tname)
        x = mdl.get_x()
        y = mdl.get_y()

        self.assertEqual(None, x, msg='No X axis for table')

        y = mdl.get_y()
        self.assertEqual(1000, y.size,
                         msg='Correct #rows for table')
        self.assertEqual(1000, y.size,
                         msg='Correct #rows for table')

        # I do not think we guarantee the type of the column, so
        # support both 32 and 64 bit data types
        self.assertIn(y.dtype, (np.float32, np.float64))

    def _test_table2(self, tname):
        """Tests for a two-column file"""
        mdl = ui.get_model_component(tname)

        x = mdl.get_x()
        y = mdl.get_y()

        self.assertEqual(1000, x.size,
                         msg='Correct #rows for table (X)')
        self.assertEqual(1000, y.size,
                         msg='Correct #rows for table (Y)')

        # I do not think we guarantee the type of the column, so
        # support both 32 and 64 bit data types. I would expect
        # both columns to have the same data type, but there
        # could be a reason that this is not true, so do not
        # enforce it here.
        self.assertIn(x.dtype, (np.float32, np.float64))
        self.assertIn(y.dtype, (np.float32, np.float64))

    def _test_yvals(self, tname1, tname2):
        """Check that the y values of the two components are similar."""

        y1 = ui.get_model_component(tname1).get_y()
        y2 = ui.get_model_component(tname2).get_y()

        # Could probably use assert_array_equal here, but use a
        # tolerance, since these are floating-point values
        assert_allclose(y1, y2, err_msg='y values for table model')

    def test_ascii_table1(self):
        """Read in a one-column ASCII file"""
        ui.load_table_model('tbl1', self.ascii_onecol)
        self._test_table1('tbl1')

    def test_ascii_table2(self):
        """Read in a two-column ASCII file"""
        ui.load_table_model('tbl2', self.ascii_twocol)
        self._test_table2('tbl2')

    def test_ascii_yvalue(self):
        """Check the Y values are as expected."""

        # this is just a consistency test, and relies on
        # the two files having the same Y values
        ui.load_table_model('tbl1', self.ascii_onecol)
        ui.load_table_model('tbl2', self.ascii_twocol)
        self._test_yvals('tbl1', 'tbl2')

    def test_fits_table1(self):
        """FITS files are unsupported"""
        self.assertRaises(IOErr, ui.load_table_model, 'ftbl1',
                          self.fits_onecol)

    def test_fits_table2(self):
        """FITS files are unsupported"""
        self.assertRaises(IOErr, ui.load_table_model, 'ftbl2',
                          self.fits_twocol)


class test_psf_ui(SherpaTestCase):

    models1d = ['gauss1d', 'delta1d', 'normgauss1d']
    models2d = ['gauss2d', 'delta2d', 'normgauss2d']

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

    import sys
    if len(sys.argv) > 1:
        datadir = sys.argv[1]
    else:
        datadir = None

    SherpaTest(ui).test(datadir=datadir)
