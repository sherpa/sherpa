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


import unittest
from sherpa.utils import SherpaTest, SherpaTestCase, test_data_missing
from sherpa.models import ArithmeticModel, Parameter
from sherpa import ui
import numpy


class UserModel(ArithmeticModel):

    def __init__(self, name='usermodel'):
        self.param1 = Parameter(name, 'param1', 1, min=0, max=100)
        self.param2 = Parameter(name, 'param2', 1, min=-100, max=100)

        ArithmeticModel.__init__(self, name, (self.param1,
                                              self.param2))

    def calc(self, p, x, *args, **kwargs):
        return p[0]*x+p[1]


class test_ui(SherpaTestCase):

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def setUp(self):
        self.ascii = self.make_path('threads/ascii_table/sim.poisson.1.dat')
        self.single = self.make_path('single.dat')
        self.double = self.make_path('double.dat')
        self.filter = self.make_path('filter_single_integer.dat')
        self.func = lambda x: x

        ui.dataspace1d(1,1000,dstype=ui.Data1D)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_ascii(self):
        ui.load_data(1, self.ascii)
        ui.load_data(1, self.ascii, 2)
        ui.load_data(1, self.ascii, 2, ("col2", "col1"))


    # Test table model
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_table_model_ascii_table(self):
        ui.load_table_model('tbl', self.single)
        ui.load_table_model('tbl', self.double)


    # Test user model
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_user_model_ascii_table(self):
        ui.load_user_model(self.func, 'mdl', self.single)
        ui.load_user_model(self.func, 'mdl', self.double)


    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_filter_ascii(self):
        ui.load_filter(self.filter)
        ui.load_filter(self.filter, ignore=True)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_add_model(self):
        ui.add_model(UserModel)
        ui.set_model('usermodel.user1')

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_set_full_model(self):
        ui.load_psf('psf1', 'gauss2d.g1')
        ui.set_full_model('psf1(gauss2d.g2)+const2d.c1')
        ui.get_model()
#        ui.get_source()

    # Bug 12644
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_source_methods_with_full_model(self):
        from sherpa.utils.err import IdentifierErr

        ui.load_data('full', self.ascii)
        ui.set_full_model('full', 'powlaw1d.p1')

        # Test Case 1
        try:
            ui.get_source('full')
        except IdentifierErr as e:
            self.assertRegexpMatches(str(e), "Convolved model\n.*\n is set for dataset full. You should use get_model instead.", str(e))
        try:
            ui.plot_source('full')
        except IdentifierErr as e:
            self.assertEquals("Convolved model\n'p1'\n is set for dataset full. You should use plot_model instead.", str(e))

        # Test Case 2
        ui.set_source('full', 'powlaw1d.p2')
        ui.get_source('full')

        # Test Case 3
        ui.load_data('not_full', self.ascii)
        try:
            ui.get_source('not_full')
        except IdentifierErr as e:
            self.assertEquals('source not_full has not been set, consider using set_source() or set_model()', str(e))


class test_psf_ui(SherpaTestCase):

    models1d = ['gauss1d', 'delta1d', 'normgauss1d']
    models2d = ['gauss2d', 'delta2d', 'normgauss2d']

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_psf_model1d(self):
        ui.dataspace1d(1, 10)
        for model in self.models1d:
            try:
                ui.load_psf('psf1d', model+'.mdl')
                ui.set_psf('psf1d')
                mdl = ui.get_model_component('mdl')
                self.assertTrue((numpy.array(mdl.get_center()) ==
                                 numpy.array([4])).all())
            except:
                print model
                raise


    def test_psf_model2d(self):
        ui.dataspace2d([216,261])
        for model in self.models2d:
            try:
                ui.load_psf('psf2d', model+'.mdl')
                ui.set_psf('psf2d')
                mdl = ui.get_model_component('mdl')
                self.assertTrue((numpy.array(mdl.get_center()) ==
                                 numpy.array([108,130])).all())
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
