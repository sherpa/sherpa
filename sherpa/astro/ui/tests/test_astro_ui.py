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

import os
import unittest
from sherpa.utils import SherpaTest, SherpaTestCase
from sherpa.utils import test_data_missing, has_fits_support
from sherpa.astro import ui
import numpy
import logging
logger = logging.getLogger("sherpa")


class test_ui(SherpaTestCase):

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def setUp(self):
        self.ascii = self.make_path('threads/ascii_table/sim.poisson.1.dat')
        self.fits = self.make_path('1838_rprofile_rmid.fits')
        self.singledat = self.make_path('single.dat')
        self.singletbl = self.make_path('single.fits')
        self.doubledat = self.make_path('double.dat')
        self.doubletbl = self.make_path('double.fits')
        self.img = self.make_path('img.fits')
        self.filter_single_int_ascii = self.make_path('filter_single_integer.dat')
        self.filter_single_int_table = self.make_path('filter_single_integer.fits')
        self.filter_single_log_table = self.make_path('filter_single_logical.fits')

        self.func = lambda x: x
        ui.dataspace1d(1,1000,dstype=ui.Data1D)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_ascii(self):
        ui.load_ascii(1, self.ascii)
        ui.load_ascii(1, self.ascii, 2)
        ui.load_ascii(1, self.ascii, 2, ("col2", "col1"))

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_table(self):
        ui.load_table(1, self.fits)
        ui.load_table(1, self.fits, 3)
        ui.load_table(1, self.fits, 3, ["RMID","SUR_BRI","SUR_BRI_ERR"])
        ui.load_table(1, self.fits, 4, ('R',"SUR_BRI",'SUR_BRI_ERR'),
                      ui.Data1DInt)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    def test_load_table_fits(self):
        # QUS: why is this not in the sherpa-test-data repository?
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ui.load_table(1, os.path.join(this_dir, 'data', 'two_column_x_y.fits.gz'))
        data = ui.get_data(1)
        self.assertEqualWithinTol(data.x, [1, 2, 3])
        self.assertEqualWithinTol(data.y, [4, 5, 6])

    # Test table model
    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_table_model_ascii_table(self):
        ui.load_table_model('tbl', self.singledat)
        ui.load_table_model('tbl', self.doubledat)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_table_model_fits_table(self):
        ui.load_table_model('tbl', self.singletbl)
        ui.load_table_model('tbl', self.doubletbl)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_table_model_fits_image(self):
        ui.load_table_model('tbl', self.img)

    # Test user model
    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_user_model_ascii_table(self):
        ui.load_user_model(self.func, 'mdl', self.singledat)
        ui.load_user_model(self.func, 'mdl', self.doubledat)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_user_model_fits_table(self):
        ui.load_user_model(self.func, 'mdl', self.singletbl)
        ui.load_user_model(self.func, 'mdl', self.doubletbl)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_filter_ascii(self):
        ui.load_filter(self.filter_single_int_ascii)
        ui.load_filter(self.filter_single_int_ascii, ignore=True)

    # Test load_filter
    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_filter_table(self):
        ui.load_filter(self.filter_single_int_table)
        ui.load_filter(self.filter_single_int_table, ignore=True)

        ui.load_filter(self.filter_single_log_table)
        ui.load_filter(self.filter_single_log_table, ignore=True)


class test_more_ui(SherpaTestCase):
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        self.img = self.make_path('img.fits')
        self.pha = self.make_path('threads/simultaneous/pi2286.fits')
        self.rmf = self.make_path('threads/simultaneous/rmf2286.fits')
        self.pha3c273 = self.make_path('ciao4.3/pha_intro/3c273.pi'

    def tearDown(self):
        logger.setLevel(self._old_logger_level)

    # bug #12732
    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_string_model_with_rmf(self):
        ui.load_pha("foo", self.pha)
        ui.load_rmf("foo", self.rmf)
        # Check that get_rmf(id)('modelexpression') works
        caught = False
        try:
            m = ui.get_rmf("foo")("powlaw1d.pl1")
        except:
            caught = True
        if caught:
            self.fail("Exception caught when it shouldn't")
        from sherpa.astro.instrument import RMFModelPHA
        self.assertTrue(isinstance(m, RMFModelPHA))

    #bug #38
    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_bug38(self):
        ui.load_pha('3c273', self.pha3c273)
        ui.notice_id('3c273', 0.3, 2)
        ui.group_counts('3c273', 30)
        ui.group_counts('3c273', 15)

class test_image_12578(SherpaTestCase):
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def setUp(self):
        self.img = self.make_path('img.fits')
        self.loggingLevel = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        ui.clean()

    def tearDown(self):
        if hasattr(self, 'loggingLevel'):
            logger.setLevel(self.loggingLevel)

    # bug #12578
    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_set_coord_bad_coord(self):
        from sherpa.utils.err import IdentifierErr, DataErr

        # Test Case #1: if the list of ids is empty, raise
        # IdentifierErr['nodatasets']
        caught = False
        try:
            ui.set_coord('image')
        except IdentifierErr:
            caught = True
        if not caught:
            self.fail("Test Case #1: IdentifierErr Exception not caught")

        # Test Case #2: check the user expected behavior. The call
        # set_coord("sky") will result in the error message
        # DataErr: unknown coordinates: 'sky'\n \
        # Valid coordinates: logical, image, physical, world, wcs
        ui.load_image(self.img)

        caught = False
        try:
            ui.set_coord("sky")
        except DataErr as e:
            okmsg = "unknown coordinates: 'sky'\nValid options: logical, image, physical, world, wcs"
            self.assertEqual(okmsg, e.message)
            caught = True
        if not caught:
            self.fail("Test Case #2: DataErr Exception not caught")


class test_psf_ui(SherpaTestCase):

    models1d = ['beta1d', 'lorentz1d', 'normbeta1d']
    models2d = ['beta2d', 'devaucouleurs2d', 'hubblereynolds', 'lorentz2d']

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @unittest.skip("TODO: failing test used to have a different name and " +
                   "was never executed")
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

    # bug #12503
    def test_psf_pars_are_frozen(self):
        ui.load_psf('psf', ui.beta2d.p1)
        self.assertEqual([], p1.thawedpars)


class test_stats_ui(SherpaTestCase):

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def setUp(self):
        self.data = self.make_path('threads/chi2/3c273.pi')
        ui.clean()

    # bugs #11400, #13297, #12365
    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits or astropy.io.fits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_chi2(self):

        # Case 1: first ds has no error, second has, chi2-derived (chi2gehrels)
        # statistic. I expect stat.name to be chi2gehrels for ds1, chi2 for
        # ds2, chi2gehrels for ds1,2
        ui.load_data(1, self.data)
        ui.load_data(2, self.data, use_errors=True)

        ui.set_source(1, "gauss1d.g1")
        ui.set_source(2, "gauss1d.g1")

        ui.set_stat("chi2gehrels")

        si = ui.get_stat_info()

        stat1 = si[0].statname
        stat2 = si[1].statname
        stat12 = si[2].statname

        self.assertEqual('chi2gehrels', stat1)
        self.assertEqual('chi2', stat2)
        self.assertEqual('chi2gehrels', stat12)

        # Case 2: first ds has errors, second has not, chi2-derived
        # (chi2gehrels) statistic. I expect stat.name to be chi2 for ds1,
        # chi2gehrels for ds2, chi2gehrels for ds1,2
        ui.load_data(2, self.data)
        ui.load_data(1, self.data, use_errors=True)

        si = ui.get_stat_info()

        stat1 = si[0].statname
        stat2 = si[1].statname
        stat12 = si[2].statname

        self.assertEqual('chi2gehrels', stat2)
        self.assertEqual('chi2', stat1)
        self.assertEqual('chi2gehrels', stat12)

        # Case 3: both datasets have errors, chi2-derived (chi2gehrels)
        # statistic. I expect stat.name to be chi2 for all of them.
        ui.load_data(2, self.data, use_errors=True)
        ui.load_data(1, self.data, use_errors=True)

        si = ui.get_stat_info()

        stat1 = si[0].statname
        stat2 = si[1].statname
        stat12 = si[2].statname

        self.assertEqual('chi2', stat2)
        self.assertEqual('chi2', stat1)
        self.assertEqual('chi2', stat12)

        # Case 4: first ds has errors, second has not, LeastSq statistic
        # I expect stat.name to be leastsq for all of them.
        ui.load_data(2, self.data)
        ui.load_data(1, self.data, use_errors=True)

        ui.set_stat("leastsq")

        si = ui.get_stat_info()

        stat1 = si[0].statname
        stat2 = si[1].statname
        stat12 = si[2].statname

        self.assertEqual('leastsq', stat2)
        self.assertEqual('leastsq', stat1)
        self.assertEqual('leastsq', stat12)

        # Case 5: both ds have errors, LeastSq statistic
        # I expect stat.name to be leastsq for all of them.
        ui.load_data(2, self.data, use_errors=True)
        ui.load_data(1, self.data, use_errors=True)

        ui.set_stat("leastsq")

        si = ui.get_stat_info()

        stat1 = si[0].statname
        stat2 = si[1].statname
        stat12 = si[2].statname

        self.assertEqual('leastsq', stat2)
        self.assertEqual('leastsq', stat1)
        self.assertEqual('leastsq', stat12)

        # Case 6: first ds has errors, second has not, CStat statistic
        # I expect stat.name to be cstat for all of them.
        ui.load_data(2, self.data)
        ui.load_data(1, self.data, use_errors=True)

        ui.set_stat("cstat")

        si = ui.get_stat_info()

        stat1 = si[0].statname
        stat2 = si[1].statname
        stat12 = si[2].statname

        self.assertEqual('cstat', stat2)
        self.assertEqual('cstat', stat1)
        self.assertEqual('cstat', stat12)

        # Case7: select chi2 as statistic. One of the ds does not provide
        # errors. I expect sherpa to raise a StatErr exception.
        ui.set_stat('chi2')

        caught = False

        from sherpa.utils.err import StatErr
        try:
            ui.get_stat_info()
        except StatErr:
            caught = True

        self.assertTrue(caught, msg='StatErr was not caught')

        # Case8: select chi2 as statistic. Both datasets provide errors
        # I expect stat to be 'chi2'
        ui.load_data(2, self.data, use_errors=True)
        si = ui.get_stat_info()

        stat1 = si[0].statname
        stat2 = si[1].statname
        stat12 = si[2].statname

        self.assertEqual('chi2', stat2)
        self.assertEqual('chi2', stat1)
        self.assertEqual('chi2', stat12)


if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        datadir = sys.argv[1]
    else:
        datadir = None

    SherpaTest(ui).test(datadir=datadir)
