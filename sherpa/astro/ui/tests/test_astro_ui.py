#
#  Copyright (C) 2012, 2015, 2016, 2017, 2018, 2019, 2020
#     Smithsonian Astrophysical Observatory
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
import re
import unittest
import tempfile

import pytest

import numpy
from numpy.testing import assert_allclose

from sherpa.utils.testing import SherpaTestCase, requires_data, \
    requires_fits, requires_group, requires_xspec
from sherpa.astro import ui
from sherpa.data import Data1D
from sherpa.astro.data import DataPHA

import logging
logger = logging.getLogger("sherpa")


@requires_fits
@requires_data
class test_ui(SherpaTestCase):

    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        self.ascii = self.make_path('sim.poisson.1.dat')
        self.fits = self.make_path('1838_rprofile_rmid.fits')
        self.singledat = self.make_path('single.dat')
        self.singletbl = self.make_path('single.fits')
        self.doubledat = self.make_path('double.dat')
        self.doubletbl = self.make_path('double.fits')
        self.img = self.make_path('img.fits')
        self.filter_single_int_ascii = self.make_path(
            'filter_single_integer.dat')
        self.filter_single_int_table = self.make_path(
            'filter_single_integer.fits')
        self.filter_single_log_table = self.make_path(
            'filter_single_logical.fits')

        self.func = lambda x: x
        ui.dataspace1d(1, 1000, dstype=ui.Data1D)

    def tearDown(self):
        if hasattr(self, '_old_logger_level'):
            logger.setLevel(self._old_logger_level)

    def test_ascii(self):
        ui.load_ascii(1, self.ascii)
        ui.load_ascii(1, self.ascii, 2)
        ui.load_ascii(1, self.ascii, 2, ("col2", "col1"))

    def test_table(self):
        ui.load_table(1, self.fits)
        ui.load_table(1, self.fits, 3)
        ui.load_table(1, self.fits, 3, ["RMID", "SUR_BRI", "SUR_BRI_ERR"])
        ui.load_table(1, self.fits, 4, ('R', "SUR_BRI", 'SUR_BRI_ERR'),
                      ui.Data1DInt)

    # Test table model
    def test_table_model_ascii_table(self):
        ui.load_table_model('tbl', self.singledat)
        ui.load_table_model('tbl', self.doubledat)

    def test_table_model_fits_table(self):
        ui.load_table_model('tbl', self.singletbl)
        ui.load_table_model('tbl', self.doubletbl)

    def test_table_model_fits_image(self):
        ui.load_table_model('tbl', self.img)

    # Test user model
    def test_user_model_ascii_table(self):
        ui.load_user_model(self.func, 'mdl', self.singledat)
        ui.load_user_model(self.func, 'mdl', self.doubledat)

    def test_user_model_fits_table(self):
        ui.load_user_model(self.func, 'mdl', self.singletbl)
        ui.load_user_model(self.func, 'mdl', self.doubletbl)

    def test_filter_ascii(self):
        ui.load_filter(self.filter_single_int_ascii)
        ui.load_filter(self.filter_single_int_ascii, ignore=True)

    # Test load_filter
    def test_filter_table(self):
        ui.load_filter(self.filter_single_int_table)
        ui.load_filter(self.filter_single_int_table, ignore=True)

        ui.load_filter(self.filter_single_log_table)
        ui.load_filter(self.filter_single_log_table, ignore=True)


@requires_data
@requires_fits
class test_more_ui(SherpaTestCase):

    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        self.img = self.make_path('img.fits')
        self.pha = self.make_path('pi2286.fits')
        self.rmf = self.make_path('rmf2286.fits')
        self.pha3c273 = self.make_path('3c273.pi')

    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    # bug #12732
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
    @requires_group
    def test_bug38(self):
        ui.load_pha('3c273', self.pha3c273)
        ui.notice_id('3c273', 0.3, 2)
        ui.group_counts('3c273', 30)
        ui.group_counts('3c273', 15)


@requires_data
@requires_fits
class test_image_12578(SherpaTestCase):
    def setUp(self):
        self.img = self.make_path('img.fits')
        self.loggingLevel = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        ui.clean()

    def tearDown(self):
        if hasattr(self, 'loggingLevel'):
            logger.setLevel(self.loggingLevel)

    # bug #12578
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
            okmsg = "unknown coordinates: 'sky'\nValid options: " + \
                    "logical, image, physical, world, wcs"
            self.assertEqual(okmsg, str(e))
            caught = True
        if not caught:
            self.fail("Test Case #2: DataErr Exception not caught")


class test_psf_ui(SherpaTestCase):

    models1d = ['beta1d', 'lorentz1d', 'normbeta1d']
    models2d = ['beta2d', 'devaucouleurs2d', 'hubblereynolds', 'lorentz2d']

    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    @unittest.skip("TODO: failing test used to have a different name and " +
                   "was never executed")
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

    # bug #12503
    def test_psf_pars_are_frozen(self):
        ui.create_model_component("beta2d", "p1")
        p1 = ui.get_model_component("p1")
        ui.load_psf('psf', p1)
        self.assertEqual([], p1.thawedpars)


@requires_data
@requires_fits
class test_stats_ui(SherpaTestCase):
    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        self.data = self.make_path('3c273.pi')
        ui.clean()

    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    # bugs #11400, #13297, #12365
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


@requires_fits
class test_save_arrays_base(SherpaTestCase):

    # _colnames:  specify column names (True) or not
    # _fits:      use FITS format for output?

    _colnames = False
    _fits = None

    def _read_func(self, filename):
        raise NotImplementedError

    def save_arrays(self):
        """Write out a small set of data using ui.save_arrays
        and then read it back in, to check it was written out
        correctly (or, at least, in a way that can be read back
        in).
        """

        # It looks like the input arrays to `write_arrays` should be numpy
        # arrays, so enforce that invariant.
        a = numpy.asarray([1, 3, 9])
        b = numpy.sqrt(numpy.asarray(a))
        c = b * 0.1

        if self._colnames:
            fields = ["x", "yy", "z"]
        else:
            fields = None

        ofh = tempfile.NamedTemporaryFile(suffix='sherpa_test')
        ui.save_arrays(ofh.name, [a, b, c], fields=fields,
                       ascii=not self._fits, clobber=True)

        out = self._read_func(ofh.name)

        rtol = 0
        atol = 1e-5
        self.assertIsInstance(out, Data1D)

        # remove potential dm syntax introduced by backend before checking for equality
        out_name = re.sub(r"\[.*\]", "", out.name)

        self.assertEqual(out_name, ofh.name, msg="file name")
        assert_allclose(out.x, a, rtol=rtol, atol=atol, err_msg="x column")
        assert_allclose(out.y, b, rtol=rtol, atol=atol, err_msg="y column")
        assert_allclose(out.staterror, c, rtol=rtol, atol=atol,
                        err_msg="staterror")
        self.assertIsNone(out.syserror, msg="syserror")


@requires_fits
class test_save_arrays_nocols_FITS(test_save_arrays_base):

    _fits = True

    def _read_func(self, filename):
        return ui.unpack_table(filename, ncols=3)

    def test_save_arrays(self):
        self.save_arrays()


@requires_fits
class test_save_arrays_cols_FITS(test_save_arrays_nocols_FITS):

    _colnames = True

    def test_save_arrays(self):
        self.save_arrays()


@requires_fits
class test_save_arrays_nocols_ASCII(test_save_arrays_base):

    _fits = False

    def _read_func(self, filename):
        return ui.unpack_ascii(filename, ncols=3)

    def test_save_arrays(self):
        self.save_arrays()


@requires_fits
class test_save_arrays_cols_ASCII(test_save_arrays_nocols_ASCII):

    _colnames = True

    def test_save_arrays(self):
        self.save_arrays()


@requires_fits
@requires_data
class test_save_pha(SherpaTestCase):
    """Write out a PHA data set as a FITS file."""

    longMessage = True

    def setUp(self):
        # hide warning messages from file I/O
        self._old_logger_level = logger.level
        logger.setLevel(logging.ERROR)

        self._id = 1
        fname = self.make_path('3c273.pi')
        ui.load_pha(self._id, fname)
        self._pha = ui.get_data(self._id)

    def tearDown(self):
        logger.setLevel(self._old_logger_level)

    def testWrite(self):
        ofh = tempfile.NamedTemporaryFile(suffix='sherpa_test')
        ui.save_pha(self._id, ofh.name, ascii=False, clobber=True)

        # limited checks
        pha = ui.unpack_pha(ofh.name)
        self.assertIsInstance(pha, DataPHA)

        for key in ["channel", "counts"]:
            newval = getattr(pha, key)
            oldval = getattr(self._pha, key)
            assert_allclose(oldval, newval, err_msg=key)

        # at present grouping and quality are not written out

        for key in ["exposure", "backscal", "areascal"]:
            newval = getattr(pha, key)
            oldval = getattr(self._pha, key)
            self.assertAlmostEqual(oldval, newval, msg=key)

        """
        Since the original file has RMF and ARF, the units
        are energy, but on write out these files are not
        created/saved, so when read back in the PHA has no
        ARF/RMF, and will have units=channel.

        for key in ["units"]:
            newval = getattr(pha, key)
            oldval = getattr(self._pha, key)
            self.assertEqual(oldval, newval, msg=key)
        """


@requires_fits
class test_basic_io(SherpaTestCase):
    def setUp(self):
        ui.clean()

    def test_load_table_fits(self):
        # QUS: why is this not in the sherpa-test-data repository?
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ui.load_table(1, os.path.join(this_dir, 'data',
                                      'two_column_x_y.fits.gz'))
        data = ui.get_data(1)
        self.assertEqualWithinTol(data.x, [1, 2, 3])
        self.assertEqualWithinTol(data.y, [4, 5, 6])


@requires_data
@requires_fits
def test_bug_276(make_data_path):
    ui.load_pha(make_data_path('3c273.pi'))
    ui.set_model('polynom1d.p1')
    ui.fit()
    ui.covar()
    scal = ui.get_covar_results().parmaxes
    ui.sample_flux(ui.get_model_component('p1'), 0.5, 1, num=5, correlated=False, scales=scal)


@requires_data
@requires_fits
@pytest.mark.usefixtures("clean_astro_ui")
def test_bug_316(make_data_path):
    """get_source_plot does not apply lo/hi directly

    This does not need a plot backend, as no plot is created,
    just the data needed to create the plot.
    """

    infile = make_data_path('3c273.pi')
    ui.load_pha('xs', infile)
    ui.set_source('xs', ui.powlaw1d.pl)

    splot = ui.get_source_plot('xs')
    xmin = splot.xlo[0]
    xmax = splot.xhi[-1]
    nbins = len(splot.xlo)

    assert xmin == pytest.approx(0.1)
    assert xmax == pytest.approx(11.0)
    assert nbins == 1090

    assert splot.mask is None

    # applying arguments doesn't filter the structure
    splot = ui.get_source_plot('xs', lo=2, hi=7)
    xmin = splot.xlo[0]
    xmax = splot.xhi[-1]
    nbins = len(splot.xlo)

    assert xmin == pytest.approx(0.1)
    assert xmax == pytest.approx(11.0)
    assert nbins == 1090

    assert splot.mask is not None
    assert splot.mask.sum() == 501

    xmin = splot.xlo[splot.mask][0]
    xmax = splot.xhi[splot.mask][-1]

    assert xmin == pytest.approx(2.0)
    assert xmax == pytest.approx(7.01)


@requires_data
@requires_fits
def test_load_multi_arfsrmfs(make_data_path, clean_astro_ui):
    """Added in #728 to ensure cache parameter is sent along by
    MultiResponseSumModel.

    This has since been simplified to switch from xsapec to
    powlaw1d as it drops the need for XSPEC and is a simpler
    model, so is less affected by changes in the model code.

    A fit of the Sherpa powerlaw-model to 3c273.pi with a
    single response in CIAO 4.11 (background subtracted,
    0.5-7 keV) returns gamma = 1.9298, ampl = 1.73862e-4
    so doubling the response should halve the amplitude but
    leave the gamma value the same when using two responses,
    as below. This is with chi2datavar.
    """

    pha_pi = make_data_path("3c273.pi")
    ui.load_pha(1, pha_pi)
    ui.load_pha(2, pha_pi)

    arf = make_data_path("3c273.arf")
    rmf = make_data_path("3c273.rmf")

    ui.load_multi_arfs(1, [arf, arf], [1, 2])
    ui.load_multi_arfs(2, [arf, arf], [1, 2])

    ui.load_multi_rmfs(1, [rmf, rmf], [1, 2])
    ui.load_multi_rmfs(2, [rmf, rmf], [1, 2])

    ui.notice(0.5, 7)
    ui.subtract(1)
    ui.subtract(2)

    src = ui.create_model_component('powlaw1d', 'src')
    ui.set_model(1, src)
    ui.set_model(2, src)

    ui.set_method('levmar')
    ui.set_stat('chi2datavar')

    ui.fit()
    fr = ui.get_fit_results()
    assert fr.succeeded
    assert fr.datasets == (1, 2)

    assert src.gamma.val == pytest.approx(1.9298, rel=1.0e-4)
    assert src.ampl.val == pytest.approx(1.73862e-4 / 2, rel=1.0e-4)
