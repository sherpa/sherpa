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

from collections import namedtuple
import os
import re
import unittest
import tempfile
import logging

import pytest

import numpy
from numpy.testing import assert_allclose

from sherpa.utils.testing import requires_data, requires_fits, \
    requires_group
from sherpa.utils.err import DataErr, IdentifierErr, StatErr
from sherpa.astro import ui
from sherpa.astro.instrument import RMFModelPHA
from sherpa.data import Data1D
from sherpa.astro.data import DataPHA

logger = logging.getLogger("sherpa")


@pytest.fixture(autouse=True)
def hide_logging():
    "hide INFO-level logging in all these tests"

    olevel = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield
    logger.setLevel(olevel)


def identity(x):
    return x


@pytest.fixture
def setup_files(make_data_path):

    out = namedtuple('store1', ['ascii', 'fits', 'img'
                                'singledat', 'singletbl',
                                'doubledat', 'doubletbl',
                                'filter_single_int_ascii',
                                'filter_single_int_table',
                                'filter_single_log_table'])

    out.ascii = make_data_path('sim.poisson.1.dat')
    out.fits = make_data_path('1838_rprofile_rmid.fits')
    out.singledat = make_data_path('single.dat')
    out.singletbl = make_data_path('single.fits')
    out.doubledat = make_data_path('double.dat')
    out.doubletbl = make_data_path('double.fits')
    out.img = make_data_path('img.fits')
    out.filter_single_int_ascii = make_data_path(
        'filter_single_integer.dat')
    out.filter_single_int_table = make_data_path(
        'filter_single_integer.fits')
    out.filter_single_log_table = make_data_path(
        'filter_single_logical.fits')

    ui.dataspace1d(1, 1000, dstype=ui.Data1D)

    return out


@requires_fits
@requires_data
def test_ui_ascii(setup_files):
    ui.load_ascii(1, setup_files.ascii)
    ui.load_ascii(1, setup_files.ascii, 2)
    ui.load_ascii(1, setup_files.ascii, 2, ("col2", "col1"))


@requires_fits
@requires_data
def test_ui_table(setup_files):
    ui.load_table(1, setup_files.fits)
    ui.load_table(1, setup_files.fits, 3)
    ui.load_table(1, setup_files.fits, 3, ["RMID", "SUR_BRI", "SUR_BRI_ERR"])
    ui.load_table(1, setup_files.fits, 4, ('R', "SUR_BRI", 'SUR_BRI_ERR'),
                  ui.Data1DInt)

# Test table model
@requires_fits
@requires_data
def test_ui_table_model_ascii_table(setup_files):
    ui.load_table_model('tbl', setup_files.singledat)
    ui.load_table_model('tbl', setup_files.doubledat)


@requires_fits
@requires_data
def test_ui_table_model_fits_table(setup_files):
    ui.load_table_model('tbl', setup_files.singletbl)
    ui.load_table_model('tbl', setup_files.doubletbl)


@requires_fits
@requires_data
def test_ui_table_model_fits_image(setup_files):
    ui.load_table_model('tbl', setup_files.img)


# Test user model
@requires_fits
@requires_data
def test_ui_user_model_ascii_table(setup_files):
    ui.load_user_model(identity, 'mdl', setup_files.singledat)
    ui.load_user_model(identity, 'mdl', setup_files.doubledat)


@requires_fits
@requires_data
def test_ui_user_model_fits_table(setup_files):
    ui.load_user_model(identity, 'mdl', setup_files.singletbl)
    ui.load_user_model(identity, 'mdl', setup_files.doubletbl)


@requires_fits
@requires_data
def test_ui_filter_ascii(setup_files):
    ui.load_filter(setup_files.filter_single_int_ascii)
    ui.load_filter(setup_files.filter_single_int_ascii, ignore=True)


# Test load_filter
@requires_fits
@requires_data
def test_ui_filter_table(setup_files):
    ui.load_filter(setup_files.filter_single_int_table)
    ui.load_filter(setup_files.filter_single_int_table, ignore=True)

    ui.load_filter(setup_files.filter_single_log_table)
    ui.load_filter(setup_files.filter_single_log_table, ignore=True)


# bug #12732
@requires_fits
@requires_data
def test_more_ui_string_model_with_rmf(make_data_path):

    ui.load_pha("foo", make_data_path('pi2286.fits'))
    ui.load_rmf("foo", make_data_path('rmf2286.fits'))

    # Check that get_rmf(id)('modelexpression') works. If it
    # raises an error the test will fail.
    #
    m = ui.get_rmf("foo")("powlaw1d.pl1")
    assert isinstance(m, RMFModelPHA)


#bug #38
@requires_fits
@requires_data
@requires_group
def test_more_ui_bug38(make_data_path):
    ui.load_pha('3c273', make_data_path('3c273.pi'))
    ui.notice_id('3c273', 0.3, 2)
    ui.group_counts('3c273', 30)
    ui.group_counts('3c273', 15)


# bug #12578
@requires_data
@requires_fits
def test_image_12578_set_coord_bad_coord(make_data_path, clean_astro_ui):

    img = make_data_path('img.fits')

    # Test Case #1: if the list of ids is empty, raise
    # IdentifierErr['nodatasets']
    #
    with pytest.raises(IdentifierErr):
        ui.set_coord('image')

    # Test Case #2: check the user expected behavior. The call
    # set_coord("sky") will result in the error message
    # DataErr: unknown coordinates: 'sky'\n \
    # Valid coordinates: logical, image, physical, world, wcs
    #
    ui.load_image(img)

    with pytest.raises(DataErr) as exc:
        ui.set_coord("sky")

    okmsg = "unknown coordinates: 'sky'\nValid options: " + \
            "logical, image, physical, world, wcs"
    assert okmsg in str(exc.value)


@unittest.skip("TODO: failing test used to have a different name and " +
               "was never executed")
@pytest.mark.parametrize("model", ['beta1d', 'lorentz1d', 'normbeta1d'])
def test_psf_model1d(model, clean_astro_ui):
    ui.dataspace1d(1, 10)
    ui.load_psf('psf1d', model + '.mdl')
    ui.set_psf('psf1d')
    mdl = ui.get_model_component('mdl')
    assert (numpy.array(mdl.get_center()) ==
            numpy.array([4])).all()


@pytest.mark.parametrize("model", ['beta2d', 'devaucouleurs2d', 'hubblereynolds', 'lorentz2d'])
def test_psf_model2d(model, clean_astro_ui):
    ui.dataspace2d([216, 261])
    ui.load_psf('psf2d', model + '.mdl')
    ui.set_psf('psf2d')
    mdl = ui.get_model_component('mdl')
    assert (numpy.array(mdl.get_center()) ==
            numpy.array([108, 130])).all()


def test_psf_pars_are_frozen(clean_astro_ui):
    "bug #12503"
    ui.create_model_component("beta2d", "p1")
    p1 = ui.get_model_component("p1")
    ui.load_psf('psf', p1)
    assert p1.thawedpars == []


@requires_data
@requires_fits
def test_chi2(make_data_path, clean_astro_ui):
    "bugs #11400, #13297, #12365"

    data = make_data_path('3c273.pi')

    # Case 1: first ds has no error, second has, chi2-derived (chi2gehrels)
    # statistic. I expect stat.name to be chi2gehrels for ds1, chi2 for
    # ds2, chi2gehrels for ds1,2
    ui.load_data(1, data)
    ui.load_data(2, data, use_errors=True)

    ui.set_source(1, "gauss1d.g1")
    ui.set_source(2, "gauss1d.g1")

    ui.set_stat("chi2gehrels")

    si = ui.get_stat_info()

    stat1 = si[0].statname
    stat2 = si[1].statname
    stat12 = si[2].statname

    assert stat1 == 'chi2gehrels'
    assert stat2 == 'chi2'
    assert stat12 == 'chi2gehrels'

    # Case 2: first ds has errors, second has not, chi2-derived
    # (chi2gehrels) statistic. I expect stat.name to be chi2 for ds1,
    # chi2gehrels for ds2, chi2gehrels for ds1,2
    ui.load_data(2, data)
    ui.load_data(1, data, use_errors=True)

    si = ui.get_stat_info()

    stat1 = si[0].statname
    stat2 = si[1].statname
    stat12 = si[2].statname

    assert stat1 == 'chi2'
    assert stat2 == 'chi2gehrels'
    assert stat12 == 'chi2gehrels'

    # Case 3: both datasets have errors, chi2-derived (chi2gehrels)
    # statistic. I expect stat.name to be chi2 for all of them.
    ui.load_data(2, data, use_errors=True)
    ui.load_data(1, data, use_errors=True)

    si = ui.get_stat_info()

    stat1 = si[0].statname
    stat2 = si[1].statname
    stat12 = si[2].statname

    assert stat1 == 'chi2'
    assert stat2 == 'chi2'
    assert stat12 == 'chi2'

    # Case 4: first ds has errors, second has not, LeastSq statistic
    # I expect stat.name to be leastsq for all of them.
    ui.load_data(2, data)
    ui.load_data(1, data, use_errors=True)

    ui.set_stat("leastsq")

    si = ui.get_stat_info()

    stat1 = si[0].statname
    stat2 = si[1].statname
    stat12 = si[2].statname

    assert stat1 == 'leastsq'
    assert stat2 == 'leastsq'
    assert stat12 == 'leastsq'

    # Case 5: both ds have errors, LeastSq statistic
    # I expect stat.name to be leastsq for all of them.
    ui.load_data(2, data, use_errors=True)
    ui.load_data(1, data, use_errors=True)

    ui.set_stat("leastsq")

    si = ui.get_stat_info()

    stat1 = si[0].statname
    stat2 = si[1].statname
    stat12 = si[2].statname

    assert stat1 == 'leastsq'
    assert stat2 == 'leastsq'
    assert stat12 == 'leastsq'

    # Case 6: first ds has errors, second has not, CStat statistic
    # I expect stat.name to be cstat for all of them.
    ui.load_data(2, data)
    ui.load_data(1, data, use_errors=True)

    ui.set_stat("cstat")

    si = ui.get_stat_info()

    stat1 = si[0].statname
    stat2 = si[1].statname
    stat12 = si[2].statname

    assert stat1 == 'cstat'
    assert stat2 == 'cstat'
    assert stat12 == 'cstat'

    # Case7: select chi2 as statistic. One of the ds does not provide
    # errors. I expect sherpa to raise a StatErr exception.
    ui.set_stat('chi2')

    with pytest.raises(StatErr):
        ui.get_stat_info()

    # Case8: select chi2 as statistic. Both datasets provide errors
    # I expect stat to be 'chi2'
    ui.load_data(2, data, use_errors=True)
    si = ui.get_stat_info()

    stat1 = si[0].statname
    stat2 = si[1].statname
    stat12 = si[2].statname

    assert stat1 == 'chi2'
    assert stat2 == 'chi2'
    assert stat12 == 'chi2'


def save_arrays(colnames, fits, read_func):
    """Write out a small set of data using ui.save_arrays
    and then read it back in, to check it was written out
    correctly (or, at least, in a way that can be read back
    in).

    Parameter
    ---------
    colnames, fits : bool
    read_func : function

    """

    # It looks like the input arrays to `write_arrays` should be numpy
    # arrays, so enforce that invariant.
    a = numpy.asarray([1, 3, 9])
    b = numpy.sqrt(numpy.asarray(a))
    c = b * 0.1

    if colnames:
        fields = ["x", "yy", "z"]
    else:
        fields = None

    ofh = tempfile.NamedTemporaryFile(suffix='sherpa_test')
    ui.save_arrays(ofh.name, [a, b, c], fields=fields,
                   ascii=not fits, clobber=True)

    out = read_func(ofh.name)
    assert isinstance(out, Data1D)

    rtol = 0
    atol = 1e-5

    # remove potential dm syntax introduced by backend before checking for equality
    out_name = re.sub(r"\[.*\]", "", out.name)

    assert ofh.name == out_name, 'file name'
    assert_allclose(out.x, a, rtol=rtol, atol=atol, err_msg="x column")
    assert_allclose(out.y, b, rtol=rtol, atol=atol, err_msg="y column")
    assert_allclose(out.staterror, c, rtol=rtol, atol=atol,
                    err_msg="staterror")
    assert out.syserror is None, 'syserror'


@requires_fits
@pytest.mark.parametrize("colnames", [False, True])
def test_save_arrays_FITS(colnames):

    def read_func(filename):
        return ui.unpack_table(filename, ncols=3)

    save_arrays(colnames=colnames, fits=True, read_func=read_func)


@requires_fits
@pytest.mark.parametrize("colnames", [False, True])
def test_save_arrays_ASCII(colnames):

    def read_func(filename):
        return ui.unpack_ascii(filename, ncols=3)

    save_arrays(colnames=colnames, fits=False, read_func=read_func)


@requires_fits
@requires_data
def testWrite(make_data_path):

    fname = make_data_path('3c273.pi')
    ui.load_pha(1, fname)
    pha_orig = ui.get_data(1)

    ofh = tempfile.NamedTemporaryFile(suffix='sherpa_test')
    ui.save_pha(1, ofh.name, ascii=False, clobber=True)

    # limited checks
    pha = ui.unpack_pha(ofh.name)
    assert isinstance(pha, DataPHA)

    for key in ["channel", "counts"]:
        newval = getattr(pha, key)
        oldval = getattr(pha_orig, key)
        assert_allclose(oldval, newval, err_msg=key)

    # at present grouping and quality are not written out

    for key in ["exposure", "backscal", "areascal"]:
        newval = getattr(pha, key)
        oldval = getattr(pha_orig, key)
        assert newval == pytest.approx(oldval), key

    """
    Since the original file has RMF and ARF, the units
    are energy, but on write out these files are not
    created/saved, so when read back in the PHA has no
    ARF/RMF, and will have units=channel.

    for key in ["units"]:
        newval = getattr(pha, key)
        oldval = getattr(self._pha, key)
        assert newval == oldval, key
    """


@requires_fits
def test_load_table_fits(clean_astro_ui):
    # QUS: why is this not in the sherpa-test-data repository?
    this_dir = os.path.dirname(os.path.abspath(__file__))
    ui.load_table(1, os.path.join(this_dir, 'data',
                                  'two_column_x_y.fits.gz'))
    data = ui.get_data(1)
    assert data.x == pytest.approx([1, 2, 3])
    assert data.y == pytest.approx([4, 5, 6])


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
def test_bug_316(make_data_path, clean_astro_ui):
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

    ui.set_stat('chi2datavar')
    ui.fit()
    fr = ui.get_fit_results()
    assert fr.succeeded
    assert fr.datasets == (1, 2)

    assert src.gamma.val == pytest.approx(1.9298, rel=1.0e-4)
    assert src.ampl.val == pytest.approx(1.73862e-4 / 2, rel=1.0e-4)
