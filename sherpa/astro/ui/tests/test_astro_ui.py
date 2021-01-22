#
#  Copyright (C) 2012, 2015, 2016, 2017, 2018, 2019, 2020, 2021
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
import tempfile
import logging

import pytest

import numpy
from numpy.testing import assert_allclose

from sherpa.astro.data import DataPHA
from sherpa.astro.instrument import RMFModelPHA
from sherpa.astro import ui
from sherpa.data import Data1D
from sherpa.utils.err import ArgumentErr, DataErr, IdentifierErr, StatErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_group, requires_xspec

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
def test_ui_ascii_noarg(setup_files, clean_astro_ui):
    """Don't give a dataset id

    It also lets us actually check the results of the load_table call
    """
    assert ui.list_data_ids() == []
    ui.load_ascii(setup_files.ascii)
    assert ui.list_data_ids() == [1]

    d = ui.get_data()
    assert d.name.endswith('sim.poisson.1.dat')
    assert isinstance(d, ui.Data1D)
    assert len(d.x) == 50
    assert d.x[0:5] == pytest.approx([0.5, 1, 1.5, 2, 2.5])
    assert d.y[0:5] == pytest.approx([27, 27, 20, 28, 27])


@requires_fits
@requires_data
def test_ui_table(setup_files):
    ui.load_table(1, setup_files.fits)
    ui.load_table(1, setup_files.fits, 3)
    ui.load_table(1, setup_files.fits, 3, ["RMID", "SUR_BRI", "SUR_BRI_ERR"])
    ui.load_table(1, setup_files.fits, 4, ('R', "SUR_BRI", 'SUR_BRI_ERR'),
                  ui.Data1DInt)


@requires_fits
@requires_data
def test_ui_table_noarg(setup_files, clean_astro_ui):
    """Don't give a dataset id

    It also lets us actually check the results of the load_table call
    """
    assert ui.list_data_ids() == []
    ui.load_table(setup_files.fits, colkeys=['RMID', 'COUNTS'])
    assert ui.list_data_ids() == [1]

    d = ui.get_data()
    assert d.name.endswith('1838_rprofile_rmid.fits')
    assert isinstance(d, ui.Data1D)
    assert len(d.x) == 38
    assert d.x[0:5] == pytest.approx([12.5, 17.5, 22.5, 27.5, 32.5])
    assert d.y[0:5] == pytest.approx([1529, 2014, 2385, 2158, 2013])


def test_dataspace1d_data1dint(clean_astro_ui):
    """Explicitly test dataspace1d for Data1DInt"""

    assert ui.list_data_ids() == []
    ui.dataspace1d(20, 30, step=2.5, id='x', dstype=ui.Data1DInt)

    assert ui.list_data_ids() == ['x']
    assert ui.get_data('x').name == 'dataspace1d'

    grid = ui.get_indep('x')
    assert len(grid) == 2

    expected = numpy.asarray([20, 22.5, 25, 27.5, 30.0])
    assert grid[0] == pytest.approx(expected[:-1])
    assert grid[1] == pytest.approx(expected[1:])

    y = ui.get_dep('x')
    assert y == pytest.approx(numpy.zeros(4))


def test_dataspace1d_datapha(clean_astro_ui):
    """Explicitly test dataspace1d for DataPHA"""

    assert ui.list_data_ids() == []

    # Note the grid is ignored, other than the number of bins
    ui.dataspace1d(20, 30, step=2.5, id='x', dstype=ui.DataPHA)

    assert ui.list_data_ids() == ['x']
    assert ui.get_data('x').name == 'dataspace1d'

    grid = ui.get_indep('x')
    assert len(grid) == 1

    expected = numpy.asarray([1, 2, 3, 4, 5])
    assert grid[0] == pytest.approx(expected)

    y = ui.get_dep('x')
    assert y == pytest.approx(numpy.zeros(5))

    assert ui.get_exposure('x') is None
    assert ui.get_grouping('x') is None
    assert ui.get_quality('x') is None

    assert ui.get_data('x').subtracted is False

    with pytest.raises(IdentifierErr):
        ui.get_bkg('x')


def test_dataspace1d_datapha_bkg_nopha(clean_astro_ui):
    """We need a PHA to create a background dataset"""

    with pytest.raises(IdentifierErr) as exc:
        ui.dataspace1d(20, 30, step=2.5, id='x', bkg_id=2, dstype=ui.DataPHA)

    assert str(exc.value) == 'data set x has not been set'


def test_dataspace1d_datapha_bkg(clean_astro_ui):
    """Explicitly test dataspace1d for DataPHA (background)"""

    # list_bkg_ids will error out until the dataset exists
    assert ui.list_data_ids() == []

    # We don't use the grid range or step size since numbins has been
    # given.
    ui.dataspace1d(20, 30, step=2.5, numbins=10, id='x', dstype=ui.DataPHA)

    assert ui.list_data_ids() == ['x']
    assert ui.list_bkg_ids('x') == []

    ui.dataspace1d(20, 30, step=2.5, numbins=10, id='x', bkg_id=2,
                   dstype=ui.DataPHA)

    assert ui.list_data_ids() == ['x']
    assert ui.list_bkg_ids('x') == [2]

    assert ui.get_data('x').name == 'dataspace1d'

    # I've explicitly not chosen the default background identifier
    with pytest.raises(IdentifierErr):
        ui.get_bkg('x')

    assert ui.get_bkg('x', 2).name == 'bkg_dataspace1d'

    grid = ui.get_indep('x', bkg_id=2)
    assert len(grid) == 1

    expected = numpy.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    assert grid[0] == pytest.approx(expected)

    y = ui.get_dep('x', bkg_id=2)
    assert y == pytest.approx(numpy.zeros(10))

    assert ui.get_exposure('x', bkg_id=2) is None
    assert ui.get_grouping('x', bkg_id=2) is None
    assert ui.get_quality('x', bkg_id=2) is None

    assert ui.get_bkg('x', bkg_id=2).subtracted is False

    # check we can subtract the dataset; as the data is all zeros
    # we don't bother checking the result.
    #
    ui.subtract('x')


@requires_fits
@requires_data
@pytest.mark.parametrize("loader", [ui.load_data, ui.load_pha])
def test_load_data(loader, make_data_path, clean_astro_ui, caplog):
    """Ensure that loading a single file to a non-integer id works.

    This is just to make sure that the support for PHA2 files in both
    load_data and load_pha does not change the single-file case.

    """
    infile = make_data_path('3c273.pi')
    bgfile = make_data_path('3c273_bg.pi')

    arf = make_data_path('3c273.arf')
    rmf = make_data_path('3c273.rmf')

    assert ui.list_data_ids() == []

    with SherpaVerbosity('INFO'):
        loader('foo', infile)

    assert ui.list_data_ids() == ['foo']

    msg1 = f"systematic errors were not found in file '{infile}'"
    msg2 = f"statistical errors were found in file '{infile}' \n" + \
        "but not used; to use them, re-read with use_errors=True"
    msg3 = f"read ARF file {arf}"
    msg4 = f"read RMF file {rmf}"

    msg5 = f"systematic errors were not found in file '{bgfile}'"
    msg6 = f"statistical errors were found in file '{bgfile}' \n" + \
        "but not used; to use them, re-read with use_errors=True"
    msg7 = f"read background file {bgfile}"

    assert caplog.record_tuples[0] == ('sherpa.astro.io', logging.WARNING, msg1)
    assert caplog.record_tuples[1] == ('sherpa.astro.io', logging.INFO, msg2)
    assert caplog.record_tuples[2] == ('sherpa.astro.io', logging.INFO, msg3)
    assert caplog.record_tuples[3] == ('sherpa.astro.io', logging.INFO, msg4)

    assert caplog.record_tuples[4] == ('sherpa.astro.io', logging.WARNING, msg5)
    assert caplog.record_tuples[5] == ('sherpa.astro.io', logging.INFO, msg6)
    assert caplog.record_tuples[6] == ('sherpa.astro.io', logging.INFO, msg7)

    assert len(caplog.records) == 7


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


@requires_fits
@requires_data
def test_ui_filter_ascii_with_id(setup_files):
    ui.load_filter(1, setup_files.filter_single_int_ascii)
    assert ui.list_data_ids() == [1]

    f = ui.get_filter()
    assert f == '2.0000,4.0000,6.0000,8.0000,10.0000:250.0000,751.0000:1000.0000'


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


@requires_fits
@requires_data
def test_bug38_filtering(make_data_path):
    """Low-level tests related to bugs #38, #917: filter"""

    from sherpa.astro.io import read_pha
    pha = read_pha(make_data_path('3c273.pi'))

    assert pha.mask is True
    pha.notice(0.3, 2)
    assert pha.mask.size == 46
    assert pha.mask.sum() == 25
    assert pha.mask[1:26].all()


@requires_fits
@requires_data
def test_bug38_filtering_grouping(make_data_path):
    """Low-level tests related to bugs #38, #917: filter+group"""

    from sherpa.astro.io import read_pha
    pha = read_pha(make_data_path('3c273.pi'))

    pha.notice(1, 6)
    pha.ignore(3, 4)

    assert pha.get_filter(group=True, format='%.4f') == '1.0147:2.7886,4.1391:6.2342'
    assert pha.get_filter(group=False, format='%.4f') == '1.0001:2.8543,4.0369:6.5627'

    pha.group_width(40)

    assert pha.get_filter(group=True, format='%.4f') == '0.8760:2.6280,3.7960:6.7160'
    assert pha.get_filter(group=False, format='%.4f') == '0.5913:2.9127,3.5113:7.0007'

    assert pha.mask.size == 26
    assert pha.mask.sum() == 10
    assert pha.mask[1:5].all()
    assert pha.mask[6:12].all()

    # get the ungrouped mask
    mask = pha.get_mask()
    assert mask.sum() == 10 * 40
    assert mask[40:200].all()
    assert mask[240:480].all()

    # check filtered bins
    elo_all, ehi_all = pha._get_ebins(group=False)
    elo, ehi = pha._get_ebins(group=True)

    assert elo[1] == elo_all[40]
    assert ehi[11] == ehi_all[479]


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


@requires_data
@requires_fits
def test_image_with_id(make_data_path, clean_astro_ui):
    """Call load_image with an identifier"""

    img = make_data_path('img.fits')

    assert ui.list_data_ids() == []
    ui.load_image('ix', img)
    assert ui.list_data_ids() == ['ix']

    d = ui.get_data('ix')
    assert isinstance(d, ui.DataIMG)
    assert d.name.endswith('img.fits')


# DJB notes (2020/02/29) that these tests used to not be run,
# because the lorentz1d model retrned 1 rather than 4. This
# has now been included in the test so that we can check
# the behavior if the PSF centering code is ever changed.
#
@pytest.mark.parametrize("model,center", [('beta1d', 4.0),
                                          ('lorentz1d', 1.0),
                                          ('normbeta1d', 4.0)])
def test_psf_model1d(model, center, clean_astro_ui):
    ui.dataspace1d(1, 10)
    ui.load_psf('psf1d', model + '.mdl')
    ui.set_psf('psf1d')
    mdl = ui.get_model_component('mdl')
    assert mdl.get_center() == (center, )


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


@pytest.mark.parametrize("loadfunc", [ui.load_multi_arfs, ui.load_multi_rmfs])
@pytest.mark.parametrize("idval", [None, 1, 'xx'])
def test_load_multi_xxx_invalid_args(loadfunc, idval):
    """Check we error out if the files and ids do not match."""

    files = ['a', 'b', 'c']
    rids = [1, 2]

    with pytest.raises(ArgumentErr) as exc:
        if idval is None:
            loadfunc(files, rids)
        else:
            loadfunc(idval, files, rids)

    assert str(exc.value) == 'A response ID is required for each file'


@requires_data
@requires_fits
def test_load_multi_arfsrmfs(make_data_path, clean_astro_ui):
    """Added in #728 to ensure cache parameter is sent along by
    MultiResponseSumModel (fix #717).

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

    # Check multiple responses have been loaded
    #
    d1 = ui.get_data(1)
    d2 = ui.get_data(2)
    assert d1.response_ids == [1, 2]
    assert d2.response_ids == [1, 2]

    # Unfortunately we load the same response so it's hard
    # to tell the difference here!
    #
    assert ui.get_arf(resp_id=1).name == arf
    assert ui.get_arf(resp_id=2).name == arf
    assert ui.get_rmf(2, resp_id=1).name == rmf
    assert ui.get_rmf(2, resp_id=2).name == rmf

    ui.notice(0.5, 7)
    ui.subtract(1)
    ui.subtract(2)

    src = ui.create_model_component('powlaw1d', 'src')
    ui.set_model(1, src)
    ui.set_model(2, src)

    # ensure the test is repeatable by running with a known
    # statistic and method
    #
    ui.set_method('levmar')
    ui.set_stat('chi2datavar')

    # Really what we care about for fixing #717 is that
    # fit does not error out, but it's useful to know that
    # the fit has changed the parameter values (which were
    # both 1 before the fit).
    #
    ui.fit()
    fr = ui.get_fit_results()
    assert fr.succeeded
    assert fr.datasets == (1, 2)

    assert src.gamma.val == pytest.approx(1.9298, rel=1.0e-4)
    assert src.ampl.val == pytest.approx(1.73862e-4 / 2, rel=1.0e-4)


@requires_xspec
@requires_fits
@requires_data
def test_pileup_model(make_data_path, clean_astro_ui):
    """Basic check of setting a pileup model.

    It is more to check we can set a pileup model, not to
    check the model works.
    """

    infile = make_data_path('3c273.pi')
    ui.load_pha('pileup', infile)
    ui.subtract('pileup')
    ui.notice(0.3, 7)

    ui.set_stat('chi2datavar')

    # pick xswabs as it is unlikely to change with XSPEC
    ui.set_source('pileup', ui.xswabs.amdl * ui.powlaw1d.pl)

    # get close to the best fit, but don't need to fit
    pl.ampl = 1.82e-4
    pl.gamma = 1.97
    amdl.nh = 0.012

    stat0 = ui.calc_stat('pileup')

    # We want to compare the data to the pileup model,
    # which should be run with the higher-energy bins included,
    # but that's not relevant here where I am just
    # checking the statistic value.

    ui.set_pileup_model('pileup', ui.jdpileup.jdp)

    # pick some values to make the model change the data
    jdp.ftime = 3.2
    jdp.fracexp = 1
    jdp.alpha = 0.95
    jdp.f = 0.91

    # Check pileup is added to get_model
    #
    mlines = str(ui.get_model('pileup')).split('\n')
    assert mlines[0] == 'apply_rmf(jdpileup.jdp((xswabs.amdl * powlaw1d.pl)))'
    assert mlines[3].strip() == 'jdp.alpha    thawed         0.95            0            1'
    assert mlines[5].strip() == 'jdp.f        thawed         0.91          0.9            1'
    assert mlines[11].strip() == 'pl.gamma     thawed         1.97          -10           10'

    # Ensure that the statistic has got worse (technically
    # it coud get better, but not for this case).
    #
    stat1 = ui.calc_stat('pileup')
    assert stat1 > stat0

    # As a test, check the actual statistic values
    # (evaluated with XSPEC 12.11.0 and Sherpa with the
    # master branch 2020/07/29, on Linux).
    #
    assert stat0 == pytest.approx(35.99899827358692)
    assert stat1 == pytest.approx(36.58791181460404)

    # Can we remove the pileup model?
    #
    ui.delete_pileup_model('pileup')

    # Check pileup is not in get_model
    #
    mlines = str(ui.get_model('pileup')).split('\n')

    # Numeric display depends on Python and/or NumPy
    # for the exposure time. I should have set the exposure
    # time to an integer value.
    #
    # assert mlines[0] == 'apply_rmf(apply_arf((38564.608926889 * (xswabs.amdl * powlaw1d.pl))))'

    toks = mlines[0].split()
    assert len(toks) == 5
    assert toks[0].startswith('apply_rmf(apply_arf((38564.608')
    assert toks[1] == '*'
    assert toks[2] == '(xswabs.amdl'
    assert toks[3] == '*'
    assert toks[4] == 'powlaw1d.pl))))'

    assert mlines[4].strip() == 'pl.gamma     thawed         1.97          -10           10'

    stat2 = ui.calc_stat('pileup')
    assert stat2 == stat0


def test_list_pileup_ids_empty(clean_astro_ui):
    assert ui.list_pileup_model_ids() == []


@pytest.mark.parametrize("id", [None, 2, "2"])
def test_list_pileup_ids_single(id, clean_astro_ui):

    # Note: do not need to set a source model
    ui.load_arrays(id, [1, 2, 3], [1, 1, 1], ui.DataPHA)
    ui.set_pileup_model(id, ui.jdpileup.jdp)

    if id is None:
        id = 1

    assert ui.list_pileup_model_ids() == [id]


def test_list_pileup_ids_multi(clean_astro_ui):

    jdp = ui.create_model_component('jdpileup', "jdp")

    for id in [1, "2"]:
        ui.load_arrays(id, [1, 2, 3], [1, 1, 1], ui.DataPHA)
        ui.set_pileup_model(id, jdp)

    ans = ui.list_pileup_model_ids()
    assert ans == [1, "2"]


def check_bad_grouping(exp_xlo,exp_xhi, exp_counts, lo1, hi1, lo2, hi2):
    """Common tests from test_grouped_pha_all_badXXX

    Sending in two ranges is a bit excessive but easiest
    thing to implement
    """

    cts = ui.get_counts()
    assert cts == pytest.approx([exp_counts])

    dplot = ui.get_data_plot()
    assert dplot.xlo == pytest.approx([exp_xlo])
    assert dplot.xhi == pytest.approx([exp_xhi])
    assert dplot.y == pytest.approx([exp_counts])

    # ignore all the data
    ui.ignore(lo1, hi1)

    # can still plot
    cts = ui.get_counts()
    assert cts == pytest.approx([exp_counts])

    cts = ui.get_counts(filter=True)
    assert len(cts) == 0

    dplot = ui.get_data_plot()
    assert len(dplot.xlo) == 0
    assert len(dplot.xhi) == 0
    assert len(dplot.y) == 0

    # ignore does not fail
    #
    ui.ignore(lo2, hi2)

    # we can restore the data
    ui.notice(None, None)

    cts = ui.get_counts()
    assert cts == pytest.approx([exp_counts])

    dplot = ui.get_data_plot()
    assert dplot.xlo == pytest.approx([exp_xlo])
    assert dplot.xhi == pytest.approx([exp_xhi])
    assert dplot.y == pytest.approx([exp_counts])

    # now ignore the bad channels (ie everything)
    #
    ui.ignore_bad()

    cts = ui.get_counts()
    assert len(cts) == 0

    dplot = ui.get_data_plot()
    assert len(dplot.xlo) == 0
    assert len(dplot.xhi) == 0
    assert len(dplot.y) == 0

    # there's nothing to notice (this line is an example of #790)
    ui.notice(lo1, hi1)

    cts = ui.get_counts()
    assert len(cts) == 0

    dplot = ui.get_data_plot()
    assert len(dplot.xlo) == 0
    assert len(dplot.xhi) == 0
    assert len(dplot.y) == 0


def test_grouped_pha_all_bad_channel(clean_astro_ui):
    """Helpdesk ticket: low-count data had no valid bins after grouping #790

    A simple PHA dataset is created, with no response, which has no
    "good" grouped data (1 group, but with a quality of 2). Several
    checks are made to ensure we can filter/notice/plot the data even
    when it is empty.  This is in channel space.

    See also test_grouped_pha_all_bad_response which is essentially the
    same test but with a response model.

    """
    chans = numpy.arange(1, 6, dtype=numpy.int16)
    counts = numpy.asarray([0, 1, 2, 0, 1], dtype=numpy.int16)
    grouping = numpy.asarray([1, -1, -1, -1, -1], dtype=numpy.int16)
    quality = numpy.asarray([2, 2, 2, 2, 2], dtype=numpy.int16)

    dset = ui.DataPHA('low', chans, counts, grouping=grouping, quality=quality)
    ui.set_data(dset)

    # Run tests
    check_bad_grouping(1, 6, 0.8, 0, 6, 2, 10)


@pytest.mark.parametrize("arf,rmf", [(True, False), (False, True), (True, True)])
@pytest.mark.parametrize("chantype,exp_counts,exp_xlo,exp_xhi,lo1,hi1,lo2,hi2",
                         [("channel", 0.8, 1.0, 6, 0, 7, 2, 6),
                          ("energy", 8.0, 0.1, 0.6, 0.05, 1.0, 0.2, 0.8),
                          ("wave", 0.03871461, 123.9841874, 20.66403123, 20, 90, 30, 85)
                         ])
def test_grouped_pha_all_bad_response(arf, rmf, chantype, exp_counts, exp_xlo, exp_xhi, lo1, hi1, lo2, hi2, clean_astro_ui):
    """Helpdesk ticket: low-count data had no valid bins after grouping #790

    A simple PHA dataset is created, which has no "good" grouped data
    (1 group, but with a quality of 2). Several checks are made to
    ensure we can filter/notice/plot the data even when it is empty.

    Checks are done for
      - arf-only
      - rmf-only
      - arf+rmf
    analysis in case there's a difference in the code paths

    """

    chans = numpy.arange(1, 6, dtype=numpy.int16)
    counts = numpy.asarray([0, 1, 2, 0, 1], dtype=numpy.int16)
    grouping = numpy.asarray([1, -1, -1, -1, -1], dtype=numpy.int16)
    quality = numpy.asarray([2, 2, 2, 2, 2], dtype=numpy.int16)

    dset = ui.DataPHA('low', chans, counts, grouping=grouping, quality=quality)
    ui.set_data(dset)

    egrid = numpy.asarray([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    elo = egrid[:-1]
    ehi = egrid[1:]

    # it is required that at least one of arf or rmf is set but this
    # is not enforced
    #
    if arf:
        ui.set_arf(ui.create_arf(elo, ehi))

    if rmf:
        # NOTE: need to set e_min/max otherwise get a 'noenergybins'
        #       error from sherpa.astro.data.DataPHA._get_ebins
        #
        ui.set_rmf(ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi))

    # plot units depend on analysis type;
    #
    ui.set_analysis(chantype)

    # Run tests
    check_bad_grouping(exp_xlo, exp_xhi, exp_counts, lo1, hi1, lo2, hi2)


@requires_fits
@requires_data
@pytest.mark.parametrize("elo,ehi,nbins",
                         [(None, 5, 41), (1, None, 33), (1, 5, 28)])
@pytest.mark.parametrize("bkg_id", [None, 1])
def test_grouped_pha_all_bad_response_bg_warning(elo, ehi, nbins, bkg_id,
                                                 caplog, make_data_path, clean_astro_ui):
    """Check we get the warning messages with background filtering"""

    ui.load_pha('check', make_data_path('3c273.pi'))

    ui.set_quality('check', 2 * numpy.ones(1024, dtype=numpy.int16), bkg_id=1)
    ui.ignore_bad('check', bkg_id=1)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.notice_id('check', elo, ehi, bkg_id=bkg_id)

    # filtering has or hasn't happened
    nsrc = ui.get_dep('check', filter=True).size
    nback = ui.get_dep('check', filter=True, bkg_id=1).size

    if bkg_id is None:
        assert nsrc == nbins
        assert nback == 0
    else:
        assert nsrc == 46  # ie no filter
        assert nback == 0

    # did we get a warning message from the background?
    assert len(caplog.records) == 1
    name, lvl, msg = caplog.record_tuples[0]
    assert name == 'sherpa.astro.data'
    assert lvl == logging.INFO
    assert msg.startswith('Skipping dataset /')
    assert msg.endswith('/3c273_bg.pi: mask excludes all data')


@pytest.mark.parametrize('dtype', [ui.Data1D, ui.Data1DInt, ui.DataPHA,
                                   ui.Data2D, ui.Data2DInt,
                                   ui.DataIMG, ui.DataIMGInt])
def test_unpack_arrays_missing(dtype):
    """What happens with the wrong number of arguments"""

    with pytest.raises(TypeError):
        ui.unpack_arrays([1, 2, 5], dtype)


def test_unpack_arrays_data1d():
    """Can we call unpack_arrays? Data1D"""

    ans = ui.unpack_arrays([1, 2, 5], [-2, 3, -1])
    assert isinstance(ans, ui.Data1D)
    assert ans.name == ''
    assert ans.x == pytest.approx([1, 2, 5])
    assert ans.y == pytest.approx([-2, 3, -1])
    assert ans.staterror is None
    assert ans.syserror is None


def test_unpack_arrays_data1d_too_many():
    """What happens with too many columns?

    It just reads in the extra columns, even if "garbage"
    (such as the negative errors). So technically this is not
    an error case.
    """

    ans = ui.unpack_arrays([1, 2, 5], [2, 3, 6], [-2, 3, -1])
    assert isinstance(ans, ui.Data1D)
    assert ans.name == ''
    assert ans.x == pytest.approx([1, 2, 5])
    assert ans.y == pytest.approx([2, 3, 6])
    assert ans.staterror == pytest.approx([-2, 3, -1])
    assert ans.syserror is None


def test_unpack_arrays_data1dint():
    """Can we call unpack_arrays? Data1DInt"""

    ans = ui.unpack_arrays([1, 2, 5], [2, 3, 10], [-2, 3, -1], ui.Data1DInt)
    assert isinstance(ans, ui.Data1DInt)
    assert ans.name == ''
    assert ans.xlo == pytest.approx([1, 2, 5])
    assert ans.xhi == pytest.approx([2, 3, 10])
    assert ans.y == pytest.approx([-2, 3, -1])
    assert ans.staterror is None
    assert ans.syserror is None


def test_unpack_arrays_datapha():
    """Can we call unpack_arrays? DataPHA"""

    ans = ui.unpack_arrays([1, 2, 3], [2, 0, 4], ui.DataPHA)
    assert isinstance(ans, ui.DataPHA)
    assert ans.name == ''
    assert ans.channel == pytest.approx([1, 2, 3])
    assert ans.counts == pytest.approx([2, 0, 4])
    assert ans.staterror is None
    assert ans.syserror is None


def test_unpack_arrays_dataimg():
    """Can we call unpack_arrays? DataIMG"""

    ivals = numpy.arange(12)
    y, x = numpy.mgrid[0:3, 0:4]
    x = x.flatten()
    y = y.flatten()
    ans = ui.unpack_arrays(x, y, ivals, (3, 4), ui.DataIMG)
    assert isinstance(ans, ui.DataIMG)
    assert ans.name == ''
    assert ans.x0 == pytest.approx(x)
    assert ans.x1 == pytest.approx(y)
    assert ans.y == pytest.approx(ivals)
    assert ans.shape == pytest.approx(3, 4)
    assert ans.staterror is None
    assert ans.syserror is None
