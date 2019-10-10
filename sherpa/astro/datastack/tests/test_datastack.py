#
#  Copyright (C) 2014, 2015, 2016, 2018, 2019
#       Smithsonian Astrophysical Observatory
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
import tempfile
import logging

import numpy as np

from sherpa.utils.testing import requires_fits, requires_stk
from sherpa.astro import ui
from sherpa.astro import datastack
from acis_bkg_model import acis_bkg_model

import pytest

logger = logging.getLogger('sherpa')

DATADIR = os.path.join(os.path.dirname(__file__), 'data')


@pytest.fixture
def ds_setup():
    """Setup and teardown code for each test."""

    # Setup
    #
    datastack.clear_stack()
    ui.clean()
    loggingLevel = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    datastack.set_stack_verbosity(logging.ERROR)
    datastack.set_template_id("__ID")

    # Run test
    #
    yield

    # Cleanup
    #
    datastack.clear_stack()
    ui.clean()
    datastack.set_template_id("__ID")
    logger.setLevel(loggingLevel)


@pytest.fixture
def ds_setup_object():
    """Setup and teardown code for each test.

    Could try and be clever and re-use ds_setup here,
    but just repeat it to be simpler.
    """

    # Setup
    #
    ds = datastack.DataStack()
    datastack.clear_stack()
    ui.clean()
    loggingLevel = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    datastack.set_stack_verbosity(logging.ERROR)
    datastack.set_template_id("__ID")

    # Run test, returning the stack object
    #
    yield ds

    # Cleanup
    #
    ds.clear_stack()
    datastack.clear_stack()
    ui.clean()
    datastack.set_template_id("__ID")
    logger.setLevel(loggingLevel)


@pytest.fixture
def ds_datadir():
    """Location of the data directory for the DataStack tests."""

    return DATADIR


@requires_fits
@requires_stk
def test_design_case_1(ds_setup, ds_datadir):

    datastack.set_template_id("ID")

    datadir = ds_datadir
    ls = '@' + '/'.join((datadir, '3c273.lis'))
    datastack.load_pha(ls, use_errors=True)

    assert 2 == len(datastack.DATASTACK.datasets)
    assert 2 == len(ui._session._data)

    datastack.load_pha("myid", '/'.join((datadir, "3c273.pi")))

    assert 2 == len(datastack.DATASTACK.datasets)
    assert 3 == len(ui._session._data)

    datastack.load_pha('/'.join((datadir, "3c273.pi")))

    assert 3 == len(datastack.DATASTACK.datasets)
    assert 4 == len(ui._session._data)

    datastack.load_pha([], '/'.join((datadir, "3c273.pi")))

    assert 4 == len(datastack.DATASTACK.datasets)
    assert 5 == len(ui._session._data)

    ds = datastack.DataStack()

    datastack.load_pha(ds, ls)

    assert 4 == len(datastack.DATASTACK.datasets)
    assert 7 == len(ui._session._data)
    assert 2 == len(ds.datasets)

    datastack.load_pha(ds, '/'.join((datadir, "3c273.pi")))

    assert 4 == len(datastack.DATASTACK.datasets)
    assert 8 == len(ui._session._data)
    assert 3 == len(ds.datasets)

    dids = datastack.DATASTACK.get_stack_ids()
    assert dids == [1, 2, 3, 4]

    sids = set(ui._session._data.keys())
    assert sids == {1, 2, 3, 4, 5, 6, 7, "myid"}

    datastack.set_source([1, 2], "powlaw1d.pID")
    datastack.set_source([3, 4], "brokenpowerlaw.bpID")

    dsids = ds.get_stack_ids()
    assert dsids == [5, 6, 7]

    p1 = ui._session._model_components['p1']
    p2 = ui._session._model_components['p2']
    bp3 = ui._session._model_components['bp3']
    bp4 = ui._session._model_components['bp4']

    assert p1 is not None
    assert p2 is not None
    assert bp3 is not None
    assert bp4 is not None

    datastack.set_source(1, "polynom1d.poly1")
    datastack.set_source([2, 3, 4], "atten.attID")

    poly1 = ui._session._model_components['poly1']
    a2 = ui._session._model_components['att2']
    a3 = ui._session._model_components['att3']
    a4 = ui._session._model_components['att4']

    assert poly1 is not None
    assert a2 is not None
    assert a3 is not None
    assert a4 is not None

    datastack.clean()

    assert 0 == len(datastack.DATASTACK.datasets)
    assert 0 == len(ui._session._data)
    assert 3 == len(ds.datasets)


def test_global_case_2(ds_setup):

    x1 = np.arange(50) + 100
    y1 = 2 * (3 * x1**2 + x1)
    x2 = np.arange(50)
    y2 = 2 * (x2**2 + 3 * x2)
    x3 = np.arange(50) + 200
    y3 = 2 * (x3**2 + 3 * x3)

    datastack.load_arrays([[x1, y1], [x2, y2], [x3, y3]])

    datastack.set_source([], 'const1d.const * polynom1d.poly__ID')

    poly1 = ui._session._get_model_component('poly1')
    poly2 = ui._session._get_model_component('poly2')
    poly3 = ui._session._get_model_component('poly3')
    const = ui._session._get_model_component('const')

    datastack.freeze([], 'poly')
    datastack.thaw([], 'poly.c0')
    datastack.thaw([], 'poly.c1')
    datastack.thaw([], 'poly.c2')
    datastack.thaw([], 'const')

    assert poly1.c0.frozen is False
    assert poly1.c1.frozen is False
    assert poly1.c2.frozen is False
    assert poly1.c3.frozen is True

    assert poly2.c0.frozen is False
    assert poly2.c1.frozen is False
    assert poly2.c2.frozen is False
    assert poly2.c3.frozen is True

    assert poly3.c0.frozen is False
    assert poly3.c1.frozen is False
    assert poly3.c2.frozen is False
    assert poly3.c3.frozen is True

    datastack.set_par([], 'poly.c1', 0.45)

    assert poly1.c1.val == 0.45
    assert poly2.c1.val == 0.45
    assert poly3.c1.val == 0.45

    datastack.set_par([1], 'poly.c1', 0.1)

    assert poly1.c1.val == 0.1
    assert poly2.c1.val == 0.45
    assert poly3.c1.val == 0.45

    datastack.set_par([], 'const.c0', 2)

    assert const.c0.val == 2

    datastack.set_par([], 'const.integrate', False)
    datastack.freeze([], 'const.c0')

    vals = datastack.get_par([], 'poly.c1.val')
    assert ([0.1, 0.45, 0.45] == vals).all()

    # QUS: pars is not checked, so is this just
    # checking that get_par doesn't fail?
    pars = datastack.get_par([], 'const.c0')

    datastack.fit([])

    assert round(poly1.c1.val) == 1
    assert round(poly1.c2.val) == 3
    assert round(poly2.c1.val) == 3
    assert round(poly2.c2.val) == 1
    assert round(poly3.c1.val) == 3
    assert round(poly3.c2.val) == 1

    datastack.clear_stack()

    x1 = np.arange(50) + 100
    y1 = 7 * (3 * x1**2 + x1)
    x2 = np.arange(50)
    y2 = 2 * (x2**2 + 3 * x2)
    x3 = np.arange(50) + 200
    y3 = 2 * (x3**2 + 3 * x3)

    datastack.load_arrays([[x1, y1], [x2, y2], [x3, y3]])

    datastack.set_template_id("foo")

    datastack.set_source([], 'const1d.constfoo * polynom1d.polyfoo')

    # const1 is not used, so should it be removed?
    const1 = ui._session._get_model_component('const1')
    const2 = ui._session._get_model_component('const2')
    const3 = ui._session._get_model_component('const3')

    datastack.link([2, 3], 'const.c0')

    datastack.set_par([2], 'const.c0', 3)
    datastack.set_par([1], 'const.c0', 7)

    datastack.freeze([1], 'const.c0')

    assert const2.c0.frozen is False

    datastack.fit([])

    assert const2.c0.val == const3.c0.val
    assert const3.c0._link is const2.c0

    datastack.unlink([], "const.c0")

    assert const3.c0._link is not const2.c0


def create_files():
    fd1, name1 = tempfile.mkstemp()
    fd2, name2 = tempfile.mkstemp()

    fdlis, lisname = tempfile.mkstemp()

    ascii1 = open(name1, 'w')
    ascii2 = open(name2, 'w')
    lis = open(lisname, 'w')

    x1 = np.arange(50) + 100
    y1 = 2 * (3 * x1**2 + x1)
    x2 = np.arange(50)
    y2 = 2 * (x2**2 + 3 * x2)

    for x, y in zip(x1, y1):
        ascii1.write("{} {}\n".format(x, y))
    ascii1.close()
    os.close(fd1)

    for x, y in zip(x2, y2):
        ascii2.write("{} {}\n".format(x, y))
    ascii2.close()
    os.close(fd2)

    lis.write("{}\n".format(name1))
    lis.write("{}\n".format(name2))
    lis.close()
    os.close(fdlis)

    return [name1, name2, lisname]


def remove_files(infiles):
    for infile in infiles:
        os.remove(infile)


@requires_fits
@requires_stk
def test_load_case_3(ds_setup, ds_datadir):

    datadir = ds_datadir
    name1, name2, lisname = create_files()

    try:
        datastack.load_ascii("@{}".format(lisname))
        assert len(ui._session._data) == 2
        assert len(datastack.DATASTACK.datasets) == 2

        datastack.load_data("@{}".format(lisname))
        assert len(ui._session._data) == 4
        assert len(datastack.DATASTACK.datasets) == 4

        datastack.load_data("@" +
                            "/".join((datadir, 'pha.lis')))
        assert len(ui._session._data) == 6
        assert len(datastack.DATASTACK.datasets) == 6

    except:
        remove_files([name1, name2, lisname])
        raise


def test_partial_oo_case_4(ds_setup_object):

    ds = ds_setup_object

    x1 = np.arange(50) + 100
    y1 = 2 * (3 * x1**2 + x1)
    x2 = np.arange(50)
    y2 = 2 * (x2**2 + 3 * x2)
    x3 = np.arange(50) + 200
    y3 = 2 * (x3**2 + 3 * x3)

    datastack.load_arrays(ds, [[x1, y1], [x2, y2], [x3, y3]])

    datastack.set_source(ds, 'const1d.const * polynom1d.poly__ID')

    poly1 = ui._session._get_model_component('poly1')
    poly2 = ui._session._get_model_component('poly2')
    poly3 = ui._session._get_model_component('poly3')
    const = ui._session._get_model_component('const')

    datastack.freeze(ds, 'poly')
    datastack.thaw(ds, 'poly.c0')
    datastack.thaw(ds, 'poly.c1')
    datastack.thaw(ds, 'poly.c2')
    datastack.thaw(ds, 'const')

    assert poly1.c0.frozen is False
    assert poly1.c1.frozen is False
    assert poly1.c2.frozen is False
    assert poly1.c3.frozen is True

    assert poly2.c0.frozen is False
    assert poly2.c1.frozen is False
    assert poly2.c2.frozen is False
    assert poly2.c3.frozen is True

    assert poly3.c0.frozen is False
    assert poly3.c1.frozen is False
    assert poly3.c2.frozen is False
    assert poly3.c3.frozen is True

    datastack.set_par(ds, 'poly.c1', 0.45)

    assert poly1.c1.val == 0.45
    assert poly2.c1.val == 0.45
    assert poly3.c1.val == 0.45

    datastack.set_par(ds[1], 'poly.c1', 0.1)

    assert poly1.c1.val == 0.1
    assert poly2.c1.val == 0.45
    assert poly3.c1.val == 0.45

    datastack.set_par(ds, 'const.c0', 2)

    assert const.c0.val == 2

    datastack.set_par(ds, 'const.integrate', False)
    datastack.freeze(ds, 'const.c0')

    vals = datastack.get_par(ds, 'poly.c1.val')
    assert ([0.1, 0.45, 0.45] == vals).all()

    # QUS: pars is not checked, so is this just
    # checking that get_par doesn't fail?
    pars = datastack.get_par(ds, 'const.c0')

    datastack.fit(ds)

    assert round(poly1.c1.val) == 1
    assert round(poly1.c2.val) == 3
    assert round(poly2.c1.val) == 3
    assert round(poly2.c2.val) == 1
    assert round(poly3.c1.val) == 3
    assert round(poly3.c2.val) == 1

    ds.clear_stack()

    x1 = np.arange(50) + 100
    y1 = 7 * (3 * x1**2 + x1)
    x2 = np.arange(50)
    y2 = 2 * (x2**2 + 3 * x2)
    x3 = np.arange(50) + 200
    y3 = 2 * (x3**2 + 3 * x3)

    datastack.load_arrays(ds, [[x1, y1], [x2, y2], [x3, y3]])

    datastack.set_template_id("foo")

    datastack.set_source(ds, 'const1d.constfoo * polynom1d.polyfoo')

    const1 = ui._session._get_model_component('const1')
    const2 = ui._session._get_model_component('const2')
    const3 = ui._session._get_model_component('const3')

    datastack.link(ds[2, 3], 'const.c0')

    datastack.set_par(ds[2], 'const.c0', 3)
    datastack.set_par(ds[1], 'const.c0', 7)

    datastack.freeze(ds[1], 'const.c0')

    assert const2.c0.frozen is False

    datastack.fit(ds)

    assert const2.c0.val == const3.c0.val
    assert const3.c0._link is const2.c0

    datastack.unlink(ds, "const.c0")

    assert const3.c0._link is not const2.c0


def test_oo_case_5(ds_setup_object):

    ds = ds_setup_object

    x1 = np.arange(50) + 100
    y1 = 2 * (3 * x1**2 + x1)
    x2 = np.arange(50)
    y2 = 2 * (x2**2 + 3 * x2)
    x3 = np.arange(50) + 200
    y3 = 2 * (x3**2 + 3 * x3)

    ds.load_arrays([[x1, y1], [x2, y2], [x3, y3]])

    ds.set_source('const1d.const * polynom1d.poly__ID')

    poly1 = ui._session._get_model_component('poly1')
    poly2 = ui._session._get_model_component('poly2')
    poly3 = ui._session._get_model_component('poly3')
    const = ui._session._get_model_component('const')

    ds.freeze('poly')
    ds.thaw('poly.c0')
    ds.thaw('poly.c1')
    ds.thaw('poly.c2')
    ds.thaw('const')

    assert poly1.c0.frozen is False
    assert poly1.c1.frozen is False
    assert poly1.c2.frozen is False
    assert poly1.c3.frozen is True

    assert poly2.c0.frozen is False
    assert poly2.c1.frozen is False
    assert poly2.c2.frozen is False
    assert poly2.c3.frozen is True

    assert poly3.c0.frozen is False
    assert poly3.c1.frozen is False
    assert poly3.c2.frozen is False
    assert poly3.c3.frozen is True

    ds.set_par('poly.c1', 0.45)

    assert poly1.c1.val == 0.45
    assert poly2.c1.val == 0.45
    assert poly3.c1.val == 0.45

    ds[1].set_par('poly.c1', 0.1)

    assert poly1.c1.val == 0.1
    assert poly2.c1.val == 0.45
    assert poly3.c1.val == 0.45

    ds.set_par('const.c0', 2)

    assert const.c0.val == 2

    ds.set_par('const.integrate', False)
    ds.freeze('const.c0')

    vals = ds.get_par('poly.c1.val')
    assert ([0.1, 0.45, 0.45] == vals).all()

    # QUS: pars is not checked, so is this just
    # checking that get_par doesn't fail?
    pars = ds.get_par('const.c0')

    ds.fit()

    assert round(poly1.c1.val) == 1
    assert round(poly1.c2.val) == 3
    assert round(poly2.c1.val) == 3
    assert round(poly2.c2.val) == 1
    assert round(poly3.c1.val) == 3
    assert round(poly3.c2.val) == 1

    ds.clear_stack()

    x1 = np.arange(50) + 100
    y1 = 7 * (3 * x1**2 + x1)
    x2 = np.arange(50)
    y2 = 2 * (x2**2 + 3 * x2)
    x3 = np.arange(50) + 200
    y3 = 2 * (x3**2 + 3 * x3)

    ds.load_arrays([[x1, y1], [x2, y2], [x3, y3]])

    datastack.set_template_id("foo")

    ds.set_source('const1d.constfoo * polynom1d.polyfoo')

    const1 = ui._session._get_model_component('const1')
    const2 = ui._session._get_model_component('const2')
    const3 = ui._session._get_model_component('const3')

    ds[2, 3].link('const.c0')

    ds[2].set_par('const.c0', 3)
    ds[1].set_par('const.c0', 7)

    ds[1].freeze('const.c0')

    assert const2.c0.frozen is False

    ds.fit()

    assert const2.c0.val == const3.c0.val
    assert const3.c0._link is const2.c0

    ds.unlink("const.c0")

    assert const3.c0._link is not const2.c0


@requires_fits
@requires_stk
def test_pha_case_6(ds_setup, ds_datadir):

    datadir = ds_datadir
    ls = '@' + '/'.join((datadir, 'pha.lis'))
    rmf1 = '/'.join((datadir, "acisf04938_000N002_r0043_rmf3.fits"))
    rmf2 = '/'.join((datadir, "acisf07867_000N001_r0002_rmf3.fits"))
    arf1 = '/'.join((datadir, "acisf04938_000N002_r0043_arf3.fits"))
    arf2 = '/'.join((datadir, "acisf07867_000N001_r0002_arf3.fits"))
    datastack.load_pha(ls)

    datastack.load_bkg_rmf([], rmf1)
    datastack.load_bkg_rmf([], rmf2)

    datastack.load_bkg_arf([], arf1)
    datastack.load_bkg_arf([], arf2)

    # Define background models
    bkg_arfs = datastack.get_bkg_arf([])
    bkg_scales = datastack.get_bkg_scale([])
    bkg_models = [ui.const1d.c1 * acis_bkg_model('acis7s'),
                  ui.const1d.c2 * acis_bkg_model('acis7s')]
    bkg_rsps = datastack.get_response([], bkg_id=1)
    for i in range(2):
        id_ = i + 1
        # Make the ARF spectral response flat.  This is required for using
        # the acis_bkg_model.
        bkg_arfs[i].specresp = bkg_arfs[i].specresp * 0 + 1.
        datastack.set_bkg_full_model(id_, bkg_rsps[i](bkg_models[i]))

    # Fit background
    datastack.notice(0.5, 8.)
    datastack.set_method("neldermead")
    datastack.set_stat("cash")

    datastack.thaw(c1.c0)
    datastack.thaw(c2.c0)
    datastack.fit_bkg()
    datastack.freeze(c1.c0)
    datastack.freeze(c2.c0)

    # Define source models
    rsps = datastack.get_response([])
    src_model = ui.powlaw1d.pow1
    src_models = [src_model,
                  src_model * ui.const1d.ratio_12]
    for i in range(2):
        id_ = i + 1
        datastack.set_full_model(id_, (rsps[i](src_models[i]) +
                                       bkg_scales[i] *
                                       bkg_rsps[i](bkg_models[i])))

    datastack.fit()


@requires_fits
@requires_stk
def test_query_case_7(ds_setup, ds_datadir):

    datadir = ds_datadir
    datastack.load_pha('@' +
                       '/'.join((datadir, 'pha.lis')))

    f = datastack.query_by_header_keyword('INSTRUME', 'ACIS')

    assert f == [1, 2]

    f = datastack.query_by_obsid(7867)

    assert f == [2]

    ds = datastack.DataStack()

    ds.load_pha('@' + '/'.join((datadir, 'pha.lis')))

    f = ds.query_by_obsid('4938')

    assert f == [3]

    f = datastack.query_by_obsid(ds, '7867')

    assert f == [4]


def test_default_instantiation(ds_setup):
    datastack.DataStack._default_instantiated = False
    ds = datastack.DataStack()
    assert ds._default_instance
    ds = datastack.DataStack()
    assert not ds._default_instance


def validate_show_stack(capsys, datadir, key1, key2):
    """Does show_stack generate the expected output?

    This is for test_show_stack and test_show_stack2.

    Parameters
    ----------
    capsys
       The capsys pytest object.
    datadir: str
       The directory containing the data files.
    key1 : str or None
       The keyword being used for the "MJD OBS" value: file 1.
       If None the value is missing.
    key2 : str or None
       The keyword being used for the "MJD OBS" value: file 2.
       If None the value is missing.

    """

    # clear out the current output
    captured = capsys.readouterr()

    datastack.show_stack()

    # Did we get the expected output?
    #
    # This depends on the serialization of numeric values, which can
    # vary with Python/system
    #
    if key1 is None:
        val1 = 'N/A'
        key1 = 'MJD-OBS'
    else:
        val1 = '53493.55477826'

    if key2 is None:
        val2 = 'N/A'
        key2 = 'MJD-OBS'
    else:
        val2 = '54374.009361043'

    l0 = '1: {}/acisf04938_000N002_r0043_pha3.fits OBS_ID: 4938 {}: {}'.format(datadir, key1, val1)
    l1 = '2: {}/acisf07867_000N001_r0002_pha3.fits OBS_ID: 7867 {}: {}'.format(datadir, key2, val2)

    captured = capsys.readouterr()
    lines = captured.out.split('\n')
    assert len(lines) == 3
    assert lines[0] == l0
    assert lines[1] == l1
    assert lines[2] == ''
    assert captured.err == ''


@requires_fits
@requires_stk
def test_show_stack(ds_setup, ds_datadir, capsys):
    """Test the show_stack handling: MJD_OBS
    """

    # These files use MJD_OBS in the header
    ls = '@' + '/'.join((ds_datadir, 'pha.lis'))
    datastack.load_pha(ls)

    validate_show_stack(capsys, ds_datadir, 'MJD_OBS', 'MJD_OBS')


@requires_fits
@requires_stk
def test_show_stack2(ds_setup, ds_datadir, capsys):
    """Test the show_stack handling: MJD-OBS

    This is test_show_stack but with the data files adjusted
    to have MJD-OBS rather than MJD_OBS keywords
    """

    # These files use MJD_OBS in the header
    ls = '@' + '/'.join((ds_datadir, 'pha.lis'))
    datastack.load_pha(ls)

    # Change to MJD-OBS
    #
    for idval in [1, 2]:
        d = datastack.get_data(idval)
        mjdobs = d.header['MJD_OBS']
        d.header['MJD-OBS'] = mjdobs
        del d.header['MJD_OBS']

    validate_show_stack(capsys, ds_datadir, 'MJD-OBS', 'MJD-OBS')


@requires_fits
@requires_stk
def test_show_stack3(ds_setup, ds_datadir, capsys):
    """Test the show_stack handling: mixed MJD_OBS and MJD-OBS
    """

    # These files use MJD_OBS in the header
    ls = '@' + '/'.join((ds_datadir, 'pha.lis'))
    datastack.load_pha(ls)

    # Change to MJD-OBS (second file only)
    #
    for idval in [2]:
        d = datastack.get_data(idval)
        mjdobs = d.header['MJD_OBS']
        d.header['MJD-OBS'] = mjdobs
        del d.header['MJD_OBS']

    validate_show_stack(capsys, ds_datadir, 'MJD_OBS', 'MJD-OBS')


@requires_fits
@requires_stk
def test_show_stack4(ds_setup, ds_datadir, capsys):
    """Test the show_stack handling: No MJD_OBS or MJD-OBS keyword
    """

    ls = '@' + '/'.join((ds_datadir, 'pha.lis'))
    datastack.load_pha(ls)

    # Remove the MJD-OBS keyword.
    #
    for idval in [1, 2]:
        d = datastack.get_data(idval)
        del d.header['MJD_OBS']

    validate_show_stack(capsys, ds_datadir, None, None)
