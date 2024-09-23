#
#  Copyright (C) 2017, 2018, 2020 - 2024
#  Smithsonian Astrophysical Observatory
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

"""
Should these tests be moved to test_astro_session.py?

This is almost a copy of sherpa/ui/tests/test_ui_unit.py, but
some of the answers are different. Is there a better way than
this duplication? Yes, we can parametrize by Session class
and run on both, to avoid duplication.
"""

import logging
import sys
import warnings

import numpy as np

import pytest

from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa.astro import io
from sherpa.astro import ui
import sherpa.models.basic
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, DataErr, \
    IdentifierErr, IOErr, ModelErr, StatErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.random import poisson_noise
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_region, requires_wcs, requires_xspec


def backend_is(name: str) -> bool:
    return io.backend.name == name


def check_table(hdu, colinfo):
    """Ensure the hdu contains the columns.

    Check that the hdu dictionary matches the expected setting.

    Parameters
    ----------
    hdu : dict
        A table HDU
    colinfo : dict
    """

    for k, v in colinfo.items():
        assert k in hdu
        assert hdu[k] == pytest.approx(v)

    assert len(hdu) == len(colinfo)


@pytest.mark.parametrize("as_string", [True, False])
def test_model_identifiers_set_globally(as_string, clean_astro_ui):
    """Check we create a global symbol for the models.

    See also the same test in
      sherpa/astro/ui/tests/test_astro_session.py
      sherpa/astro/ui/tests/test_astro_ui_import.py

    """

    # The "global" symbol table depends on what has been run before. We
    # could try and make sure that we are "clean", but this makes checking
    # what this test is doing hard to do, so we remove the symbols just
    # in case.
    #
    for name in ["mdl1", "mdl2"]:
        try:
            del sys.modules["__main__"].__dict__[name]
        except KeyError:
            pass

    ui.dataspace1d(1, 10, 1)

    for store in [globals(), locals(), sys.modules["__main__"].__dict__]:
        assert "mdl1" not in store
        assert "mdl2" not in store

    if as_string:
        ui.set_source("const1d.mdl1 + gauss1d.mdl2")
    else:
        # unlike the direct result in test_astro_session.py we can
        # call set_source without a string.
        #
        ui.set_source(ui.const1d.mdl1 + ui.gauss1d.mdl2)

    for store in [globals(), locals()]:
        assert "mdl1" not in store
        assert "mdl2" not in store

    assert "mdl1" in sys.modules["__main__"].__dict__
    assert "mdl2" in sys.modules["__main__"].__dict__

    assert isinstance(sys.modules["__main__"].__dict__["mdl1"],
                      sherpa.models.basic.Const1D)
    assert isinstance(sys.modules["__main__"].__dict__["mdl2"],
                      sherpa.models.basic.Gauss1D)


# This is part of #397
#
def test_list_samplers():
    """Ensure list_samplers returns a list."""

    samplers = ui.list_samplers()

    assert isinstance(samplers, list)
    assert len(samplers) > 0


def test_list_samplers_contents():
    """Are the expected values included"""

    # Test that the expected values exist in this list,
    # but do not enforce these are the only values. This is
    # a slightly-different return list to the non-astro version.
    #
    samplers = ui.list_samplers()
    for expected in ["mh", "metropolismh", "pragbayes", "fullbayes"]:
        assert expected in samplers


def test_all_has_no_repeated_elements():
    """Ensure __all__ does not contain repeated elements.

    It is not actually a problem if this happens, but it does
    indicate the possibility of confusion over what functionality
    the repeated symbol has (originally noticed with erf, which
    could be either sherpa.utils.erf or the ModelWrapper version
    of sherpa.models.basic.Erf). See
    https://github.com/sherpa/sherpa/issues/502
    """

    n1 = len(ui.__all__)
    n2 = len(set(ui.__all__))
    assert n1 == n2


def setup_data1d_fit():
    """Create a 1D dataset for fitting a gaussian-like profile."""

    x = np.linspace(2300, 2400, 51)

    cpt = ui.voigt1d("cpt")
    cpt.pos = 2345
    cpt.fwhm_g = 20
    cpt.fwhm_l = 15
    cpt.ampl = 480

    y = poisson_noise(cpt(x), rng=np.random.RandomState(72383))

    ui.load_arrays(1, x, y, ui.Data1D)
    ui.set_stat("leastsq")
    ui.set_method("simplex")
    ui.delete_model_component("cpt")


@pytest.mark.parametrize("model,stat,pars",
                         [(ui.gauss1d, 168.05369410194248,
                           [33.39014224934166, 2343.993186258643, 10.837604379839043]),
                          (ui.normgauss1d, 168.05369410194248,
                           [33.39014225250392, 2343.9931862492167, 385.19777742500384]),
                          (ui.lorentz1d, 175.07991080617057,
                           [27.355631135552674, 2342.94368338116, 504.57588886948827]),
                          (ui.pseudovoigt1d, 163.17653966109435,
                           [0.5253350907054826, 31.35910472670331, 2343.6786455810307, 440.7219934649222]),
                          (ui.voigt1d, 162.4321009471359,
                           [24.290859063445527, 12.35931466539915, 2343.751436403753, 437.0342780289447])])
def test_fit_profile(model, stat, pars, clean_astro_ui):
    """Regression test simple 1D fits"""

    setup_data1d_fit()
    ui.set_source(model("mdl"))
    ui.guess()
    ui.fit()

    assert ui.calc_stat() == pytest.approx(stat)
    assert np.asarray(mdl.thawedpars) == pytest.approx(np.asarray(pars))


@pytest.mark.parametrize("func", [ui.notice_id, ui.ignore_id])
def test_check_ids_not_none(func):
    """Check they error out when id is None"""

    with pytest.raises(ArgumentTypeErr,
                       match="'ids' must be an identifier or list of identifiers"):
        func(None)


@pytest.mark.parametrize("f", [[False, True], np.asarray([False, True])])
@pytest.mark.parametrize("bid", [None, 1])
def test_set_filter_mismatch(f, bid, clean_astro_ui):
    """Does set_filter error when there's a mismatch?
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)
    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and mask: 3 vs 2"):
        ui.set_filter(f, bkg_id=bid)


@pytest.mark.parametrize("f", [[False, True], np.asarray([False, True])])
@pytest.mark.parametrize("bid", [None, 1])
def test_set_filter_mismatch_with_filter(f, bid, clean_astro_ui):
    """Does set_filter error when there's a mismatch after a filter?

    test_set_filter_mismatch checks when .mask is a scalar,
    so now check when it's a NumPy array.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    ui.ignore(3, None)  # set the .mask attribute to an array

    with pytest.raises(DataErr,
                       match="size mismatch between data and filter: 3 vs 2"):
        ui.set_filter(f, bkg_id=bid)


@pytest.mark.parametrize("bid", [None, 1])
def test_get_syserror_missing(bid, clean_astro_ui):
    """Does get_syserror error out?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)
    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror(bkg_id=bid)


@requires_fits
@pytest.mark.parametrize("ascii,reader",
                         [(True, ui.unpack_ascii), (False, ui.unpack_table)])
@pytest.mark.parametrize("bid", [None, 1])
def test_save_filter_pha(ascii, reader, bid, tmp_path, clean_astro_ui):
    """Does save_filter work? [PHA]

    Does background copy the source filter by default?  No, hence we
    get different filters.

    """

    x = np.arange(1, 11, dtype=np.int16)
    ui.load_arrays(1, x, x, ui.DataPHA)
    bkg = ui.DataPHA("bkg", x, x)
    ui.set_bkg(bkg)

    ui.notice(2, 4)
    ui.notice(6, 8)

    outfile = tmp_path / "filter.dat"
    ui.save_filter(str(outfile), bkg_id=bid, ascii=ascii)

    expected = [0, 1, 1, 1, 0, 1, 1, 1, 0, 0]

    d = reader(str(outfile), colkeys=["X", "FILTER"])
    assert isinstance(d, ui.Data1D)
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(expected)
    assert d.staterror is None
    assert d.syserror is None


@requires_fits
@pytest.mark.parametrize("ascii,reader",
                         [(True, ui.unpack_ascii), (False, ui.unpack_table)])
@pytest.mark.parametrize("idval", [None, 2])
def test_save_filter_data1d(ascii, reader, idval, tmp_path, clean_astro_ui):
    """Does save_filter work? [Data1D]"""

    x = np.arange(1, 11, dtype=np.int16)
    ui.load_arrays(1, x, x)

    ui.notice(2, 4)
    ui.notice(6, 8)

    outfile = tmp_path / "filter.dat"
    ui.save_filter(str(outfile), ascii=ascii)

    expected = [0, 1, 1, 1, 0, 1, 1, 1, 0, 0]

    d = reader(str(outfile), colkeys=["X", "FILTER"])
    assert isinstance(d, ui.Data1D)
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(expected)
    assert d.staterror is None
    assert d.syserror is None


@requires_fits
@pytest.mark.parametrize("ascii,reader",
                         [(True, ui.unpack_ascii), (False, ui.unpack_table)])
@pytest.mark.parametrize("bid,expected",
                         [(None, [1, -1, 1, 1, 1, -1, -1, 1, 1, -1]),
                          (1, [1] * 10)])
def test_save_grouping_pha(ascii, reader, bid, expected, tmp_path, clean_astro_ui):
    """Does save_grouping work? [PHA]"""

    x = np.arange(1, 11, dtype=np.int16)
    ui.load_arrays(1, x, x, ui.DataPHA)
    bkg = ui.DataPHA("bkg", x, x)
    ui.set_bkg(bkg)

    ui.set_grouping([1, -1, 1, 1, 1, -1, -1, 1, 1, -1])
    ui.set_grouping([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], bkg_id=1)

    outfile = tmp_path / "grouping.dat"
    ui.save_grouping(str(outfile), bkg_id=bid, ascii=ascii)

    d = reader(str(outfile), colkeys=["CHANNEL", "GROUPS"])
    assert isinstance(d, ui.Data1D)
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(expected)
    assert d.staterror is None
    assert d.syserror is None


@requires_fits
@pytest.mark.parametrize("ascii,reader",
                         [(True, ui.unpack_ascii), (False, ui.unpack_table)])
@pytest.mark.parametrize("bid,expected",
                         [(None, [0] * 10),
                          (1, [5, 0, 0, 0, 0, 0, 0, 0, 0, 2])])
def test_save_quality_pha(ascii, reader, bid, expected, tmp_path, clean_astro_ui):
    """Does save_quality work? [PHA]"""

    x = np.arange(1, 11, dtype=np.int16)
    ui.load_arrays(1, x, x, ui.DataPHA)
    bkg = ui.DataPHA("bkg", x, x)
    ui.set_bkg(bkg)

    ui.set_quality([0] * 10)
    ui.set_quality([5, 0, 0, 0, 0, 0, 0, 0, 0, 2], bkg_id=1)

    outfile = tmp_path / "quality.dat"
    ui.save_quality(str(outfile), bkg_id=bid, ascii=ascii)

    d = reader(str(outfile), colkeys=["CHANNEL", "QUALITY"])
    assert isinstance(d, ui.Data1D)
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(expected)
    assert d.staterror is None
    assert d.syserror is None


@requires_fits
@pytest.mark.parametrize("ascii,reader",
                         [(True, ui.unpack_ascii), (False, ui.unpack_table)])
def test_save_table_data1d(ascii, reader, tmp_path, clean_astro_ui):
    """Does save_table work? [Data1D]"""

    x = np.asarray([1, 2, 5])
    y = np.asarray([2, -3, 5])
    ui.load_arrays(1, x, y)

    outfile = tmp_path / "table.dat"
    ui.save_table(str(outfile), ascii=ascii)

    d = reader(str(outfile))
    assert isinstance(d, ui.Data1D)
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(y)
    assert d.staterror is None
    assert d.syserror is None


@requires_fits
@pytest.mark.parametrize("ascii,reader",
                         [(True, ui.unpack_ascii), (False, ui.unpack_table)])
def test_save_table_pha(ascii, reader, tmp_path, clean_astro_ui):
    """Does save_table work? [PHA] See issue #47"""

    x = np.arange(1, 11, dtype=np.int16)
    ui.load_arrays(1, x, x, ui.DataPHA)
    ui.set_quality(1, [0] * 10)

    # Note that the background is ignored on output
    #
    bkg = ui.DataPHA("bkg", x, x)
    ui.set_bkg(bkg)
    ui.set_quality(1, [5, 0, 0, 0, 0, 0, 0, 0, 0, 2], bkg_id=1)

    outfile = tmp_path / "table.dat"
    ui.save_table(str(outfile), ascii=ascii)

    # The channel and counts arrays are written out. Things
    # we check
    #  - column names (CHANNEL and COUNTS)
    #  - the data can be read back in
    #  - we don't read in statistical errors
    #
    d = reader(str(outfile), colkeys=["CHANNEL", "COUNTS"])
    assert isinstance(d, ui.Data1D)
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(x)
    assert d.staterror is None
    assert d.syserror is None


@pytest.mark.parametrize("bid", [None, 1])
def test_save_filter_ignored(bid, tmp_path, clean_astro_ui):
    """Does save_filter error out if everything is masked?

    We should be able to write out the filter in this case,
    as it's easy (all False).
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    ui.ignore(None, None)

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(DataErr,
                       match="mask excludes all data"):
        ui.save_filter(str(outfile), bkg_id=bid)


@pytest.mark.parametrize("func,etype,emsg",
                         [(ui.save_filter, DataErr,
                           "data set '1' has no filter"),
                          (ui.save_grouping, DataErr,
                           "data set '1' does not specify grouping flags"),
                          (ui.save_quality,  DataErr,
                           "data set '1' does not specify quality flags"),
                          (ui.save_staterror,  StatErr,
                           "If you select chi2 as the statistic, all datasets must provide a staterror column"),
                          (ui.save_syserror,  DataErr,
                           "data set '1' does not specify systematic errors"),
                          (ui.save_error,  StatErr,
                           "If you select chi2 as the statistic, all datasets must provide a staterror column"),
                          (ui.save_source,  IdentifierErr,
                           r"source 1 has not been set, consider using set_source\(\) or set_model\(\)"),
                          (ui.save_model,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_resid,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_delchi,  IdentifierErr,
                           "model 1 has not been set")])
def test_save_xxx_nodata(func, etype, emsg, tmp_path, clean_astro_ui):
    """Does save_xxx error out if there's no data to save? DataPHA
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    # Without this we can't check save_staterror or save_error
    ui.set_stat("chi2")

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(etype,
                       match=emsg):
        func(str(outfile))


@requires_fits
def test_save_image_nodata(tmp_path, clean_astro_ui):
    """Does save_image error out if there's no data to save? DataPHA

    Unlike the other calls, this requires a FITS library
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(IOErr,
                       match="data set '' does not contain an image"):
        ui.save_image(str(outfile))


@pytest.mark.parametrize("func,etype,emsg",
                         [(ui.save_filter, DataErr,
                          "data set '1' has no filter"),
                          (ui.save_grouping, DataErr,
                           "data set '1' does not specify grouping flags"),
                          (ui.save_quality,  DataErr,
                           "data set '1' does not specify quality flags"),
                          (ui.save_source,  ModelErr,
                           "background model 1 for data set 1 has not been set"),
                          (ui.save_staterror,  StatErr,
                           "If you select chi2 as the statistic, all datasets must provide a staterror column"),
                          (ui.save_syserror,  DataErr,
                           "data set '1' does not specify systematic errors"),
                          (ui.save_error,  StatErr,
                           "If you select chi2 as the statistic, all datasets must provide a staterror column"),
                          (ui.save_model,  ModelErr,
                           "background model 1 for data set 1 has not been set"),
                          (ui.save_resid,  ModelErr,
                           "background model 1 for data set 1 has not been set"),
                          (ui.save_delchi,  ModelErr,
                           "background model 1 for data set 1 has not been set")])
def test_save_xxx_bkg_nodata(func, etype, emsg, tmp_path, clean_astro_ui):
    """Does save_xxx error out if there's no data to save? DataPHA + bkg

    Note that save_image does not support a bkg_id parameter so is
    not tested here.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    # Without this we can't check save_staterror or save_error
    ui.set_stat("chi2")

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(etype,
                       match=emsg):
        func(str(outfile), bkg_id=1)


@pytest.mark.parametrize("func,etype,emsg",
                         [(ui.save_filter, DataErr,
                           "data set '1' has no filter"),
                          (ui.save_grouping, ArgumentErr,
                           "data set 1 does not contain PHA data"),
                          (ui.save_quality,  ArgumentErr,
                           "data set 1 does not contain PHA data"),
                          (ui.save_source,  IdentifierErr,
                           r"source 1 has not been set, consider using set_source\(\) or set_model\(\)"),
                          (ui.save_model,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_resid,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_delchi,  IdentifierErr,
                           "model 1 has not been set")])
def test_save_xxx_data1d_nodata(func, etype, emsg, tmp_path, clean_astro_ui):
    """Does save_xxx error out if there's no data to save? Data1D
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(etype,
                       match=emsg):
        func(str(outfile))


@requires_fits
def test_save_image_data1d_nodata(tmp_path, clean_astro_ui):
    """Does save_image error out if there's no data to save? Data1D

    Unlike the other cases we need a FITS backend.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(IOErr,
                       match="data set '' does not contain an image"):
        ui.save_image(str(outfile))


@pytest.mark.parametrize("savefunc", [ui.save_data,
                                      ui.save_filter,
                                      ui.save_grouping,
                                      ui.save_quality,
                                      ui.save_image,
                                      ui.save_resid])
def test_save_xxx_not_a_string(savefunc, tmp_path, clean_astro_ui):
    """Does save_xxx error out if not given a filename?

    There are times it would be nice to give a non-string filename
    (pathlib.Path or StringIO), but this is currently not supported.

    At the moment this check is generic for some save commands,
    so just worry about a restricted set of commands.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)
    out = tmp_path / "data.dat"

    with pytest.raises(ArgumentTypeErr,
                       match="'filename' must be a string"):
        savefunc(out)


def test_save_delchi_image_fails(tmp_path, clean_astro_ui):
    """save_delchi doesn't work for image datasets

    Only test out DataIMG
    """

    y, x = np.mgrid[10:12, 20:23]
    x = x.flatten()
    y = y.flatten()
    z = (x - 11)**2 + (y - 21)**2
    ui.load_arrays(1, x, y, z, (2, 3), ui.DataIMG)

    ui.set_source(ui.const2d.cmdl)

    out = tmp_path / "does-not-exist"
    with pytest.raises(AttributeError,
                       match=r"^save_delchi\(\) does not apply for images$"):
        ui.save_delchi(str(out))


def check_clobber(outpath, func):
    """Does calling func raise a clobber error and not change the contents

    Parameters
    ----------
    outpath : pathlib.Path instance
    func : function reference
        Called with a single string argument.
    """

    old = outpath.read_text()

    # check it clobbers
    with pytest.raises(IOErr,
                       match="^file '.*' exists and clobber is not set$"):
        func(str(outpath))

    # Check the file hasn't changed
    new = outpath.read_text()
    assert new == old


@requires_fits
def test_save_data_data1d_no_clobber(tmp_path, clean_astro_ui):
    """save_data: does clobber=False work? Data1D"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)

    out.write_text("some text")
    check_clobber(out, ui.save_data)


@requires_fits
def test_save_data_datapha_no_clobber(tmp_path, clean_astro_ui):
    """save_data: does clobber=False work? DataPHA"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    out = tmp_path / "data.dat"
    out.write_text("some text")
    check_clobber(out, ui.save_data)


@requires_fits
def test_save_pha_no_clobber(tmp_path, clean_astro_ui):
    """save_pha: does clobber=False work?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    out = tmp_path / "data.dat"

    out.write_text("some text")
    check_clobber(out, ui.save_data)


@requires_fits
@pytest.mark.parametrize("writer", [ui.save_data, ui.save_image])
def test_save_image_no_clobber(writer, tmp_path, clean_astro_ui):
    """save_image: does clobber=False work?"""

    y, x = np.mgrid[10:12, 20:23]
    x = x.flatten()
    y = y.flatten()
    z = (x - 11)**2 + (y - 21)**2
    ui.load_arrays(1, x, y, z, (2, 3), ui.DataIMG)

    out = tmp_path / "data.dat"

    out.write_text("some text")
    check_clobber(out, writer)


def check_output(out, colnames, rows):
    """Is the output as expected?

    The format of the ASCII output depends on the backend, but
    does so in a 'principled' manner, so it can be handled
    here.

    Parameters
    ----------
    out : str
        The actual output
    colnames : sequence of str
        The column names (assumed to be upper case).
    rows : sequence of sequence of numbers
        The row data. We convert the rows from the input file to
        numeric values to avoid issues with the formatting (e.g. is 5
        or 5.0 written out). There is no check that each row has the
        same number of columns.

    """

    lines = out.split("\n")
    assert len(lines) > 2

    cols = " ".join(colnames)
    if backend_is("crates"):
        assert lines[0] == "#TEXT/SIMPLE"
        assert lines[1] == f"# {cols}"
        lines = lines[2:]
    elif backend_is("pyfits"):
        assert lines[0] == f"# {cols}"
        lines = lines[1:]
    else:
        raise RuntimeError(f"UNKNOWN I/O BACKEND: {io.backend}")

    assert lines[-1] == ""
    lines = lines[:-1]

    assert len(lines) == len(rows)
    for l, r in zip(lines, rows):
        got = [float(t) for t in l.split()]
        assert got == pytest.approx(r)


@requires_fits
def test_save_data_data1d_clobber(tmp_path, clean_astro_ui):
    """save_data: does clobber=True work?"""

    ui.load_arrays(1, [1], [5], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)

    out.write_text("some text")

    ui.save_data(outfile, clobber=True)
    cts = out.read_text()
    check_output(cts, ["X", "Y"], [[1, 5]])


@requires_fits
def test_save_data_data1d(tmp_path, clean_astro_ui):
    """Does save_data work for Data1D?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ["X", "Y"], [[1, 5], [2, 4], [3, 3]])


@requires_fits
def test_save_data_data1d_fits(tmp_path, clean_astro_ui):
    """Does save_data work for Data1D? FITS"""

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile, ascii=False)

    _, blocks, _ = read_table_blocks(outfile)
    assert len(blocks) == 2
    check_table(blocks[2],
                {"X": [1, 2, 3],
                 "Y": [5, 4, 3]})


@requires_fits
def test_save_data_data1dint(tmp_path, clean_astro_ui):
    """Does save_data work for Data1DInt?"""

    ui.load_arrays(1, [1, 2, 4], [2, 3, 5], [5, 4, 3], ui.Data1DInt)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ["XLO", "XHI", "Y"],
                 [[1, 2, 5], [2, 3, 4], [4, 5, 3]])


@requires_fits
def test_save_data_data1dint_fits(tmp_path, clean_astro_ui):
    """Does save_data work for Data1DInt? FITS"""

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [1, 2, 4], [2, 3, 5], [5, 4, 3], ui.Data1DInt)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile, ascii=False)

    _, blocks, _ = read_table_blocks(outfile)
    assert len(blocks) == 2
    check_table(blocks[2],
                {"XLO": [1, 2, 4],
                 "XHI": [2, 3, 5],
                 "Y": [5, 4, 3]})


@requires_fits
def test_save_data_data2d(tmp_path, clean_astro_ui):
    """Does save_data work for Data2D?"""

    y, x = np.mgrid[20:22, 10:13]
    x = x.flatten()
    y = y.flatten()
    z = (x - 15)**2 + (y - 12)**2
    ui.load_arrays(1, x, y, z, (2, 3), ui.Data2D)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()

    # the output depends on the backend, and neither seems ideal
    #
    if backend_is("crates"):
        expected = ["#TEXT/SIMPLE", "# X0 X1 Y SHAPE"]
        s = [2, 3, 0, 0, 0, 0]
        for xi, yi, zi, si in zip(x, y, z, s):
            expected.append(f"{xi} {yi} {zi} {si}")

        expected = "\n".join(expected) + "\n"

    elif backend_is("pyfits"):
        expected = "\n".join([str(zz) for zz in z]) + "\n"
    else:
        raise RuntimeError(f"UNKNOWN I/O BACKEND: {io.backend}")

    assert cts == expected


@requires_fits
def test_save_data_data2d_fits(tmp_path, clean_astro_ui):
    """Does save_data work for Data2D? FITS"""

    from sherpa.astro.io import read_image

    y, x = np.mgrid[10:12, 20:23]
    x = x.flatten()
    y = y.flatten()
    z = (x - 12)**2 + (y - 22)**2
    ui.load_arrays(1, x, y, z, (2, 3), ui.Data2D)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile, ascii=False)

    ans = read_image(outfile)
    assert ans.shape == (2, 3)  # Is this correct?

    yl, xl = np.mgrid[1:3, 1:4]
    xl = xl.flatten()
    yl = yl.flatten()
    assert ans.x0 == pytest.approx(xl)
    assert ans.x1 == pytest.approx(yl)
    assert ans.y == pytest.approx(z)


@requires_fits
def test_save_data_dataimg(tmp_path, clean_astro_ui):
    """Does save_data work for DataIMG? ASCII"""

    # Can not write out an ASCII image with crates
    if backend_is("crates"):
        pytest.skip("ASCII not supported for images with pycrates")

    y, x = np.mgrid[0:2, 0:3]
    x = x.flatten()
    y = y.flatten()
    z = (x - 1)**2 + y**2
    ui.load_arrays(1, x, y, z, (2, 3), ui.DataIMG)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    expected = "\n".join([str(zz) for zz in z]) + "\n"
    assert cts == expected


@requires_fits
def test_save_data_dataimg_fits(tmp_path, clean_astro_ui):
    """Does save_data work for DataIMG? FITS"""

    from sherpa.astro.io import read_image

    y, x = np.mgrid[10:12, 20:23]
    x = x.flatten()
    y = y.flatten()
    z = (x - 11)**2 + (y - 21)**2
    ui.load_arrays(1, x, y, z, (2, 3), ui.DataIMG)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile, ascii=False)

    ans = read_image(outfile)
    assert ans.shape == (2, 3)  # Is this correct?

    yl, xl = np.mgrid[1:3, 1:4]
    xl = xl.flatten()
    yl = yl.flatten()
    assert ans.x0 == pytest.approx(xl)
    assert ans.x1 == pytest.approx(yl)
    assert ans.y == pytest.approx(z)


@requires_data
@requires_fits
@requires_wcs
@pytest.mark.parametrize("writer", [ui.save_image, ui.save_data])
def test_save_image_dataimg_fits_wcs(writer, make_data_path, tmp_path, clean_astro_ui):
    """Does save_image work for a FITS image with WCS

    We also check the input file, is read in, for fun.

    dmlist reports that this file has

    Block    1: EVENTS_IMAGE                   Image      Int2(316x313)

    Physical Axis Transforms for Image Block EVENTS_IMAGE

    Group# Axis#
       1   1,2    sky(x) = (+2986.8401) +(+1.0)* ((#1)-(+0.50))
                     (y)   (+4362.6602)  (+1.0)  ((#2) (+0.50))

    World Coordinate Axis Transforms for Image Block EVENTS_IMAGE

    Group# Axis#
       1   1,2    EQPOS(RA ) = (+149.8853)[deg] +TAN[(-0.000136667)* (sky(x)-(+4096.50))]
                       (DEC)   (+2.6079  )           (+0.000136667)  (   (y) (+4096.50))

    """

    from sherpa.astro.io import read_image
    from sherpa.astro.io.wcs import WCS

    # It looks like we don't write out the header info, at least for
    # pyfits. Probably a bug.
    #
    def check_data(d, header=True):
        assert isinstance(d, ui.DataIMG)
        assert d.shape == (313, 316)
        assert d.y.shape == (313 * 316, )
        assert d.y.sum() == 965
        assert d.y.max() == 3

        hdr = d.header
        assert isinstance(hdr, dict)
        if header:
            assert hdr["OBJECT"] == "CSC"
            assert hdr["ONTIME"] == pytest.approx(18220.799932122)
            assert hdr["TSTOP"] == pytest.approx(280770182.77747)

        assert isinstance(d.sky, WCS)
        assert d.sky.name == "physical"
        assert d.sky.type == "LINEAR"
        assert d.sky.crval == pytest.approx([2986.84008789, 4362.66015625])
        assert d.sky.crpix == pytest.approx([0.5, 0.5])
        assert d.sky.cdelt == pytest.approx([1, 1])

        assert isinstance(d.eqpos, WCS)
        assert d.eqpos.name == "world"
        assert d.eqpos.type == "WCS"
        assert d.eqpos.epoch == pytest.approx(2000)
        assert d.eqpos.equinox == pytest.approx(2000)
        assert d.eqpos.crval == pytest.approx([149.88533198,   2.60794887])
        assert d.eqpos.crpix == pytest.approx([4096.5, 4096.5])
        assert d.eqpos.cdelt * 3600 == pytest.approx([-0.492, 0.492])


    infile = make_data_path("acisf08478_000N001_r0043_regevt3_srcimg.fits")
    ui.load_data(2, infile)

    dorig = ui.get_data(2)
    check_data(dorig)

    out = tmp_path / "data.dat"
    outfile = str(out)
    writer(2, outfile, ascii=False)

    ans = read_image(outfile)
    check_data(ans, header=False)


@requires_fits
def test_save_data_datapha(tmp_path, clean_astro_ui):
    """Does save_data work for DataPHA?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ["CHANNEL", "COUNTS"],
                 [[1, 5], [2, 4], [3, 3]])


@requires_fits
@pytest.mark.parametrize("idval", [None, 2])
def test_save_data_datapha_bkg(idval, tmp_path, clean_astro_ui):
    """Does save_data work for DataPHA [bkg]?"""

    idarg = 1 if idval is None else idval
    ui.load_arrays(idarg, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    ui.set_bkg(idarg, ui.DataPHA("bg", [1, 2, 3], [0, 2, 1]))

    out = tmp_path / "data.dat"
    outfile = str(out)
    if idval is None:
        ui.save_data(outfile, bkg_id=1)
    else:
        ui.save_data(idval, outfile, bkg_id=1)

    cts = out.read_text()
    check_output(cts, ["CHANNEL", "COUNTS"],
                 [[1, 0], [2, 2], [3, 1]])


@requires_fits
def test_save_pha(tmp_path, clean_astro_ui):
    """Does save_pha work for DataPHA?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_pha(outfile)

    # just check the first line; the output may depend on the FITS backend
    cts = out.read_text()[:80].split()
    assert cts[0] == "SIMPLE"
    assert cts[1] == "="
    assert cts[2] == "T"


@requires_fits
@pytest.mark.parametrize("idval", [None, 2])
def test_save_pha_bkg(idval, tmp_path, clean_astro_ui):
    """Does save_pha work for DataPHA (bkg, ascii)?"""

    idarg = 1 if idval is None else idval
    ui.load_arrays(idarg, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    ui.set_exposure(idarg, 10)

    bkg1 = ui.DataPHA("bkg1", [1, 2, 3], [2, 1, 4])
    ui.set_bkg(idarg, bkg1)
    ui.set_exposure(idarg, 20, bkg_id=1)

    bkg2 = ui.DataPHA("bkg2", [1, 2, 3], [1, 0, 3])
    ui.set_bkg(idarg, bkg2, bkg_id=2)
    ui.set_exposure(idarg, 40, bkg_id=2)

    out = tmp_path / "data.dat"
    outfile = str(out)
    if idval is None:
        ui.save_pha(outfile, bkg_id=2, ascii=True)
    else:
        ui.save_pha(idval, outfile, bkg_id=2, ascii=True)

    # The exact format depends on the I/O backend so just check we can
    # read it back in using Sherpa I/O routines.
    #
    data = io.read_ascii(outfile, colkeys=["CHANNEL", "COUNTS"])
    assert data.x == pytest.approx([1, 2, 3])
    assert data.y == pytest.approx([1, 0, 3])


@requires_fits
@pytest.mark.parametrize("savefunc,mtype", [(ui.save_source, "SOURCE"),
                                            (ui.save_model, "MODEL")])
def test_save_model_ascii(savefunc, mtype, clean_astro_ui, tmp_path):
    """Can we write out data for save_source/model? Data1D and ASCII

    As this is not a PHA dataset, the two should be the same bar the
    header line.
    """

    ui.load_arrays(1, [1, 2], [5, 10])
    ui.set_source(ui.polynom1d.cmdl)
    cmdl.c0 = -5
    cmdl.c1 = 10

    out = tmp_path / "model.dat"
    savefunc(str(out), ascii=True)

    cts = out.read_text()
    check_output(cts, ["X", mtype], [[1, 5], [2, 15]])


@requires_fits
def test_save_source_pha_ascii(clean_astro_ui, tmp_path):
    """Can we write out data for save_source? DataPHA and ASCII"""

    ui.load_arrays(1, [1, 2], [5, 10], ui.DataPHA)

    # we need a response
    egrid = np.asarray([0.1, 0.2, 0.4])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf(rmf)

    yarf = np.asarray([10, 20])
    arf = create_arf(elo, ehi, yarf)
    ui.set_arf(arf)

    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 2

    out = tmp_path / "model.dat"
    ui.save_source(str(out), ascii=True)

    cts = out.read_text()
    check_output(cts, ["XLO", "XHI", "SOURCE"],
                 [[0.1, 0.2, 2], [0.2, 0.4, 2]])


@requires_fits
def test_save_model_pha_ascii(clean_astro_ui, tmp_path):
    """Can we write out data for save_model? DataPHA and ASCII"""

    ui.load_arrays(1, [1, 2], [5, 10], ui.DataPHA)

    # we need a response
    egrid = np.asarray([0.1, 0.2, 0.4])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf(rmf)

    yarf = np.asarray([10, 20])
    arf = create_arf(elo, ehi, yarf)
    ui.set_arf(arf)

    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 2

    out = tmp_path / "model.dat"
    ui.save_model(str(out), ascii=True)

    cts = out.read_text()
    check_output(cts, ["XLO", "XHI", "MODEL"],
                 [[0.1, 0.2, 20], [0.2, 0.4, 40]])


@requires_fits
def test_save_model_pha_bkg_ascii(clean_astro_ui, tmp_path):
    """Can we write out data for save_model (background)? DataPHA and ASCII"""

    data = ui.DataPHA("ex", [1, 2], [5, 10])

    # we need a response
    egrid = np.asarray([0.1, 0.2, 0.4])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    data.set_rmf(rmf)

    yarf = np.asarray([10, 20])
    arf = create_arf(elo, ehi, yarf)
    data.set_arf(arf)

    # It is not 100% clear when the code picks up the response
    # for the background data from the source region. I had
    # thought this would do it, but the plot output below
    # suggests the RMF has not been used.
    #
    bkg = ui.DataPHA("bex", [1, 2], [2, 3])
    data.set_background(bkg)

    ui.set_data(1, data)

    ui.set_source(ui.const1d.cmdl)
    ui.set_bkg_source(ui.const1d.bmdl)
    cmdl.c0 = 2
    bmdl.c0 = 10

    out = tmp_path / "model.dat"
    ui.save_model(str(out), ascii=True, bkg_id=1)

    # As mentioned above, I had expected this to be in energy units
    # but it is not. So for now treat this as a regression test.
    #
    cts = out.read_text()
    print(cts)
    check_output(cts, ["XLO", "XHI", "MODEL"],
                 [[1, 2, 10], [2, 3, 40]])  # note: channels!


@requires_fits
@pytest.mark.parametrize("savefunc,mtype", [(ui.save_source, "SOURCE"),
                                            (ui.save_model, "MODEL")])
def test_save_model_fits(savefunc, mtype, clean_astro_ui, tmp_path):
    """Can we write out data for save_source/model? Data1D and FITS

    As this is not a PHA dataset, the two should be the same bar the
    header line.
    """

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [1, 2], [5, 10])
    ui.set_source(ui.polynom1d.cmdl)
    cmdl.c0 = -5
    cmdl.c1 = 10

    out = tmp_path / "model.dat"
    outfile = str(out)
    savefunc(outfile)

    _, blocks, _ = read_table_blocks(outfile)
    assert len(blocks) == 2
    check_table(blocks[2],
                {"X": [1, 2],
                 mtype: [5, 15]})


@requires_fits
def test_save_source_pha_fits(clean_astro_ui, tmp_path):
    """Can we write out data for save_source? DataPHA and FITS
    """

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [1, 2], [5, 10], ui.DataPHA)

    # we need a response
    egrid = np.asarray([0.1, 0.2, 0.4])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf(rmf)

    yarf = np.asarray([10, 20])
    arf = create_arf(elo, ehi, yarf)
    ui.set_arf(arf)

    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 2

    out = tmp_path / "model.dat"
    outfile = str(out)
    ui.save_source(outfile)

    _, blocks, _ = read_table_blocks(outfile)
    assert len(blocks) == 2
    check_table(blocks[2],
                {"XLO": [0.1, 0.2],
                 "XHI": [0.2, 0.4],
                 "SOURCE": [2, 2]})


@requires_fits
def test_save_model_pha_fits(clean_astro_ui, tmp_path):
    """Can we write out data for save_model? DataPHA and FITS
    """

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [1, 2], [5, 10], ui.DataPHA)

    # we need a response
    egrid = np.asarray([0.1, 0.2, 0.4])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf(rmf)

    yarf = np.asarray([10, 20])
    arf = create_arf(elo, ehi, yarf)
    ui.set_arf(arf)

    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 2

    out = tmp_path / "model.dat"
    outfile = str(out)
    ui.save_model(outfile)

    _, blocks, _ = read_table_blocks(outfile)
    assert len(blocks) == 2
    check_table(blocks[2],
                {"XLO": [0.1, 0.2],
                 "XHI": [0.2, 0.4],
                 "MODEL": [20, 40]})


@requires_fits
def test_save_resid_data1d(clean_astro_ui, tmp_path):
    """Residual, Data1D, ASCII"""

    ui.load_arrays(1, [100, 200], [20, 230], ui.Data1D)
    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 220

    out = tmp_path / "resid.out"
    outfile = str(out)
    ui.save_resid(outfile, ascii=True)

    cts = out.read_text()
    check_output(cts, ["X", "RESID"], [[100, -200], [200, 10]])


@requires_fits
def test_save_resid_data1d_fits(clean_astro_ui, tmp_path):
    """Residual, Data1D, FITS"""

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [100, 200], [20, 230], ui.Data1D)
    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 220

    out = tmp_path / "resid.out"
    outfile = str(out)
    ui.save_resid(outfile)

    _, blocks, _ = read_table_blocks(outfile)
    assert len(blocks) == 2
    check_table(blocks[2],
                {"X": [100, 200],
                 "RESID": [-200, 10]})


@requires_fits
def test_save_resid_datapha(clean_astro_ui, tmp_path):
    """Residual, DataPHA, ASCII"""

    ui.load_arrays(1, [1, 2], [5, 10], ui.DataPHA)

    # we need a response
    egrid = np.asarray([0.1, 0.2, 0.4])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf(rmf)

    yarf = np.asarray([10, 20])
    arf = create_arf(elo, ehi, yarf)
    ui.set_arf(arf)

    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 2

    out = tmp_path / "resid.out"
    outfile = str(out)
    ui.save_resid(outfile, ascii=True)

    cts = out.read_text()
    check_output(cts, ["X", "RESID"], [[0.15, 30], [0.3, 10]])


@requires_fits
def test_save_resid_datapha_fits(clean_astro_ui, tmp_path):
    """Residual, DataPHA, FITS"""

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [1, 2], [5, 10], ui.DataPHA)

    # we need a response
    egrid = np.asarray([0.1, 0.2, 0.4])
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf(rmf)

    yarf = np.asarray([10, 20])
    arf = create_arf(elo, ehi, yarf)
    ui.set_arf(arf)

    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 2

    out = tmp_path / "resid.out"
    outfile = str(out)
    ui.save_resid(outfile)

    _, blocks, _ = read_table_blocks(outfile)
    assert len(blocks) == 2
    check_table(blocks[2],
                {"X": [0.15, 0.3],
                 "RESID": [30, 10]})


@requires_fits
def test_save_resid_dataimg(clean_astro_ui, tmp_path):
    """Residual, DataIMG, ASCII"""

    # Can not write out an ASCII image with crates
    if backend_is("crates"):
        pytest.skip("ASCII not supported for images with pycrates")

    y, x = np.mgrid[10:12, 20:23]
    x = x.flatten()
    y = y.flatten()
    z = (x - 11)**2 + (y - 21)**2
    ui.load_arrays(1, x, y, z, (2, 3), ui.DataIMG)

    ui.set_source(1, ui.const2d.cmdl)
    cmdl.c0 = 10

    out = tmp_path / "resid"
    outfile = str(out)
    ui.save_resid(outfile, ascii=True)

    cts = out.read_text()
    expected = "\n".join([str(zz - 10) for zz in z]) + "\n"
    assert cts == expected


@requires_fits
def test_save_resid_dataimg_fits(clean_astro_ui, tmp_path):
    """Residual, DataIMG, FITS"""

    from sherpa.astro.io import read_image

    y, x = np.mgrid[10:12, 20:23]
    x = x.flatten()
    y = y.flatten()
    z = (x - 11)**2 + (y - 21)**2
    ui.load_arrays(1, x, y, z, (2, 3), ui.DataIMG)

    ui.set_source(1, ui.const2d.cmdl)
    cmdl.c0 = 100

    out = tmp_path / "resid"
    outfile = str(out)
    ui.save_resid(outfile)

    ans = read_image(outfile)
    assert ans.shape == (2, 3)  # Is this correct?

    yl, xl = np.mgrid[1:3, 1:4]
    xl = xl.flatten()
    yl = yl.flatten()
    assert ans.x0 == pytest.approx(xl)
    assert ans.x1 == pytest.approx(yl)
    assert ans.y == pytest.approx(z - 100)


def test_delete_bkg_model(clean_astro_ui):
    """Check we can delete a background model"""

    channels = np.arange(1, 5)
    counts = np.zeros(channels.size)
    d = ui.DataPHA("src", channels, counts)
    b = ui.DataPHA("bkg", channels, counts)
    d.set_background(b, id=2)
    ui.set_data(d)

    ui.set_bkg_source(ui.const1d.bmdl + ui.gauss1d.gmdl, bkg_id=2)
    assert ui.get_bkg_source(bkg_id=2) is not None

    ui.delete_bkg_model(bkg_id=2)

    # Expression has been removed
    #
    with pytest.raises(ModelErr,
                       match="background model 2 for data set 1 has not been set"):
        ui.get_bkg_source(bkg_id=2)

    # components still exist
    mdls = ui.list_model_components()
    assert set(mdls) == set(["bmdl", "gmdl"])


def test_default_background_issue(clean_astro_ui):
    """Test issue #943"""

    ui.set_default_id("x")

    # use least-square as we don't really care about the fit
    ui.set_stat("leastsq")

    ui.load_arrays("x", [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    arf = ui.create_arf(np.asarray([0.1, 0.2, 0.3]), np.asarray([0.2, 0.3, 0.4]))
    bkg.set_arf(arf)
    ui.set_bkg(bkg)

    ui.set_bkg_source(ui.const1d.mdl2)

    # Ensure we can fit the background model. Prior to #943 being
    # fixed the fit_bkg call would error out.
    #
    ui.fit_bkg()
    assert mdl2.c0.val == pytest.approx(2 / 3 / 0.1)


def test_show_bkg_model_issue943(clean_astro_ui):
    """Test issue #943

    We do not check that show_bkg_model is creating anything
    useful, just that it can be called.

    See https://github.com/sherpa/sherpa/issues/943#issuecomment-696119982
    """

    ui.set_default_id("x")

    ui.load_arrays("x", [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    arf = ui.create_arf(np.asarray([0.1, 0.2, 0.3]), np.asarray([0.2, 0.3, 0.4]))
    bkg.set_arf(arf)
    ui.set_bkg(bkg)

    ui.set_bkg_source(ui.const1d.mdl2)
    ui.show_bkg_model()


def test_default_background_issue_fit(clean_astro_ui):
    """Test issue #943 with fit

    See https://github.com/sherpa/sherpa/issues/943#issuecomment-696119982
    """

    ui.set_default_id("x")

    # use least-square as we don't really care about the fit
    ui.set_stat("leastsq")

    ui.load_arrays("x", [1, 2, 3, 4], [5, 4, 3, 4], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3, 4]), [1, 1, 0, 1])
    arf = ui.create_arf(np.asarray([0.1, 0.2, 0.3, 0.4]),
                        np.asarray([0.2, 0.3, 0.4, 0.5]))
    ui.set_arf(arf)
    bkg.set_arf(arf)
    ui.set_bkg(bkg)

    # The model being fitted is a constant to 1,1,0,1 for
    # the background, so that should be 0.75 / 0.1 (as the
    # bin width is constant), and for the source it is
    # 5,4,3,4 - <0.75> [here ignoring the bin-width],
    # so [4.25,3.25,2.25,3.25] -> 13 / 4 -> 3.25
    #
    ui.set_source(ui.const1d.mdl1)
    ui.set_bkg_source(ui.const1d.mdl2)

    # Prior to #943 this would give a confusing error.
    #
    ui.fit()
    assert mdl1.c0.val == pytest.approx(3.25 / 0.1)
    assert mdl2.c0.val == pytest.approx(0.75 / 0.1)


def test_bkg_id_get_bkg_source(clean_astro_ui):
    """Check the error message when the background model has not been set (issue #943)"""

    ui.set_default_id("x")

    ui.load_arrays("x", [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    with pytest.raises(ModelErr,
                       match="background model 1 for data set x has not been set"):
        ui.get_bkg_source()


def test_fix_background_id_error_checks1(clean_astro_ui):
    """Check error handling of background id"""

    ui.load_arrays(2, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(2, bkg)

    with pytest.raises(ArgumentTypeErr,
                       match="identifiers must be integers or strings"):
        ui.get_bkg_source(id=2, bkg_id=bkg)


def test_fix_background_id_error_checks2(clean_astro_ui):
    """Check error handling of background id"""

    ui.load_arrays(2, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(2, bkg)

    with pytest.raises(IdentifierErr,
                       match="identifier 'bkg' is a reserved word"):
        ui.get_bkg_source(id=2, bkg_id="bkg")


@pytest.mark.parametrize("idval", [1, "x"])
def test_delete_bkg_model_with_bkgid(idval, clean_astro_ui):
    """Check we call delete_bkg_model with non-default bkg_id"""

    ui.load_arrays(idval, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(idval, bkg, bkg_id=2)

    ui.set_bkg_source(idval, ui.const1d.bmdl, bkg_id=2)
    assert ui.list_model_components() == ["bmdl"]
    assert ui.get_bkg_source(idval, 2).name == "const1d.bmdl"

    ui.delete_bkg_model(idval, bkg_id=2)
    assert ui.list_model_components() == ["bmdl"]

    with pytest.raises(ModelErr,
                       match=f"background model 2 for data set {idval} has not been set"):
        ui.get_bkg_source(idval, 2)


@requires_fits
@pytest.mark.parametrize("loadfunc", [ui.load_grouping, ui.load_quality])
def test_load_xxx_no_data(loadfunc, clean_astro_ui, tmp_path):
    """What happens when there's no data??"""

    path = tmp_path / "data"
    path.write_text("1\n0\n0\n")

    with pytest.raises(IdentifierErr,
                       match="data set 2 has not been set"):
        loadfunc(2, str(path))


@requires_fits
@pytest.mark.parametrize("loadfunc", [ui.load_grouping, ui.load_quality])
def test_load_xxx_not_pha(loadfunc, clean_astro_ui, tmp_path):
    """What happens with a non-PHA dataset?"""

    ui.load_arrays(2, [1, 2, 3], [2, 4, 9])

    path = tmp_path / "data"
    path.write_text("1\n0\n0\n")

    with pytest.raises(ArgumentErr,
                       match="data set 2 does not contain PHA data"):
        loadfunc(2, str(path))


@requires_fits
@pytest.mark.parametrize("idval", [None, 1, "xx"])
def test_load_grouping(idval, clean_astro_ui, tmp_path, caplog):
    """Simple grouping check"""

    x = [1, 2, 3]
    y = [0, 4, 3]
    if idval is None:
        ui.load_arrays(1, x, y, ui.DataPHA)
        idstr = "1"
    else:
        ui.load_arrays(idval, x, y, ui.DataPHA)
        idstr = str(idval)

    path = tmp_path / "group.dat"
    path.write_text("1\n-1\n1")

    data = ui.get_data(idval)
    assert data.grouping is None

    assert len(caplog.records) == 0
    if idval is None:
        ui.load_grouping(str(path))
    else:
        ui.load_grouping(idval, str(path))

    assert len(caplog.records) == 1
    r = caplog.record_tuples[-1]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == f"dataset {idstr}: 1:3 Channel (unchanged)"

    assert not data.grouped
    assert data.grouping is not None

    assert len(caplog.records) == 1
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.group(idval)

    assert data.grouped
    assert len(caplog.records) == 2
    r = caplog.record_tuples[-1]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == f"dataset {idstr}: 1:3 Channel (unchanged)"

    grps = ui.get_grouping(idval)
    assert grps.shape == (3, )

    # It's not clear what requirements load_grouping makes of the
    # data, so do not enforce a data type. At a minimum there
    # would be potential backend differences.
    #
    # assert grps.dtype == np.int16

    assert grps == pytest.approx([1, -1, 1])

    # Note that get_dep is returning the sum per group / channel width
    # (since we have no instrument response).
    #
    y = ui.get_dep(idval)
    assert y.shape == (2, )
    assert y == pytest.approx([2, 3])

    # Just to check no more logging has been done.
    #
    assert len(caplog.records) == 2


@requires_fits
@pytest.mark.parametrize("idval", [None, 1, "xx"])
def test_load_quality(idval, clean_astro_ui, tmp_path, caplog):
    """Simple quality check"""

    x = [1, 2, 3]
    y = [0, 4, 3]
    qual = [0, 2, 0]
    if idval is None:
        ui.load_arrays(1, x, y, ui.DataPHA)
        idstr = "1"
    else:
        ui.load_arrays(idval, x, y, ui.DataPHA)
        idstr = str(idval)

    path = tmp_path / "qual.dat"
    path.write_text("0\n2\n0")

    data = ui.get_data(idval)
    assert data.quality is None

    if idval is None:
        ui.load_quality(str(path))
    else:
        ui.load_quality(idval, str(path))

    assert not data.grouped
    assert data.quality == pytest.approx(qual)

    qual = ui.get_quality(idval)

    assert qual == pytest.approx(qual)

    # Because the data is not grouped the quality filter is reported
    # as part of the filter string.
    #
    assert ui.get_filter(idval) == "1:3"
    ui.ignore_bad(idval)
    assert ui.get_filter(idval) == "1,3"

    y = ui.get_dep(idval, filter=True)
    assert y == pytest.approx([0, 3])

    y = ui.get_dep(idval, filter=False)
    assert y == pytest.approx([0, 4, 3])


@pytest.mark.parametrize("idval", [None, 1, "xx"])
def test_group_already_grouped(idval, clean_astro_ui, caplog):
    """Does group still work if the data is already grouped?"""

    x = [1, 2, 3]
    y = [0, 4, 3]
    if idval is None:
        ui.load_arrays(1, x, y, ui.DataPHA)
        assert len(caplog.records) == 0
        ui.set_grouping([1, -1, 1])
        idstr = "1"
    else:
        ui.load_arrays(idval, x, y, ui.DataPHA)
        assert len(caplog.records) == 0
        ui.set_grouping(idval, [1, -1, 1])
        idstr = str(idval)

    assert len(caplog.records) == 1
    r = caplog.record_tuples[-1]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == f"dataset {idstr}: 1:3 Channel (unchanged)"

    data = ui.get_data(idval)
    assert not data.grouped

    assert len(caplog.records) == 1
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.group(idval)

    assert data.grouped
    assert ui.get_dep(idval) == pytest.approx([2, 3])
    assert len(caplog.records) == 2
    r = caplog.record_tuples[-1]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == f"dataset {idstr}: 1:3 Channel (unchanged)"

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.group(idval)

    assert ui.get_dep(idval) == pytest.approx([2, 3])
    assert data.grouped
    assert len(caplog.records) == 3
    r = caplog.record_tuples[-1]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == f"dataset {idstr}: 1:3 Channel (unchanged)"


@pytest.mark.parametrize("idval", [None, 1, "xx"])
def test_subtract_already_subtracted(idval, clean_astro_ui):
    """Does subtract still work if the data is already subtracted?"""

    x = [1, 2, 3]
    y = [0, 4, 3]
    ui.load_arrays("bgnd", x, y, ui.DataPHA)
    bkg = ui.get_data("bgnd")

    if idval is None:
        ui.load_arrays(1, x, y, ui.DataPHA)
        ui.set_bkg(bkg)
    else:
        ui.load_arrays(idval, x, y, ui.DataPHA)
        ui.set_bkg(idval, bkg)

    data = ui.get_data(idval)
    assert not data.subtracted

    ui.subtract(idval)
    assert data.subtracted

    ui.subtract(idval)
    assert ui.get_dep(idval) == pytest.approx([0, 0, 0])


@pytest.mark.parametrize("callfunc", [ui.group, ui.subtract])
def test_xxx_not_pha(callfunc, clean_astro_ui):
    """Just check that these commands require a PHA dataset"""

    ui.load_arrays(1, [1, 2, 3], [4, 5, 6])

    with pytest.raises(ArgumentErr,
                       match="data set 1 does not contain PHA data"):
        callfunc()


def test_get_axes_data1d(clean_astro_ui):
    ui.load_arrays(1, [2, 10, 20], [1, 2, 3])
    ax = ui.get_axes()
    assert len(ax) == 1
    assert ax[0] == pytest.approx([2, 10, 20])


def test_get_axes_data1dint(clean_astro_ui):
    ui.load_arrays(1, [2, 10, 20], [10, 12, 22], [1, 2, 3], ui.Data1DInt)
    ax = ui.get_axes()
    assert len(ax) == 2
    assert ax[0] == pytest.approx([2, 10, 20])
    assert ax[1] == pytest.approx([10, 12, 22])


def test_get_axes_datapha_no_response(clean_astro_ui):
    """Since we have no response this is a bit odd.

    I have noted that it's unclear whether the bin edges are
    1-2, 2-3, 3-4 or 0.5-1.5, 1.5-2.5, 2.5-3.5. Let's test
    the status quo here. It used to return the latter (bin edges
    at half-integer values) but in the Sherpa 4.13 development
    it was changed to use integer values (the former case).
    """

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3], ui.DataPHA)
    ax = ui.get_axes()
    assert len(ax) == 2
    assert ax[0] == pytest.approx([1, 2, 3])
    assert ax[1] == pytest.approx([2, 3, 4])


def test_get_axes_datapha_rmf(clean_astro_ui):
    """RMF only"""

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3], ui.DataPHA)

    ebins = np.asarray([0.1, 0.2, 0.4, 0.8])
    elo = ebins[:-1]
    ehi = ebins[1:]
    ui.set_rmf(ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi))

    ax = ui.get_axes()
    assert len(ax) == 2
    assert ax[0] == pytest.approx([0.1, 0.2, 0.4])
    assert ax[1] == pytest.approx([0.2, 0.4, 0.8])


def test_get_axes_datapha_arf(clean_astro_ui):
    """ARF only"""

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3], ui.DataPHA)

    ebins = np.asarray([0.1, 0.2, 0.4, 0.8])
    elo = ebins[:-1]
    ehi = ebins[1:]
    ui.set_arf(ui.create_arf(elo, ehi))

    ax = ui.get_axes()
    assert len(ax) == 2
    assert ax[0] == pytest.approx([0.1, 0.2, 0.4])
    assert ax[1] == pytest.approx([0.2, 0.4, 0.8])


def test_get_axes_datapha_rsp(clean_astro_ui):
    """Let's have a RMF and ARF for fun"""

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3], ui.DataPHA)

    ebins = np.asarray([0.1, 0.2, 0.4, 0.8])
    elo = ebins[:-1]
    ehi = ebins[1:]
    ui.set_arf(ui.create_arf(elo, ehi))
    ui.set_rmf(ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi))

    ax = ui.get_axes()
    assert len(ax) == 2
    assert ax[0] == pytest.approx([0.1, 0.2, 0.4])
    assert ax[1] == pytest.approx([0.2, 0.4, 0.8])


def test_get_axes_datapha_rsp_bkg(clean_astro_ui):
    """Let's have a RMF and ARF for fun (but only for source)

    It is not clear whether we expect the response to be automatically
    applied to the background here, so treat this as a regression test
    (because, at present, we do not have units=energy for the
    background).
    """

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3], ui.DataPHA)
    ui.set_bkg(1, ui.DataPHA("x", [1, 2, 3], [9, 0, 7]))

    ebins = np.asarray([0.1, 0.2, 0.4, 0.8])
    elo = ebins[:-1]
    ehi = ebins[1:]
    ui.set_arf(ui.create_arf(elo, ehi))
    ui.set_rmf(ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi))

    ax = ui.get_axes(bkg_id=1)
    assert len(ax) == 2
    assert ax[0] == pytest.approx([1, 2, 3])
    assert ax[1] == pytest.approx([2, 3, 4])


def test_get_axes_dataimg_logical(clean_astro_ui):
    """We don't set up a coordinate system so we just get logical back"""

    y, x = np.mgrid[-5:-1, 5:8]
    x = x.flatten()
    y = y.flatten()
    z = np.zeros(x.size)

    # Not sure what ordering the shape array is meant to have
    # and whether there is an ordering to x, y, z values.
    #
    ui.load_arrays(1, x, y, z, (4, 3), ui.DataIMG)

    ax = ui.get_axes()
    assert len(ax) == 2
    assert ax[0] == pytest.approx([1, 2, 3])
    assert ax[1] == pytest.approx([1, 2, 3, 4])


def test_get_rate_data1d(clean_astro_ui):
    """Check this falls over"""

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
    with pytest.raises(ArgumentErr,
                       match="data set 1 does not contain PHA data"):
        ui.get_rate()


def test_get_rate_datapha_rsp(clean_astro_ui):
    """Let's have a RMF and ARF for fun"""

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3], ui.DataPHA)
    ui.set_exposure(100)

    ebins = np.asarray([0.1, 0.2, 0.4, 0.8])
    elo = ebins[:-1]
    ehi = ebins[1:]
    ui.set_arf(ui.create_arf(elo, ehi, np.asarray([10, 20, 15])))
    ui.set_rmf(ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi))

    # check we do not have rate selected
    ui.set_analysis("energy", type="counts")
    assert not ui.get_data().rate

    # counts
    assert ui.get_dep() == pytest.approx([1, 2, 3])

    # counts / bin width / exposure time
    expected = np.asarray([1, 2, 3]) / np.asarray([0.1, 0.2, 0.4]) / 100
    assert ui.get_rate() == pytest.approx(expected)

    # and now switch to rate
    ui.set_analysis("energy", type="rate")
    assert ui.get_data().rate

    assert ui.get_dep() == pytest.approx([1, 2, 3])
    assert ui.get_rate() == pytest.approx(expected)


@requires_fits
@requires_data
def test_wstat_errors_data1d(clean_astro_ui, make_data_path):
    """Check we error out with a mixture of data.

    This was hit during some test that needed a clean_astro_ui
    fixture, so I just wanted to make sure we ran a similar
    test.
    """

    infile = make_data_path("3c273.pi")
    ui.load_pha(infile)

    x = np.asarray([1, 2, 3])
    y = np.asarray([2, 0, 4])
    ui.load_arrays(2, x, y)

    ui.set_source(ui.powlaw1d.m1)
    ui.set_source(2, ui.polynom1d.m2)

    ui.set_stat("wstat")

    with pytest.raises(StatErr,
                       match="No background data has been supplied. Use cstat"):
        ui.get_stat_info()


@pytest.mark.parametrize("label,vals",
                         [("grouping", [1, 1, -1, -1, 1]),
                          ("quality", [0, 0, 0, 0, 2])])
def test_pha_column_set_none(label, vals, clean_astro_ui):
    """Can clear with label=None."""

    setfunc = getattr(ui, f"set_{label}")
    getfunc = getattr(ui, f"get_{label}")

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.DataPHA)
    setfunc(vals)
    assert getfunc() == pytest.approx(vals)

    setfunc(None)
    assert getfunc() is None


@pytest.mark.parametrize("label", ["grouping", "quality"])
@pytest.mark.parametrize("vals", [True, False])
def test_pha_column_has_correct_size_scalar(label, vals, clean_astro_ui):
    """Check that the label column is the right size: scalar

    This is related to issue #1160.
    """

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.DataPHA)

    with pytest.raises(DataErr,
                       match="Array must be a sequence or None"):
        getattr(ui, f"set_{label}")(vals)


@pytest.mark.parametrize("label", ["grouping", "quality"])
@pytest.mark.parametrize("vals", [[1, 2, 3], np.ones(20)])
def test_pha_column_has_correct_size_sequence(label, vals, clean_astro_ui):
    """Check that the label column is the right size: sequence

    This is related to issue #1160.
    """

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.DataPHA)

    with pytest.raises(DataErr,
                       match=f"size mismatch between independent axis and {label}: 5 vs {len(vals)}"):
        getattr(ui, f"set_{label}")(vals)


@pytest.mark.parametrize("label", ["grouping", "quality"])
def test_pha_column_checks_pha(label, clean_astro_ui):
    """Check we error out with a non-PHA class"""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.Data1D)

    with pytest.raises(ArgumentErr,
                       match="data set 1 does not contain PHA data"):
        getattr(ui, f"get_{label}")()


@pytest.mark.parametrize("label", ["grouping", "quality"])
def test_pha_column_converts_to_int(label, clean_astro_ui):
    """Check that the label column is converted to an int"""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.DataPHA)

    getfunc = getattr(ui, f"get_{label}")
    setfunc = getattr(ui, f"set_{label}")
    assert getfunc() is None

    vals = [2.0, -2.0, 1.5, 1.9, 1.2]
    setfunc(vals)

    got = getfunc()
    assert len(got) == 5
    assert isinstance(got, np.ndarray)
    assert issubclass(got.dtype.type, np.integer)
    assert got == pytest.approx([2, -2, 1, 1, 1])


def test_pha_does_set_grouping_set_grouped(clean_astro_ui):
    """If the grouping column is set is the data automatically grouped?"""

    counts = [5, 4, 2, 3, 7]
    ui.load_arrays(1, [1, 2, 3, 4, 5], counts, ui.DataPHA)
    pha = ui.get_data()

    assert pha.grouping is None
    assert not pha.grouped
    assert ui.get_dep() == pytest.approx(counts)

    ui.set_grouping([1, -1, 1, 1, 1])

    # ui.get_dep returns something when grouped
    assert pha.grouping is not None
    assert not pha.grouped
    assert ui.get_dep() == pytest.approx(counts)


def test_group_when_background_has_no_grouping(clean_astro_ui, caplog):
    """How does group behave when the background has no grouping?

    This is something of a corner case, so we act as a regression
    test.
    """

    src = ui.DataPHA('src', [1, 2, 3], [10, 12, 2])
    bkg = ui.DataPHA('bkg', [1, 2, 3], [2, 1, 1])
    src.set_background(bkg)
    ui.set_data(src)

    assert len(caplog.records) == 0

    # Done this way, the data set is not grouped
    ui.set_grouping([1, 1, 1])
    assert len(caplog.records) == 1
    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset 1: 1:3 Channel (unchanged)"

    assert not ui.get_data().grouped
    assert ui.get_data().get_background().grouping is None

    assert len(caplog.records) == 1
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.group()

    assert len(caplog.records) == 3
    r = caplog.record_tuples[1]
    assert r[0] == "sherpa.astro.data"
    assert r[1] == logging.INFO
    assert r[2] == "data set 'bkg' does not specify grouping flags"

    r = caplog.record_tuples[2]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset 1: 1:3 Channel (unchanged)"



@pytest.mark.xfail
def test_set_grouping_1636_indirect(clean_astro_ui, caplog):
    """See issue #1636. Also #1635.

    I assume this was broken by #1477. Issue #1636 points out this is
    broken but #1635 asks what the behavior should be, so this is
    a regression test.

    See also test_set_grouping_1636_direct.

    """

    src = ui.DataPHA('src', [1, 2, 3], [10, 12, 2], grouping=[1, 1, 1])
    bkg = ui.DataPHA('bkg', [1, 2, 3], [2, 1, 1])
    src.set_background(bkg)
    ui.set_data(src)

    # Done this way, the data set is grouped, unlike
    # test_group_when_background_has_no_grouping
    assert ui.get_data().grouped
    assert ui.get_data().get_background().grouping == pytest.approx([1, 1, 1])

    # For some reason, sending both the id and the value, when val is None,
    # is a problem.
    #
    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        # FAILS with DataErr: Array must be a sequence or None
        ui.set_grouping(1, None)

    # For the moment there's no logging
    assert len(caplog.records) == 0

    assert ui.get_data().grouped
    assert ui.get_data().grouping is None


@pytest.mark.xfail
def test_set_grouping_1636_direct(clean_astro_ui, caplog):
    """See issue #1636. Also #1635.

    See also test_set_grouping_1636_indirect

    """

    src = ui.DataPHA('src', [1, 2, 3], [10, 12, 2], grouping=[1, 1, 1])
    bkg = ui.DataPHA('bkg', [1, 2, 3], [2, 1, 1])
    src.set_background(bkg)
    ui.set_data(src)

    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        # FAILS with DataErr: Array must be a sequence or None
        ui.set_grouping(id=1, val=None)

    # For the moment there's no logging
    assert len(caplog.records) == 0

    assert ui.get_data().grouped
    assert ui.get_data().grouping is None


def test_pha_what_does_get_dep_return_when_grouped(clean_astro_ui, caplog):
    """Regression test for get_dep with grouped data"""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.DataPHA)
    assert len(caplog.records) == 0

    ui.set_grouping([1, -1, 1, -1, -1])
    assert len(caplog.records) == 1
    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset 1: 1:5 Channel (unchanged)"

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.group()

    # Looks like it's returning mean of channel values in group
    assert ui.get_dep() == pytest.approx([4.5, 4])

    assert len(caplog.records) == 2
    r = caplog.record_tuples[-1]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset 1: 1:5 Channel (unchanged)"


@requires_fits
@requires_data
@requires_region
@requires_wcs
def test_image_filter_coord_change_same(make_data_path, clean_astro_ui, caplog):
    """What happens to the mask after a coordinate change? NO CHANGE

    This is really just a way to test the DataIMG class without having
    to set up all the values. There was some interesting behavior with
    calc_data_sum2d which was down to whether changing the coord would
    change the mask, so this is an explicit test of that behavior.

    """

    ui.load_image("foo", make_data_path("image2.fits"))
    assert ui.get_filter("foo") == ""

    ui.set_coord("foo", "physical")

    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.notice2d_id("foo", "rect(4000, 4200, 4100 , 4300 ) ")

    assert ui.get_filter("foo") == "Rectangle(4000,4200,4100,4300)"
    assert len(caplog.records) == 1

    # Is there no way to get the mask data via the UI interface,
    # without calling save_filter?
    #
    expected = 2500
    d = ui.get_data("foo")
    assert d.mask.sum() == expected

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.set_coord("foo", "physical")

    assert ui.get_filter("foo") == "Rectangle(4000,4200,4100,4300)"
    assert d.mask.sum() == expected

    assert len(caplog.records) == 1

    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset foo: Field() -> Rectangle(4000,4200,4100,4300)"


@requires_fits
@requires_data
@requires_region
@requires_wcs
def test_image_filter_coord_change(make_data_path, clean_astro_ui, caplog):
    """What happens to the mask after a coordinate change?

    This is the more-general case of test_image_filter_coord_change_same
    but it's not actually clear what should be done here when the
    coordinate setting is changed - should the mask be removed as
    we don't store the mask expression with a way to update it to
    match the new coordianate setting.

    """

    ui.load_image("foo", make_data_path("image2.fits"))

    ui.set_coord("foo", "physical")

    # We test out the "send an array if ids" here, just because
    # we can. Because it is slightly different to
    # test_image_filter_coord_change_same_negate we re-run some
    # of the tests.
    #
    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.notice2d_id(["foo"], "rect(4000, 4200, 4100 , 4300 ) ")

    assert ui.get_filter("foo") == "Rectangle(4000,4200,4100,4300)"
    assert len(caplog.records) == 1

    expected = 2500
    d = ui.get_data("foo")
    assert d.mask.sum() == expected

    assert len(caplog.records) == 1
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.set_coord("foo", "logical")

    assert len(caplog.records) == 2

    # The region filter has been removed
    assert ui.get_filter("foo") == ""
    assert d.mask is True

    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset foo: Field() -> Rectangle(4000,4200,4100,4300)"

    r = caplog.record_tuples[1]
    assert r[0] == "sherpa.astro.data"
    assert r[1] == logging.WARN
    assert r[2].startswith("Region filter has been removed from '")


@requires_fits
@requires_data
@requires_region
@requires_wcs
def test_image_filter_coord_change_same_negate(make_data_path, clean_astro_ui, caplog):
    """A ignore2d case of test_image_filter_coord_change_same

    """

    ui.load_image("foo", make_data_path("image2.fits"))

    ui.set_coord("foo", "physical")
    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.ignore2d_id("foo", "rect(4000, 4200, 4100 , 4300 ) ")

    assert ui.get_filter("foo") == "Field()&!Rectangle(4000,4200,4100,4300)"
    assert len(caplog.records) == 1

    # The image has 56376 pixels and we know from
    # test_image_filter_coord_change that that filter selects 2500
    # pixels, so we should have the following (technically thanks to
    # edge effects this may not hold, but it does here).
    #
    expected = 56376 - 2500
    d = ui.get_data("foo")
    assert d.mask.sum() == expected

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.set_coord("foo", "physical")

    assert ui.get_filter("foo") == "Field()&!Rectangle(4000,4200,4100,4300)"
    assert d.mask.sum() == expected

    assert len(caplog.records) == 1

    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset foo: Field() -> Field()&!Rectangle(4000,4200,4100,4300)"


@requires_fits
@requires_data
@requires_region
@requires_wcs
def test_image_filter_coord_change_negate(make_data_path, clean_astro_ui, caplog):
    """negate the filter used in test_image_filter_coord_change

    This was created as we didn't have sensible tests of ignore2d_id
    and it seemed easiest to do it this way.
    """

    ui.load_image("foo", make_data_path("image2.fits"))

    ui.set_coord("foo", "physical")

    # We test out the "send an array if ids" here, just because
    # we can. Because it is slightly different to
    # test_image_filter_coord_change_same_negate we re-run some
    # of the tests.
    #
    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.ignore2d_id(["foo"], "rect(4000, 4200, 4100 , 4300 ) ")

    assert len(caplog.records) == 1
    assert ui.get_filter("foo") == "Field()&!Rectangle(4000,4200,4100,4300)"

    expected = 56376 - 2500
    d = ui.get_data("foo")
    assert d.mask.sum() == expected

    assert len(caplog.records) == 1
    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.set_coord("foo", "logical")

    assert len(caplog.records) == 2

    assert ui.get_filter("foo") == ""
    assert d.mask is True

    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset foo: Field() -> Field()&!Rectangle(4000,4200,4100,4300)"

    r = caplog.record_tuples[1]
    assert r[0] == "sherpa.astro.data"
    assert r[1] == logging.WARN
    assert r[2].startswith("Region filter has been removed from '")


def test_set_source_checks_dimensionality_1d(clean_astro_ui):
    """Check that you can not set a model of the wrong dimensionality"""

    chans = np.arange(1, 5)
    ui.load_arrays(1, chans, chans * 0, ui.DataPHA)

    with pytest.raises(DataErr,
                       match="Data and model dimensionality do not match: 1D and 2D"):
        ui.set_source(ui.const2d.mdl)


def test_set_source_checks_dimensionality_2d(clean_astro_ui):
    """Check that you can not set a model of the wrong dimensionality"""

    x0 = np.asarray([1, 2, 1, 2])
    x1 = np.asarray([1, 1, 1, 1])
    y = np.asarray([1] * 4)
    ui.load_arrays(1, x0, x1, y, (2, 2), ui.DataIMG)

    with pytest.raises(DataErr,
                       match="Data and model dimensionality do not match: 2D and 1D"):
        ui.set_source(ui.const1d.mdl)


def test_model_dimensionality_check_is_not_triggered_calc_stat(clean_astro_ui):
    """We can still end up where a model and dataset dimensionality do not match.

    So check that other parts of the system catch this error: calc_stat
    """

    ui.load_arrays(1, np.asarray([1, 2, 1, 1]), np.asarray([1, 1, 2, 2]),
                   np.asarray([1, 1, 1, 1]), (2, 2), ui.DataIMG)
    ui.set_source(ui.const2d.mdl)

    chans = np.arange(1, 5)
    ui.load_arrays(1, chans, chans * 0)
    assert ui.get_data().ndim == 1

    # Changing the data set currently does not reset the source model
    # so check this.
    src = ui.get_source()
    assert src.ndim == 2

    # However, this creates a Fit object and so does error out
    #
    with pytest.raises(DataErr,
                       match="Data and model dimensionality do not match: 1D and 2D"):
        ui.calc_stat()


def test_model_dimensionality_check_is_not_triggered_plot_model(clean_astro_ui):
    """We can still end up where a model and dataset dimensionality do not match.

    So check that other parts of the system catch this error: plot_model
    """

    ui.load_arrays(1, np.asarray([1, 2, 1, 1]), np.asarray([1, 1, 2, 2]),
                   np.asarray([1, 1, 1, 1]), (2, 2), ui.DataIMG)
    ui.set_source(ui.const2d.mdl)

    chans = np.arange(1, 5)
    ui.load_arrays(1, chans, chans * 0)
    assert ui.get_data().ndim == 1

    src = ui.get_source()
    assert src.ndim == 2

    with pytest.raises(DataErr,
                       match="Data and model dimensionality do not match: 1D and 2D"):
        ui.plot_model()


def simple_data1dint(idval):
    """Create a simple dataset for the #1444 checks.

    Usw Data1DInt rather than Data1D as the fact it's integrated
    gives a chance to check the model evaluation.
    """

    xlo = np.asarray([1, 3, 7])
    xhi = np.asarray([2, 7, 9])
    y = np.asarray([5, 17, 7])

    ui.load_arrays(idval, xlo, xhi, y, ui.Data1DInt)
    ui.set_source(idval, ui.const1d.mdl)

    mdl.c0 = 100


def simple_dataimg(idval):
    """Create a simple dataset for the #1444 checks.

    We want to use DataIMG and not Data2D to ensure we have
    notice2d.
    """

    x1, x0 = np.mgrid[1:4, 1:5]
    y = x1 + 10 * x0
    ui.load_arrays(idval, x0.flatten(), x1.flatten(), y.flatten(),
                   x0.shape, ui.DataIMG)

    ui.set_source(idval, ui.polynom2d.mdl)
    mdl.c = 100
    mdl.cx1 = 10
    mdl.cy1 = 1


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_1d_data_1d(idval, clean_astro_ui):
    """Check issue #1444 for a 1D data set: calc_data_sum

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_data1dint(idval)
    assert ui.calc_data_sum(id=idval) == pytest.approx(29)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_1d_data_2d(idval, clean_astro_ui):
    """Check issue #1444 for a 1D data set: calc_data_sum2d

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_data1dint(idval)
    with pytest.raises(DataErr,
                       match="data set '' does not contain 2-D data"):
        ui.calc_data_sum2d(id=idval)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_1d_model_1d(idval, clean_astro_ui):
    """Check issue #1444 for a 1D data set: calc_model_sum

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_data1dint(idval)
    assert ui.calc_model_sum(id=idval) == pytest.approx(700)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_1d_model_2d(idval, clean_astro_ui):
    """Check issue #1444 for a 1D data set: calc_model_sum2d

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_data1dint(idval)
    with pytest.raises(DataErr,
                       match="data set '' does not contain 2-D data"):
        ui.calc_model_sum2d(id=idval)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_1d_source_1d(idval, clean_astro_ui):
    """Check issue #1444 for a 1D data set: calc_source_sum

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_data1dint(idval)
    # TODO: what is calc_source_sum calculating here?
    assert ui.calc_source_sum(id=idval) == pytest.approx(300)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_1d_source_2d(idval, clean_astro_ui):
    """Check issue #1444 for a 1D data set: calc_source_sum2d

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_data1dint(idval)
    with pytest.raises(DataErr,
                       match="data set '' does not contain 2-D data"):
        ui.calc_source_sum2d(id=idval)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_2d_data_1d(idval, clean_astro_ui):
    """Check issue #1444 for a 2D data set: calc_data_sum

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_dataimg(idval)
    with pytest.raises(DataErr,
                       match="data set '' does not contain 1-D data"):
        ui.calc_data_sum(id=idval)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_2d_data_2d(idval, clean_astro_ui):
    """Check issue #1444 for a 2D data set: calc_data_sum2d

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_dataimg(idval)
    assert ui.calc_data_sum2d(id=idval) == pytest.approx(324)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_2d_model_1d(idval, clean_astro_ui):
    """Check issue #1444 for a 2D data set: calc_model_sum

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_dataimg(idval)
    with pytest.raises(DataErr,
                       match="data set '' does not contain 1-D data"):
        ui.calc_model_sum(id=idval)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_2d_model_2d(idval, clean_astro_ui):
    """Check issue #1444 for a 2D data set: calc_model_sum2d

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_dataimg(idval)
    assert ui.calc_model_sum2d(id=idval) == pytest.approx(1524)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_2d_source_1d(idval, clean_astro_ui):
    """Check issue #1444 for a 2D data set: calc_source_sum

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_dataimg(idval)
    with pytest.raises(DataErr,
                       match="data set '' does not contain 1-D data"):
        ui.calc_source_sum(id=idval)


@pytest.mark.parametrize("idval", [None, 1, 9, "foo"])
def test_1444_2d_source_2d(idval, clean_astro_ui):
    """Check issue #1444 for a 2D data set: calc_source_sum2d

    This is a regression test, in that it just tests the existing behavior.
    """

    simple_dataimg(idval)
    assert ui.calc_source_sum2d(id=idval) == pytest.approx(1524)


@pytest.fixture
def simple_pha(clean_astro_ui):
    """A simple grouped PHA dataset."""

    chans = np.arange(1, 10)
    counts = chans * 2
    grouping = np.asarray([1, -1, -1, 1, 1, -1, -1, 1, 1])
    data = ui.DataPHA('x', chans, counts)
    data.grouping = grouping

    assert not data.grouped
    data.group()

    ui.set_data(data)
    ui.set_stat("chi2datavar")

    # return the expected counts
    return 2 * np.asarray([1 + 2 + 3, 4, 5 + 6 + 7, 8, 9])


@pytest.fixture
def background_pha(clean_astro_ui):
    """A simple PHA dataset with background and response.

    The source and background ARFs are different but the
    RMF is the same.
    """

    chans = np.arange(1, 10)
    counts = chans * 2
    grouping = np.asarray([1, -1, -1, 1, 1, -1, -1, 1, 1])
    data = ui.DataPHA('x', chans, counts, grouping=grouping)
    data.grouping = grouping
    assert data.grouped

    bkg = ui.DataPHA('y', chans, np.ones_like(chans))

    ui.set_data(data)
    ui.set_bkg(bkg)
    ui.set_stat("chi2datavar")

    egrid = np.arange(1, 11) * 0.1
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

    ui.set_rmf(rmf)
    ui.set_rmf(rmf, bkg_id=1)

    sarf = create_arf(elo, ehi, np.ones_like(elo) * 0.8)
    barf = create_arf(elo, ehi, np.ones_like(elo) * 0.4)

    ui.set_arf(sarf)
    ui.set_arf(barf, bkg_id=1)


def test_pha_set_dep_none(simple_pha):
    """What happens if set_dep is called with None?"""

    ui.set_dep(None)

    # We get [nan, nan, nan, nan, nan]
    assert np.isnan(ui.get_dep()) == pytest.approx([1, 1, 1, 1, 1])


def test_pha_set_dep_scalar(simple_pha):
    """What happens if set_dep is called with a scalar?"""

    ui.set_dep(3)
    assert ui.get_dep() == pytest.approx(3 * np.ones(5))


def test_pha_set_dep_array(simple_pha):
    """What happens if set_dep is called with an array?"""

    # The grouping is 1-3, 4, 5-7, 8, 9
    ui.set_dep([2, 4, 5, 19, 11, 12, 13, 0, 9])

    # I do not understand why get_dep returns the "average" value
    # rather than the sum, but let's test the current behavior.
    #
    assert ui.get_dep() == pytest.approx([11 / 3, 19, 36 / 3, 0, 9])


def test_pha_set_dep_array_wrong(simple_pha):
    """What happens if set_dep is called with an array with the wrong length?"""

    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and y: 9 vs 3"):
        ui.set_dep([2, 4, 5])


def test_pha_get_staterror(simple_pha):
    """Check get_staterror with a PHA dataset"""

    error = np.sqrt(simple_pha)
    assert ui.get_staterror() == pytest.approx(error)
    assert ui.get_error() == pytest.approx(error)


def test_pha_get_syserror(simple_pha):
    """Check get_syserror with a PHA dataset"""

    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror()


def test_pha_set_staterror_scalar_no_fractional(simple_pha):
    """What happens when we set the staterror to a scalar fractional=False?"""

    # The errors are calculated from
    #    (3,3,3), 3, (3, 3, 3), 3, 3
    # thanks to the grouping. Those bins with multiple values have
    # the error values combined in quadrature.
    #
    error = 3 * np.sqrt([3, 1, 3, 1, 1])

    ui.set_staterror(3)
    assert ui.get_staterror() == pytest.approx(error)
    assert ui.get_error() == pytest.approx(error)


def test_pha_set_syserror_scalar_no_fractional(simple_pha):
    """What happens when we set the syserror to a scalar fractional=False?"""

    staterror = np.sqrt(simple_pha)
    syserror = 2 * np.sqrt([3, 1, 3, 1, 1])
    combo = np.sqrt(simple_pha + syserror * syserror)

    ui.set_syserror(2)
    assert ui.get_staterror() == pytest.approx(staterror)
    assert ui.get_syserror() == pytest.approx(syserror)
    assert ui.get_error() == pytest.approx(combo)


def test_pha_set_staterror_scalar_fractional(simple_pha):
    """What happens when we set the staterror to a scalar fractional=True?"""

    vals = 0.4 * 2 * np.arange(1, 10)
    vals2 = vals * vals
    error = np.sqrt([vals2[0] + vals2[1] + vals2[2], vals2[3],
                     vals2[4] + vals2[5] + vals2[6], vals2[7],
                     vals2[8]])

    ui.set_staterror(0.4, fractional=True)
    assert ui.get_staterror() == pytest.approx(error)
    assert ui.get_error() == pytest.approx(error)


def test_pha_set_syserror_scalar_fractional(simple_pha):
    """What happens when we set the syserror to a scalar fractional=True?"""

    staterror = np.sqrt(simple_pha)

    vals = 0.5 * 2 * np.arange(1, 10)
    vals2 = vals * vals
    syserror = np.sqrt([vals2[0] + vals2[1] + vals2[2], vals2[3],
                        vals2[4] + vals2[5] + vals2[6], vals2[7],
                        vals2[8]])

    combo = np.sqrt(simple_pha + syserror * syserror)

    ui.set_syserror(0.5, fractional=True)
    assert ui.get_staterror() == pytest.approx(staterror)
    assert ui.get_syserror() == pytest.approx(syserror)
    assert ui.get_error() == pytest.approx(combo)


def test_pha_set_staterror_array(simple_pha):
    """What happens when we set the staterror to an array?"""

    # The errors are given per channel, but the response is
    # grouped.
    #
    vals = 0.1 * np.arange(1, 10)
    vals2 = vals * vals
    error = np.sqrt([vals2[0] + vals2[1] + vals2[2], vals2[3],
                     vals2[4] + vals2[5] + vals2[6], vals2[7],
                     vals2[8]])

    ui.set_staterror(vals)
    assert ui.get_staterror() == pytest.approx(error)
    assert ui.get_error() == pytest.approx(error)


def test_pha_set_syserror_array(simple_pha):
    """What happens when we set the syserror to an array?"""

    vals = 0.1 * np.arange(1, 10)
    vals2 = vals * vals
    syserror = np.sqrt([vals2[0] + vals2[1] + vals2[2], vals2[3],
                        vals2[4] + vals2[5] + vals2[6], vals2[7],
                        vals2[8]])

    combo = np.sqrt(simple_pha + syserror * syserror)

    ui.set_syserror(vals)
    assert ui.get_syserror() == pytest.approx(syserror)
    assert ui.get_error() == pytest.approx(combo)


@pytest.mark.parametrize("field", ["staterror", "syserror"])
def test_pha_set_error_array_wrong(field, simple_pha):
    """What happens when we set the stat/syserror to an array of the wrong length?"""

    setfunc = getattr(ui, f"set_{field}")

    # NOTE: this error is confusing, because the length should be 5 not 9
    #       ie the filter does not match the grouped data
    #
    with pytest.raises(DataErr,
                       match=f"size mismatch between independent axis and {field}: 9 vs 4"):
        setfunc(np.asarray([1, 2, 3, 4]))


def test_pha_set_filter_unmasked(simple_pha):
    """What happens when we call set_filter to an unfiltered dataset?"""

    data = ui.get_data()
    assert data.mask

    expected = [True, True, False, True, False]
    ui.set_filter(expected)

    assert data.mask == pytest.approx(expected)


def test_pha_set_filter_unmasked_wrong(simple_pha):
    """What happens when we call set_filter to an unfiltered dataset with the wrong size?"""

    with pytest.raises(DataErr,
                       match="size mismatch between grouped data and mask: 5 vs 2"):
        ui.set_filter(np.asarray([True, False]))


def test_pha_set_filter_masked(simple_pha):
    """What happens when we call set_filter to a filtered dataset?"""

    data = ui.get_data()

    ui.ignore(4, 8)
    assert data.mask == pytest.approx([True, False, False, False, True])

    ui.set_filter(np.asarray([True, False, True, False, False]))
    assert data.mask == pytest.approx([True, False, True, False, True])


def test_pha_set_filter_masked_wrong(simple_pha):
    """What happens when we call set_filter to a filtered dataset with the wrong size?"""

    ui.ignore(4, 8)

    with pytest.raises(DataErr,
                       match="size mismatch between data and filter: 5 vs 2"):
        ui.set_filter(np.asarray([True, False]))


def test_pha_get_specresp_src(background_pha):

    arf = np.ones(9) * 0.8
    assert ui.get_specresp() == pytest.approx(arf)


def test_pha_get_specresp_bkg(background_pha):

    arf = np.ones(9) * 0.4
    assert ui.get_specresp(bkg_id=1) == pytest.approx(arf)


def test_get_arf_not_pha(clean_astro_ui):
    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    ui.load_arrays(2, [1, 2], [1, 2])
    with pytest.raises(ArgumentErr,
                       match="data set 2 does not contain PHA data"):
        ui.get_arf_plot(2)


def test_get_arf_no_pha(clean_astro_ui):
    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    with pytest.raises(DataErr,
                       match="data set '1' does not have an associated ARF"):
        ui.get_arf_plot()


def test_get_order_not_pha(clean_astro_ui):
    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    ui.load_arrays(2, [1, 2], [1, 2])
    with pytest.raises(ArgumentErr,
                       match="data set 2 does not contain PHA data"):
        ui.get_order_plot(2)


@pytest.mark.parametrize("label", ["area", "back"])
def test_set_xxxscal_no_data(label, clean_astro_ui):
    """What happens when the dataset does not exist?"""

    func = getattr(ui, f"get_{label}scal")

    with pytest.raises(IdentifierErr,
                       match="data set 1 has not been set"):
        func()


@pytest.mark.parametrize("label", ["area", "back"])
def test_set_xxxscal_no_data_id(label, clean_astro_ui):
    """What happens when the dataset does not exist?"""

    func = getattr(ui, f"get_{label}scal")

    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    with pytest.raises(IdentifierErr,
                       match="data set 2 has not been set"):
        func(2)


@pytest.mark.parametrize("label", ["area", "back"])
def test_set_xxxscal_can_be_none(label, clean_astro_ui):
    """What happens when xxxscal is not set?

    This is a regression test. We may decide to handle this
    differently in the future.

    """
    func = getattr(ui, f"get_{label}scal")

    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    assert func() is None


@pytest.mark.parametrize("label", ["area", "back"])
def test_set_xxxscal_invalid(label, clean_astro_ui):
    """What happens when xxxscal is obviously wrong?

    At the moment there's no check, so it is allowed.
    """
    getfunc = getattr(ui, f"get_{label}scal")
    setfunc = getattr(ui, f"set_{label}scal")

    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    setfunc(0)
    assert getfunc() == pytest.approx(0)


@pytest.mark.parametrize("label", ["area", "back"])
def test_set_xxxscal_integer(label, clean_astro_ui):
    """What happens when xxxscal is set to an integer?"""
    getfunc = getattr(ui, f"get_{label}scal")
    setfunc = getattr(ui, f"set_{label}scal")

    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    setfunc(12)
    got = getfunc()
    assert got == pytest.approx(12)
    assert isinstance(got, np.float64)


@pytest.mark.parametrize("label", ["area", "back"])
def test_set_xxxscal_float(label, clean_astro_ui):
    """What happens when xxxscal is set to a float?"""
    getfunc = getattr(ui, f"get_{label}scal")
    setfunc = getattr(ui, f"set_{label}scal")

    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    setfunc(1.2e-3)
    got = getfunc()
    assert got == pytest.approx(1.2e-3)
    assert isinstance(got, np.float64)


@pytest.mark.parametrize("label", ["area", "back"])
def test_set_xxxscal_bkg(label, clean_astro_ui):
    getfunc = getattr(ui, f"get_{label}scal")
    setfunc = getattr(ui, f"set_{label}scal")

    ui.load_arrays(1, [1, 2], [1, 2], ui.DataPHA)
    bkg = ui.DataPHA("bkg", [1, 2], [1, 2])
    ui.set_bkg(bkg)

    setfunc(1.2e-2)
    setfunc(2.4e-2, bkg_id=1)
    got1 = getfunc()
    got2 = getfunc(bkg_id=1)

    assert got2 / got1 == pytest.approx(2)
    assert got2 == pytest.approx(2.4e-2)


@pytest.mark.parametrize("label", ["area", "back"])
def test_set_xxxscal_id(label, clean_astro_ui):
    getfunc = getattr(ui, f"get_{label}scal")
    setfunc = getattr(ui, f"set_{label}scal")

    ui.load_arrays("foo", [1, 2], [1, 2], ui.DataPHA)
    bkg = ui.DataPHA("bkg", [1, 2], [1, 2])
    ui.set_bkg("foo", bkg)

    setfunc("foo", 1.2e-2)
    setfunc("foo", 2.4e-2, bkg_id=1)
    got1 = getfunc("foo")
    got2 = getfunc("foo", bkg_id=1)

    assert got2 / got1 == pytest.approx(2)
    assert got2 == pytest.approx(2.4e-2)


@pytest.mark.parametrize("idval", [None, 1, "up"])
@pytest.mark.parametrize("bkg_id", [1, "up"])
def test_dataspace1d_is_a_background(idval, bkg_id, clean_astro_ui):
    """What happens if we try to make a PHA background?"""

    ui.dataspace1d(1, 10, id=idval, dstype=ui.DataPHA)
    ui.dataspace1d(1, 10, id=idval, bkg_id=bkg_id, dstype=ui.DataPHA)
    expected = [1] if idval is None else [idval]
    assert ui.list_data_ids() == pytest.approx(expected)
    assert ui.list_bkg_ids(idval) == pytest.approx([bkg_id])


def test_dataspace1d_not_a_background(clean_astro_ui):
    """What happens if we try to make a non-PHA background?"""

    ui.dataspace1d(1, 10, dstype=ui.DataPHA)
    with pytest.raises(ArgumentTypeErr,
                       match="^'dstype' must be set to DataPHA$"):
        ui.dataspace1d(1, 10, bkg_id=1)


def test_set_bkg_not_a_background(clean_astro_ui):
    """What happens if try to make add a non-PHA background?"""

    ui.dataspace1d(1, 10, dstype=ui.DataPHA)
    with pytest.raises(ArgumentTypeErr,
                       match="^'bkg' must be a PHA data set$"):
        ui.set_bkg(ui.Data1D("x", [1, 2], [1, 2]))


@requires_fits
@requires_data
@requires_wcs  # only needed for physical / world
@pytest.mark.parametrize("incoord,outcoord,x0,x1",
                         [("logical", None, 1, 1),
                          ("image", "logical", 1, 1),
                          ("physical", None, 2987.3400878906, 4363.16015625),
                          ("world", None, 150.03707837630427, 2.644383154504334),
                          ("wcs", "world", 150.03707837630427, 2.644383154504334)])
def test_1762_ui(incoord, outcoord, x0, x1, clean_astro_ui, make_data_path):
    """A version of sherpa/astro/ui/tests/test_io_img.py:test_1762"""

    infile = make_data_path("acisf08478_000N001_r0043_regevt3_srcimg.fits")
    ui.load_image("img", infile, coord=incoord)

    # Not really needed but added in anyway
    assert ui.list_data_ids() == ["img"]

    got = ui.get_coord("img")
    if outcoord is None:
        assert got == incoord
    else:
        assert got == outcoord

    # Just check the first element
    i0, i1 = ui.get_indep("img")
    assert i0[0] == pytest.approx(x0)
    assert i1[0] == pytest.approx(x1)


@requires_xspec
@requires_fits
@requires_data
@pytest.mark.parametrize("idval", [None, 2])
def test_guess_with_response_and_multiple_models(idval, clean_astro_ui, caplog, make_data_path):
    """Check we can call guess on XSPEC models.

    This should really be with made-up responses and data to
    check they work, but for now read in the data.
    """

    infile = make_data_path("3c273.pi")

    with SherpaVerbosity("ERROR"):
        ui.load_pha(id=idval, arg=infile)
        ui.notice(lo=0.3, hi=6)
        ui.subtract(id=idval)

    cpt1 = ui.create_model_component("xsphabs", "gal")
    cpt2 = ui.create_model_component("xsapec", "src1")
    cpt3 = ui.create_model_component("xsgaussian", "src2")
    ui.set_model(id=idval, model=cpt1 * (cpt2 + cpt3))

    # A limited check on the parameter values - e.g. this does not
    # check the other parameters.
    #
    assert cpt1.nH.val == pytest.approx(1)
    assert cpt1.nH.min == pytest.approx(0)
    assert cpt1.nH.max == pytest.approx(1e6)

    assert cpt2.norm.val == pytest.approx(1)
    assert cpt2.norm.min == pytest.approx(0)
    assert cpt2.norm.max == pytest.approx(1e24)

    assert cpt3.norm.val == pytest.approx(1)
    assert cpt3.norm.min == pytest.approx(0)
    assert cpt3.norm.max == pytest.approx(1e24)

    ui.guess(id=idval)

    assert cpt1.nH.val == pytest.approx(1)
    assert cpt1.nH.min == pytest.approx(0)
    assert cpt1.nH.max == pytest.approx(1e6)

    expected = 3.539017671368409
    assert cpt2.norm.val == pytest.approx(expected / 1000)
    assert cpt2.norm.min == pytest.approx(expected / 1000 / 1000)
    assert cpt2.norm.max == pytest.approx(expected)

    expected = 9.910362995161027
    assert cpt3.norm.val == pytest.approx(expected / 1000)
    assert cpt3.norm.min == pytest.approx(expected / 1000 / 1000)
    assert cpt3.norm.max == pytest.approx(expected)

    # The only message is from xsphabs
    assert len(caplog.records) == 1
    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.models.model"
    assert r[1] == logging.WARN
    assert r[2] == "No guess found for xsphabs.gal"


def test_pha_ignore_bad_no_quality(clean_astro_ui):
    """Check error case."""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 2, 3, 4, 55], ui.DataPHA)
    with pytest.raises(DataErr,
                       match="^data set '' does not specify quality flags$"):
        ui.ignore_bad()


def test_pha_ignore_bad(caplog, clean_astro_ui):
    """Regression test."""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 2, 3, 4, 55], ui.DataPHA)
    ui.set_quality([0, 0, 2, 0, 5])

    assert len(caplog.records) == 0
    ui.ignore_bad()
    assert len(caplog.records) == 1

    assert ui.get_dep(filter=False) == pytest.approx([5, 2, 3, 4, 55])
    assert ui.get_dep(filter=True) == pytest.approx([5, 2, 4])

    r = caplog.records[0]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:5 -> 1:2,4 Channel"


def test_pha_filter_ignore_bad(caplog, clean_astro_ui):
    """Regression test."""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 2, 3, 4, 55], ui.DataPHA)
    ui.set_quality([0, 0, 2, 0, 5])

    ui.notice(hi=4)
    assert len(caplog.records) == 1

    assert ui.get_dep(filter=False) == pytest.approx([5, 2, 3, 4, 55])
    assert ui.get_dep(filter=True) == pytest.approx([5, 2, 3, 4])

    ui.ignore_bad()
    assert len(caplog.records) == 2

    assert ui.get_dep(filter=False) == pytest.approx([5, 2, 3, 4, 55])
    assert ui.get_dep(filter=True) == pytest.approx([5, 2, 4])

    r = caplog.records[0]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:5 -> 1:4 Channel"

    r = caplog.records[1]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:4 -> 1:2,4 Channel"


def test_pha_group_filter_ignore_bad(caplog, clean_astro_ui):
    """Regression test."""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 2, 3, 4, 55], ui.DataPHA)
    ui.set_quality([0, 0, 2, 0, 5])
    ui.set_grouping([1, -1, 1, 1, 1])
    ui.group()
    assert len(caplog.records) == 2

    ui.notice(hi=4)
    assert len(caplog.records) == 3

    assert ui.get_dep(filter=False) == pytest.approx([3.5, 3, 4, 55])
    assert ui.get_dep(filter=True) == pytest.approx([3.5, 3, 4])

    ui.ignore_bad()
    assert len(caplog.records) == 5

    assert ui.get_dep(filter=False) == pytest.approx([3.5, 4])
    assert ui.get_dep(filter=True) == pytest.approx([3.5, 4])

    r = caplog.records[0]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:5 Channel (unchanged)"

    r = caplog.records[1]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:5 Channel (unchanged)"

    r = caplog.records[2]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:5 -> 1:4 Channel"

    r = caplog.records[3]
    assert r.name == "sherpa.astro.data"
    assert r.levelname == "WARNING"
    assert r.getMessage() == "filtering grouped data with quality flags, previous filters deleted"

    r = caplog.records[4]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:4 -> 1:5 Channel"


def test_pha_group_ignore_bad_group(caplog, clean_astro_ui):
    """Regression test.

    Check what happens after ignore_bad(); group_counts()
    """

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 2, 3, 4, 55], ui.DataPHA)
    ui.set_quality([0, 0, 2, 0, 5])
    ui.set_grouping([1, -1, 1, 1, 1])
    ui.group()
    ui.ignore_bad()
    assert len(caplog.records) == 3

    assert ui.get_dep(filter=False) == pytest.approx([3.5, 4])
    assert ui.get_dep(filter=True) == pytest.approx([3.5, 4])

    ui.group_counts(3)
    assert len(caplog.records) == 4

    assert ui.get_dep(filter=False) == pytest.approx([5, 2, 4])
    assert ui.get_dep(filter=True) == pytest.approx([5, 2, 4])

    r = caplog.records[3]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:5 Channel (unchanged)"


def test_pha_group_filter_ignore_bad_group(caplog, clean_astro_ui):
    """Regression test.

    Check what happens after ignore_bad(); notice(..); notice()
    """

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 2, 3, 4, 55], ui.DataPHA)
    ui.set_quality([0, 0, 2, 0, 5])
    ui.set_grouping([1, -1, 1, 1, 1])
    ui.group()
    ui.ignore_bad()
    assert len(caplog.records) == 3

    assert ui.get_dep(filter=False) == pytest.approx([3.5, 4])
    assert ui.get_dep(filter=True) == pytest.approx([3.5, 4])

    ui.notice(lo=2, hi=4)
    assert len(caplog.records) == 4

    assert ui.get_dep(filter=False) == pytest.approx([3.5, 4])
    assert ui.get_dep(filter=True) == pytest.approx([3.5, 4])

    ui.notice()
    assert len(caplog.records) == 5

    assert ui.get_dep(filter=False) == pytest.approx([3.5, 3, 4, 55])
    assert ui.get_dep(filter=True) == pytest.approx([3.5, 3, 4, 55])

    r = caplog.records[4]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    assert r.getMessage() == "dataset 1: 1:4 -> 1:5 Channel"
