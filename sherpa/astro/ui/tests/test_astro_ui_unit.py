#
#  Copyright (C) 2017, 2018, 2020, 2021, 2022
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

import numpy as np

import pytest

from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa.astro import io
from sherpa.astro import ui
from sherpa.utils import poisson_noise
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, DataErr, \
    IdentifierErr, IOErr, ModelErr, StatErr
from sherpa.utils.testing import requires_data, requires_fits


def backend_is(name):
    return io.backend.__name__ == f"sherpa.astro.io.{name}_backend"


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
    """Create a 1D dataset for fitting a gaussian-like profile.

    This sets the random seed to a fixed value, so use the
    reset_seed fixture.
    """

    x = np.linspace(2300, 2400, 51)

    cpt = ui.voigt1d("cpt")
    cpt.pos = 2345
    cpt.fwhm_g = 20
    cpt.fwhm_l = 15
    cpt.ampl = 480

    np.random.seed(72383)
    y = poisson_noise(cpt(x))

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
def test_fit_profile(model, stat, pars, reset_seed, clean_astro_ui):
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

    with pytest.raises(ArgumentTypeErr) as exc:
        func(None)

    assert str(exc.value) == "'ids' must be an identifier or list of identifiers"


@pytest.mark.parametrize("f", [[False, True], np.asarray([False, True])])
@pytest.mark.parametrize("bid", [None, 1])
def test_set_filter_mismatch(f, bid):
    """Does set_filter error when there's a mis-match?
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)
    with pytest.raises(DataErr) as exc:
        ui.set_filter(f, bkg_id=bid)

    assert str(exc.value) == "size mismatch between 3 and 2"


@pytest.mark.parametrize("f", [[False, True], np.asarray([False, True])])
@pytest.mark.parametrize("bid", [None, 1])
def test_set_filter_mismatch_with_filter(f, bid):
    """Does set_filter error when there's a mis-match after a filter?

    test_set_filter_mismatch checks when .mask is a scalar,
    so now check when it's a NumPy array.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    ui.ignore(3, None)  # set the .mask attribute to an array

    with pytest.raises(DataErr) as exc:
        ui.set_filter(f, bkg_id=bid)

    assert str(exc.value) == "size mismatch between 3 and 2"


@pytest.mark.parametrize("bid", [None, 1])
def test_get_syserror_missing(bid):
    """Does get_syserror error out?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)
    with pytest.raises(DataErr) as exc:
        ui.get_syserror(bkg_id=bid)

    assert str(exc.value) == "data set '1' does not specify systematic errors"


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
def test_save_filter_ignored(bid, tmp_path):
    """Does save_filter error out if everything is masked?

    We should be able to write out the filter in this case,
    as it's easy (all False).
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    ui.ignore(None, None)

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(DataErr) as exc:
        ui.save_filter(str(outfile), bkg_id=bid)

    assert str(exc.value) == "mask excludes all data"


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
                           "source 1 has not been set, consider using set_source() or set_model()"),
                          (ui.save_model,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_resid,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_delchi,  IdentifierErr,
                           "model 1 has not been set")])
def test_save_xxx_nodata(func, etype, emsg, tmp_path):
    """Does save_xxx error out if there's no data to save? DataPHA
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    # Without this we can't check save_staterror or save_error
    ui.set_stat("chi2")

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(etype) as exc:
        func(str(outfile))

    assert str(exc.value) == emsg


@requires_fits
def test_save_image_nodata(tmp_path):
    """Does save_image error out if there's no data to save? DataPHA

    Unlike the other calls, this requires a FITS library
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(IOErr) as exc:
        ui.save_image(str(outfile))

    assert str(exc.value) == "data set '' does not contain an image"


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
def test_save_xxx_bkg_nodata(func, etype, emsg, tmp_path):
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
    with pytest.raises(etype) as exc:
        func(str(outfile), bkg_id=1)

    assert str(exc.value) == emsg


@pytest.mark.parametrize("func,etype,emsg",
                         [(ui.save_filter, DataErr,
                           "data set '1' has no filter"),
                          (ui.save_grouping, ArgumentErr,
                           "data set 1 does not contain PHA data"),
                          (ui.save_quality,  ArgumentErr,
                           "data set 1 does not contain PHA data"),
                          (ui.save_source,  IdentifierErr,
                           "source 1 has not been set, consider using set_source() or set_model()"),
                          (ui.save_model,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_resid,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_delchi,  IdentifierErr,
                           "model 1 has not been set")])
def test_save_xxx_data1d_nodata(func, etype, emsg, tmp_path):
    """Does save_xxx error out if there's no data to save? Data1D
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(etype) as exc:
        func(str(outfile))

    assert str(exc.value) == emsg


@requires_fits
def test_save_image_data1d_nodata(tmp_path):
    """Does save_image error out if there's no data to save? Data1D

    Unlike the other cases we need a FITS backend.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    outfile = tmp_path / "temp-file-that-should-not-be-created"
    with pytest.raises(IOErr) as exc:
        ui.save_image(str(outfile))

    assert str(exc.value) == "data set '' does not contain an image"


@pytest.mark.parametrize("savefunc", [ui.save_data,
                                      ui.save_filter,
                                      ui.save_grouping,
                                      ui.save_quality,
                                      ui.save_image,
                                      ui.save_resid])
def test_save_xxx_not_a_string(savefunc, tmp_path):
    """Does save_xxx error out if not given a filename?

    There are times it would be nice to give a non-string filename
    (pathlib.Path or StringIO), but this is currently not supported.

    At the moment this check is generic for some save commands,
    so just worry about a restricted set of commands.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)
    out = tmp_path / "data.dat"

    with pytest.raises(ArgumentTypeErr) as exc:
        savefunc(out)

    assert str(exc.value) == "'filename' must be a string"


def test_save_delchi_image_fails(tmp_path):
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
    with pytest.raises(AttributeError) as exc:
        ui.save_delchi(str(out))

    assert str(exc.value) == "save_delchi() does not apply for images"


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
    with pytest.raises(IOErr) as exc:
        func(str(outpath))

    emsg = str(exc.value)
    assert emsg.startswith("file '")
    assert emsg.endswith("' exists and clobber is not set")

    # Check the file hasn't changed
    new = outpath.read_text()
    assert new == old


@requires_fits
def test_save_data_data1d_no_clobber(tmp_path):
    """save_data: does clobber=False work? Data1D"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)

    out.write_text("some text")
    check_clobber(out, ui.save_data)


@requires_fits
def test_save_data_datapha_no_clobber(tmp_path):
    """save_data: does clobber=False work? DataPHA"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    out = tmp_path / "data.dat"
    out.write_text("some text")
    check_clobber(out, ui.save_data)


@requires_fits
def test_save_pha_no_clobber(tmp_path):
    """save_pha: does clobber=False work?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    out = tmp_path / "data.dat"

    out.write_text("some text")
    check_clobber(out, ui.save_data)


@requires_fits
@pytest.mark.parametrize("writer", [ui.save_data, ui.save_image])
def test_save_image_no_clobber(writer, tmp_path):
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
        assert lines[0] == f"#{cols}"
        lines = lines[1:]
    else:
        raise RuntimeError(f"UNKNOWN I/O BACKEND: {io.backend.__name__}")

    assert lines[-1] == ""
    lines = lines[:-1]

    assert len(lines) == len(rows)
    for l, r in zip(lines, rows):
        got = [float(t) for t in l.split()]
        assert got == pytest.approx(r)


@requires_fits
def test_save_data_data1d_clobber(tmp_path):
    """save_data: does clobber=True work?"""

    ui.load_arrays(1, [1], [5], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)

    out.write_text("some text")

    ui.save_data(outfile, clobber=True)
    cts = out.read_text()
    check_output(cts, ["X", "Y"], [[1, 5]])


@requires_fits
def test_save_data_data1d(tmp_path):
    """Does save_data work for Data1D?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ["X", "Y"], [[1, 5], [2, 4], [3, 3]])


@requires_fits
def test_save_data_data1d_fits(tmp_path):
    """Does save_data work for Data1D? FITS"""

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile, ascii=False)

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {"X": [1, 2, 3],
                 "Y": [5, 4, 3]})


@requires_fits
def test_save_data_data1dint(tmp_path):
    """Does save_data work for Data1DInt?"""

    ui.load_arrays(1, [1, 2, 4], [2, 3, 5], [5, 4, 3], ui.Data1DInt)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ["XLO", "XHI", "Y"],
                 [[1, 2, 5], [2, 3, 4], [4, 5, 3]])


@requires_fits
def test_save_data_data1dint_fits(tmp_path):
    """Does save_data work for Data1DInt? FITS"""

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [1, 2, 4], [2, 3, 5], [5, 4, 3], ui.Data1DInt)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile, ascii=False)

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {"XLO": [1, 2, 4],
                 "XHI": [2, 3, 5],
                 "Y": [5, 4, 3]})


@requires_fits
def test_save_data_data2d(tmp_path):
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
        raise RuntimeError(f"UNKNOWN I/O BACKEND: {io.backend.__name__}")

    assert cts == expected


@requires_fits
def test_save_data_data2d_fits(tmp_path):
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
def test_save_data_dataimg(tmp_path):
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
def test_save_data_dataimg_fits(tmp_path):
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
@requires_fits  # assume this means we have WCS too
@pytest.mark.parametrize("writer", [ui.save_image, ui.save_data])
def test_save_image_dataimg_fits_wcs(writer, make_data_path, tmp_path):
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
def test_save_data_datapha(tmp_path):
    """Does save_data work for DataPHA?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ["CHANNEL", "COUNTS"],
                 [[1, 5], [2, 4], [3, 3]])


@requires_fits
def test_save_pha(tmp_path):
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
@pytest.mark.parametrize("savefunc,mtype", [(ui.save_source, "SOURCE"),
                                            (ui.save_model, "MODEL")])
def test_save_model_fits(savefunc, mtype, clean_astro_ui, tmp_path):
    """Can we write out data for save_source/model? Data1D and FITS

    As this is not a PHA dataset, the two shouldbe the same bar the
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

    ans = read_table_blocks(outfile)
    blocks = ans[1]
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

    ans = read_table_blocks(outfile)
    blocks = ans[1]
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

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {"XLO": [0.1, 0.2],
                 "XHI": [0.2, 0.4],
                 "MODEL": [20, 40]})


@requires_fits
def test_save_resid_data1d(tmp_path):
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
def test_save_resid_data1d_fits(tmp_path):
    """Residual, Data1D, FITS"""

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [100, 200], [20, 230], ui.Data1D)
    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 220

    out = tmp_path / "resid.out"
    outfile = str(out)
    ui.save_resid(outfile)

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {"X": [100, 200],
                 "RESID": [-200, 10]})


@requires_fits
def test_save_resid_datapha(tmp_path):
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
def test_save_resid_datapha_fits(tmp_path):
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

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {"X": [0.15, 0.3],
                 "RESID": [30, 10]})


@requires_fits
def test_save_resid_dataimg(tmp_path):
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
def test_save_resid_dataimg_fits(tmp_path):
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
    with pytest.raises(ModelErr) as exc:
        ui.get_bkg_source(bkg_id=2)

    assert str(exc.value) == "background model 2 for data set 1 has not been set"

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

    with pytest.raises(ModelErr) as exc:
        ui.get_bkg_source()

    assert str(exc.value) == "background model 1 for data set x has not been set"


def test_fix_background_id_error_checks1():
    """Check error handling of background id"""

    ui.load_arrays(2, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(2, bkg)

    with pytest.raises(ArgumentTypeErr) as exc:
        ui.get_bkg_source(id=2, bkg_id=bkg)

    assert str(exc.value) == "identifiers must be integers or strings"


def test_fix_background_id_error_checks2():
    """Check error handling of background id"""

    ui.load_arrays(2, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(2, bkg)

    with pytest.raises(IdentifierErr) as exc:
        ui.get_bkg_source(id=2, bkg_id="bkg")

    assert str(exc.value) == "identifier 'bkg' is a reserved word"


@pytest.mark.parametrize("id", [1, "x"])
def test_delete_bkg_model_with_bkgid(id, clean_astro_ui):
    """Check we call delete_bkg_model with non-default bkg_id"""

    ui.load_arrays(id, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA("bkg", np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(id, bkg, bkg_id=2)

    ui.set_bkg_source(id, ui.const1d.bmdl, bkg_id=2)
    assert ui.list_model_components() == ["bmdl"]
    assert ui.get_bkg_source(id, 2).name == "const1d.bmdl"

    ui.delete_bkg_model(id, bkg_id=2)
    assert ui.list_model_components() == ["bmdl"]

    with pytest.raises(ModelErr) as exc:
        ui.get_bkg_source(id, 2)

    emsg = f"background model 2 for data set {id} has not been set"
    assert str(exc.value) == emsg


@requires_fits
@pytest.mark.parametrize("loadfunc", [ui.load_grouping, ui.load_quality])
def test_load_xxx_no_data(loadfunc, clean_astro_ui, tmp_path):
    """What happens when there's no data??"""

    path = tmp_path / "data"
    path.write_text("1\n0\n0\n")

    with pytest.raises(IdentifierErr) as exc:
        loadfunc(2, str(path))

    assert str(exc.value) == "data set 2 has not been set"


@requires_fits
@pytest.mark.parametrize("loadfunc", [ui.load_grouping, ui.load_quality])
def test_load_xxx_not_pha(loadfunc, clean_astro_ui, tmp_path):
    """What happens with a non-PHA dataset?"""

    ui.load_arrays(2, [1, 2, 3], [2, 4, 9])

    path = tmp_path / "data"
    path.write_text("1\n0\n0\n")

    with pytest.raises(ArgumentErr) as exc:
        loadfunc(2, str(path))

    assert str(exc.value) == "data set 2 does not contain PHA data"


@requires_fits
@pytest.mark.parametrize("idval", [None, 1, "xx"])
def test_load_grouping(idval, clean_astro_ui, tmp_path):
    """Simple grouping check"""

    x = [1, 2, 3]
    y = [0, 4, 3]
    if idval is None:
        ui.load_arrays(1, x, y, ui.DataPHA)
    else:
        ui.load_arrays(idval, x, y, ui.DataPHA)

    path = tmp_path / "group.dat"
    path.write_text("1\n-1\n1")

    data = ui.get_data(idval)
    assert data.grouping is None

    if idval is None:
        ui.load_grouping(str(path))
    else:
        ui.load_grouping(idval, str(path))

    assert not data.grouped
    assert data.grouping is not None

    ui.group(idval)

    assert data.grouped

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


@pytest.mark.parametrize("idval", [None, 1, "xx"])
def test_group_already_grouped(idval):
    """Does group still work if the data is already grouped?"""

    x = [1, 2, 3]
    y = [0, 4, 3]
    if idval is None:
        ui.load_arrays(1, x, y, ui.DataPHA)
        ui.set_grouping([1, -1, 1])
    else:
        ui.load_arrays(idval, x, y, ui.DataPHA)
        ui.set_grouping(idval, [1, -1, 1])

    data = ui.get_data(idval)
    assert not data.grouped

    ui.group(idval)
    assert data.grouped
    assert ui.get_dep(idval) == pytest.approx([2, 3])

    ui.group(idval)
    assert ui.get_dep(idval) == pytest.approx([2, 3])
    assert data.grouped


@pytest.mark.parametrize("idval", [None, 1, "xx"])
def test_subtract_already_subtracted(idval):
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
def test_xxx_not_pha(callfunc):
    """Just check that these commands require a PHA dataset"""

    ui.load_arrays(1, [1, 2, 3], [4, 5, 6])

    with pytest.raises(ArgumentErr) as exc:
        callfunc()

    assert str(exc.value) == "data set 1 does not contain PHA data"


def test_get_axes_data1d():
    ui.load_arrays(1, [2, 10, 20], [1, 2, 3])
    ax = ui.get_axes()
    assert len(ax) == 1
    assert ax[0] == pytest.approx([2, 10, 20])


def test_get_axes_data1dint():
    ui.load_arrays(1, [2, 10, 20], [10, 12, 22], [1, 2, 3], ui.Data1DInt)
    ax = ui.get_axes()
    assert len(ax) == 2
    assert ax[0] == pytest.approx([2, 10, 20])
    assert ax[1] == pytest.approx([10, 12, 22])


def test_get_axes_datapha_no_response():
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


def test_get_axes_datapha_rmf():
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


def test_get_axes_datapha_arf():
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


def test_get_axes_datapha_rsp():
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


def test_get_axes_dataimg_logical():
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

    with pytest.raises(StatErr) as exc:
        ui.get_stat_info()

    assert str(exc.value) == "No background data has been supplied. Use cstat"


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

    # This does not throw an error
    getattr(ui, f"set_{label}")(vals)


@pytest.mark.parametrize("label", ["grouping", "quality"])
@pytest.mark.parametrize("vals", [[1, 2, 3], np.ones(20)])
def test_pha_column_has_correct_size_sequence(label, vals, clean_astro_ui):
    """Check that the label column is the right size: sequence

    This is related to issue #1160.
    """

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.DataPHA)

    # This does not throw an error
    getattr(ui, f"set_{label}")(vals)


@pytest.mark.parametrize("label", ["grouping", "quality"])
def test_pha_column_checks_pha(label, clean_astro_ui):
    """Check we error out with a non-PHA class"""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.Data1D)

    with pytest.raises(ArgumentErr) as err:
        getattr(ui, f"get_{label}")()

    assert str(err.value) == "data set 1 does not contain PHA data"


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



def test_pha_what_does_get_dep_return_when_grouped(clean_astro_ui):
    """Regression test for get_dep with grouped data"""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [5, 4, 2, 3, 7], ui.DataPHA)
    ui.set_grouping([1, -1, 1, -1, -1])
    ui.group()

    # Looks like it's returning mean of channel values in group
    assert ui.get_dep() == pytest.approx([4.5, 4])


@requires_fits
@requires_data
def test_image_filter_coord_change_same(make_data_path, clean_astro_ui):
    """What happens to the mask after a coordinate change? NO CHANGE

    This is really just a way to test the DataIMG class without having
    to set up all the values. There was some interesting behavior with
    calc_data_sum2d which was down to whether changing the coord would
    change the mask, so this is an explicit test of that behavior.

    """

    ui.load_image("foo", make_data_path("image2.fits"))
    assert ui.get_filter("foo") == ""

    ui.set_coord("foo", "physical")
    ui.notice2d_id("foo", "rect(4000, 4200, 4100 , 4300 ) ")
    assert ui.get_filter("foo") == "Rectangle(4000,4200,4100,4300)"

    # Is there no way to get the mask data via the UI interface,
    # without calling save_filter?
    #
    d = ui.get_data("foo")
    assert d.mask.sum() == 2500

    ui.set_coord("foo", "physical")
    assert ui.get_filter("foo") == "Rectangle(4000,4200,4100,4300)"
    assert d.mask.sum() == 2500


@requires_fits
@requires_data
def test_image_filter_coord_change(make_data_path, clean_astro_ui):
    """What happens to the mask after a coordinate change?

    This is the more-general case of test_image_filter_coord_change_same
    but it's not actually clear what should be done here when the
    coordinate setting is changed - should the mask be removed as
    we don't store the mask expression with a way to update it to
    match the new coordianate setting.

    """

    ui.load_image("foo", make_data_path("image2.fits"))
    assert ui.get_filter("foo") == ""

    ui.set_coord("foo", "physical")
    ui.notice2d_id("foo", "rect(4000, 4200, 4100 , 4300 ) ")
    assert ui.get_filter("foo") == "Rectangle(4000,4200,4100,4300)"

    # Is there no way to get the mask data via the UI interface,
    # without calling save_filter?
    #
    d = ui.get_data("foo")
    assert d.mask.sum() == 2500

    ui.set_coord("foo", "logical")

    # This expression is no-longer technically valid, but we do
    # not change it yet.
    #
    assert ui.get_filter("foo") == "Rectangle(4000,4200,4100,4300)"

    # We know what the mask was, so assume it hasn't changed.
    # We could also just clear the mask so d.mask would be True.
    #
    assert d.mask.sum() == 2500
