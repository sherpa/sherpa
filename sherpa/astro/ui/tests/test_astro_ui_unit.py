#
#  Copyright (C) 2017, 2018, 2020, 2021  Smithsonian Astrophysical Observatory
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
from sherpa.astro import ui
from sherpa.utils import poisson_noise
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, DataErr, \
    IdentifierErr, IOErr, ModelErr
from sherpa.utils.testing import requires_fits


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
    for expected in ['mh', 'metropolismh', 'pragbayes', 'fullbayes']:
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

    cpt = ui.voigt1d('cpt')
    cpt.pos = 2345
    cpt.fwhm_g = 20
    cpt.fwhm_l = 15
    cpt.ampl = 480

    np.random.seed(72383)
    y = poisson_noise(cpt(x))

    ui.load_arrays(1, x, y, ui.Data1D)
    ui.set_stat('leastsq')
    ui.set_method('simplex')
    ui.delete_model_component('cpt')


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
                           [24.290859063445527, 12.35931466539915, 2343.751436403753, 437.0342780289447])
                         ])
def test_fit_profile(model, stat, pars, reset_seed, clean_astro_ui):
    """Regression test simple 1D fits"""

    setup_data1d_fit()
    ui.set_source(model('mdl'))
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
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)
    with pytest.raises(DataErr) as exc:
        ui.set_filter(f, bkg_id=bid)

    assert str(exc.value) == 'size mismatch between 3 and 2'


@pytest.mark.parametrize("f", [[False, True], np.asarray([False, True])])
@pytest.mark.parametrize("bid", [None, 1])
def test_set_filter_mismatch_with_filter(f, bid):
    """Does set_filter error when there's a mis-match after a filter?

    test_set_filter_mismatch checks when .mask is a scalar,
    so now check when it's a NumPy array.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    ui.ignore(3, None)  # set the .mask attribute to an array

    with pytest.raises(DataErr) as exc:
        ui.set_filter(f, bkg_id=bid)

    assert str(exc.value) == 'size mismatch between 3 and 2'


@pytest.mark.parametrize("bid", [None, 1])
def test_get_syserror_missing(bid):
    """Does get_syserror error out?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)
    with pytest.raises(DataErr) as exc:
        ui.get_syserror(bkg_id=bid)

    assert str(exc.value) == "data set '1' does not specify systematic errors"


@pytest.mark.parametrize("bid", [None, 1])
def test_save_filter_ignored(bid):
    """Does save_filter error out if everything is masked?

    We should be able to write out the filter in this case,
    as it's easy (all False).
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    ui.ignore(None, None)

    with pytest.raises(DataErr) as exc:
        ui.save_filter("temp-file-that-should-not-be-created",
                       bkg_id=bid)

    assert str(exc.value) == "mask excludes all data"


@pytest.mark.parametrize("func,etype,emsg",
                         [(ui.save_filter, DataErr,
                           "data set '1' has no filter"),
                          (ui.save_grouping, DataErr,
                           "data set '1' does not specify grouping flags"),
                          (ui.save_quality,  DataErr,
                           "data set '1' does not specify quality flags"),
                          (ui.save_source,  IdentifierErr,
                           "source 1 has not been set, consider using set_source() or set_model()"),
                          (ui.save_model,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_resid,  IdentifierErr,
                           "model 1 has not been set"),
                          (ui.save_delchi,  IdentifierErr,
                           "model 1 has not been set")])
def test_save_xxx_nodata(func, etype, emsg):
    """Does save_xxx error out if there's no data to save? DataPHA
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    with pytest.raises(etype) as exc:
        func("temp-file-that-should-not-be-created")

    assert str(exc.value) == emsg


@requires_fits
def test_save_image_nodata():
    """Does save_image error out if there's no data to save? DataPHA

    Unlike the other calls, this requires a FITS library
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    with pytest.raises(IOErr) as exc:
        ui.save_image("temp-file-that-should-not-be-created")

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
                          (ui.save_model,  ModelErr,
                           "background model 1 for data set 1 has not been set"),
                          (ui.save_resid,  ModelErr,
                           "background model 1 for data set 1 has not been set"),
                          (ui.save_delchi,  ModelErr,
                           "background model 1 for data set 1 has not been set")])
def test_save_xxx_bkg_nodata(func, etype, emsg):
    """Does save_xxx error out if there's no data to save? DataPHA + bkg

    Note that save_image does not support a bkg_id parameter so is
    not tested here.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    with pytest.raises(etype) as exc:
        func("temp-file-that-should-not-be-created", bkg_id=1)

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
def test_save_xxx_data1d_nodata(func, etype, emsg):
    """Does save_xxx error out if there's no data to save? Data1D
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    with pytest.raises(etype) as exc:
        func("temp-file-that-should-not-be-created")

    assert str(exc.value) == emsg


@requires_fits
def test_save_image_data1d_nodata():
    """Does save_image error out if there's no data to save? Data1D

    Unlike the other cases we need a FITS backend.
    """

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    with pytest.raises(IOErr) as exc:
        ui.save_image("temp-file-that-should-not-be-created")

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

    out = tmp_path / 'does-not-exist'
    with pytest.raises(AttributeError) as exc:
        ui.save_delchi(str(out))

    assert str(exc.value) == 'save_delchi() does not apply for images'


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


@pytest.mark.xfail(reason='fall through does not work')
def test_save_data_data1d_no_clobber(tmp_path):
    """save_data: does clobber=False work? Data1D"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)

    out.write_text('some text')
    check_clobber(out, ui.save_data)


@pytest.mark.xfail(reason='fall through does not work')
@requires_fits
def test_save_data_datapha_no_clobber(tmp_path):
    """save_data: does clobber=False work? DataPHA"""

    # This import requires an I/O backend hence it is done here
    #
    from sherpa.astro.io.meta import Meta

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)

    # The code requires DataPHA to have a "valid" header so
    # fake one. Ideally we would not require it.
    #
    hdr = Meta()
    ui.get_data().header = hdr

    out = tmp_path / "data.dat"

    out.write_text('some text')
    check_clobber(out, ui.save_data)


@pytest.mark.xfail(reason='fall through does not work')
@requires_fits
def test_save_pha_no_clobber(tmp_path):
    """save_pha: does clobber=False work?"""

    # This import requires an I/O backend hence it is done here
    #
    from sherpa.astro.io.meta import Meta

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)

    # The code requires DataPHA to have a "valid" header so
    # fake one. Ideally we would not require it.
    #
    hdr = Meta()
    ui.get_data().header = hdr

    out = tmp_path / "data.dat"

    out.write_text('some text')
    check_clobber(out, ui.save_data)


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

    from sherpa.astro.io import backend

    lines = out.split('\n')
    assert len(lines) > 2

    cols = ' '.join(colnames)
    if backend.__name__ == 'sherpa.astro.io.crates_backend':
        assert lines[0] == "#TEXT/SIMPLE"
        assert lines[1] == f"# {cols}"
        lines = lines[2:]
    elif backend.__name__ == 'sherpa.astro.io.pyfits_backend':
        assert lines[0] == f"#{cols}"
        lines = lines[1:]
    else:
        raise RuntimeError(f"UNKNOWN I/O BACKEND: {backend.__name__}")

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

    out.write_text('some text')

    ui.save_data(outfile, clobber=True)
    cts = out.read_text()
    check_output(cts, ['X', 'Y'], [[1, 5]])


@requires_fits
def test_save_data_data1d(tmp_path):
    """Does save_data work for Data1D?"""

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.Data1D)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ['X', 'Y'], [[1, 5], [2, 4], [3, 3]])


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
                {'X': [1, 2, 3],
                 'Y': [5, 4, 3]})


@requires_fits
def test_save_data_data1dint(tmp_path):
    """Does save_data work for Data1DInt?"""

    ui.load_arrays(1, [1, 2, 4], [2, 3, 5], [5, 4, 3], ui.Data1DInt)

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ['XLO', 'XHI', 'Y'],
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
                {'XLO': [1, 2, 4],
                 'XHI': [2, 3, 5],
                 'Y': [5, 4, 3]})


@requires_fits
def test_save_data_data2d(tmp_path):
    """Does save_data work for Data2D?"""

    from sherpa.astro.io import backend

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
    if backend.__name__ == 'sherpa.astro.io.crates_backend':
        expected = ["#TEXT/SIMPLE", "# X0 X1 Y SHAPE"]
        s = [2, 3, 0, 0, 0, 0]
        for xi, yi, zi, si in zip(x, y, z, s):
            expected.append(f"{xi} {yi} {zi} {si}")

        expected = "\n".join(expected) + "\n"

    elif backend.__name__ == 'sherpa.astro.io.pyfits_backend':
        expected = "\n".join([str(zz) for zz in z]) + "\n"
    else:
        raise RuntimeError(f"UNKNOWN I/O BACKEND: {backend.__name__}")

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
    from sherpa.astro.io import backend
    if backend.__name__ == 'sherpa.astro.io.crates_backend':
        pytest.skip('ASCII not supported for images with pycrates')

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


@requires_fits
def test_save_data_datapha(tmp_path):
    """Does save_data work for DataPHA?"""

    # This import requires an I/O backend hence it is done here
    #
    from sherpa.astro.io.meta import Meta

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)

    # The code requires DataPHA to have a "valid" header so
    # fake one. Ideally we would not require it.
    #
    hdr = Meta()
    ui.get_data().header = hdr

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_data(outfile)

    cts = out.read_text()
    check_output(cts, ['CHANNEL', 'COUNTS'],
                 [[1, 5], [2, 4], [3, 3]])


@requires_fits
def test_save_pha(tmp_path):
    """Does save_pha work for DataPHA?"""

    # This import requires an I/O backend hence it is done here
    #
    from sherpa.astro.io.meta import Meta

    ui.load_arrays(1, [1, 2, 3], [5, 4, 3], ui.DataPHA)

    # The code requires DataPHA to have a "valid" header so
    # fake one. Ideally we would not require it.
    #
    hdr = Meta()
    ui.get_data().header = hdr

    out = tmp_path / "data.dat"
    outfile = str(out)
    ui.save_pha(outfile)

    # just check the first line; the output may depend on the FITS backend
    cts = out.read_text()[:80].split()
    assert cts[0] == 'SIMPLE'
    assert cts[1] == '='
    assert cts[2] == 'T'


@requires_fits
@pytest.mark.parametrize("savefunc,mtype", [(ui.save_source, 'SOURCE'),
                                            (ui.save_model, 'MODEL')])
def test_save_model_ascii(savefunc, mtype, clean_astro_ui, tmp_path):
    """Can we write out data for save_source/model? Data1D and ASCII

    As this is not a PHA dataset, the two should be the same bar the
    header line.
    """

    ui.load_arrays(1, [1, 2], [5, 10])
    ui.set_source(ui.polynom1d.cmdl)
    cmdl.c0 = -5
    cmdl.c1 = 10

    out = tmp_path / 'model.dat'
    savefunc(str(out), ascii=True)

    cts = out.read_text()
    check_output(cts, ['X', mtype], [[1, 5], [2, 15]])


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

    out = tmp_path / 'model.dat'
    ui.save_source(str(out), ascii=True)

    cts = out.read_text()
    check_output(cts, ['XLO', 'XHI', 'SOURCE'],
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

    out = tmp_path / 'model.dat'
    ui.save_model(str(out), ascii=True)

    cts = out.read_text()
    check_output(cts, ['XLO', 'XHI', 'MODEL'],
                 [[0.1, 0.2, 20], [0.2, 0.4, 40]])


@requires_fits
@pytest.mark.parametrize("savefunc,mtype", [(ui.save_source, 'SOURCE'),
                                            (ui.save_model, 'MODEL')])
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

    out = tmp_path / 'model.dat'
    outfile = str(out)
    savefunc(outfile)

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {'X': [1, 2],
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

    out = tmp_path / 'model.dat'
    outfile = str(out)
    ui.save_source(outfile)

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {'XLO': [0.1, 0.2],
                 'XHI': [0.2, 0.4],
                 'SOURCE': [2, 2]})


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

    out = tmp_path / 'model.dat'
    outfile = str(out)
    ui.save_model(outfile)

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {'XLO': [0.1, 0.2],
                 'XHI': [0.2, 0.4],
                 'MODEL': [20, 40]})


@requires_fits
def test_save_resid_data1d(tmp_path):
    """Residual, Data1D, ASCII"""

    ui.load_arrays(1, [100, 200], [20, 230], ui.Data1D)
    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 220

    out = tmp_path / 'resid.out'
    outfile = str(out)
    ui.save_resid(outfile, ascii=True)

    cts = out.read_text()
    check_output(cts, ['X', 'RESID'], [[100, -200], [200, 10]])


@requires_fits
def test_save_resid_data1d_fits(tmp_path):
    """Residual, Data1D, FITS"""

    from sherpa.astro.io import read_table_blocks

    ui.load_arrays(1, [100, 200], [20, 230], ui.Data1D)
    ui.set_source(ui.const1d.cmdl)
    cmdl.c0 = 220

    out = tmp_path / 'resid.out'
    outfile = str(out)
    ui.save_resid(outfile)

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {'X': [100, 200],
                 'RESID': [-200, 10]})


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

    out = tmp_path / 'resid.out'
    outfile = str(out)
    ui.save_resid(outfile, ascii=True)

    cts = out.read_text()
    check_output(cts, ['X', 'RESID'], [[0.15, 30], [0.3, 10]])


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

    out = tmp_path / 'resid.out'
    outfile = str(out)
    ui.save_resid(outfile)

    ans = read_table_blocks(outfile)
    blocks = ans[1]
    assert len(blocks) == 2
    check_table(blocks[2],
                {'X': [0.15, 0.3],
                 'RESID': [30, 10]})


@requires_fits
def test_save_resid_dataimg(tmp_path):
    """Residual, DataIMG, ASCII"""

    # Can not write out an ASCII image with crates
    from sherpa.astro.io import backend
    if backend.__name__ == 'sherpa.astro.io.crates_backend':
        pytest.skip('ASCII not supported for images with pycrates')

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
    d = ui.DataPHA('src', channels, counts)
    b = ui.DataPHA('bkg', channels, counts)
    d.set_background(b, id=2)
    ui.set_data(d)

    ui.set_bkg_source(ui.const1d.bmdl + ui.gauss1d.gmdl, bkg_id=2)
    assert ui.get_bkg_source(bkg_id=2) is not None

    ui.delete_bkg_model(bkg_id=2)

    # Expression has been removed
    #
    with pytest.raises(ModelErr) as exc:
        ui.get_bkg_source(bkg_id=2)

    assert str(exc.value) == 'background model 2 for data set 1 has not been set'

    # components still exist
    mdls = ui.list_model_components()
    assert set(mdls) == set(['bmdl', 'gmdl'])


def test_default_background_issue(clean_astro_ui):
    """Test issue #943"""

    ui.set_default_id('x')

    # use least-square as we don't really care about the fit
    ui.set_stat('leastsq')

    ui.load_arrays('x', [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
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

    ui.set_default_id('x')

    ui.load_arrays('x', [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    arf = ui.create_arf(np.asarray([0.1, 0.2, 0.3]), np.asarray([0.2, 0.3, 0.4]))
    bkg.set_arf(arf)
    ui.set_bkg(bkg)

    ui.set_bkg_source(ui.const1d.mdl2)
    ui.show_bkg_model()


def test_default_background_issue_fit(clean_astro_ui):
    """Test issue #943 with fit

    See https://github.com/sherpa/sherpa/issues/943#issuecomment-696119982
    """

    ui.set_default_id('x')

    # use least-square as we don't really care about the fit
    ui.set_stat('leastsq')

    ui.load_arrays('x', [1, 2, 3, 4], [5, 4, 3, 4], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3, 4]), [1, 1, 0, 1])
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

    ui.set_default_id('x')

    ui.load_arrays('x', [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(bkg)

    with pytest.raises(ModelErr) as exc:
        ui.get_bkg_source()

    assert str(exc.value) == 'background model 1 for data set x has not been set'


def test_fix_background_id_error_checks1():
    """Check error handling of background id"""

    ui.load_arrays(2, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(2, bkg)

    with pytest.raises(ArgumentTypeErr) as exc:
        ui.get_bkg_source(id=2, bkg_id=bkg)

    assert str(exc.value) == 'identifiers must be integers or strings'


def test_fix_background_id_error_checks2():
    """Check error handling of background id"""

    ui.load_arrays(2, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(2, bkg)

    with pytest.raises(IdentifierErr) as exc:
        ui.get_bkg_source(id=2, bkg_id='bkg')

    assert str(exc.value) == "identifier 'bkg' is a reserved word"


@pytest.mark.parametrize("id", [1, 'x'])
def test_delete_bkg_model_with_bkgid(id, clean_astro_ui):
    """Check we call delete_bkg_model with non-default bkg_id"""

    ui.load_arrays(id, [1, 2, 3], [5, 4, 3], ui.DataPHA)
    bkg = ui.DataPHA('bkg', np.asarray([1, 2, 3]), [1, 1, 0])
    ui.set_bkg(id, bkg, bkg_id=2)

    ui.set_bkg_source(id, ui.const1d.bmdl, bkg_id=2)
    assert ui.list_model_components() == ['bmdl']
    assert ui.get_bkg_source(id, 2).name == 'const1d.bmdl'

    ui.delete_bkg_model(id, bkg_id=2)
    assert ui.list_model_components() == ['bmdl']

    with pytest.raises(ModelErr) as exc:
        ui.get_bkg_source(id, 2)

    emsg = 'background model 2 for data set {} has not been set'.format(id)
    assert str(exc.value) == emsg


@pytest.mark.parametrize("loadfunc", [ui.load_grouping, ui.load_quality])
@requires_fits
def test_load_xxx_no_data(loadfunc, clean_astro_ui, tmp_path):
    """What happens when there's no data??"""

    path = tmp_path / 'data'
    path.write_text('1\n0\n0\n')

    with pytest.raises(IdentifierErr) as exc:
        loadfunc(2, str(path))

    assert str(exc.value) == 'data set 2 has not been set'


@pytest.mark.parametrize("loadfunc", [ui.load_grouping, ui.load_quality])
@requires_fits
def test_load_xxx_not_pha(loadfunc, clean_astro_ui, tmp_path):
    """What happens with a non-PHA dataset?"""

    ui.load_arrays(2, [1, 2, 3], [2, 4, 9])

    path = tmp_path / 'data'
    path.write_text('1\n0\n0\n')

    with pytest.raises(ArgumentErr) as exc:
        loadfunc(2, str(path))

    assert str(exc.value) == 'data set 2 does not contain PHA data'


@pytest.mark.parametrize("idval", [None, 1, 'xx'])
def test_load_grouping(idval, clean_astro_ui, tmp_path):
    """Simple grouping check"""

    x = [1, 2, 3]
    y = [0, 4, 3]
    if idval is None:
        ui.load_arrays(1, x, y, ui.DataPHA)
    else:
        ui.load_arrays(idval, x, y, ui.DataPHA)

    path = tmp_path / 'group.dat'
    path.write_text('1\n-1\n1')

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
