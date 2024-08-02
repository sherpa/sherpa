#
#  Copyright (C) 2016 - 2018, 2020, 2021, 2023, 2024
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

import struct
import warnings

import numpy as np

import pytest

from sherpa.astro import io
from sherpa.astro import ui
from sherpa.data import Data1D, Data2DInt
from sherpa.models.basic import Box1D, Const1D
from sherpa.utils.err import IOErr
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_group, requires_xspec


@requires_data
@requires_fits
@requires_xspec
def test_mod_fits(make_data_path, clean_astro_ui, caplog):
    """Can we read in an XSPEC table model with load_table_model.

    This approach is deprecated. Use load_xstable_model instead.
    """

    tablemodelfile = make_data_path("xspec-tablemodel-RCS.mod")
    with warnings.catch_warnings(record=True) as warn:
        ui.load_table_model("tmod", tablemodelfile)

    msg = "Use load_xstable_model to load XSPEC table models"
    assert len(warn) == 1
    assert warn[0].category == DeprecationWarning
    assert str(warn[0].message) == msg

    tmod = ui.get_model_component("tmod")
    assert tmod.name == "xstablemodel.tmod"

    assert len(caplog.records) == 1
    assert caplog.records[0].name == "sherpa.astro.ui.utils"
    assert caplog.records[0].levelname == "WARNING"
    assert caplog.records[0].getMessage() == msg


@requires_data
@requires_fits
@requires_xspec
def test_xsmod_fits(make_data_path, clean_astro_ui, caplog):
    """Can we read in an XSPEC table model with load_xstable_model."""

    tablemodelfile = make_data_path("xspec-tablemodel-RCS.mod")
    with warnings.catch_warnings(record=True) as warn:
        ui.load_xstable_model("tmod", tablemodelfile)

    assert len(warn) == 0

    tmod = ui.get_model_component("tmod")
    assert tmod.name == "xstablemodel.tmod"

    # just check the log output is empty
    assert len(caplog.records) == 0


@requires_fits
@pytest.mark.parametrize("flag", [True, False])
def test_warnings_are_gone_arrays(tmp_path, flag):
    ui.load_arrays(1, [1, 2, 3], [4, 5, 6])
    #  We now have logic in conftest.py to catch white-listed warnings and fail on unexpected ones.
    #  We just need to make any warnings bubble up, here and in the following test.
    outpath = tmp_path / "test.dat"
    ui.save_data(1, str(outpath), ascii=flag, clobber=True)


@requires_fits
@requires_data
def test_warnings_are_gone_pha(make_data_path, tmp_path):
    pha = make_data_path("3c273.pi")
    ui.load_pha(pha)
    outpath = tmp_path / "test.pi"
    ui.save_data(1, str(outpath), ascii=False, clobber=True)


def assert_staterr(use_errors):
    assert np.all(ui.get_data("phacounts").counts ==
                  pytest.approx(ui.get_data("pharate").counts))
    if use_errors is True:
        assert np.all(ui.get_data("phacounts").staterror ==
                      pytest.approx(ui.get_data("pharate").staterror))
    else:
        assert ui.get_data("phacounts").staterror is None
        assert ui.get_data("pharate").staterror is None

    for n in ['phacounts', 'pharate']:
        ui.group_bins(n, 16)
    ui.set_analysis('energy')
    ui.ignore(None, 3.)
    if use_errors is True:
        assert np.all(ui.get_data("phacounts").get_error() ==
                      pytest.approx(ui.get_data("pharate").get_error()))
    else:
        assert ui.get_data("phacounts").get_error() is None
        assert ui.get_data("pharate").get_error() is None


@requires_fits
@requires_data
@requires_group
@pytest.mark.parametrize("use_errors", [True, False])
def test_scaling_staterr(make_data_path, use_errors):
    '''Regression test for https://github.com/sherpa/sherpa/issues/800

    Notes
    -----
    Test files are made with (in sherpa-test-data/sherpatest)
    dmtcalc source.pi source_edit.pi expression="stat_err=stat_err/exposure"
    dmcopy "source_edit.pi[cols channel,pi,RATE=count_rate,stat_err]" source_rate.pi clobber=yes
    # Then fix units to counts/s manually in the header
    '''
    ui.load_pha("phacounts", make_data_path("source.pi"),
                use_errors=use_errors)
    ui.load_pha("pharate", make_data_path("source_rate.pi"),
                use_errors=use_errors)
    for n in ['phacounts', 'pharate']:
        ui.load_arf(n, make_data_path("source.arf"))
        ui.load_rmf(n, make_data_path("source.rmf"))

    assert_staterr(use_errors)


@requires_fits
@requires_data
@requires_group
@pytest.mark.parametrize("use_errors", [True, False])
def test_scaling_staterr_pha2(make_data_path, use_errors):
    '''Regression test for https://github.com/sherpa/sherpa/issues/800

    Notes
    -----
    Test files are made with (in sherpa-test-data/sherpatest)
    dmtcalc 3c120_pha2.gz 3c120_pha2_edit.gz expression="RATE=(float)counts/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit.gz 3c120_pha2_edit2.gz expression="stat_err=stat_err/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit2.gz 3c120_pha2_edit3.gz expression="backgroun_up=(float)background_up/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit3.gz 3c120_pha2_edit4.gz expression="backgroun_down=(float)background_down/exposure" clobber=yes

    dmcopy "3c120_pha2_edit4.gz[cols -COUNTS,-BACKGROUND_UP,-BACKGROUND_DOWN][cols *,BACKGROUND_UP=backgroun_up,BACKGROUND_DOWN=backgroun_down][col -COUNTS]" 3c120_pha2_rate.gz clobber=yes
    Then fix units manually.
    '''
    ui.load_pha(make_data_path("3c120_pha2_rate.gz"),
                use_errors=use_errors)
    ui.copy_data(9, "phacounts")
    ui.load_pha(make_data_path("3c120_pha2_rate.gz"),
                use_errors=use_errors)
    ui.copy_data(9, "pharate")
    for n in ['phacounts', 'pharate']:
        ui.load_arf(n, make_data_path("3c120_meg_-1.arf.gz"))
        ui.load_rmf(n, make_data_path("3c120_meg_-1.rmf.gz"))
    assert_staterr(use_errors)


def fake_rmf(outfile):
    """Create a "perfect" RMF with a scalar MATRIX column

    We do this rather than add another file to the test
    data directory as it is easier.
    """

    ebins = np.arange(0.15, 0.2, 0.01, dtype=np.float32)
    elo = ebins[:-1]
    ehi = ebins[1:]
    nchan = elo.size
    chans = np.arange(1, nchan + 1, dtype=np.int16)

    def hdr(key, value):
        if isinstance(value, str):
            value = "'{}'".format(value)
        else:
            if isinstance(value, bool):
                value = 'T' if value else 'F'

            # add spacing to make FVERIFY happy
            value = f"{str(value):>20s}"

        out = "{:8s}= {}".format(key, value)
        return out.ljust(80)

    hdr_img = [hdr('SIMPLE', True),
               hdr('BITPIX', 8),
               hdr('NAXIS', 0),
               hdr('EXTEND', True)]

    hdr_ebounds = [hdr('XTENSION', 'BINTABLE'),
                   hdr('BITPIX', 8),
                   hdr('NAXIS', 2),
                   hdr('NAXIS1', 10),
                   hdr('NAXIS2', nchan),
                   hdr('PCOUNT', 0),
                   hdr('GCOUNT', 1),
                   hdr('TFIELDS', 3),
                   hdr('TTYPE1', 'CHANNEL'),
                   hdr('TFORM1', 'I'),
                   hdr('TTYPE2', 'E_MIN'),
                   hdr('TFORM2', 'E'),
                   hdr('TTYPE3', 'E_MAX'),
                   hdr('TFORM3', 'E'),
                   hdr('EXTNAME', 'EBOUNDS'),
                   hdr('HDUCLASS', 'OGIP'),
                   hdr('HDUCLAS1', 'RESPONSE'),
                   hdr('HDUCLAS2', 'EBOUNDS'),
                   hdr('HDUVERS', '1.3.0'),
                   hdr('CHANTYPE', 'PI'),
                   hdr('DETCHANS', nchan)]

    hdr_matrix = [hdr('XTENSION', 'BINTABLE'),
                  hdr('BITPIX', 8),
                  hdr('NAXIS', 2),
                  hdr('NAXIS1', 18),
                  hdr('NAXIS2', nchan),
                  hdr('PCOUNT', 0),
                  hdr('GCOUNT', 1),
                  hdr('TFIELDS', 6),
                  hdr('TTYPE1', 'ENERG_LO'),
                  hdr('TFORM1', 'E'),
                  hdr('TTYPE2', 'ENERG_HI'),
                  hdr('TFORM2', 'E'),
                  hdr('TTYPE3', 'N_GRP'),
                  hdr('TFORM3', 'I'),
                  hdr('TTYPE4', 'F_CHAN'),
                  hdr('TFORM4', 'I'),
                  hdr('TTYPE5', 'N_CHAN'),
                  hdr('TFORM5', 'I'),
                  hdr('TTYPE6', 'MATRIX'),
                  hdr('TFORM6', 'E'),
                  hdr('EXTNAME', 'MATRIX'),
                  hdr('HDUCLASS', 'OGIP'),
                  hdr('HDUCLAS1', 'RESPONSE'),
                  hdr('HDUCLAS2', 'RSP_MATRIX'),
                  hdr('HDUCLAS3', 'REDIST'),
                  hdr('HDUVERS', '1.3.0'),
                  hdr('CHANTYPE', 'PI'),
                  hdr('TLMIN4', 1),
                  hdr('DETCHANS', nchan),
                  hdr('NUMGRP', nchan),
                  hdr('NUMELT', nchan),
                  hdr('LO_THRES', '1E-06')]

    ngrps = np.ones(nchan, dtype=np.int16)
    fchans = np.arange(1, nchan + 1, dtype=np.int16)
    nchans = ngrps
    matrix = np.ones(nchan, dtype=np.float32)

    def row1(chan, e1, e2):
        return struct.pack('>hff', chan, e1, e2)

    def row2(e1, e2, ngrp, fchan, nchan, mat):
        return struct.pack('>ffhhhf', e1, e2, ngrp, fchan, nchan, mat)

    with open(outfile, 'wb') as fh:

        def extend():
            """extend to next 2880 block"""
            p = fh.tell()
            n = p % 2880
            if n != 0:
                fh.write(b' ' * (2880 - n))

        def mkhdr(cards):
            for card in cards:
                fh.write(card.encode('ascii'))

            fh.write(b'END' + b' ' * 77)
            extend()

        mkhdr(hdr_img)

        mkhdr(hdr_ebounds)
        for row in zip(chans, elo, ehi):
            fh.write(row1(*row))

        extend()

        mkhdr(hdr_matrix)
        for row in zip(elo, ehi, ngrps, fchans, nchans, matrix):
            fh.write(row2(*row))

        extend()


@requires_fits
def test_read_ideal_rmf(tmp_path):
    """Can a RMF similar to issue #862 be read in?

    The MATRIX column in this file is a scalar rather than array,
    and let's do EBOUNDS then MATRIX blocks.
    """

    ebins = np.arange(0.15, 0.2, 0.01)
    elo = ebins[:-1]
    ehi = ebins[1:]

    outpath = tmp_path / "test.rmf"
    outfile = str(outpath)
    fake_rmf(outfile)
    r = io.read_rmf(outfile)

    # Can we read in the data
    #
    assert r.detchans == 5
    assert r.energ_lo == pytest.approx(elo)
    assert r.energ_hi == pytest.approx(ehi)
    assert (r.n_grp == [1, 1, 1, 1, 1]).all()
    assert (r.f_chan == [1, 2, 3, 4, 5]).all()
    assert (r.n_chan == [1, 1, 1, 1, 1]).all()
    assert r.offset == 1
    assert r.e_min == pytest.approx(elo)
    assert r.e_max == pytest.approx(ehi)
    assert r.ethresh == 1e-10

    # Can we apply it?
    #
    # The cmdl evaluates to a value of 2 * bin width
    # The bmdl evalates to the bin width * x
    # where x = [0, 1, 0.5, 0, 0]
    #
    cmdl = Const1D()
    cmdl.c0 = 2

    bmdl = Box1D()
    bmdl.xlow = 0.16
    bmdl.xhi = 0.175

    mdl = bmdl + cmdl

    # Multiply by 100 so numbers are close to unity
    expected = 100 * 0.01 * np.asarray([2, 3, 2.5, 2, 2])
    y = 100 * r.eval_model(mdl)
    assert y == pytest.approx(expected, rel=2e-6)


@requires_fits
@requires_data
def test_fits_file_lower_case(make_data_path):
    """Caused issue #143

    The file contains

        MTYPE1       = sky                  / DM Keyword: Descriptor name.
        MFORM1       = X,Y                  / [pixel]
        MTYPE2       = EQPOS                / DM Keyword: Descriptor name.
        MFORM2       = RA,Dec               / [deg]

    so it has - for the transformed case - a column name that
    is not upper case.

    """

    infile = make_data_path("1838_rprofile_rmid.fits")
    tbl = ui.unpack_table(infile, colkeys=["RA", "Dec"])

    assert isinstance(tbl, Data1D)
    assert len(tbl.x) == len(tbl.y)
    assert len(tbl.x) == 38
    assert tbl.staterror is None
    assert tbl.syserror is None

    # Every point is the same, which makes it easy to check
    #
    assert (tbl.x == tbl.x[0]).all()
    assert (tbl.y == tbl.y[0]).all()

    assert tbl.x[0] == pytest.approx(278.3897960639)
    assert tbl.y[0] == pytest.approx(-10.5690222237)


@requires_fits
@requires_data
def test_fits_file_missing_column(make_data_path):
    """Follow on from #143

    Ensure we try to access a missing column. This is low-level
    (i.e. calls a sherpa.astro.io routine) so that we aren't bothered
    with the cascading error fall through of the ui layer code.

    """

    infile = make_data_path("1838_rprofile_rmid.fits")

    # The error message depends on the backend
    # - crates lists the available columns
    # - pyfits lists the filename
    # so just check the common part. Unfortunately the name of the
    # column depends on the backend too; crates uses the user value
    # whereas pyfits converts it to a capital. The test could have used
    # "FOO" as the column name to ignore this but let's keep this.
    #
    emsg = "^Required column 'F(oo|OO)' not found in "
    with pytest.raises(IOErr, match=emsg):
        io.read_table(infile, colkeys=["ra", "Foo"])


@requires_fits
@requires_data
def test_get_header_data_missing_key(make_data_path):
    """What happens if a requested key is missing?

    TODO: If get_header_data is useful should we export if from
    the io, and not backend, level?
    """

    infile = make_data_path("1838_rprofile_rmid.fits")
    with pytest.raises(IOErr,
                       match=" does not have a 'NOTAKEYWORD' keyword$"):
        io.backend.get_header_data(infile, hdrkeys=["NOTAKEYWORD"])


@requires_fits
def test_set_arrays_not_sequence_of_seqence(tmp_path):
    """Check this error condition.

    TODO: Should set_arrays be part of the backend API?
    """

    outpath = tmp_path / 'do-not-edit.dat'
    with pytest.raises(IOErr,
                       match=r"^please supply array\(s\) to write to file$"):
        io.backend.set_arrays(str(outpath),
                              args=[np.arange(3), True],
                              clobber=True)


def test_read_arrays_no_data():
    """This can run even with the dummy backend"""

    emsg = "^no arrays found to be loaded$"
    with pytest.raises(IOErr, match=emsg):
        io.read_arrays()


@requires_fits
def test_read_arrays_no_data_but_dstype():
    """This is a slightly-different error path than read_arrays()"""

    with pytest.raises(IOErr,
                       match="^no arrays found to be loaded$"):
        io.read_arrays(Data2DInt)


@requires_fits
@pytest.mark.parametrize("dstype,dname,nargs",
                         [(Data1D, 'Data1D', 2),
                          (Data2DInt, 'Data2DInt', 5)])
def test_read_arrays_not_enough_data(dstype, dname, nargs):

    emsg = f"^data set '{dname}' takes at least {nargs} args$"
    with pytest.raises(TypeError, match=emsg):
        io.read_arrays([1, 2, 3], dstype)


@requires_fits
def test_read_arrays_not_an_array_type():

    emsg = "^'foo' must be a Numpy array, list, or tuple$"
    with pytest.raises(IOErr, match=emsg):
        io.read_arrays("foo")


@requires_fits
def test_read_arrays_data1d():

    dset = io.read_arrays([1, 2, 3], (4, 5, 6), np.asarray([0.1, 0.2, 0.1]))
    assert isinstance(dset, Data1D)
    assert dset.name == ''
    assert dset.x == pytest.approx(np.asarray([1, 2, 3]))
    assert dset.y == pytest.approx(np.asarray([4, 5, 6]))
    assert dset.staterror == pytest.approx(np.asarray([0.1, 0.2, 0.1]))
    assert dset.syserror is None

    # Check it creates NumPy arrays
    assert isinstance(dset.x, np.ndarray)
    assert isinstance(dset.y, np.ndarray)
    assert isinstance(dset.staterror, np.ndarray)


@requires_fits
def test_read_arrays_data1d_combined():

    arg = np.asarray([[1, 4, 0.2, 0.1],
                      [2, 5, 0.3, 0.05]])
    dset = io.read_arrays(arg)
    assert isinstance(dset, Data1D)
    assert dset.name == ''
    assert dset.x == pytest.approx(np.asarray([1, 2]))
    assert dset.y == pytest.approx(np.asarray([4, 5]))
    assert dset.staterror == pytest.approx(np.asarray([0.2, 0.3]))
    assert dset.syserror == pytest.approx(np.asarray([0.1, 0.05]))

    # Check it creates NumPy arrays
    assert isinstance(dset.x, np.ndarray)
    assert isinstance(dset.y, np.ndarray)
    assert isinstance(dset.staterror, np.ndarray)
    assert isinstance(dset.syserror, np.ndarray)


@requires_fits
@pytest.mark.parametrize("arg", [[], None, (None, None)])
def test_write_arrays_no_data(arg, tmp_path):

    tmpfile = tmp_path / 'test.dat'
    emsg = r"^please supply array\(s\) to write to file$"
    with pytest.raises(IOErr, match=emsg):
        io.write_arrays(str(tmpfile), arg, clobber=True)


@requires_fits
def test_write_arrays_wrong_lengths(tmp_path):

    tmpfile = tmp_path / 'test.dat'
    emsg = "^not all arrays are of equal length$"
    with pytest.raises(IOErr, match=emsg):
        io.write_arrays(str(tmpfile), ([1, 2], [1, 2, 3]), clobber=True)


@requires_fits
def test_write_arrays_wrong_field_length(tmp_path):

    tmpfile = tmp_path / 'test.dat'
    emsg = "^Expected 2 columns but found 3$"
    with pytest.raises(IOErr, match=emsg):
        io.write_arrays(str(tmpfile), ([1, 2], [1, 2]), fields=['a', 'b', 'c'], clobber=True)


@requires_fits
def test_write_arrays_clobber(tmp_path):

    outfile = tmp_path / 'test.dat'
    outfile.write_text('x')

    emsg = "^file '.*test.dat' exists and clobber is not set$"
    with pytest.raises(IOErr, match=emsg):
        io.write_arrays(str(outfile), None, clobber=False)


@requires_data
@requires_fits
def test_read_table_pha(make_data_path):
    """Do we pick up CHANNEL and COUNTS by default?"""

    # This file has columns:
    #   CHANNEL
    #   PI
    #   COUNTS
    #   STAT_ERR
    #   COUNT_RATE
    # Check we pick up CHANNEL and COUNTS
    #
    infile = make_data_path('source.pi')
    tbl = io.read_table(infile)
    assert isinstance(tbl, Data1D)
    assert len(tbl.x) == 1024
    assert len(tbl.y) == 1024
    assert tbl.staterror is None
    assert tbl.syserror is None

    # These are meant to be integer values so we can test
    # for equality. Except that it depends on the backend:
    # pyfits returns Int32 whereas crates returns Float64.
    #
    assert tbl.x == pytest.approx(np.arange(1, 1025))
    assert tbl.y.min() == pytest.approx(0)
    assert tbl.y[0] == pytest.approx(0)
    assert tbl.y[-1] == pytest.approx(128)
    assert tbl.y[:-1].max() == pytest.approx(78)
    assert np.argmax(tbl.y[:-1]) == 74


@requires_data
@requires_fits
def test_read_ascii_3_col(make_data_path):
    """Found when working on #1921 so explicitly test this case."""

    infile = make_data_path("data1.dat")
    tbl = io.read_ascii(infile, 3)
    assert isinstance(tbl, Data1D)
    assert tbl.name.find("/data1.dat")
    assert tbl.x == pytest.approx(np.arange(0.5, 11.5))
    assert tbl.y == pytest.approx([1.6454, 1.7236, 1.9472,
                                   2.2348, 2.6187, 2.8642,
                                   3.1263, 3.2073, 3.2852,
                                   3.3092, 3.4496])
    assert tbl.staterror == pytest.approx(0.04114 * np.ones(11))
    assert tbl.syserror is None


@requires_fits
@requires_data
def test_read_table_object(make_data_path):
    """Check we can send in a 'object'.

    This support is not well advertised so it could perhaps be
    removed, but let's add a basic test at least.

    """

    # Could create the "object" manually, but let's just read one in
    #
    infile = make_data_path("1838_rprofile_rmid.fits")
    close = False

    if io.backend.name == "crates":
        import pycrates  # type: ignore
        arg = pycrates.read_file(infile)

    elif io.backend.name == "pyfits":
        from astropy.io import fits  # type: ignore
        arg = fits.open(infile)
        close = True

    else:
        assert False, f"unknown backend: {io.backend}"

    try:
        # this implicitly checks case insensitivity on column names
        tbl = io.read_table(arg, colkeys=["area", "counts"])

    finally:
        if close:
            arg.close()

    assert isinstance(tbl, Data1D)
    assert tbl.x.sum() == pytest.approx(125349.55)
    assert tbl.y.sum() == pytest.approx(21962)
