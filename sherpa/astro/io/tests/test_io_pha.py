#
#  Copyright (C) 2020 - 2024
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

"""Check that we can write out and read in PHA files.

Check that the PHA files match OGIP standards where appropriate and
that we can do things like read in the file we've just created.

The Sherpa I/O routines are deliberately not general purpose, so in
order to check the output files we use the backend-specific
functionality.  There may be a way to abstract this somewhat, but wait
until we have enough tests that make it worthwhile.

"""

import logging

import numpy as np

import pytest

from sherpa.astro.data import DataPHA
from sherpa.astro import io
from sherpa.utils.err import DataErr, IOErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, requires_fits


def backend_is(name: str) -> bool:
    """Are we using the specified backend?"""
    return io.backend.name == name


@requires_fits
@requires_data
def test_pha_write_basic(make_data_path, tmp_path):
    """Check we can write out a PHA file as a FITS file.

    This uses an existing PHA file rather than creating one manually.

    """

    def check_header(obj):
        # check header for a selected set of keywords
        hdr = obj.header

        assert "HDUNAME" not in hdr

        # selected OGIP keywords
        assert hdr["HDUCLASS"] == "OGIP"
        assert hdr["HDUCLAS1"] == "SPECTRUM"
        assert hdr["HDUCLAS2"] == "TOTAL"
        assert hdr["HDUCLAS3"] == "TYPE:I"
        assert hdr["HDUCLAS4"] == "COUNT"
        assert hdr["HDUVERS"] == "1.1.0"
        assert hdr["HDUVERS1"] == "1.1.0"

        # a few header keywords to check they are handled,
        # including data types (string, float, integer,
        # logical).
        #
        assert hdr["OBJECT"] == "3C 273"
        assert hdr["INSTRUME"] == "ACIS"
        assert hdr["GRATING"] == "HETG"
        assert hdr["DETNAM"] == "ACIS-56789"
        assert hdr["RA_NOM"] == pytest.approx(187.28186566)
        assert hdr["DEC_NOM"] == pytest.approx(2.05610034)
        assert hdr["SIM_X"] == pytest.approx(-0.68282252)
        assert hdr["DATE-OBS"] == "2000-01-10T06:47:15"
        assert isinstance(hdr["CLOCKAPP"], (bool, np.bool_))
        assert hdr["CLOCKAPP"]
        assert hdr["TIMEZERO"] == 0
        assert hdr["DETCHANS"] == 1024

        assert hdr["SYS_ERR"] == pytest.approx(0.0)
        assert hdr["GROUPING"] == 0
        assert hdr["QUALITY"] == 0

    def check_data(obj, roundtrip=False):
        """Basic checks of the data"""

        assert obj.staterror is None
        assert obj.syserror is None
        assert obj.bin_lo is None
        assert obj.bin_hi is None

        assert obj.exposure == pytest.approx(38564.608926889)
        assert np.log10(obj.backscal) == pytest.approx(-5.597491618115439)
        assert obj.areascal == pytest.approx(1.0)
        for f in ["grouped", "subtracted", "rate"]:
            assert isinstance(getattr(obj, f), bool)

        assert obj.grouped
        assert not obj.subtracted
        assert obj.rate

        assert obj.plot_fac == 0

        assert obj.channel.dtype == np.dtype("float64")
        assert obj.counts.dtype == np.dtype("float64")

        assert obj.channel == pytest.approx(np.arange(1, 1025))
        assert len(obj.counts) == 1024
        assert obj.counts[0:11] == pytest.approx(np.zeros(11))
        cvals = [1, 3, 2, 3, 7, 1, 6, 4, 4, 0]
        assert obj.counts[12:22] == pytest.approx(cvals)

        assert len(obj.grouping) == 1024
        if roundtrip:
            assert obj.quality is None
        else:
            assert len(obj.quality) == 1024
            assert obj.quality == pytest.approx(np.zeros(1024))

        assert obj.grouping.dtype == np.dtype("int16")
        if not roundtrip:
            assert obj.quality.dtype == np.dtype("int16")

        one, = np.where(obj.grouping == 1)
        expected = [0,  17,  21,  32,  39,  44,  48,  51,  54,  56,  59,  61,  65,
                    68,  71,  75,  78,  82,  88,  96, 101, 110, 116, 124, 130, 133,
                    139, 143, 150, 156, 164, 177, 186, 196, 211, 232, 244, 260, 276,
                    291, 323, 344, 368, 404, 450, 676]
        assert one == pytest.approx(expected)
        assert obj.grouping.sum() == -932
        assert set(obj.grouping) == {-1, 1}

    infile = make_data_path("3c273.pi")
    indata = io.read_pha(infile)
    assert isinstance(indata, DataPHA)
    assert indata.name.endswith("/3c273.pi")
    check_header(indata)
    check_data(indata)

    # The responses and background should be read in
    #
    assert indata.response_ids == pytest.approx([1])
    assert indata.background_ids == pytest.approx([1])

    assert indata.get_arf().name.endswith("/3c273.arf")
    assert indata.get_rmf().name.endswith("/3c273.rmf")
    assert indata.get_background().name.endswith("/3c273_bg.pi")

    # check we can write it out - be explicit with all options
    #
    outfile = tmp_path / "test.pha"
    outfile = str(outfile)  # our IO routines don't recognize paths
    io.write_pha(outfile, indata, ascii=False, clobber=False)

    outdata = io.read_pha(outfile)
    assert isinstance(outdata, DataPHA)
    assert outdata.name.endswith("/test.pha")
    check_header(outdata)
    check_data(outdata, roundtrip=True)

    # The responses and background should NOT be read in
    # (we are in a different directory to infile so we can't
    # find these files).
    #
    assert outdata.response_ids == []
    assert outdata.background_ids == []

    assert outdata.get_arf() is None
    assert outdata.get_rmf() is None
    assert outdata.get_background() is None


@requires_fits
@requires_data
def test_pha_write_basic_errors(make_data_path, tmp_path):
    """test_pha_write_basic but with use_errors=True.
    """

    def check_header(obj):
        # check header for a selected set of keywords
        hdr = obj.header

        assert "HDUNAME" not in hdr

        # selected OGIP keywords
        assert hdr["HDUCLASS"] == "OGIP"
        assert hdr["HDUCLAS1"] == "SPECTRUM"
        assert hdr["HDUCLAS2"] == "TOTAL"
        assert hdr["HDUCLAS3"] == "TYPE:I"
        assert hdr["HDUCLAS4"] == "COUNT"
        assert hdr["HDUVERS"] == "1.1.0"
        assert hdr["HDUVERS1"] == "1.1.0"

        # a few header keywords to check they are handled,
        # including data types (string, float, integer,
        # logical).
        #
        assert hdr["OBJECT"] == "3C 273"
        assert hdr["INSTRUME"] == "ACIS"
        assert hdr["GRATING"] == "HETG"
        assert hdr["DETNAM"] == "ACIS-56789"
        assert hdr["RA_NOM"] == pytest.approx(187.28186566)
        assert hdr["DEC_NOM"] == pytest.approx(2.05610034)
        assert hdr["SIM_X"] == pytest.approx(-0.68282252)
        assert hdr["DATE-OBS"] == "2000-01-10T06:47:15"
        assert isinstance(hdr["CLOCKAPP"], (bool, np.bool_))
        assert hdr["CLOCKAPP"]
        assert hdr["TIMEZERO"] == 0
        assert hdr["DETCHANS"] == 1024

        assert hdr["SYS_ERR"] == pytest.approx(0.0)
        assert hdr["GROUPING"] == 0
        assert hdr["QUALITY"] == 0

    def check_data(obj, roundtrip=False):
        """Basic checks of the data"""

        assert (obj.staterror == 0).sum() == 714
        assert obj.staterror.sum() == pytest.approx(449.71593831)
        assert obj.staterror[0] == pytest.approx(0.0)
        assert obj.staterror[13] == pytest.approx(1.73205081)
        assert obj.staterror[51] == pytest.approx(3.1622776602)
        assert obj.staterror[1023] == pytest.approx(2.82842712)
        assert np.argmax(obj.staterror) == 51

        assert obj.syserror is None
        assert obj.bin_lo is None
        assert obj.bin_hi is None

        assert obj.exposure == pytest.approx(38564.608926889)
        assert np.log10(obj.backscal) == pytest.approx(-5.597491618115439)
        assert obj.areascal == pytest.approx(1.0)
        for f in ["grouped", "subtracted", "rate"]:
            assert isinstance(getattr(obj, f), bool)

        assert obj.grouped
        assert not obj.subtracted
        assert obj.rate

        assert obj.plot_fac == 0

        assert obj.channel.dtype == np.dtype("float64")
        assert obj.counts.dtype == np.dtype("float64")

        assert obj.channel == pytest.approx(np.arange(1, 1025))
        assert len(obj.counts) == 1024
        assert obj.counts[0:11] == pytest.approx(np.zeros(11))
        cvals = [1, 3, 2, 3, 7, 1, 6, 4, 4, 0]
        assert obj.counts[12:22] == pytest.approx(cvals)

        assert len(obj.grouping) == 1024
        if roundtrip:
            assert obj.quality is None
        else:
            assert len(obj.quality) == 1024
            assert obj.quality == pytest.approx(np.zeros(1024))

        assert obj.grouping.dtype == np.dtype("int16")
        if not roundtrip:
            assert obj.quality.dtype == np.dtype("int16")

        one, = np.where(obj.grouping == 1)
        expected = [0,  17,  21,  32,  39,  44,  48,  51,  54,  56,  59,  61,  65,
                    68,  71,  75,  78,  82,  88,  96, 101, 110, 116, 124, 130, 133,
                    139, 143, 150, 156, 164, 177, 186, 196, 211, 232, 244, 260, 276,
                    291, 323, 344, 368, 404, 450, 676]
        assert one == pytest.approx(expected)
        assert obj.grouping.sum() == -932
        assert set(obj.grouping) == {-1, 1}

    infile = make_data_path("3c273.pi")
    indata = io.read_pha(infile, use_errors=True)
    assert isinstance(indata, DataPHA)
    assert indata.name.endswith("/3c273.pi")
    check_header(indata)
    check_data(indata)

    # The responses and background should be read in
    #
    assert indata.response_ids == pytest.approx([1])
    assert indata.background_ids == pytest.approx([1])

    assert indata.get_arf().name.endswith("/3c273.arf")
    assert indata.get_rmf().name.endswith("/3c273.rmf")
    assert indata.get_background().name.endswith("/3c273_bg.pi")

    # check we can write it out - be explicit with all options
    #
    outfile = tmp_path / "test.pha"
    outfile = str(outfile)  # our IO routines don't recognize paths
    io.write_pha(outfile, indata, ascii=False, clobber=False)

    outdata = io.read_pha(outfile, use_errors=True)
    assert isinstance(outdata, DataPHA)
    assert outdata.name.endswith("/test.pha")
    check_header(outdata)
    check_data(outdata, roundtrip=True)

    # The responses and background should NOT be read in
    # (we are in a different directory to infile so we can't
    # find these files).
    #
    assert outdata.response_ids == []
    assert outdata.background_ids == []

    assert outdata.get_arf() is None
    assert outdata.get_rmf() is None
    assert outdata.get_background() is None


@requires_fits
@requires_data
def test_pha_write_xmm_grating(make_data_path, tmp_path):
    """Check we can handle an XMM grating.

    This has an AREASCAL column and WCS attached to the CHANNEL
    column. We don't guarantee the WCS info will be retained or
    propagated.

    """

    def check_header(obj, roundtrip=False):
        # check header for a selected set of keywords
        hdr = obj.header

        # selected OGIP keywords
        assert "HDUNAME" not in hdr

        assert hdr["HDUCLASS"] == "OGIP"
        assert hdr["HDUCLAS1"] == "SPECTRUM"
        assert hdr["HDUCLAS2"] == "NET"
        assert hdr["HDUCLAS3"] == "COUNT"
        if roundtrip:
            assert hdr["HDUCLAS4"] == "TYPE:I"
        else:
            assert "HDUCLAS4" not in hdr
        assert hdr["HDUVERS"] == "1.2.0"
        assert "HDUVERS1" not in hdr

        # a few header keywords to check they are handled,
        # including data types (string, float, integer,
        # logical).
        #
        assert hdr["OBJECT"] == "TW Hya"
        assert hdr["TELESCOP"] == "XMM"
        assert hdr["INSTRUME"] == "RGS1"
        assert hdr["FILTER"] == "UNKNOWN"

        assert not hdr["POISSERR"]
        assert hdr["CHANTYPE"] == "PI"
        assert hdr["DETCHANS"] == 3600

        assert "EXPOSIRE" not in hdr

        assert hdr["SYS_ERR"] == pytest.approx(0.0)
        assert hdr["GROUPING"] == 0
        assert "QUALITY" not in hdr
        assert hdr["CORRSCAL"] == 1
        assert "AREASCAL" not in hdr
        assert "BACKSCAL" not in hdr

        assert hdr["CORRFILE"] == "none"
        for key in ["RESPFILE", "ANCRFILE", "BACKFILE"]:
            assert key not in hdr

        # check WCS that was attached to the CHANNEL column
        assert "TCRYP1" not in hdr
        assert "TCRVL1" not in hdr
        assert "TCUNI1" not in hdr
        assert "TCRPX1" not in hdr
        assert "TCDLT1" not in hdr
        assert "TLMIN1" not in hdr
        assert "TLMAX1" not in hdr

    def check_data(obj, roundtrip=False):
        """Basic checks of the data"""

        assert obj.staterror.min() == pytest.approx(1.0)
        assert obj.staterror.sum() == pytest.approx(4496.9490242)
        assert obj.staterror[1498] == pytest.approx(7.937253952)
        assert np.argmax(obj.staterror) == 1498

        assert obj.syserror is None
        assert obj.bin_lo is None
        assert obj.bin_hi is None

        assert obj.exposure == pytest.approx(28965.6406250)
        assert obj.backscal == pytest.approx(1.0)
        assert len(obj.areascal) == 3600
        assert obj.areascal.min() == pytest.approx(0.0)
        assert obj.areascal.max() == pytest.approx(1.0)
        assert obj.areascal.sum() == pytest.approx(2950.6953666)

        for f in ["grouped", "subtracted", "rate"]:
            assert isinstance(getattr(obj, f), bool)

        assert not obj.grouped
        assert not obj.subtracted
        assert obj.rate

        assert obj.plot_fac == 0

        assert obj.channel.dtype == np.dtype("float64")
        assert obj.counts.dtype == np.dtype("float64")

        assert obj.channel == pytest.approx(np.arange(1, 3601))
        assert len(obj.counts) == 3600

        expected = [0, 0, -1, -3, 0, -1, -1, 0, -2, -1]
        assert obj.counts[112:122] == pytest.approx(expected)

        assert obj.grouping is None
        assert len(obj.quality) == 3600
        assert obj.quality.min() == 0
        assert obj.quality.max() == 1
        assert obj.quality.sum() == 741

    infile = make_data_path("xmmrgs/P0112880201R1S004SRSPEC1003.FTZ")
    indata = io.read_pha(infile, use_errors=True)
    assert isinstance(indata, DataPHA)
    assert indata.name.endswith("/P0112880201R1S004SRSPEC1003.FTZ")
    check_header(indata)
    check_data(indata)

    assert indata.response_ids == []
    assert indata.background_ids == []

    assert indata.get_arf() is None
    assert indata.get_rmf() is None
    assert indata.get_background() is None

    # check we can write it out - be explicit with all options
    #
    outfile = tmp_path / "test.pha"
    outfile = str(outfile)  # our IO routines don't recognize paths

    io.write_pha(outfile, indata, ascii=False, clobber=False)

    outdata = io.read_pha(outfile, use_errors=True)
    assert isinstance(outdata, DataPHA)
    assert outdata.name.endswith("/test.pha")
    check_header(outdata, roundtrip=True)
    check_data(outdata, roundtrip=True)

    # The responses and background should NOT be read in
    # (we are in a different directory to infile so we can't
    # find these files).
    #
    assert outdata.response_ids == []
    assert outdata.background_ids == []

    assert outdata.get_arf() is None
    assert outdata.get_rmf() is None
    assert outdata.get_background() is None


def check_write_pha_fits_basic_roundtrip_crates(path):
    import pycrates  # type: ignore
    ds = pycrates.CrateDataset(str(path), mode="r")

    assert ds.get_ncrates() == 2
    cr = ds.get_crate(2)

    assert cr.name == "SPECTRUM"
    assert cr.get_colnames() == ["CHANNEL", "COUNTS"]

    # unfortunately crates auto-converts int32 to int64
    # (this is actually the underlying cxcdm module).
    #
    c0 = cr.get_column(0)
    assert c0.name == "CHANNEL"
    assert c0.values.dtype == np.int64
    assert c0.get_tlmin() == 1
    assert c0.get_tlmax() == 4

    c1 = cr.get_column(1)
    assert c1.name == "COUNTS"
    assert c1.values.dtype == np.int64
    # crates circa CIAO 4.14 doesn't return the correct values here;
    # it's returning int32 values and has got the negative value
    # wrong by 1
    assert c1.get_tlmin() <= -2147483647
    assert c1.get_tlmax() == 2147483647

    assert cr.get_key_value("HDUCLASS") == "OGIP"
    assert cr.get_key_value("HDUCLAS1") == "SPECTRUM"
    assert cr.get_key_value("HDUCLAS2") == "TOTAL"
    assert cr.get_key_value("HDUCLAS3") == "TYPE:I"
    assert cr.get_key_value("HDUCLAS4") == "COUNT"
    assert cr.get_key_value("HDUVERS") == "1.2.1"
    assert cr.get_key_value("TELESCOP") == "none"
    assert cr.get_key_value("INSTRUME") == "none"
    assert cr.get_key_value("FILTER") == "none"
    assert cr.get_key_value("POISSERR")

    assert cr.get_key_value("CHANTYPE") == "PI"
    assert cr.get_key_value("DETCHANS") == 4

    assert cr.get_key_value("SYS_ERR") == 0
    assert cr.get_key_value("QUALITY") == 0
    assert cr.get_key_value("GROUPING") == 0

    assert cr.get_key_value("AREASCAL") == pytest.approx(1.0)
    assert cr.get_key_value("BACKSCAL") == pytest.approx(1.0)
    assert cr.get_key_value("CORRSCAL") == 0

    for key in ["BACKFILE", "CORRFILE", "RESPFILE", "ANCRFILE"]:
        assert cr.get_key_value(key) == "none"

    assert cr.get_key_value("EXPOSURE") == pytest.approx(1.0)


def check_write_pha_fits_basic_roundtrip_pyfits(path):
    from astropy.io import fits  # type: ignore
    hdus = fits.open(str(path))
    try:
        assert len(hdus) == 2
        hdu = hdus[1]
        assert hdu.name == "SPECTRUM"
        assert hdu.level == 1
        assert hdu.ver == 1
        assert len(hdu.columns) == 2

        assert hdu.columns[0].name == "CHANNEL"
        assert hdu.columns[0].format == "J"
        assert hdu.columns[1].name == "COUNTS"
        assert hdu.columns[1].format == "J"

        assert hdu.header["HDUCLASS"] == "OGIP"
        assert hdu.header["HDUCLAS1"] == "SPECTRUM"
        assert hdu.header["HDUCLAS2"] == "TOTAL"
        assert hdu.header["HDUCLAS3"] == "TYPE:I"
        assert hdu.header["HDUCLAS4"] == "COUNT"
        assert hdu.header["HDUVERS"] == "1.2.1"
        assert hdu.header["TELESCOP"] == "none"
        assert hdu.header["INSTRUME"] == "none"
        assert hdu.header["FILTER"] == "none"

        assert hdu.header["POISSERR"]

        assert hdu.header["CHANTYPE"] == "PI"
        assert hdu.header["DETCHANS"] == 4

        assert hdu.header['SYS_ERR'] == 0
        assert hdu.header['QUALITY'] == 0
        assert hdu.header['GROUPING'] == 0

        assert hdu.header["AREASCAL"] == pytest.approx(1.0)
        assert hdu.header["BACKSCAL"] == pytest.approx(1.0)
        assert hdu.header["CORRSCAL"] == 0

        for key in ["BACKFILE", "CORRFILE", "RESPFILE", "ANCRFILE"]:
            assert hdu.header[key] == "none"

        assert hdu.header["TLMIN1"] == 1
        assert hdu.header["TLMAX1"] == 4

        # check some keywords we don't expect
        #
        assert "TLMIN2" not in hdu.header
        assert "TLMAX2" not in hdu.header

        assert hdu.header["EXPOSURE"] == pytest.approx(1.0)

    finally:
        hdus.close()


@requires_fits
def test_write_pha_fits_basic_roundtrip(tmp_path):
    """A very-basic PHA output

    No ancillary information and no extra header information.

    """

    chans = np.arange(1, 5, dtype=np.int16)
    counts = np.asarray([1, 0, 3, 2], dtype=np.int16)
    pha = DataPHA("testy", chans, counts)

    outfile = tmp_path / "out.pi"
    io.write_pha(str(outfile), pha, ascii=False, clobber=False)
    pha = None

    inpha = io.read_pha(str(outfile))
    assert isinstance(inpha, DataPHA)
    assert inpha.channel == pytest.approx(chans)
    assert inpha.counts == pytest.approx(counts)
    for field in ["staterror", "syserror", "bin_lo", "bin_hi",
                  "grouping", "quality"]:
        assert getattr(inpha, field) is None

    assert inpha.exposure == pytest.approx(1.0)
    assert inpha.backscal == pytest.approx(1.0)
    assert inpha.areascal == pytest.approx(1.0)

    assert not inpha.grouped
    assert not inpha.subtracted
    assert inpha.units == "channel"
    assert inpha.rate
    assert inpha.plot_fac == 0
    assert inpha.response_ids == []
    assert inpha.background_ids == []

    if backend_is("crates"):
        check_write_pha_fits_basic_roundtrip_crates(outfile)

    elif backend_is("pyfits"):
        check_write_pha_fits_basic_roundtrip_pyfits(outfile)

    else:
        # Technically this could be dummy_backend but it would have
        # failed earlier.
        #
        raise RuntimeError(f"Unknown io backend: {io.backend}")


def check_write_pha_fits_with_extras_roundtrip_crates(path, etime, bscal):
    import pycrates  # type: ignore
    ds = pycrates.CrateDataset(str(path), mode="r")

    assert ds.get_ncrates() == 2
    cr = ds.get_crate(2)

    assert cr.name == "SPECTRUM"
    assert cr.get_colnames() == ["CHANNEL", "COUNTS", "GROUPING", "QUALITY", "AREASCAL"]

    c0 = cr.get_column(0)
    assert c0.name == "CHANNEL"
    assert c0.values.dtype == np.int64

    # The CHANNEL data range is 1-4 but the header contains
    # DETCHANS=10, so the data does not really make sense. What should
    # be written out?
    #
    assert c0.get_tlmin() == 1
    assert c0.get_tlmax() == 4

    c1 = cr.get_column(1)
    assert c1.name == "COUNTS"
    assert c1.values.dtype == np.float32
    assert c1.get_tlmin() == pytest.approx(-3.4028235e+38)
    assert c1.get_tlmax() == pytest.approx(3.4028235e+38)

    c2 = cr.get_column(2)
    assert c2.name == "GROUPING"
    assert c2.values.dtype == np.int16
    assert c2.get_tlmin() == -32768
    assert c2.get_tlmax() == 32767

    c3 = cr.get_column(3)
    assert c3.name == "QUALITY"
    assert c3.values.dtype == np.int16
    assert c3.get_tlmin() == -32768
    assert c3.get_tlmax() == 32767

    c4 = cr.get_column(4)
    assert c4.name == "AREASCAL"
    # FITS standard has as this 4-byte real
    assert c4.values.dtype == np.float32
    assert c4.get_tlmin() < -3e38
    assert c4.get_tlmax() > 3e38

    assert cr.get_key_value("HDUCLASS") == "OGIP"
    assert cr.get_key_value("HDUCLAS1") == "SPECTRUM"
    assert cr.get_key_value("HDUCLAS2") == "TOTAL"
    assert cr.get_key_value("HDUCLAS3") == "COUNT"
    assert cr.get_key_value("HDUCLAS4") == "TYPE:I"
    assert cr.get_key_value("HDUVERS") == "1.2.1"
    assert cr.get_key_value("HDUVERS1") is None

    assert cr.get_key_value("POISSERR")

    assert cr.get_key_value("TELESCOP") == "CHANDRA"
    assert cr.get_key_value("INSTRUME") == "ACIS"
    assert cr.get_key_value("FILTER") == "NONE"

    assert cr.get_key_value("EXPOSURE") == pytest.approx(etime)
    assert cr.get_key_value("BACKSCAL") == pytest.approx(bscal)
    assert cr.get_key_value("CORRFILE") == "None"
    assert cr.get_key_value("ANCRFILE") == "made-up-ancrfile.fits"

    assert cr.get_key_value("CHANTYPE") == "PI"
    assert cr.get_key_value("DETCHANS") == 10

    assert cr.get_key_value("CORRSCAL") == 0

    for key in ["BACKFILE", "RESPFILE"]:
        assert cr.get_key_value(key) == "none"

    assert cr.get_key_value("SYS_ERR") == 0

    # We do not have these keywords as they are stored as columns
    for key in ["AREASCAL", "QUALITY", "GROUPING"]:
        assert cr.get_key_value(key) is None


def check_write_pha_fits_with_extras_roundtrip_pyfits(path, etime, bscal):
    from astropy.io import fits  # type: ignore
    hdus = fits.open(str(path))
    try:
        assert len(hdus) == 2
        hdu = hdus[1]
        assert hdu.name == "SPECTRUM"
        assert hdu.level == 1
        assert hdu.ver == 1
        assert len(hdu.columns) == 5

        assert hdu.columns[0].name == "CHANNEL"
        assert hdu.columns[0].format == "J"
        assert hdu.columns[1].name == "COUNTS"
        assert hdu.columns[1].format == "E"
        assert hdu.columns[2].name == "GROUPING"
        assert hdu.columns[2].format == "I"
        assert hdu.columns[3].name == "QUALITY"
        assert hdu.columns[3].format == "I"
        assert hdu.columns[4].name == "AREASCAL"
        # FITS standard has as this 4-byte real
        assert hdu.columns[4].format == "E"

        assert hdu.header["HDUCLASS"] == "OGIP"
        assert hdu.header["HDUCLAS1"] == "SPECTRUM"
        assert hdu.header["HDUCLAS2"] == "TOTAL"
        assert hdu.header["HDUCLAS3"] == "COUNT"  # WRONG I think
        assert hdu.header["HDUCLAS4"] == "TYPE:I"  # WRONG I think
        assert hdu.header["HDUVERS"] == "1.2.1"  # should this be 1.1.0?
        assert "HDUVERS1" not in hdu.header

        assert hdu.header["POISSERR"]

        assert hdu.header["TELESCOP"] == "CHANDRA"
        assert hdu.header["INSTRUME"] == "ACIS"
        assert hdu.header["FILTER"] == "NONE"

        assert hdu.header["EXPOSURE"] == pytest.approx(etime)
        assert hdu.header["BACKSCAL"] == pytest.approx(bscal)
        assert hdu.header["CORRFILE"] == "None"
        assert hdu.header["ANCRFILE"] == "made-up-ancrfile.fits"

        assert hdu.header["CHANTYPE"] == "PI"
        assert hdu.header["DETCHANS"] == 10

        assert hdu.header["CORRSCAL"] == 0

        for key in ["BACKFILE", "RESPFILE"]:
            assert hdu.header[key] == "none"

        assert hdu.header["SYS_ERR"] == 0

        # We do not have these keywords as they are stored as columns
        for key in ["AREASCAL", "QUALITY", "GROUPING"]:
            assert key not in hdu.header

        # The CHANNEL data range is 1-4 but the header contains
        # DETCHANS=10, so the data does not really make sense. What
        # should be written out?
        #
        assert hdu.header["TLMIN1"] == 1
        assert hdu.header["TLMAX1"] == 4

    finally:
        hdus.close()


@requires_fits
def test_write_pha_fits_with_extras_roundtrip(tmp_path, caplog):
    """PHA-I with grouping/quality/errors/header

    This covers issue #488 - that is, should the output use the correct
    data type for columns? At present the code does not.

    """

    chans = np.arange(1, 5, dtype=np.int32)
    counts = np.asarray([1, 0, 3, 2], dtype=np.int32)
    grouping = np.asarray([1, -1, 1, 1], dtype=np.int16)
    quality = np.asarray([0, 0, 0, 2], dtype=np.int16)
    etime = 1023.4
    bscal = 0.05
    ascal = np.asarray([1, 1, 0.9, 0.9])

    hdr = {"TELESCOP": "CHANDRA", "INSTRUME": "ACIS", "FILTER": "NONE",
           "CHANTYPE": "PI",
           # DETCHANS intentionally does not match the data, so what gets written out?
           "DETCHANS": 10,
           "OBJECT": "Made up source",
           "CORRFILE": "None",
           # This will cause a warning when reading in the file
           "ANCRFILE": "made-up-ancrfile.fits",
           # "structural keywords" which match DETCHANS
           "TLMIN1": 1, "TLMAX1": 10}

    pha = DataPHA("testy",
                  chans.astype(np.float64),
                  counts.astype(np.float32),
                  grouping=grouping.astype(np.float32),
                  quality=quality.astype(np.float64),
                  exposure=etime,
                  backscal=bscal,
                  areascal=ascal,
                  header=hdr)

    outfile = tmp_path / "out.pi"
    io.write_pha(str(outfile), pha, ascii=False, clobber=False)
    pha = None

    assert len(caplog.record_tuples) == 0

    with SherpaVerbosity("INFO"):
        inpha = io.read_pha(str(outfile))

    assert len(caplog.record_tuples) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.astro.io"
    assert lvl == logging.WARNING

    # message depends on the backend.
    start = "unable to read ARF: "
    if backend_is("crates"):
        assert msg.startswith(f"{start}File")
        assert msg.endswith("/made-up-ancrfile.fits does not exist.")
    elif backend_is("pyfits"):
        assert msg.startswith(f"{start}file '")
        assert msg.endswith("/made-up-ancrfile.fits' not found")
    else:
        raise RuntimeError(f"Unknown io backend: {io.backend}")

    assert isinstance(inpha, DataPHA)
    assert inpha.channel == pytest.approx(chans)
    assert inpha.counts == pytest.approx(counts)
    assert inpha.grouping == pytest.approx(grouping)
    assert inpha.quality == pytest.approx(quality)
    assert inpha.exposure == pytest.approx(etime)
    assert inpha.backscal == pytest.approx(bscal)
    assert inpha.areascal == pytest.approx(ascal)
    for field in ["staterror", "syserror", "bin_lo", "bin_hi"]:
        assert getattr(inpha, field) is None

    assert inpha.grouped
    assert not inpha.subtracted
    assert inpha.units == "channel"
    assert inpha.rate
    assert inpha.plot_fac == 0
    assert inpha.response_ids == []
    assert inpha.background_ids == []

    if backend_is("crates"):
        check_write_pha_fits_with_extras_roundtrip_crates(outfile, etime, bscal)

    elif backend_is("pyfits"):
        check_write_pha_fits_with_extras_roundtrip_pyfits(outfile, etime, bscal)

    else:
        raise RuntimeError(f"Unknown io backend: {io.backend}")


@requires_fits
def test_pha_missing_backfile(tmp_path, caplog):
    """What happens if BACKFILE does not exist.

    This goes through a slightly-different path to the missing
    ANCRFILE case tested in test_write_pha_fits_with_extras_roundtrip.

    """

    chans = np.arange(1, 5, dtype=np.int32)
    counts = np.asarray([1, 0, 3, 2], dtype=np.int32)
    etime = 1023.4
    bscal = 0.05

    hdr = {"TELESCOP": "CHANDRA", "INSTRUME": "ACIS", "FILTER": "NONE",
           "CHANTYPE": "PI",
           "DETCHANS": 5,
           "OBJECT": "Made up source",
           "CORRFILE": "None",
           "BACKFILE": "made-up-backfile.fits",
           # "structural keywords" which match DETCHANS
           "TLMIN1": 1, "TLMAX1": 5}

    pha = DataPHA("testy",
                  chans.astype(np.float64),
                  counts.astype(np.float32),
                  exposure=etime,
                  backscal=bscal,
                  header=hdr)

    outfile = tmp_path / "out.pi"
    io.write_pha(str(outfile), pha, ascii=False, clobber=False)
    pha = None

    assert len(caplog.record_tuples) == 0

    with SherpaVerbosity("INFO"):
        inpha = io.read_pha(str(outfile))

    assert len(caplog.record_tuples) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.astro.io"
    assert lvl == logging.WARNING

    # message depends on the backend.
    start = "unable to read background: "
    if backend_is("crates"):
        assert msg.startswith(f"{start}File")
        assert msg.endswith("/made-up-backfile.fits does not exist.")
    elif backend_is("pyfits"):
        assert msg.startswith(f"{start}file '")
        assert msg.endswith("/made-up-backfile.fits' not found")
    else:
        raise RuntimeError(f"Unknown io backend: {io.backend}")

    assert inpha.channel == pytest.approx(chans)
    assert inpha.counts == pytest.approx(counts)

    assert inpha.response_ids == []
    assert inpha.background_ids == []

    # Checks related to issue #1885 (so this is partly a regression
    # test).
    #
    for key in ["ANCRFILE", "BACKFILE", "RESPFILE"]:
        assert key not in inpha.header

    assert inpha.header["CORRFILE"] == "None"

    # DJB is interested to see when this gets set, so treat this
    # as a regression test.
    #
    assert inpha.header["SYS_ERR"] == 0


@requires_fits
@requires_data
def test_chandra_phaII_roundtrip(make_data_path, tmp_path):
    """Can we read in/write out/read in a PHA-II dataset.

    There is no abiility to write out a PHA-II file,
    only PHA-I.
    """

    def check_header(pha, roundtrip=False):
        hdr = pha.header
        assert hdr["TG_M"] == 2
        assert hdr["TG_PART"] == 1

        assert "HDUNAME" not in hdr

        assert hdr["HDUCLASS"] == "OGIP"
        assert hdr["HDUCLAS1"] == "SPECTRUM"
        assert hdr["HDUCLAS2"] == "TOTAL"
        assert hdr["HDUCLAS3"] == "COUNT"
        assert hdr["HDUCLAS4"] == "TYPE:II"
        assert hdr["HDUVERS"] == "1.0.0"
        assert hdr["HDUVERS1"] == "1.0.0"
        assert "HDUVERS2" not in hdr

        assert hdr["ORIGIN"] == "ASC"
        assert hdr["CREATOR"] == "tgextract - Version CIAO 4.8"

        assert hdr["OBJECT"] == "3C 120"
        assert hdr["MISSION"] == "AXAF"
        assert hdr["TELESCOP"] == "CHANDRA"
        assert hdr["INSTRUME"] == "ACIS"
        assert hdr["GRATING"] == "HETG"
        assert hdr["DETNAM"] == "ACIS-56789"

        assert hdr["DETCHANS"] == 8192
        assert hdr["CHANTYPE"] == "PI"
        assert "TOTCTS" not in hdr

        assert not hdr["POISSERR"]
        if roundtrip:
            assert hdr["SYS_ERR"] == 0
        else:
            assert not "SYS_ERR" in hdr

        assert hdr["QUALITY"] == 0
        assert hdr["GROUPING"] == 0

        assert hdr["CORRFILE"] == "none"
        assert hdr["CORRSCAL"] == pytest.approx(1.0)

        for key in ["ANCRFILE", "BACKFILE", "RESPFILE"]:
            assert key not in hdr

    def check_data(pha, bkg=False, roundtrip=False):
        assert len(pha.channel) == 8192
        assert len(pha.counts) == 8192
        assert pha.staterror is None
        assert pha.syserror is None
        assert len(pha.channel) == 8192
        assert len(pha.channel) == 8192
        assert pha.grouping is None
        assert pha.quality is None
        assert pha.exposure == pytest.approx(77716.294300039)
        if bkg:
            assert pha.backscal == pytest.approx(4.0188284)
            assert pha.areascal is None
        else:
            assert pha.backscal == pytest.approx(1.0)
            assert pha.areascal == pytest.approx(1.0)

        assert pha.rate
        assert pha.response_ids == []
        expected = [] if bkg or roundtrip else [1, 2]
        assert pha.background_ids == expected

    infile = make_data_path("3c120_pha2.gz")
    phas = io.read_pha(infile)
    assert len(phas) == 12

    # pick the fifth element and a few quick checks
    pha = phas[4]
    assert isinstance(pha, DataPHA)
    assert pha.name.endswith("/3c120_pha2.gz")
    check_header(pha)
    check_data(pha)

    bkg1 = pha.get_background(1)
    bkg2 = pha.get_background(2)
    check_header(bkg1)
    check_header(bkg1)
    check_data(bkg1, bkg=True)
    check_data(bkg2, bkg=True)

    outfile = tmp_path / "test.pha"
    io.write_pha(str(outfile), pha, ascii=False, clobber=False)

    outpha = io.read_pha(str(outfile))
    assert isinstance(outpha, DataPHA)
    assert outpha.name.endswith("/test.pha")
    check_header(outpha, roundtrip=True)
    check_data(outpha, roundtrip=True)


def check_csc_pha_roundtrip_crates(path):
    import pycrates  # type: ignore
    ds = pycrates.CrateDataset(str(path), mode="r")

    assert ds.get_ncrates() == 2
    cr = ds.get_crate(2)

    assert cr.name == "SPECTRUM"
    assert cr.get_colnames() == ["CHANNEL", "COUNTS"]

    c0 = cr.get_column(0)
    assert c0.name == "CHANNEL"
    assert c0.values.dtype == np.int64
    assert c0.get_tlmin() == 1
    assert c0.get_tlmax() == 1024

    c1 = cr.get_column(1)
    assert c1.name == "COUNTS"
    assert c1.values.dtype == np.float32
    assert c1.get_tlmin() == pytest.approx(-3.4028235e+38)
    assert c1.get_tlmax() == pytest.approx(3.4028235e+38)

    assert cr.get_key_value("HDUCLASS") == "OGIP"
    assert cr.get_key_value("HDUCLAS1") == "SPECTRUM"
    assert cr.get_key_value("HDUCLAS2") == "TOTAL"
    assert cr.get_key_value("HDUCLAS3") == "COUNT"
    assert cr.get_key_value("HDUCLAS4") == "TYPE:I"
    assert cr.get_key_value("HDUVERS") == "1.1.0"
    assert cr.get_key_value("HDUVERS1") == "1.1.0"

    assert cr.get_key_value("POISSERR")

    assert cr.get_key_value("TELESCOP") == "CHANDRA"
    assert cr.get_key_value("INSTRUME") == "ACIS"
    assert cr.get_key_value("FILTER") == "none"

    assert cr.get_key_value("EXPOSURE") == pytest.approx(37664.157219191)
    assert cr.get_key_value("BACKSCAL") == pytest.approx(2.2426552620567e-06)
    assert cr.get_key_value("CORRFILE") == "none"
    assert cr.get_key_value("ANCRFILE") == "acisf01575_001N001_r0085_arf3.fits"
    assert cr.get_key_value("BACKFILE") == "acisf01575_001N001_r0085_pha3.fits"
    assert cr.get_key_value("RESPFILE") == "acisf01575_001N001_r0085_rmf3.fits"

    assert cr.get_key_value("CHANTYPE") == "PI"
    assert cr.get_key_value("DETCHANS") == 1024

    assert cr.get_key_value("CORRSCAL") == 0
    assert cr.get_key_value("SYS_ERR") == 0

    assert cr.get_key_value("AREASCAL") == pytest.approx(1.0)
    assert cr.get_key_value("QUALITY") == 0
    assert cr.get_key_value("GROUPING") == 0


def check_csc_pha_roundtrip_pyfits(path):
    from astropy.io import fits  # type: ignore
    hdus = fits.open(str(path))
    try:
        assert len(hdus) == 2
        hdu = hdus[1]
        assert hdu.name == "SPECTRUM"
        assert hdu.level == 1
        assert hdu.ver == 1
        assert len(hdu.columns) == 2

        assert hdu.columns[0].name == "CHANNEL"
        assert hdu.columns[0].format == "J"
        assert hdu.columns[1].name == "COUNTS"
        assert hdu.columns[1].format == "E"

        assert hdu.header["HDUCLASS"] == "OGIP"
        assert hdu.header["HDUCLAS1"] == "SPECTRUM"
        assert hdu.header["HDUCLAS2"] == "TOTAL"
        assert hdu.header["HDUCLAS3"] == "COUNT"
        assert hdu.header["HDUCLAS4"] == "TYPE:I"
        assert hdu.header["HDUVERS"] == "1.1.0"
        assert hdu.header["HDUVERS1"] == "1.1.0"

        assert hdu.header["POISSERR"]

        assert hdu.header["TELESCOP"] == "CHANDRA"
        assert hdu.header["INSTRUME"] == "ACIS"
        assert hdu.header["FILTER"] == "none"

        assert hdu.header["EXPOSURE"] == pytest.approx(37664.157219191)
        assert hdu.header["BACKSCAL"] == pytest.approx(2.2426552620567e-06)
        assert hdu.header["CORRFILE"] == "none"
        assert hdu.header["ANCRFILE"] == "acisf01575_001N001_r0085_arf3.fits"
        assert hdu.header["BACKFILE"] == "acisf01575_001N001_r0085_pha3.fits"
        assert hdu.header["RESPFILE"] == "acisf01575_001N001_r0085_rmf3.fits"

        assert hdu.header["CHANTYPE"] == "PI"
        assert hdu.header["DETCHANS"] == 1024

        assert hdu.header["CORRSCAL"] == 0
        assert hdu.header["SYS_ERR"] == 0

        assert hdu.header["AREASCAL"] == pytest.approx(1.0)
        assert hdu.header["QUALITY"] == 0
        assert hdu.header["GROUPING"] == 0

        assert hdu.header["TLMIN1"] == 1
        assert hdu.header["TLMAX1"] == 1024

    finally:
        hdus.close()


@requires_fits
@requires_data
def test_csc_pha_roundtrip(make_data_path, tmp_path):
    """Can we read in/write out/read in a CSC file.

    These files contain both source and background spectra.
    """

    def check_header(obj, background=False):
        # check header for a selected set of keywords
        hdr = obj.header

        assert "HDUNAME" not in hdr

        # selected OGIP keywords
        assert hdr["HDUCLASS"] == "OGIP"
        assert hdr["HDUCLAS1"] == "SPECTRUM"
        assert hdr["HDUCLAS2"] == "BKG" if background else "TOTAL"
        assert hdr["HDUCLAS3"] == "COUNT"
        assert hdr["HDUVERS"] == "1.1.0"
        assert hdr["HDUVERS1"] == "1.1.0"

        assert hdr["ORIGIN"] == "ASC"
        assert hdr["CREATOR"] == "dmextract - Version CAT 2.7"
        assert hdr["REGIONID"] == "0085"

        assert hdr["OBJECT"] == "CSC"
        assert hdr["MISSION"] == "AXAF"
        assert hdr["TELESCOP"] == "CHANDRA"
        assert hdr["INSTRUME"] == "ACIS"
        assert hdr["GRATING"] == "NONE"
        assert hdr["DETNAM"] == "ACIS-67"

        assert hdr["DETCHANS"] == 1024
        assert hdr["CHANTYPE"] == "PI"
        assert hdr["TOTCTS"] == 303 if background else 855

        assert hdr["POISSERR"]
        assert hdr["SYS_ERR"] == pytest.approx(0)
        assert hdr["QUALITY"] == 0
        assert hdr["GROUPING"] == 0

        assert hdr["CORRFILE"] == "none"
        assert hdr["CORRSCAL"] == pytest.approx(0)

        for key in ["ANCRFILE", "BACKFILE", "RESPFILE"]:
            assert key not in hdr

    def check_data(obj, background=False):
        """Basic checks of the data"""

        lbscal = -4.269110831682076 if background else -5.649237480448729

        assert obj.staterror is None
        assert obj.syserror is None
        assert obj.bin_lo is None
        assert obj.bin_hi is None

        assert obj.exposure == pytest.approx(37664.157219191)
        assert np.log10(obj.backscal) == pytest.approx(lbscal)
        assert obj.areascal == pytest.approx(1.0)
        for f in ["grouped", "subtracted", "rate"]:
            assert isinstance(getattr(obj, f), bool)

        assert not obj.grouped
        assert not obj.subtracted
        assert obj.rate

        assert obj.plot_fac == 0

        assert obj.channel.dtype == np.dtype("float64")
        assert obj.counts.dtype == np.dtype("float64")

        assert obj.channel == pytest.approx(np.arange(1, 1025))
        assert len(obj.counts) == 1024
        assert obj.counts[0:11] == pytest.approx(np.zeros(11))
        assert obj.counts.sum() == 303 if background else 855
        assert obj.counts.max() == 59 if background else 15
        assert np.argmax(obj.counts) == 1023 if background else 60

        assert obj.grouping is None
        assert obj.quality is None

    infile = make_data_path("acisf01575_001N001_r0085_pha3.fits.gz")
    inpha = io.read_pha(infile)
    assert isinstance(inpha, DataPHA)
    assert inpha.name.endswith("/acisf01575_001N001_r0085_pha3.fits.gz")
    check_header(inpha)
    assert "HDUCLAS4" not in inpha.header
    check_data(inpha)

    # We do read in a response
    assert inpha.get_arf().name.endswith("/acisf01575_001N001_r0085_arf3.fits")
    assert inpha.get_rmf().name.endswith("/acisf01575_001N001_r0085_rmf3.fits")

    # Did we read in the background?
    assert inpha.background_ids == [1]
    inbkg = inpha.get_background()
    assert isinstance(inbkg, DataPHA)
    # Note that the name is slightly different to inpha
    assert inbkg.name.endswith("/acisf01575_001N001_r0085_pha3.fits")
    check_header(inbkg, background=True)
    assert "HDUCLAS4" not in inbkg.header
    check_data(inbkg, background=True)

    # We do read in a response
    assert inbkg.get_arf().name.endswith("/acisf01575_001N001_r0085_arf3.fits")
    assert inbkg.get_rmf().name.endswith("/acisf01575_001N001_r0085_rmf3.fits")

    # but not a background
    assert inbkg.get_background() is None

    outfile = tmp_path / "test.pha"
    io.write_pha(str(outfile), inpha, ascii=False, clobber=False)

    outpha = io.read_pha(str(outfile))
    assert isinstance(outpha, DataPHA)
    assert outpha.name.endswith("/test.pha")
    check_header(outpha)
    assert outpha.header["HDUCLAS4"] == "TYPE:I"
    check_data(outpha)

    # The responses and background should NOT be read in
    # (we are in a different directory to infile so we can't
    # find these files).
    #
    assert outpha.response_ids == []
    assert outpha.background_ids == []

    assert outpha.get_arf() is None
    assert outpha.get_rmf() is None
    assert outpha.get_background() is None

    if backend_is("crates"):
        check_csc_pha_roundtrip_crates(outfile)

    elif backend_is("pyfits"):
        check_csc_pha_roundtrip_pyfits(outfile)

    else:
        raise RuntimeError(f"Unknown io backend: {io.backend}")


@pytest.mark.parametrize("column", ["channel", "counts"])
@pytest.mark.parametrize("is_ascii", [True, False])
@requires_fits
def test_write_pha_missing_column(column, is_ascii, tmp_path):
    """What happens when there's no "required" column.

    At the moment we write out the file.
    """

    if column == "channel":
        chans = None
        counts = np.arange(1, 4, dtype=np.int32)
    elif column == "counts":
        chans = np.arange(1, 4, dtype=np.int32)
        counts = None
    else:
        assert False, "unknown column"

    pha = DataPHA("tmp", chans, counts)

    tmpfile = tmp_path / 'test.out'
    assert not tmpfile.exists()
    io.write_pha(str(tmpfile), pha, ascii=is_ascii)

    # We want to check we've created a single-column file, but
    # it depends on the backend and whether it was written
    # out in ASCII or FITS format. As we don't really care here,
    # all we do is check the output file exists and is not empty.
    #
    assert tmpfile.exists()
    assert tmpfile.is_file()
    assert tmpfile.stat().st_size > 0


def test_pack_pha_invalid_counts():
    """Use pack_pha rather than write_pha for a change

    This should work even with the dummy backend.
    """

    pha = DataPHA("dummy",
                  np.asarray([1, 2], dtype=np.int16),
                  np.asarray([complex(1), complex(2)]))

    emsg = "The PHA dataset 'dummy' contains an unsupported COUNTS column"
    with pytest.raises(DataErr, match=f"^{emsg}$"):
        io.pack_pha(pha)


@requires_data
@requires_fits
def test_read_hrci_rmf(make_data_path):
    """Check out issue #1830

    Can the pyfits backend read in the HRC-I RMF?
    """

    rmffile = make_data_path("chandra_hrci/hrcf24564_000N030_r0001"
                             "_rmf3.fits.gz")
    rmf = io.read_rmf(rmffile)

    NROWS = 16921
    NGRP = 16960
    NMATRIX = 8525597
    NCHAN = 1024

    assert len(rmf.energ_lo) == NROWS
    assert rmf.energ_lo[0] == pytest.approx(0.06)
    assert rmf.energ_hi[-1] == pytest.approx(9.9902925491)

    assert len(rmf.n_grp) == NROWS
    assert rmf.n_grp.sum() == NGRP

    assert len(rmf.f_chan) == NGRP
    assert len(rmf.n_chan) == NGRP
    assert rmf.n_chan.sum() == NMATRIX

    # The value checks are taken from a crates run and so
    # act as a regression test.
    #
    assert len(rmf.matrix) == NMATRIX
    assert rmf.matrix[0] == pytest.approx(0.00013532221782952547)
    assert rmf.matrix[-100] == pytest.approx(0.00029959273524582386)

    assert rmf.detchans == NCHAN
    assert len(rmf.e_min) == NCHAN
    assert len(rmf.e_max) == NCHAN

    assert rmf.e_min[0] == pytest.approx(0.06)
    assert rmf.e_max[-1] == pytest.approx(10)


@requires_fits
@requires_data
def test_read_pha_object(make_data_path):
    """Check we can send in a 'object'.

    This support is not well advertised so it could perhaps be
    removed, but let's add a basic test at least.

    """

    # Could create the "object" manually, but let's just read one in
    #
    infile = make_data_path("acisf01575_001N001_r0085_pha3.fits.gz")
    close = False

    if backend_is("crates"):
        import pycrates  # type: ignore
        arg = pycrates.PHACrateDataset(infile, mode="r")

    elif backend_is("pyfits"):
        from astropy.io import fits  # type: ignore
        arg = fits.open(infile)
        close = True

    else:
        assert False, f"unknown backend: {io.backend}"

    try:
        pha = io.read_pha(arg)

    finally:
        if close:
            arg.close()

    assert isinstance(pha, DataPHA)

    # basic test that the data is read in (do not do a complete
    # test)
    assert pha.counts.max() == pytest.approx(15.0)
    assert pha.exposure == pytest.approx(37664.157219191)

    # Check we've loaded in the response and background
    #
    # For 4.16.0 and earlier, the crates backend "fails" because it
    # does not load in the response or background datasets, instead
    # saying:
    #
    # WARNING: File acisf01575_001N001_r0085_arf3.fits does not exist.
    # WARNING: File acisf01575_001N001_r0085_rmf3.fits does not exist.
    # WARNING: File acisf01575_001N001_r0085_pha3.fits does not exist.
    #
    # So there is an issue in handling compressed files, since the
    # responses are saved as .gz files, although perhaps its more
    # that the paths aren't being set up correctly. The background
    # should be taken from the object itself, but perhaps it's
    # easiest just to re-read it, and then we hit the same issue.
    #
    # For now treat this as a regression test so we can find out
    # when it's fixed.
    #
    if backend_is("crates"):
        assert pha.response_ids == []
        assert pha.background_ids == []
        return

    assert len(pha.response_ids) == 1
    assert len(pha.background_ids) == 1

    assert pha.get_arf().specresp.max() == pytest.approx(711.6605834960938)
    assert pha.get_rmf().matrix.max() == pytest.approx(0.1611403226852417)

    assert pha.get_background().counts[-1] == pytest.approx(59)


@requires_fits
def test_roundtrip_channel0(tmp_path):
    """Check what happens if we write out and read back in channel 0.

    This is a regression test.
    """

    pha = DataPHA("chan0", [0, 1, 2], [2, 5, 9])
    pha.exposure = 12.0
    pha.backscal = 1.2
    pha.areascal = 0.4

    outfile = str(tmp_path / "out.pha")
    io.write_pha(outfile, pha, ascii=False)
    pha = None

    pha2 = io.read_pha(outfile)
    assert pha2.channel == pytest.approx([0, 1, 2])
    assert pha2.counts == pytest.approx([2, 5, 9])
    assert pha2.exposure == pytest.approx(12)
    assert pha2.backscal == pytest.approx(1.2)
    assert pha2.areascal == pytest.approx(0.4)


@requires_fits
def test_write_pha_with_bad_quality(tmp_path):
    """If ignore_bad has been called what is written out?"""

    # This gives groups of
    #
    # | channel | counts | filtered |
    # | ------- | ------ | -------- |
    # |  1 - 3  |   12   |     8    |
    # |  4 - 5  |   18   |    18    |
    # |    6    |   12   |    12    |
    # |    7    |   14   |    14    |
    # |  8 - 9  |   34   |     -    |
    #
    chans = np.arange(1, 10)
    counts = chans * 2
    group = [1, -1, -1, 1, -1, 1, 1, 1, -1]
    quality = [0, 5, 0, 0, 0, 0, 0, 2, 2]
    qfilt = [True, False] + [True] * 5 + [False] * 2

    pha0 = DataPHA("qual", chans, counts, grouping=group,
                   quality=quality)
    pha0.ignore_bad()
    assert pha0.quality_filter == pytest.approx(qfilt)

    outpath = tmp_path / "test.pha"
    outfile = str(outpath)
    io.write_pha(outfile, pha0, ascii=False, clobber=True)

    # Check that the quality filter field hasn't changed
    assert pha0.quality_filter == pytest.approx(qfilt)

    # Ensure all rows are written out
    pha1 = io.read_pha(outfile)
    assert pha1.channel == pytest.approx(chans)
    assert pha1.counts == pytest.approx(counts)
    assert pha1.grouping == pytest.approx(group)
    assert pha1.quality == pytest.approx(quality)
    assert pha1.quality_filter is None
