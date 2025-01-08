#
#  Copyright (C) 2024
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

"""Test handling of Chandra data.

Most of the Chandra testing is done in other files.

"""

import numpy as np

import pytest

from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.astro import io
from sherpa.astro import ui
from sherpa.utils.err import IOErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, requires_fits


# We potentially have different instruments to test.
#
PHAFILE = "obs3.pi"
ARFFILE = "obs3.arf"
RMFFILE = "obs3.rmf"
PHA_EXPOSURE = 68936.36552833

PHAFILE_GRATING = "3c120_heg_-1.pha.gz"
ARFFILE_GRATING = "3c120_heg_-1.arf.gz"
RMFFILE_GRATING = "3c120_heg_-1.rmf.gz"
PHA_GRATING_EXPOSURE = 77716.294300039
ARF_GRATING_EXPOSURE = 77717.297715256


def check_pha(pha, responses=True):
    """Does the data match expectation?"""

    assert isinstance(pha, DataPHA)

    assert pha.exposure == pytest.approx(PHA_EXPOSURE)
    assert pha.areascal == pytest.approx(1.0)
    assert pha.backscal == pytest.approx(1.63383e-7)

    nchan = 1024
    assert pha.channel == pytest.approx(np.arange(1, nchan + 1))

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(pha.counts) == nchan
    assert pha.counts.min() == pytest.approx(0.0)
    assert pha.counts.max() == pytest.approx(2.0)
    assert pha.counts.sum() == pytest.approx(77)
    assert np.argmax(pha.counts) == 10

    for field in ['staterror', 'syserror', 'quality', 'grouping',
                  'bin_lo', 'bin_hi']:
        assert getattr(pha, field) is None

    assert pha.grouped is False
    assert pha.subtracted is False
    assert pha.rate
    assert pha.plot_fac == 0

    if responses:
        assert pha.units == 'energy'
        assert pha.response_ids == [1]
        assert pha.background_ids == [1]

        barf, brmf = pha.get_response()
        assert barf.name.endswith("/obs3.arf")
        assert brmf.name.endswith("/obs3.rmf")

        b1 = pha.get_background(1)
        assert b1.name.endswith("obs3_bg.pi")
        assert b1.areascal == pytest.approx(1.0)
        assert b1.backscal == pytest.approx(3.26767e-7)
        assert b1.exposure == pytest.approx(PHA_EXPOSURE)

        assert b1.counts.min() == pytest.approx(0.0)
        assert b1.counts.max() == pytest.approx(4.0)
        assert b1.counts.sum() == pytest.approx(76.0)

    else:
        assert pha.units == 'channel'
        assert pha.response_ids == []
        assert pha.background_ids == []

        assert pha.get_background() is None

    # Select a few header keywords. Some are OGIP related,
    # some are "we expect them", and some are there just to
    # check we encode the type / value correctly.
    #
    assert pha.header["HDUCLASS"] == "OGIP"
    assert pha.header["HDUCLAS1"] == "SPECTRUM"
    assert pha.header["HDUCLAS2"] == "TOTAL"
    assert pha.header["HDUCLAS3"] == "TYPE:I"
    assert pha.header["HDUCLAS4"] == "COUNT"
    assert pha.header["DETCHANS"] == nchan
    assert pha.header["CHANTYPE"] == "PI"

    # This keyword is removed from the header.
    assert "POISERR" not in pha.header

    assert pha.header["TELESCOP"] == "CHANDRA"
    assert pha.header["INSTRUME"] == "ACIS"
    assert pha.header["GRATING"] == "NONE"
    assert "FILTER" not in pha.header["FILTER"]
    assert pha.header["OBJECT"] == "4C19.44"

    assert pha.header["RA_TARG"] == pytest.approx(209.26875)
    assert pha.header["DEC_TARG"] == pytest.approx(19.316944)
    assert "RA_OBJ" not in pha.header
    assert "DEC_OBJ" not in pha.header


def check_arf(arf):
    """Does the ARF contain the expected values?"""

    assert isinstance(arf, DataARF)

    nbins = 1078

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert len(arf.energ_lo) == nbins
    assert (arf.energ_lo[1:] == arf.energ_hi[:-1]).all()
    assert arf.energ_lo[0] == pytest.approx(0.22)
    assert arf.energ_hi[0] == pytest.approx(0.23)
    assert arf.energ_lo[-1] == pytest.approx(10.99)
    assert arf.energ_hi[-1] == pytest.approx(11.0)

    assert arf.bin_lo is None
    assert arf.bin_hi is None

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(arf.specresp) == nbins
    assert arf.specresp.min() == pytest.approx(0.5785124898)
    assert arf.specresp.max() == pytest.approx(686.1730957)
    assert arf.specresp.sum() == pytest.approx(254262.18811)
    assert np.argmax(arf.specresp) == 132

    assert arf.exposure == pytest.approx(PHA_EXPOSURE)


def check_rmf(rmf):
    """Does the RMF contain the expected values?"""

    assert isinstance(rmf, DataRMF)

    nchans = 1024
    nbins = 1078
    nsize1 = 1484
    nsize2 = 429463

    assert rmf.detchans == nchans
    assert rmf.offset == 1

    for field in ['energ_lo', 'energ_hi']:
        assert getattr(rmf, field).size == nbins

    assert rmf.n_grp.size == nbins
    assert rmf.n_grp.min() == 1
    assert rmf.n_grp.max() == 3
    assert rmf.n_grp.sum() == nsize1

    assert rmf.f_chan.size == nsize1
    assert rmf.f_chan.min() == 1
    assert rmf.f_chan.max() == 728
    assert rmf.f_chan.sum() == 267802

    assert rmf.n_chan.size == nsize1
    assert rmf.n_chan.min() == 6
    assert rmf.n_chan.max() == 1023
    assert rmf.n_chan.sum() == nsize2

    assert rmf.matrix.size == nsize2
    assert rmf.matrix.min() == pytest.approx(8.832556908089373e-09)
    assert rmf.matrix.max() == pytest.approx(0.13701042532920837)
    assert rmf.matrix.sum() == pytest.approx(1078)
    assert np.argmax(rmf.matrix) == 5086

    for field in ['e_min', 'e_max']:
        assert getattr(rmf, field).size == nchans

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert rmf.energ_lo[1:] == pytest.approx(rmf.energ_hi[:-1])

    assert rmf.energ_lo[0] == pytest.approx(0.22)
    assert rmf.energ_hi[0] == pytest.approx(0.23)
    assert rmf.energ_lo[-1] == pytest.approx(10.99)
    assert rmf.energ_hi[-1] == pytest.approx(11.0)

    assert len(rmf.e_min) == 1024
    assert rmf.e_min[1:] == pytest.approx(rmf.e_max[:-1])

    assert rmf.e_min[0] == pytest.approx(0.0014600000577047467)
    assert rmf.e_max[0] == pytest.approx(0.014600000344216824)
    assert rmf.e_min[-1] == pytest.approx(14.935799598693848)
    assert rmf.e_max[-1] == pytest.approx(14.950400352478027)


def check_pha_grating(pha, errors=False, responses=True):
    """Does the data match expectation?"""

    assert isinstance(pha, DataPHA)

    assert pha.exposure == pytest.approx(PHA_GRATING_EXPOSURE)
    assert pha.areascal == pytest.approx(1.0)
    assert pha.backscal == pytest.approx(1.0)

    nchan = 8192
    assert pha.channel == pytest.approx(np.arange(1, nchan + 1))

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(pha.counts) == nchan
    assert pha.counts.min() == pytest.approx(0.0)
    assert pha.counts.max() == pytest.approx(15.0)
    assert pha.counts.sum() == pytest.approx(16298.0)
    assert np.argmax(pha.counts) == 5871

    assert len(pha.bin_lo) == nchan
    assert len(pha.bin_hi) == nchan
    assert pha.bin_hi[1:] == pytest.approx(pha.bin_lo[:-1])
    assert pha.bin_lo[0] == pytest.approx(21.4775)
    assert pha.bin_hi[0] == pytest.approx(21.48)
    assert pha.bin_lo[-1] == pytest.approx(1.0)
    assert pha.bin_hi[-1] == pytest.approx(1.0025)

    for field in ['syserror', 'quality', 'grouping']:
        assert getattr(pha, field) is None

    if errors:
        assert len(pha.staterror) == nchan
        assert pha.staterror.min() == pytest.approx(1.8660254)
        assert pha.staterror.max() == pytest.approx(4.968627)
        assert pha.staterror.sum() == pytest.approx(20342.717)
    else:
        assert pha.staterror is None

    assert pha.grouped is False
    assert pha.subtracted is False
    assert pha.rate
    assert pha.plot_fac == 0

    if responses:
        assert pha.units == 'energy'
        assert pha.response_ids == [1]
        assert pha.background_ids == [1, 2]

        barf, brmf = pha.get_response()
        assert barf.name.endswith("/3c120_heg_-1.arf")
        assert brmf.name.endswith("/3c120_heg_-1.rmf")

        b1 = pha.get_background(1)
        b2 = pha.get_background(2)
        assert b1.name.endswith("3c120_heg_-1.pha.gz")
        assert b2.name.endswith("3c120_heg_-1.pha.gz")

        for b in [b1, b2]:
            assert getattr(b, "backscal") == pytest.approx(4.0188284)
            assert getattr(b, "exposure") == pytest.approx(PHA_GRATING_EXPOSURE)
            assert getattr(b, "staterror") is None

        assert b1.counts.min() == pytest.approx(0.0)
        assert b1.counts.max() == pytest.approx(3.0)
        assert b1.counts.sum() == pytest.approx(513.0)

        assert b2.counts.min() == pytest.approx(0.0)
        assert b2.counts.max() == pytest.approx(2.0)
        assert b2.counts.sum() == pytest.approx(405.0)

    else:
        assert pha.units == 'channel'
        assert pha.response_ids == []
        assert pha.background_ids == []

        assert pha.get_background() is None

    # Select a few header keywords. Some are OGIP related,
    # some are "we expect them", and some are there just to
    # check we encode the type / value correctly.
    #
    assert pha.header["HDUCLASS"] == "OGIP"
    assert pha.header["HDUCLAS1"] == "SPECTRUM"
    assert pha.header["HDUCLAS2"] == "TOTAL"
    assert pha.header["HDUCLAS3"] == "COUNT"
    assert pha.header["HDUCLAS4"] == "TYPE:I"
    assert pha.header["DETCHANS"] == nchan
    assert pha.header["CHANTYPE"] == "PI"

    # This keyword is removed from the header.
    assert "POISERR" not in pha.header

    assert pha.header["TELESCOP"] == "CHANDRA"
    assert pha.header["INSTRUME"] == "ACIS"
    assert pha.header["GRATING"] == "HETG"
    assert pha.header["FILTER"] == ""
    assert pha.header["OBJECT"] == "3C 120"

    assert pha.header["RA_TARG"] == pytest.approx(68.29625)
    assert pha.header["DEC_TARG"] == pytest.approx(5.354333)
    assert "RA_OBJ" not in pha.header
    assert "DEC_OBJ" not in pha.header

    assert pha.header["SPEC_NUM"] == 3
    assert pha.header["TG_M"] == -1
    assert pha.header["TG_PART"] == 1
    assert pha.header["TG_SRCID"] == 1


def check_arf_grating(arf):
    """Does the ARF contain the expected values?"""

    assert isinstance(arf, DataARF)

    nbins = 8192

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert len(arf.energ_lo) == nbins
    assert (arf.energ_lo[1:] == arf.energ_hi[:-1]).all()
    assert arf.energ_lo[0] == pytest.approx(0.57720757)
    assert arf.energ_hi[0] == pytest.approx(0.57727474)
    assert arf.energ_lo[-1] == pytest.approx(12.36749935)
    assert arf.energ_hi[-1] == pytest.approx(12.39841843)

    assert len(arf.bin_lo) == nbins
    assert len(arf.bin_hi) == nbins
    assert arf.bin_hi[1:] == pytest.approx(arf.bin_lo[:-1])
    assert arf.bin_lo[0] == pytest.approx(21.4775)
    assert arf.bin_hi[0] == pytest.approx(21.48)
    assert arf.bin_lo[-1] == pytest.approx(1.0)
    assert arf.bin_hi[-1] == pytest.approx(1.0025)

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(arf.specresp) == nbins
    assert arf.specresp.min() == pytest.approx(0.0)
    assert arf.specresp.max() == pytest.approx(28.142807006835938)
    assert arf.specresp.sum() == pytest.approx(61266.833162411756)
    assert np.argmax(arf.specresp) == 5849

    assert arf.exposure == pytest.approx(ARF_GRATING_EXPOSURE)


def check_rmf_grating(rmf):
    """Does the RMF contain the expected values?"""

    assert isinstance(rmf, DataRMF)

    nchans = 8192
    nbins = 8192
    nsize1 = 8192
    nsize2 = 841124

    assert rmf.detchans == nchans
    assert rmf.offset == 1

    for field in ['energ_lo', 'energ_hi']:
        assert getattr(rmf, field).size == nbins

    assert rmf.n_grp.size == nbins
    assert rmf.n_grp.min() == 1
    assert rmf.n_grp.max() == 1
    assert rmf.n_grp.sum() == nsize1

    assert rmf.f_chan.size == nsize1
    assert rmf.f_chan.min() == 1
    assert rmf.f_chan.max() == 8141
    assert rmf.f_chan.sum() == 33142062

    assert rmf.n_chan.size == nsize1
    assert rmf.n_chan.min() == 52
    assert rmf.n_chan.max() == 103
    assert rmf.n_chan.sum() == nsize2

    assert rmf.matrix.size == nsize2
    assert rmf.matrix.min() == pytest.approx(1.1740531890333088e-06)
    assert rmf.matrix.max() == pytest.approx(0.31429963429969676)
    assert rmf.matrix.sum() == pytest.approx(7860.054242052272)
    assert np.argmax(rmf.matrix) == 0

    for field in ['e_min', 'e_max']:
        assert getattr(rmf, field).size == nchans

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert rmf.energ_lo[1:] == pytest.approx(rmf.energ_hi[:-1])

    assert rmf.energ_lo[0] == pytest.approx(0.57720757)
    assert rmf.energ_hi[0] == pytest.approx(0.57727474)
    assert rmf.energ_lo[-1] == pytest.approx(12.36749935)
    assert rmf.energ_hi[-1] == pytest.approx(12.39841843)

    assert rmf.e_min == pytest.approx(rmf.energ_lo)
    assert rmf.e_max == pytest.approx(rmf.energ_hi)


def record_starts_with(r, start, level="INFO"):
    """This only checks the start of the message, as it contains
    test-specific information (the location of the file name).
    """

    assert r.name == "sherpa.astro.io"
    assert r.levelname == level

    # We use this, rather than msg.startswith(start), as the error
    # message is more useful if there is a difference.
    #
    nchar = len(start)
    msg = r.getMessage()
    assert msg[:nchar] == start


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_read_pha(errors, make_data_path, caplog):
    """Can we read in a Chandra PHA file."""

    infile = make_data_path(PHAFILE)
    with SherpaVerbosity("INFO"):
        pha = io.read_pha(infile, use_errors=errors)

    check_pha(pha, responses=True)

    assert len(caplog.records) == 3

    record_starts_with(caplog.records[0],
                       "read ARF file /")
    record_starts_with(caplog.records[1],
                       "read RMF file /")
    record_starts_with(caplog.records[2],
                       "read background file /")


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_roundtrip_pha(errors, make_data_path, tmp_path, caplog):
    """Can we write out a Chandra PHA file and then read it back in?"""

    infile = make_data_path(PHAFILE)
    with SherpaVerbosity("ERROR"):
        pha1 = io.read_pha(infile, use_errors=errors)

    assert len(caplog.records) == 0

    outpath = tmp_path / "test.pha"
    outfile = str(outpath)
    io.write_pha(outfile, pha1, ascii=False, clobber=True)

    with SherpaVerbosity("INFO"):
        pha2 = io.read_pha(outfile, use_errors=errors)

    check_pha(pha2, responses=False)

    assert len(caplog.records) == 3

    record_starts_with(caplog.records[0],
                       "unable to read ARF: ",
                       level="WARNING")
    record_starts_with(caplog.records[1],
                       "unable to read RMF: ",
                       level="WARNING")
    record_starts_with(caplog.records[2],
                       "unable to read background: ",
                       level="WARNING")


@requires_data
@requires_fits
def test_read_arf(make_data_path):
    """Can we read in a Chandra ARF."""

    infile = make_data_path(ARFFILE)
    arf = io.read_arf(infile)
    check_arf(arf)


@requires_data
@requires_fits
def test_roundtrip_arf(make_data_path, tmp_path):
    """Can we write out a Chandra ARF file and then read it back in?"""

    infile = make_data_path(ARFFILE)
    arf1 = io.read_arf(infile)

    outpath = tmp_path / "test.arf"
    outfile = str(outpath)
    io.write_arf(outfile, arf1, ascii=False, clobber=True)

    arf2 = io.read_arf(outfile)
    check_arf(arf2)


@requires_data
@requires_fits
def test_read_rmf(make_data_path):
    """Can we read in a Chandra RMF."""

    infile = make_data_path(RMFFILE)
    rmf = io.read_rmf(infile)
    check_rmf(rmf)


@requires_data
@requires_fits
def test_roundtrip_rmf(make_data_path, tmp_path):
    """Can we write out a Chandra RMF file and then read it back in?"""

    infile = make_data_path(RMFFILE)
    rmf1 = io.read_rmf(infile)

    outpath = tmp_path / "test.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, rmf1, clobber=True)

    rmf2 = io.read_rmf(outfile)
    check_rmf(rmf2)


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_read_pha_grating(errors, make_data_path, caplog):
    """Can we read in a Chandra grating PHA file."""

    infile = make_data_path(PHAFILE_GRATING)
    with SherpaVerbosity("INFO"):
        pha = io.read_pha(infile, use_errors=errors)

    check_pha_grating(pha, errors=errors, responses=True)

    if errors:
        assert len(caplog.records) == 4
        idx = 0
    else:
        assert len(caplog.records) == 6

        record_starts_with(caplog.records[0],
                           "systematic errors were not found in file '/",
                           level="WARNING")
        record_starts_with(caplog.records[1],
                           "statistical errors were found in file '/")

        idx = 2

    record_starts_with(caplog.records[idx],
                       "read ARF file /")
    idx += 1
    record_starts_with(caplog.records[idx],
                       "read RMF file /")
    idx += 1
    record_starts_with(caplog.records[idx],
                       "read background_up into a dataset from file /")
    idx += 1
    record_starts_with(caplog.records[idx],
                       "read background_down into a dataset from file /")


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_roundtrip_pha_grating(errors, make_data_path, tmp_path, caplog):
    """Can we write out a Chandra grating PHA file and then read it back in?"""

    infile = make_data_path(PHAFILE_GRATING)
    with SherpaVerbosity("ERROR"):
        pha1 = io.read_pha(infile, use_errors=errors)

    assert len(caplog.records) == 0

    outpath = tmp_path / "test.pha"
    outfile = str(outpath)
    io.write_pha(outfile, pha1, ascii=False, clobber=True)

    with SherpaVerbosity("INFO"):
        pha2 = io.read_pha(outfile, use_errors=errors)

    check_pha_grating(pha2, errors=errors, responses=False)

    assert len(caplog.records) == 3

    record_starts_with(caplog.records[0],
                       "unable to read ARF: ",
                       level="WARNING")
    record_starts_with(caplog.records[1],
                       "unable to read RMF: ",
                       level="WARNING")
    record_starts_with(caplog.records[2],
                       "unable to read background: ",
                       level="WARNING")


@requires_data
@requires_fits
def test_read_arf_grating(make_data_path):
    """Can we read in a Chandra grating ARF."""

    infile = make_data_path(ARFFILE_GRATING)
    arf = io.read_arf(infile)
    check_arf_grating(arf)


@requires_data
@requires_fits
def test_roundtrip_arf_grating(make_data_path, tmp_path):
    """Can we write out a Chandra grating ARF file and then read it back in?"""

    infile = make_data_path(ARFFILE_GRATING)
    arf1 = io.read_arf(infile)

    outpath = tmp_path / "test.arf"
    outfile = str(outpath)
    io.write_arf(outfile, arf1, ascii=False, clobber=True)

    arf2 = io.read_arf(outfile)
    check_arf_grating(arf2)


@requires_data
@requires_fits
def test_read_rmf_grating(make_data_path):
    """Can we read in a Chandra grating RMF."""

    infile = make_data_path(RMFFILE_GRATING)
    rmf = io.read_rmf(infile)
    check_rmf_grating(rmf)


@requires_data
@requires_fits
def test_roundtrip_rmf_grating(make_data_path, tmp_path):
    """Can we write out a Chandra grating RMF file and then read it back in?"""

    infile = make_data_path(RMFFILE_GRATING)
    rmf1 = io.read_rmf(infile)

    outpath = tmp_path / "test.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, rmf1, clobber=True)

    rmf2 = io.read_rmf(outfile)
    check_rmf_grating(rmf2)
