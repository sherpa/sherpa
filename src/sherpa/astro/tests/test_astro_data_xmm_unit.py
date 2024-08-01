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

"""Test handling of XMM data.

Most of the XMM testing is done in other files.

"""

import warnings

import numpy as np

import pytest

from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.astro import io
from sherpa.astro import ui
from sherpa.utils.err import IOErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, requires_fits, requires_xspec


# We potentially have different instruments to test.
#
PHAFILE = "MNLup_2138_0670580101_EMOS1_S001_spec.15grp"
ARFFILE = "MNLup_2138_0670580101_EMOS1_S001_spec.arf"
RMFFILE = "MNLup_2138_0670580101_EMOS1_S001_spec.rmf"
PHA_EXPOSURE = 42589.380

PHAFILE_GRATING = "xmmrgs/P0112880201R1S004SRSPEC1003.FTZ"
RMFFILE_GRATING = "xmmrgs/P0112880201R1S004RSPMAT1003.FTZ"

EMIN = 1.0e-10


def check_pha(pha, responses=True):
    """Does the data match expectation?"""

    assert isinstance(pha, DataPHA)

    assert pha.exposure == pytest.approx(PHA_EXPOSURE)
    assert pha.areascal == pytest.approx(1.0)
    assert pha.backscal == pytest.approx(1136000.)

    nchan = 800
    assert pha.channel == pytest.approx(np.arange(0, nchan))

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(pha.counts) == nchan
    assert pha.counts.min() == pytest.approx(0.0)
    assert pha.counts.max() == pytest.approx(40.0)
    assert pha.counts.sum() == pytest.approx(1834)
    assert np.argmax(pha.counts) == 61

    assert len(pha.grouping) == nchan
    assert (pha.grouping == 1).sum() == 223
    assert (pha.grouping == -1).sum() == (nchan - 223)

    assert len(pha.quality) == nchan
    assert pha.quality.min() == 0
    assert pha.quality.max() == 2
    assert pha.quality.sum() == 276
    assert pha.quality.argmax() == 662

    for field in ['staterror', 'syserror',
                  'bin_lo', 'bin_hi']:
        assert getattr(pha, field) is None

    assert pha.grouped is True
    assert pha.subtracted is False
    assert pha.rate
    assert pha.plot_fac == 0

    if responses:
        assert pha.units == 'energy'
        assert pha.response_ids == [1]
        assert pha.background_ids == [1]

        barf, brmf = pha.get_response()
        assert barf.name.endswith("/MNLup_2138_0670580101_EMOS1_S001_spec.arf")
        assert brmf.name.endswith("/MNLup_2138_0670580101_EMOS1_S001_spec.rmf")

        b1 = pha.get_background(1)
        assert b1.name.endswith("MNLup_2138_0670580101_EMOS1_S001_specbg.fits")
        assert b1.areascal == pytest.approx(1.0)
        assert b1.backscal == pytest.approx(19487580)
        assert b1.exposure == pytest.approx(PHA_EXPOSURE)

        assert b1.counts.min() == pytest.approx(0.0)
        assert b1.counts.max() == pytest.approx(42.0)
        assert b1.counts.sum() == pytest.approx(2880.0)

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
    assert pha.header["DETCHANS"] == nchan
    assert pha.header["CHANTYPE"] == "PI"

    # This keyword is removed from the header.
    assert "POISERR" not in pha.header

    assert pha.header["TELESCOP"] == "XMM"
    assert pha.header["INSTRUME"] == "EMOS1"
    assert "GRATING" not in pha.header
    assert pha.header["FILTER"] == "Medium"
    assert "OBJECT" not in pha.header


def check_arf(arf):
    """Does the ARF contain the expected values?"""

    assert isinstance(arf, DataARF)

    nbins = 2400

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert len(arf.energ_lo) == nbins
    assert (arf.energ_lo[1:] == arf.energ_hi[:-1]).all()
    assert arf.energ_lo[0] == pytest.approx(EMIN)
    assert arf.energ_hi[0] == pytest.approx(0.005)
    assert arf.energ_lo[-1] == pytest.approx(11.9949998856)
    assert arf.energ_hi[-1] == pytest.approx(12.0)

    assert arf.bin_lo is None
    assert arf.bin_hi is None

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(arf.specresp) == nbins
    assert arf.specresp.min() == pytest.approx(0.0)
    assert arf.specresp.max() == pytest.approx(449.26184082)
    assert arf.specresp.sum() == pytest.approx(396003.62812)
    assert np.argmax(arf.specresp) == 367

    assert arf.exposure is None


def check_rmf(rmf):
    """Does the RMF contain the expected values?"""

    assert isinstance(rmf, DataRMF)

    nchans = 800
    nbins = 2400
    nsize1 = 2394
    nsize2 = 1281216

    assert rmf.detchans == nchans
    assert rmf.offset == 0

    for field in ['energ_lo', 'energ_hi']:
        assert getattr(rmf, field).size == nbins

    assert rmf.n_grp.size == nbins
    assert rmf.n_grp.min() == 0
    assert rmf.n_grp.max() == 1
    assert rmf.n_grp.sum() == nsize1

    assert rmf.f_chan.size == nsize1
    assert rmf.f_chan.min() == 0
    assert rmf.f_chan.max() == 0

    assert rmf.n_chan.size == nsize1
    assert rmf.n_chan.min() == 24
    assert rmf.n_chan.max() == 800
    assert rmf.n_chan.sum() == nsize2

    assert rmf.matrix.size == nsize2
    assert rmf.matrix.min() == pytest.approx(1.0056578503281344e-06)
    assert rmf.matrix.max() == pytest.approx(0.5105319619178772)
    assert rmf.matrix.sum() == pytest.approx(2394.0000001526996)
    assert np.argmax(rmf.matrix) == 0

    for field in ['e_min', 'e_max']:
        assert getattr(rmf, field).size == nchans

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert rmf.energ_lo[1:] == pytest.approx(rmf.energ_hi[:-1])

    assert rmf.energ_lo[0] == pytest.approx(EMIN)
    assert rmf.energ_hi[0] == pytest.approx(0.005)
    assert rmf.energ_lo[-1] == pytest.approx(11.9949998856)
    assert rmf.energ_hi[-1] == pytest.approx(12.0)

    assert len(rmf.e_min) == nchans
    assert rmf.e_min[1:] == pytest.approx(rmf.e_max[:-1])

    assert rmf.e_min[0] == pytest.approx(0.0)
    assert rmf.e_max[0] == pytest.approx(0.014999999665)
    assert rmf.e_min[-1] == pytest.approx(11.984999657)
    assert rmf.e_max[-1] == pytest.approx(12.0)


def check_pha_grating(pha, errors=False):
    """Does the data match expectation?"""

    assert isinstance(pha, DataPHA)

    assert pha.exposure == pytest.approx(28965.640625)
    assert pha.backscal == pytest.approx(1.0)

    nchan = 3600
    assert pha.channel == pytest.approx(np.arange(1, nchan + 1))

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(pha.areascal) == nchan
    assert pha.areascal.min() == pytest.approx(0.0)
    assert pha.areascal.max() == pytest.approx(1.0)
    assert pha.areascal.sum() == pytest.approx(2950.6953666)
    assert pha.areascal.min() == pytest.approx(0.0)

    non_zero = pha.areascal > 0
    assert non_zero.sum() == 2964
    assert pha.areascal[non_zero].min() == pytest.approx(0.01542056445)

    assert len(pha.counts) == nchan
    assert pha.counts.min() == pytest.approx(-6.0)
    assert pha.counts.max() == pytest.approx(63.0)
    assert pha.counts.sum() == pytest.approx(2724.0)
    assert np.argmax(pha.counts) == 1498

    for field in ['syserror', 'bin_lo', 'bin_hi',
                  'grouping']:
        assert getattr(pha, field) is None

    if errors:
        assert len(pha.staterror) == nchan
        assert pha.staterror.min() == pytest.approx(1.0)
        assert pha.staterror.max() == pytest.approx(7.937253952)
        assert pha.staterror.sum() == pytest.approx(4496.9490242)
    else:
        assert pha.staterror is None

    assert len(pha.quality) == nchan
    assert pha.quality.min() == pytest.approx(0.0)
    assert pha.quality.max() == pytest.approx(1.0)
    assert pha.quality.sum() == pytest.approx(741.0)

    assert pha.grouped is False
    assert pha.subtracted is False
    assert pha.units == 'channel'
    assert pha.rate
    assert pha.plot_fac == 0
    assert pha.response_ids == []
    assert pha.background_ids == []
    assert pha.get_background() is None

    # Select a few header keywords. Some are OGIP related,
    # some are "we expect them", and some are there just to
    # check we encode the type / value correctly.
    #
    assert pha.header["HDUCLASS"] == "OGIP"
    assert pha.header["HDUCLAS1"] == "SPECTRUM"
    assert pha.header["HDUCLAS2"] == "NET"
    assert pha.header["HDUCLAS3"] == "COUNT"
    assert pha.header["DETCHANS"] == nchan
    assert pha.header["CHANTYPE"] == "PI"

    # This keyword is removed from the header.
    assert "POISERR" not in pha.header

    assert pha.header["TELESCOP"] == "XMM"
    assert pha.header["INSTRUME"] == "RGS1"
    assert pha.header["FILTER"] == "UNKNOWN"
    assert pha.header["OBJECT"] == "TW Hya"

    assert pha.header["RA_OBJ"] == pytest.approx(165.4667511)
    assert pha.header["DEC_OBJ"] == pytest.approx(-34.70463943)

    assert pha.header["SOURCEID"] == 3


def check_rmf_grating(rmf):
    """Does the RMF contain the expected values?"""

    assert isinstance(rmf, DataRMF)

    nchans = 3600
    nbins = 4000
    nsize1 = 82878
    nsize2 = 8475822

    assert rmf.detchans == nchans
    assert rmf.offset == 1

    for field in ['energ_lo', 'energ_hi']:
        assert getattr(rmf, field).size == nbins

    assert rmf.n_grp.size == nbins
    assert rmf.n_grp.min() == 12
    assert rmf.n_grp.max() == 48
    assert rmf.n_grp.sum() == nsize1

    assert rmf.f_chan.size == nsize1
    assert rmf.f_chan.min() == 114
    assert rmf.f_chan.max() == 3417
    assert rmf.f_chan.sum() == 173869713

    assert rmf.n_chan.size == nsize1
    assert rmf.n_chan.min() == 1
    assert rmf.n_chan.max() == 281
    assert rmf.n_chan.sum() == nsize2

    assert rmf.matrix.size == nsize2
    assert rmf.matrix.min() == pytest.approx(1.0000043459967856e-07)
    assert rmf.matrix.max() == pytest.approx(5.669764518737793)
    assert rmf.matrix.sum() == pytest.approx(102430.35533461263)
    assert np.argmax(rmf.matrix) == 5859948

    for field in ['e_min', 'e_max']:
        assert getattr(rmf, field).size == nchans

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1. Note that we check that the low
    # bin has been replaced (i.e. is not 0.0).
    #
    assert rmf.energ_lo[1:] == pytest.approx(rmf.energ_hi[:-1])

    assert rmf.energ_lo[0] == pytest.approx(0.30000001)
    assert rmf.energ_hi[0] == pytest.approx(0.30006698)
    assert rmf.energ_lo[-1] == pytest.approx(2.79407907)
    assert rmf.energ_hi[-1] == pytest.approx(2.79989982)

    # This is in reverse!
    assert rmf.e_min[:-1] == pytest.approx(rmf.e_max[1:])
    assert (rmf.e_min < rmf.e_max).all()

    assert rmf.e_min[0] == pytest.approx(3.09187627)
    assert rmf.e_max[0] == pytest.approx(3.09960628)
    assert rmf.e_min[-1] == pytest.approx(0.3099606)
    assert rmf.e_max[-1] == pytest.approx(0.31003815)


def check_warning(ftype, fname, w):

    assert w.category == UserWarning

    # Can use startswith/endswith but if there is an error then it is
    # harder to see the difference.
    #
    msg = str(w.message)
    assert msg[:33] == f"The minimum ENERG_LO in the {ftype} '"
    epos = len(fname)
    assert msg[-(39 + epos):] == f"/{fname}' was 0 and has been replaced by 1e-10"


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
def test_read_pha(errors, make_data_path,caplog):
    """Can we read in a XMM PHA file."""

    infile = make_data_path(PHAFILE)
    with warnings.catch_warnings(record=True) as ws:
        pha = io.read_pha(infile, use_errors=errors)

    check_pha(pha, responses=True)

    assert len(ws) == 2
    check_warning("ARF", ARFFILE, ws[0])
    check_warning("RMF", RMFFILE, ws[1])

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
    """Can we write out a XMM PHA file and then read it back in?"""

    infile = make_data_path(PHAFILE)
    with SherpaVerbosity("ERROR"):
        with warnings.catch_warnings(record=True) as ws:
            pha1 = io.read_pha(infile, use_errors=errors)

    assert len(ws) == 2

    outpath = tmp_path / "test.pha"
    outfile = str(outpath)
    io.write_pha(outfile, pha1, ascii=False, clobber=True)

    pha2 = io.read_pha(outfile, use_errors=errors)
    with SherpaVerbosity("INFO"):
        with warnings.catch_warnings(record=True) as ws2:
            check_pha(pha2, responses=False)

    assert len(ws2) == 0

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
    """Can we read in a XMM ARF."""

    infile = make_data_path(ARFFILE)
    with warnings.catch_warnings(record=True) as ws:
        arf = io.read_arf(infile)

    check_arf(arf)

    assert len(ws) == 1
    check_warning("ARF", ARFFILE, ws[0])


@requires_data
@requires_fits
def test_roundtrip_arf(make_data_path, tmp_path):
    """Can we write out a XMM ARF file and then read it back in?"""

    infile = make_data_path(ARFFILE)
    with warnings.catch_warnings(record=True) as ws:
        arf1 = io.read_arf(infile)

    assert len(ws) == 1

    outpath = tmp_path / "test.arf"
    outfile = str(outpath)
    io.write_arf(outfile, arf1, ascii=False, clobber=True)

    with warnings.catch_warnings(record=True) as ws2:
        arf2 = io.read_arf(outfile)

    check_arf(arf2)

    assert len(ws2) == 0


@requires_data
@requires_fits
def test_read_rmf(make_data_path):
    """Can we read in a XMM RMF."""

    infile = make_data_path(RMFFILE)
    with warnings.catch_warnings(record=True) as ws:
        rmf = io.read_rmf(infile)

    check_rmf(rmf)

    assert len(ws) == 1
    check_warning("RMF", RMFFILE, ws[0])


@requires_data
@requires_fits
def test_roundtrip_rmf(make_data_path, tmp_path):
    """Can we write out a XMM RMF file and then read it back in?"""

    infile = make_data_path(RMFFILE)
    with warnings.catch_warnings(record=True) as ws:
        rmf1 = io.read_rmf(infile)

    assert len(ws) == 1

    outpath = tmp_path / "test.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, rmf1, clobber=True)

    with warnings.catch_warnings(record=True) as ws2:
        rmf2 = io.read_rmf(outfile)

    check_rmf(rmf2)

    assert len(ws2) == 1


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_read_pha_grating(errors, make_data_path, caplog):
    """Can we read in a XMM grating PHA file."""

    infile = make_data_path(PHAFILE_GRATING)
    with SherpaVerbosity("INFO"):
        pha = io.read_pha(infile, use_errors=errors)

    check_pha_grating(pha, errors=errors)

    if errors:
        assert len(caplog.records) == 0
        return

    assert len(caplog.records) == 2

    record_starts_with(caplog.records[0],
                       "systematic errors were not found in file '",
                       level="WARNING")
    record_starts_with(caplog.records[1],
                       "statistical errors were found in file '")


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_roundtrip_pha_grating(errors, make_data_path, tmp_path, caplog):
    """Can we write out a XMM grating PHA file and then read it back in?"""

    infile = make_data_path(PHAFILE_GRATING)
    with SherpaVerbosity("ERROR"):
        pha1 = io.read_pha(infile, use_errors=errors)

    assert len(caplog.records) == 0

    outpath = tmp_path / "test.pha"
    outfile = str(outpath)
    io.write_pha(outfile, pha1, ascii=False, clobber=True)

    pha2 = io.read_pha(outfile, use_errors=errors)
    check_pha_grating(pha2, errors=errors)

    # Unlike the original file there is no output when errors=False
    #
    assert len(caplog.records) == 0


@requires_data
@requires_fits
def test_read_rmf_grating(make_data_path):
    """Can we read in a XMM grating RMF."""

    infile = make_data_path(RMFFILE_GRATING)
    rmf = io.read_rmf(infile)
    check_rmf_grating(rmf)


@requires_data
@requires_fits
def test_roundtrip_rmf(make_data_path, tmp_path):
    """Can we write out a XMM grating RMF file and then read it back in?"""

    infile = make_data_path(RMFFILE_GRATING)
    rmf1 = io.read_rmf(infile)

    outpath = tmp_path / "test.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, rmf1, clobber=True)

    rmf2 = io.read_rmf(outfile)
    check_rmf_grating(rmf2)


@requires_xspec
@requires_data
@requires_fits
def test_can_use_pha_grating(make_data_path, clean_astro_ui):
    """A basic check that we can use the XMM grating data.

    Unlike the previous tests, that directly access the io module,
    this uses the ui interface.
    """

    ui.set_stat("chi2xspecvar")

    ui.load_pha(make_data_path(PHAFILE_GRATING), use_errors=True)
    assert ui.get_analysis() == "channel"

    ui.load_rmf(make_data_path(RMFFILE_GRATING))
    assert ui.get_analysis() == "energy"

    d = ui.get_data()
    assert d.get_filter(format='%.4f') == "0.3100:3.0919"

    # XSPEC 12.14.1 has the energy ranges for a plot
    #
    # 2.40980291
    # 2.405128
    # ...
    # 0.324863762
    # 0.324778676
    #
    # and there are 2859 elements.
    #
    # The EBOUNDS block has
    #
    #   #  CHANNEL    E_MIN                E_MAX
    #             1         3.0918762684         3.0996062756
    #             2         3.0841853619         3.0918762684
    #             3         3.0765318871         3.0841853619
    #   ...
    #          3598     0.31011566519737     0.31019327044487
    #          3599     0.31003814935684     0.31011566519737
    #          3600     0.30996060371399     0.31003814935684
    #
    dp = ui.get_data_plot()
    assert len(dp.x) == 3600
    assert dp.x[0] == pytest.approx((3.0918762684 + 3.0996062756) / 2)
    assert dp.x[1] == pytest.approx((3.0841853619 + 3.0918762684) / 2)
    assert dp.x[-2] == pytest.approx((0.31003814935684 + 0.31011566519737) / 2)
    assert dp.x[-1] == pytest.approx((0.30996060371399 + 0.31003814935684) / 2)

    # If we ignore_bad we should match XSPEC.
    #
    ui.ignore_bad()

    dp = ui.get_data_plot()
    assert len(dp.x) == 2859
    assert dp.x[0] == pytest.approx(2.40980291)
    assert dp.x[1] == pytest.approx(2.405128)
    assert dp.x[-2] == pytest.approx(0.324863762)
    assert dp.x[-1] == pytest.approx(0.324778676)

    # XSPEC 12.14.1 gives
    #
    #   XSPEC12>data xmmrgs/P0112880201R1S004SRSPEC1003.FTZ
    #
    #   1 spectrum  in use
    #
    #   Spectral Data File: xmmrgs/P0112880201R1S004SRSPEC1003.FTZ  Spectrum 1
    #   Net count rate (cts/s) for Spectrum:1  9.342e-02 +/- 2.623e-03
    #    Assigned to Data Group 1 and Plot Group 1
    #     Noticed Channels:  1-2859
    #     Telescope: XMM Instrument: RGS1  Channel Type: PI
    #     Exposure Time: 2.897e+04 sec
    #    Using fit statistic: chi
    #    No response loaded.
    #
    #   ***Warning!  One or more spectra are missing responses,
    #                  and are not suitable for fit.
    #   XSPEC12>response xmmrgs/P0112880201R1S004RSPMAT1003.FTZ
    #   Response successfully loaded.
    #   XSPEC12>ignore **-0.64
    #     1792 channels (1068-2859) ignored in spectrum #     1
    #
    #   XSPEC12>ignore 0.68-**
    #      957 channels (1-957) ignored in spectrum #     1
    #
    #   XSPEC12>mo gauss
    #
    #   Input parameter value, delta, min, bot, top, and max values for ...
    #
    #   Input parameter value, delta, min, bot, top, and max values for ...
    #               6.5       0.05(     0.065)          0          0      1e+06      1e+06
    #   1:gaussian:LineE>0.6530
    #               0.1       0.05(     0.001)          0          0         10         20
    #   2:gaussian:Sigma>1e-4,-1
    #                 1       0.01(      0.01)          0          0      1e+20      1e+24
    #   3:gaussian:norm>3.63e-4
    #
    #   ========================================================================
    #   Model gaussian<1> Source No.: 1   Active/On
    #   Model Model Component  Parameter  Unit     Value
    #    par  comp
    #      1    1   gaussian   LineE      keV      0.653000     +/-  0.0
    #      2    1   gaussian   Sigma      keV      1.00000E-04  frozen
    #      3    1   gaussian   norm                3.63000E-04  +/-  0.0
    #   ________________________________________________________________________
    #
    #
    #   Fit statistic  : Chi-Squared                  118.23     using 110 bins.
    #
    #   Test statistic : Chi-Squared                  118.23     using 110 bins.
    #    Null hypothesis probability of 2.36e-01 with 108 degrees of freedom
    #    Current data and model not fit yet.
    #

    mdl = ui.create_model_component("xsgaussian", "gline")
    mdl.linee = 0.6530
    mdl.sigma.set(1e-4, frozen=True)
    mdl.norm = 3.63e-4

    ui.set_source("gline")

    # This energy range contains two "bad-channel" ranges
    # at ~ 0.655 keV and ~ 0.675 keV.
    #
    ui.ignore(hi=0.64)
    ui.ignore(lo=0.68)

    stats = ui.get_stat_info()
    assert len(stats) == 1
    stat = stats[0]

    assert stat.numpoints == 110
    assert stat.dof == 108

    assert stat.statval == pytest.approx(118.23, rel=0, abs=0.005)
