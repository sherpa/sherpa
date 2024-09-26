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

"""Test handling of ASCA data.

This has the advantage of having, in our test data

- CHANNEL=0 and CHANNEL=1
- different DETCHAN values

"""

import numpy as np

import pytest

from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.astro import io
from sherpa.astro import ui
from sherpa.utils.err import IdentifierErr, IOErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, requires_fits, requires_xspec


# We have an ASCA SIS file which uses offset=0 and another with
# offset=1, so check both.
#
PHA0FILE = "s0_mar24_bin.pha"
ARF0FILE = "s0_mar24.arf"
RMF0FILE = "s0_mar24.rmf"

PHA1FILE = "sis0.pha"
ARF1FILE = "sis0.arf"
RMF1FILE = "sis0.rmf"


def check_pha0(pha, errors=True, responses=True):
    """Does the data match expectation?"""

    assert isinstance(pha, DataPHA)

    assert pha.exposure == pytest.approx(33483.25)
    assert pha.backscal == pytest.approx(4.4189453125e-2)
    assert pha.areascal == pytest.approx(1.0)

    nchan = 512
    assert pha.channel == pytest.approx(np.arange(0, nchan))

    assert len(pha.counts) == nchan
    assert pha.counts.min() == pytest.approx(0.0)
    assert pha.counts.max() == pytest.approx(35.0)
    assert pha.counts.sum() == pytest.approx(1688.0)
    assert np.argmax(pha.counts) == 36

    for field in ['staterror', 'syserror', 'bin_lo', 'bin_hi']:
        assert getattr(pha, field) is None

    assert len(pha.grouping) == nchan
    assert len(pha.quality) == nchan

    assert sum(pha.grouping == 1) == 79
    assert sum(pha.grouping == -1) == (nchan - 79)

    assert sum(pha.quality == 0) == 496
    assert pha.quality[0] == 0
    assert pha.quality[17] == 0
    vals = pha.quality[np.arange(1, 17)]
    assert vals[0] == 5
    assert np.ptp(vals) == 0

    assert pha.grouped
    assert pha.subtracted is False
    assert pha.rate
    assert pha.plot_fac == 0
    assert pha.background_ids == []
    assert pha.get_background() is None

    if responses:
        assert pha.units == 'energy'
        assert pha.response_ids == [1]

        barf, brmf = pha.get_response()
        assert barf.name.endswith("/s0_mar24.arf")
        assert brmf.name.endswith("/s0_mar24.rmf")

    else:
        assert pha.units == 'channel'
        assert pha.response_ids == []

    # Select a few header keywords. Some are OGIP related,
    # some are "we expect them", and some are there just to
    # check we encode the type / value correctly.
    #
    assert pha.header["HDUCLASS"] == "OGIP"
    assert pha.header["HDUCLAS1"] == "SPECTRUM"
    assert pha.header["HDUCLAS2"] == ""
    assert pha.header["HDUCLAS3"] == "COUNT"
    assert pha.header["DETCHANS"] == nchan
    assert pha.header["CHANTYPE"] == "PI"

    # This keyword is removed from the header.
    assert "POISERR" not in pha.header

    assert pha.header["TELESCOP"] == "ASCA"
    assert pha.header["INSTRUME"] == "SIS0"
    assert pha.header["FILTER"] == "NONE"
    assert pha.header["OBJECT"] == "HE_1104-1805"


def check_pha1(pha, errors=True, responses=True):
    """Does the data match expectation?

    For this file neither option does anything,
    """

    assert isinstance(pha, DataPHA)

    assert pha.exposure == pytest.approx(28677.96586243063)
    assert pha.backscal == pytest.approx(2.550781332e-2)
    assert pha.areascal == pytest.approx(1.0)

    nchan = 1024
    assert pha.channel == pytest.approx(np.arange(1, nchan + 1))

    assert len(pha.counts) == nchan
    assert pha.counts.min() == pytest.approx(0.0)
    assert pha.counts.max() == pytest.approx(487.0)
    assert pha.counts.sum() == pytest.approx(50389.0)
    assert np.argmax(pha.counts) == 85

    for field in ['staterror', 'syserror',
                  'bin_lo', 'bin_hi',
                  'grouping', 'quality'
                  ]:
        assert getattr(pha, field) is None

    assert pha.grouped is False
    assert pha.subtracted is False
    assert pha.units == 'channel'
    assert pha.rate
    assert pha.plot_fac == 0
    assert pha.background_ids == []
    assert pha.response_ids == []
    assert pha.get_background() is None

    # Select a few header keywords. Some are OGIP related,
    # some are "we expect them", and some are there just to
    # check we encode the type / value correctly.
    #
    assert pha.header["HDUCLASS"] == "OGIP"
    assert pha.header["HDUCLAS1"] == "SPECTRUM"
    assert pha.header["HDUCLAS2"] == ""
    assert pha.header["HDUCLAS3"] == "COUNT"
    assert pha.header["DETCHANS"] == nchan
    assert pha.header["CHANTYPE"] == "PI"

    # This keyword is removed from the header.
    assert "POISERR" not in pha.header

    assert pha.header["TELESCOP"] == "ASCA"
    assert pha.header["INSTRUME"] == "SIS0"
    assert pha.header["FILTER"] == "none"
    assert pha.header["OBJECT"] == "A2029"


def check_arf0(arf):
    """Does the ARF contain the expected values?"""

    assert isinstance(arf, DataARF)

    nbins = 1180

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert len(arf.energ_lo) == nbins
    assert (arf.energ_lo[1:] == arf.energ_hi[:-1]).all()
    assert arf.energ_lo[0] == pytest.approx(0.20000000298023)
    assert arf.energ_hi[0] == pytest.approx(0.21000000834465)
    assert arf.energ_lo[-1] == pytest.approx(11.9901790619)
    assert arf.energ_hi[-1] == pytest.approx(12.0)

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(arf.specresp) == nbins
    assert arf.specresp.min() == pytest.approx(5.23834753036499)
    assert arf.specresp.max() == pytest.approx(221.3942108154297)
    assert arf.specresp.sum() == pytest.approx(115339.06293487549)
    assert np.argmax(arf.specresp) == 193

    for field in ['bin_lo', 'bin_hi', 'exposure']:
        assert getattr(arf, field) is None


def check_arf1(arf):
    """Does the ARF contain the expected values?"""

    assert isinstance(arf, DataARF)

    nbins = 1180

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert len(arf.energ_lo) == nbins
    assert (arf.energ_lo[1:] == arf.energ_hi[:-1]).all()
    assert arf.energ_lo[0] == pytest.approx(0.20000000298023)
    assert arf.energ_hi[0] == pytest.approx(0.21000000834465)
    assert arf.energ_lo[-1] == pytest.approx(11.9901790619)
    assert arf.energ_hi[-1] == pytest.approx(12.0)

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(arf.specresp) == nbins
    assert arf.specresp.min() == pytest.approx(5.3960146904)
    assert arf.specresp.max() == pytest.approx(235.13063049)
    assert arf.specresp.sum() == pytest.approx(125779.58236)
    assert np.argmax(arf.specresp) == 194

    for field in ['bin_lo', 'bin_hi', 'exposure']:
        assert getattr(arf, field) is None


def check_rmf0(rmf):
    """Does the RMF contain the expected values?"""

    assert isinstance(rmf, DataRMF)

    nchans = 512
    nbins = 1180
    nsize1 = 1189
    nsize2 = 260775

    assert rmf.detchans == nchans
    assert rmf.offset == 0

    for field in ['energ_lo', 'energ_hi']:
        assert getattr(rmf, field).size == nbins

    assert rmf.n_grp.size == nbins
    assert rmf.n_grp.min() == 0
    assert rmf.n_grp.max() == 2
    assert rmf.n_grp.sum() == nsize1

    assert rmf.f_chan.size == nsize1
    assert rmf.f_chan.min() == 14
    assert rmf.f_chan.max() == 19
    assert rmf.f_chan.sum() == 16711

    assert rmf.n_chan.size == nsize1
    assert rmf.n_chan.min() == 1
    assert rmf.n_chan.max() == 426
    assert rmf.n_chan.sum() == nsize2

    assert rmf.matrix.size == nsize2
    assert rmf.matrix.min() == pytest.approx(1.0003149952808599e-07)
    assert rmf.matrix.max() == pytest.approx(0.11968599259853363)
    assert rmf.matrix.sum() == pytest.approx(544.8837793875416)
    assert np.argmax(rmf.matrix) == 4307

    for field in ['e_min', 'e_max']:
        assert getattr(rmf, field).size == nchans

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1. Note that we check that the low
    # bin has been replaced (i.e. is not 0.0).
    #
    assert rmf.energ_lo[1:] == pytest.approx(rmf.energ_hi[:-1])
    assert rmf.energ_lo[0] == pytest.approx(0.20000000298023)
    assert rmf.energ_hi[0] == pytest.approx(0.21000000834465)
    assert rmf.energ_lo[-1] == pytest.approx(11.9901790619)
    assert rmf.energ_hi[-1] == pytest.approx(12.0)

    assert rmf.e_min[1:] == pytest.approx(rmf.e_max[:-1])
    assert rmf.e_min[0] == pytest.approx(0.013292968273162842)
    assert rmf.e_max[0] == pytest.approx(0.04193559288978577)
    assert rmf.e_min[-1] == pytest.approx(14.649674415588379)
    assert rmf.e_max[-1] == pytest.approx(14.678317070007324)


def check_rmf1(rmf):

    assert isinstance(rmf, DataRMF)

    nchans = 1024
    nbins = 1180
    nsize1 = 1171
    nsize2 = 505483

    assert rmf.detchans == nchans
    assert rmf.offset == 1

    for field in ['energ_lo', 'energ_hi']:
        assert getattr(rmf, field).size == nbins

    assert rmf.n_grp.size == nbins
    assert rmf.n_grp.min() == 0
    assert rmf.n_grp.max() == 1
    assert rmf.n_grp.sum() == nsize1

    assert rmf.f_chan.size == nsize1
    assert rmf.f_chan.min() == 26
    assert rmf.f_chan.max() == 26
    assert rmf.f_chan.sum() == 30446

    assert rmf.n_chan.size == nsize1
    assert rmf.n_chan.min() == 1
    assert rmf.n_chan.max() == 851
    assert rmf.n_chan.sum() == nsize2

    assert rmf.matrix.size == nsize2
    # TODO: asking for the min value causes this test to hang on my
    # machine. I can read in the file outside of the test and get
    # the minimum value. Why is it just here?
    #
    # assert rmf.matrix.min() == pytest.approx(1.0003149952808599e-7)
    assert (rmf.matrix > 1e-7).all()  # approximate the above
    assert rmf.matrix.max() == pytest.approx(0.11497905850410461)
    assert rmf.matrix.sum() == pytest.approx(552.932040504181)
    assert np.argmax(rmf.matrix) == 7324

    for field in ['e_min', 'e_max']:
        assert getattr(rmf, field).size == nchans

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1. Note that we check that the low
    # bin has been replaced (i.e. is not 0.0).
    #
    assert rmf.energ_lo[1:] == pytest.approx(rmf.energ_hi[:-1])
    assert rmf.energ_lo[0] == pytest.approx(0.20000000298023)
    assert rmf.energ_hi[0] == pytest.approx(0.21000000834465)
    assert rmf.energ_lo[-1] == pytest.approx(11.9901790619)
    assert rmf.energ_hi[-1] == pytest.approx(12.0)

    assert rmf.e_min[1:] == pytest.approx(rmf.e_max[:-1])
    assert rmf.e_min[0] == pytest.approx(0.0)
    assert rmf.e_max[0] == pytest.approx(0.013666952028870583)
    assert rmf.e_min[-1] == pytest.approx(14.228230476379395)
    assert rmf.e_max[-1] == pytest.approx(14.240230560302734)


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


# Since the log messages we get for the PHA0 and PHA1 cases are
# different, it's easiest just to repeat the tests rather than to try
# and write a single test.
#
@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_read_pha0(errors, make_data_path, caplog):
    """Can we read in an ASCS PHA file."""

    with SherpaVerbosity("INFO"):
        infile = make_data_path(PHA0FILE)

    pha = io.read_pha(infile, use_errors=errors)
    check_pha0(pha, errors=errors, responses=True)

    # What messages were created?
    #
    assert len(caplog.records) == 3

    rec = caplog.records[0]
    record_starts_with(rec, "read ARF file /")

    rec = caplog.records[1]
    record_starts_with(rec, "read RMF file /")

    # The end of the message comes from the I/O backend and we do not
    # bother testing that here.
    #
    rec = caplog.records[2]
    record_starts_with(rec, "unable to read background: ",
                       level="WARNING")


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_read_pha1(errors, make_data_path, caplog):
    """Can we read in an ASCS PHA file."""

    with SherpaVerbosity("INFO"):
        infile = make_data_path(PHA1FILE)

    pha = io.read_pha(infile, use_errors=errors)
    check_pha1(pha, errors=errors, responses=True)

    assert len(caplog.records) == 0


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_roundtrip_pha0(errors, make_data_path, tmp_path, caplog):
    """Can we write out an ASCA PHA file and then read it back in?"""

    infile = make_data_path(PHA0FILE)
    with SherpaVerbosity("ERROR"):
        pha1 = io.read_pha(infile, use_errors=errors)

    assert len(caplog.records) == 0

    outpath = tmp_path / "test.pha"
    outfile = str(outpath)
    io.write_pha(outfile, pha1, ascii=False, clobber=True)

    with SherpaVerbosity("INFO"):
        pha2 = io.read_pha(outfile, use_errors=errors)

    check_pha0(pha2, responses=False, errors=errors)

    assert len(caplog.records) == 2

    # The end of the message comes from the I/O backend and we do not
    # bother testing that here.
    #
    rec = caplog.records[0]
    record_starts_with(rec, "unable to read ARF: ",
                       level="WARNING")

    rec = caplog.records[1]
    record_starts_with(rec, "unable to read RMF: ",
                       level="WARNING")

    # TODO: why does this no-longer mention being unable to read in
    # the background?


@requires_data
@requires_fits
@pytest.mark.parametrize("errors", [False, True])
def test_roundtrip_pha1(errors, make_data_path, tmp_path, caplog):
    """Can we write out an ASCA PHA file and then read it back in?"""

    infile = make_data_path(PHA1FILE)
    with SherpaVerbosity("ERROR"):
        pha1 = io.read_pha(infile, use_errors=errors)

    assert len(caplog.records) == 0

    outpath = tmp_path / "test.pha"
    outfile = str(outpath)
    io.write_pha(outfile, pha1, ascii=False, clobber=True)

    with SherpaVerbosity("INFO"):
        pha2 = io.read_pha(outfile, use_errors=errors)

    check_pha1(pha2, responses=False, errors=errors)


@requires_data
@requires_fits
@pytest.mark.parametrize("filename,check",
                         [(ARF0FILE, check_arf0),
                          (ARF1FILE, check_arf1)])
def test_read_arf(filename, check, make_data_path):
    """Can we read in an ASCA ARF."""

    infile = make_data_path(filename)
    arf = io.read_arf(infile)
    check(arf)


@requires_data
@requires_fits
@pytest.mark.parametrize("filename,check",
                         [(ARF0FILE, check_arf0),
                          (ARF1FILE, check_arf1)])
def test_roundtrip_arf(filename, check, make_data_path, tmp_path):
    """Can we write out an ASCA ARF file and then read it back in?"""

    infile = make_data_path(filename)
    arf1 = io.read_arf(infile)

    outpath = tmp_path / "test.arf"
    outfile = str(outpath)
    io.write_arf(outfile, arf1, ascii=False, clobber=True)

    arf2 = io.read_arf(outfile)
    check(arf2)


def validate_record(r):
    """Check it represents the expected log output."""

    assert r.name == "sherpa.astro.io"
    assert r.levelname == "ERROR"
    msg = r.getMessage()

    # Can use startswith/endswith, but then the error message when the
    # test fails is less useful. This way pytest reports the location
    # of the difference.
    #
    assert msg[:62] == "Failed to locate TLMIN keyword for F_CHAN column in RMF file '"
    assert msg[-102:] == "sis0.rmf'; Update the offset value in the RMF data set to the appropriate TLMIN value prior to fitting"


@requires_data
@requires_fits
@pytest.mark.parametrize("filename,check,tlmin",
                         [(RMF0FILE, check_rmf0, False),
                          (RMF1FILE, check_rmf1, True)])
def test_read_rmf(filename, check, tlmin, make_data_path, caplog):
    """Can we read in an ASCA RMF."""

    infile = make_data_path(filename)
    with SherpaVerbosity("INFO"):
        rmf = io.read_rmf(infile)

    if tlmin:
        assert len(caplog.records) == 1
        validate_record(caplog.records[0])
    else:
        assert len(caplog.records) == 0

    check(rmf)


@requires_data
@requires_fits
@pytest.mark.parametrize("filename,check,tlmin",
                         [(RMF0FILE, check_rmf0, False),
                          (RMF1FILE, check_rmf1, True)])
def test_roundtrip_rmf(filename, check, tlmin, make_data_path, tmp_path, caplog):
    """Can we write out an ASCA RMF file and then read it back in?"""

    infile = make_data_path(filename)
    with SherpaVerbosity("INFO"):
        rmf1 = io.read_rmf(infile)

    if tlmin:
        assert len(caplog.records) == 1
        validate_record(caplog.records[0])
        ntot = 1
    else:
        assert len(caplog.records) == 0
        ntot = 0

    outpath = tmp_path / "test.rmf"
    outfile = str(outpath)
    with SherpaVerbosity("INFO"):
        io.write_rmf(outfile, rmf1, clobber=True)

    # Note that the TLMIN value is written out so there is no more
    # message about the value being missing.
    #
    with SherpaVerbosity("INFO"):
        rmf2 = io.read_rmf(outfile)

    check(rmf2)
    assert len(caplog.records) == ntot


@requires_xspec
@requires_data
@requires_fits
def test_can_use_pha0(make_data_path, clean_astro_ui):
    """A basic check that we can use the ASCA data.

    Unlike the previous tests, that directly access the io module,
    this uses the ui interface.
    """

    ui.set_stat("chi2xspecvar")

    # Load the errors to match XSPEC.
    #
    ui.load_pha(make_data_path(PHA0FILE), use_errors=True)
    assert ui.get_analysis() == 'energy'

    assert ui.get_arf().name.endswith("/s0_mar24.arf")
    assert ui.get_rmf().name.endswith("/s0_mar24.rmf")
    with pytest.raises(IdentifierErr):
        assert ui.get_bkg()

    # As noted below, this drops the first bin that XSPEC uses. XSPEC
    # reports mid-energy values for its first 5 groups as:
    #
    #    2.76142806E-2
    #    0.27107656
    #    0.514538884
    #    0.557502866
    #    0.614788055
    #
    # In CIAO 4.16 we have "collapsed" the first three groups as our
    # first couple of bins are
    #
    #        mid     |     lo      |      hi
    #    0.27107659    0.01329297     0.52886021
    #    0.55750284    0.52886021     0.58614546
    #    0.61478809    0.58614546     0.64343071
    #
    # This is because some of the channels are labelled "bad"
    #
    # % dmlist s0_mar24_bin.pha"[cols -sys_err][#row=:24]" data,clean
    # #  CHANNEL COUNTS     QUALITY    GROUPING
    #       0          0          0          1
    #       1          0          5         -1
    #       2          0          5         -1
    #       3          0          5         -1
    #       4          0          5         -1
    #       5          0          5         -1
    #       6          0          5         -1
    #       7          0          5         -1
    #       8          0          5         -1
    #       9          0          5         -1
    #      10          0          5         -1
    #      11          0          5         -1
    #      12          0          5         -1
    #      13          0          5         -1
    #      14         13          5         -1
    #      15         15          5         -1
    #      16          8          5         -1
    #      17         17          0         -1
    #      18         12          0          1
    #      19          9          0         -1
    #      20          7          0          1
    #      21         12          0         -1
    #      22         18          0          1
    #      23         15          0          1
    #
    # So, does XSPEC ignore the QUALITY=5 rows so it treats the
    # first group as channels 0 and 17,
    # second group as channels 18 and 19 ?
    #
    ui.ignore(hi=0.5)
    ui.ignore(lo=7.0)

    ui.set_source(ui.xsphabs.gal * ui.powlaw1d.pl)
    ui.set_par('gal.nh', 0.05)
    ui.set_par('pl.gamma', 1.7)
    ui.set_par('pl.ampl', 2.7e-4)

    # Value obtained from XSPEC 12.14.1
    #
    #    XSPEC12> data s0_mar24_bin.pha
    #    - complains about FILTER none/NONE mis-match
    #    - complains about missing background
    #    XSPEC12> setplot energy
    #    XSPEC12> ignore **-0.5,7.0-**
    #    - reports ignores "channels" 1-2 and 74-81
    #    XSPEC12> mo phabs * powerlaw
    #    XSPEC12> newpar 1 0.05
    #    XSPEC12> newpar 2 1.7
    #    XSPEC12> newpar 3 2.7e-4
    #    - reports a statistic of 69.34
    #

    # XSPEC reports 71, we get 70 (CIAO 4.16)
    #
    NBINS = 70

    # XSPEC reports 69.34
    #
    STATVAL = 65.06460092188946

    s = ui.get_stat_info()[0]
    assert s.numpoints == NBINS
    assert s.dof == (NBINS - 3)

    assert s.statval == pytest.approx(STATVAL, rel=0, abs=0.005)

    # The plot data starts
    #
    #   0.514538884 1.43213272E-2 1.77258905E-2 4.29915963E-3 8.83223955E-3
    #   0.557502866 2.86426246E-2 1.09483553E-2 2.38912692E-3 7.78075447E-3
    #
    # and ends with
    #
    #  5.79910278 0.229140997 9.77531774E-4 2.52397615E-4 8.68383213E-4
    #  6.35763454 0.329390287 6.80021825E-4 1.75580892E-4 5.45510848E-4
    #
    fp = ui.get_fit_plot()
    dp = fp.dataplot
    mp = fp.modelplot

    assert len(dp.x) == NBINS

    assert dp.xlo[0] == pytest.approx(0.52886021)
    assert dp.xhi[0] == pytest.approx(0.58614546)
    assert dp.x[0] == pytest.approx(0.55750284)

    assert dp.xlo[-1] == pytest.approx(6.02824402)
    assert dp.xhi[-1] == pytest.approx(6.68702459)
    assert dp.x[-1] == pytest.approx(6.35763431)

    # Check the model values, as they get convolved through the
    # response, and then filtered and grouped (since using the
    # modelplot from a fit rather than from get_model_plot, which is
    # ungrouped).
    #
    assert len(mp.y) == NBINS

    # Values from XSPEC:
    #
    assert mp.y[0] == pytest.approx(7.78075447E-3)
    assert mp.y[-1] == pytest.approx(5.45510848e-4)


@requires_xspec
@requires_data
@requires_fits
def test_can_use_pha1(make_data_path, clean_astro_ui):
    """A basic check that we can use the ASCA data.

    Unlike the previous tests, that directly access the io module,
    this uses the ui interface.
    """

    ui.set_stat("chi2xspecvar")

    # Load the errors to match XSPEC.
    #
    ui.load_pha(make_data_path(PHA1FILE), use_errors=True)
    assert ui.get_analysis() == 'channel'

    ui.load_rmf(make_data_path("sis0.rmf"))
    ui.load_arf(make_data_path("sis0.arf"))

    assert ui.get_analysis() == 'energy'

    with pytest.raises(IdentifierErr):
        assert ui.get_bkg()

    ui.ignore(hi=0.4)
    ui.ignore(lo=7.0)

    # This is not a great fit. Use mekal rather than apec in the
    # assumption it isn't going to change significantly over time.
    #
    ui.set_source(ui.xsphabs.gal * ui.xsmekal.clus)
    ui.set_par('gal.nh', 0.09)
    ui.set_par('clus.kt', 5.9)
    ui.set_par('clus.abundanc', 0.3)
    ui.set_par('clus.redshift', 0.0767)
    ui.set_par('clus.norm', 4.65e-2)

    # Value obtained from XSPEC 12.14.1
    #
    #    XSPEC12> data s0_mar24_bin.pha
    #    - complains about FILTER none/NONE mis-match
    #    - complains about missing background
    #    XSPEC12> setplot energy
    #    XSPEC12> ignore **-0.5,7.0-**
    #    - reports ignores "channels" 1-2 and 74-81
    #    XSPEC12> mo phabs * mekal
    #    XSPEC12> newpar 1 9e-2
    #    XSPEC12> newpar 2 5.9
    #    XSPEC12> newpar 4 0.3
    #    XSPEC12> newpar 7 4.65e-2
    #    - reports a statistic of 544.23
    #

    NBINS = 450
    STATVAL = 544.23

    s = ui.get_stat_info()[0]
    assert s.numpoints == NBINS
    assert s.dof == (NBINS - 3)

    assert s.statval == pytest.approx(STATVAL, rel=0, abs=0.005)

    fp = ui.get_fit_plot()
    dp = fp.dataplot
    mp = fp.modelplot

    assert len(dp.x) == NBINS
    assert len(mp.y) == NBINS

    # Values taken from
    #
    #   XSPEC12> iplot
    #   PLT> wdata foo.txt
    #
    assert dp.x[0] == pytest.approx(0.416999876)
    assert dp.x[-1] == pytest.approx(6.98360157)

    assert mp.y[0] == pytest.approx(0.153414682)
    assert mp.y[-1] == pytest.approx(1.06004886e-2)
