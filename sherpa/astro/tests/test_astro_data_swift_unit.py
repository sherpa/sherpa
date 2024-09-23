#
#  Copyright (C) 2017, 2018, 2021, 2023, 2024
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

"""Test handling of Swift [1]_ data.

At the moment the tests are of routines in the sherpa.astro.io module,
since the high-level routines like sherpa.astro.ui.unpack_pha are a
very-thin wrapper around the io routines.

References
----------

.. [1] The Swift Gamma Ray-Burst Mission https://swift.gsfc.nasa.gov/

"""

import warnings

import numpy as np

import pytest

from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.astro import io
from sherpa.astro import ui
from sherpa.utils.err import IOErr
from sherpa.utils.testing import requires_data, requires_fits


def backend_is(name: str) -> bool:
    """Are we using the specified backend?"""
    return io.backend.name == name


# PHA and ARF are stored uncompressed while RMF is gzip-compressed.
# Use the same name (no .gz suffix) to access all files.
#
PHAFILE = 'target_sr.pha'
ARFFILE = 'target_sr.arf'
RMFFILE = 'swxpc0to12s6_20130101v014.rmf'

# The minimum energy value for ENERG_LO bins. This should match the
# value in the [ogip] minimum_energy field of the configuration file.
#
EMIN = 1.0e-10


def check_pha(pha):
    """Does the data match expectation?"""

    assert isinstance(pha, DataPHA)

    assert pha.exposure == pytest.approx(4812.26699895)
    assert pha.backscal == pytest.approx(1.148e-3)
    assert pha.areascal == pytest.approx(1.0)

    # Although channel is an integer, it's treated as a
    # float. Note that the channels start at 0 in the file,
    # but they are promoted to 1 by the generic backend.
    # XSPEC 12.9.1b, on reading this file, reports channel numbers
    # between 1 and 1024 inclusive.
    #
    nchan = 1024
    assert pha.channel == pytest.approx(np.arange(0, nchan))

    # Rather than check each element, use some simple summary
    # statistics.
    assert len(pha.counts) == nchan
    assert pha.counts.min() == pytest.approx(0.0)
    assert pha.counts.max() == pytest.approx(3.0)
    assert pha.counts.sum() == pytest.approx(58.0)
    assert np.argmax(pha.counts) == 110

    for field in ['staterror', 'syserror', 'bin_lo', 'bin_hi',
                  'grouping', 'quality']:
        assert getattr(pha, field) is None

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
    assert pha.header["HDUCLAS2"] == "TOTAL"
    assert pha.header["HDUCLAS3"] == "COUNT"
    assert pha.header["DETCHANS"] == 1024
    assert pha.header["CHANTYPE"] == "PI"

    # This keyword is removed from the header.
    assert "POISERR" not in pha.header

    assert pha.header["TELESCOP"] == "SWIFT"
    assert pha.header["INSTRUME"] == "XRT"
    assert pha.header["FILTER"] == "NONE"
    assert pha.header["OBJECT"] == "2M02594633+3855363"

    assert pha.header["RA_OBJ"] == pytest.approx(44.9445)
    assert pha.header["DEC_OBJ"] == pytest.approx(38.9263)
    assert pha.header["UTCFINIT"] == pytest.approx(-15.64118)

    assert not pha.header["XBADFRAM"]


def check_arf(arf):
    """Does the ARF contain the expected values?"""

    assert isinstance(arf, DataARF)

    nbins = 2400

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1. Note that we check that the low
    # bin has been replaced (i.e. is not 0.0).
    #
    assert len(arf.energ_lo) == nbins
    assert (arf.energ_lo[1:] == arf.energ_hi[:-1]).all()
    assert arf.energ_lo[0] == pytest.approx(EMIN)
    assert arf.energ_hi[0] == pytest.approx(0.005)
    assert arf.energ_lo[-1] == pytest.approx(11.995)
    assert arf.energ_hi[-1] == pytest.approx(12.0)

    # Rather than check each element, use some simple summary
    # statistics.
    #
    assert len(arf.specresp) == nbins
    assert arf.specresp.min() == pytest.approx(0.1506345272064209)
    assert arf.specresp.max() == pytest.approx(145.38259887695312)
    assert arf.specresp.sum() == pytest.approx(178696.1007297188)
    assert np.argmax(arf.specresp) == 309

    for field in ['bin_lo', 'bin_hi', 'exposure']:
        assert getattr(arf, field) is None


def check_rmf(rmf):
    """Does the RMF contain the expected values?"""

    assert isinstance(rmf, DataRMF)

    nchans = 1024
    nbins = 2400

    assert rmf.detchans == nchans
    assert rmf.offset == 0

    for field in ['energ_lo', 'energ_hi', 'n_grp', 'f_chan',
                  'n_chan']:
        assert getattr(rmf, field).size == nbins

    assert rmf.matrix.size == nchans * nbins

    for field in ['e_min', 'e_max']:
        assert getattr(rmf, field).size == nchans

    assert (rmf.n_grp == 1).all()
    assert (rmf.f_chan == 0).all()
    assert (rmf.n_chan == nchans).all()

    # Rather than check each element, use some simple summary
    # statistics.
    assert rmf.matrix.min() == pytest.approx(0.0)
    assert rmf.matrix.max() == pytest.approx(0.05834391713142395)
    assert rmf.matrix.sum() == pytest.approx(1133.8128197300257)
    assert np.argmax(rmf.matrix) == 269443

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1. Note that we check that the low
    # bin has been replaced (i.e. is not 0.0).
    #
    assert rmf.energ_lo[1:] == pytest.approx(rmf.energ_hi[:-1])
    assert rmf.e_min[1:] == pytest.approx(rmf.e_max[:-1])

    assert rmf.energ_lo[0] == pytest.approx(EMIN)
    assert rmf.energ_hi[0] == pytest.approx(0.005)
    assert rmf.energ_lo[-1] == pytest.approx(11.995)
    assert rmf.energ_hi[-1] == pytest.approx(12.0)

    assert rmf.e_min[0] == pytest.approx(0.0)
    assert rmf.e_max[0] == pytest.approx(0.01)
    assert rmf.e_min[-1] == pytest.approx(10.23)
    assert rmf.e_max[-1] == pytest.approx(10.24)


@requires_data
@requires_fits
def test_read_pha(make_data_path):
    """Can we read in a Swift PHA file."""

    infile = make_data_path(PHAFILE)
    pha = io.read_pha(infile)
    check_pha(pha)


@requires_data
@requires_fits
def test_roundtrip_pha(make_data_path, tmp_path):
    """Can we write out a Swift PHA file and then read it back in?"""

    infile = make_data_path(PHAFILE)
    pha1 = io.read_pha(infile)

    outpath = tmp_path / "test.pha"
    outfile = str(outpath)
    io.write_pha(outfile, pha1, ascii=False, clobber=True)

    pha2 = io.read_pha(outfile)
    check_pha(pha2)


# In the error checks (i.e. checking that invalid files fail),
# the error messages and types depend on the backend. An attempt
# is made to be specific when testing these, primarily to point
# out when any refactoring changes this (since in many ways it
# would be nice to have the same error messages between the
# backends when possible).
#
@requires_data
@requires_fits
def test_read_pha_fails_arf(make_data_path):
    """Just check in we can't read in an ARF as a PHA file."""

    if backend_is("pyfits"):
        emsg = " does not appear to be a PHA spectrum"
        etype = IOErr
    elif backend_is("crates"):
        emsg = "File must be a PHA file."
        etype = TypeError
    else:
        assert False, f"Internal error: unknown backend {io.backend}"

    infile = make_data_path(ARFFILE)
    with pytest.raises(etype) as excinfo:
        io.read_pha(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_pha_fails_rmf(make_data_path):
    """Just check in we can't read in a RMF as a PHA file."""

    if backend_is("pyfits"):
        emsg = " does not appear to be a PHA spectrum"
        etype = IOErr
    elif backend_is("crates"):
        emsg = "File must be a PHA file."
        etype = TypeError
    else:
        assert False, f"Internal error: unknown backend {io.backend}"

    infile = make_data_path(RMFFILE)
    with pytest.raises(etype) as excinfo:
        io.read_pha(infile)

    assert emsg in str(excinfo.value)


def validate_replacement_warning(ws, rtype, label):
    """Check that there is one expected warning."""

    assert len(ws) == 1
    w = ws[0]
    assert w.category == UserWarning

    emsg = f"The minimum ENERG_LO in the {rtype} '{label}' " + \
           f"was 0 and has been replaced by {EMIN}"
    assert str(w.message) == emsg


@requires_data
@requires_fits
def test_read_arf(make_data_path):
    """Can we read in a Swift ARF."""

    infile = make_data_path(ARFFILE)

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        arf = io.read_arf(infile)

    validate_replacement_warning(ws, 'ARF', infile)
    check_arf(arf)


@requires_data
@requires_fits
def test_roundtrip_arf(make_data_path, tmp_path):
    """Can we write out a Swift ARF file and then read it back in?"""

    infile = make_data_path(ARFFILE)

    # Skip the warning check here.
    #
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")
        arf1 = io.read_arf(infile)

    outpath = tmp_path / "test.arf"
    outfile = str(outpath)
    io.write_arf(outfile, arf1, ascii=False, clobber=True)

    # This should have no warnings as the ENERG_LO values were updated
    # when creating arf1.
    #
    arf2 = io.read_arf(outfile)
    check_arf(arf2)


@requires_data
@requires_fits
def test_read_arf_fails_pha(make_data_path):
    """Just check in we can't read in a PHA file as an ARF."""

    if backend_is("pyfits"):
        emsg = " does not appear to be an ARF"
    elif backend_is("crates"):
        emsg = "Required column 'ENERG_LO' not found in "
    else:
        assert False, f"Internal error: unknown backend {io.backend}"

    infile = make_data_path(PHAFILE)
    with pytest.raises(IOErr) as excinfo:
        io.read_arf(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_arf_fails_rmf(make_data_path):
    """Just check in we can't read in a RNF as an ARF."""

    if backend_is("pyfits"):
        emsg = " does not appear to be an ARF"
    elif backend_is("crates"):
        emsg = "Required column 'SPECRESP' not found in "
    else:
        assert False, f"Internal error: unknown backend {io.backend}"

    infile = make_data_path(RMFFILE)
    with pytest.raises(IOErr) as excinfo:
        io.read_arf(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_rmf(make_data_path):
    """Can we read in a Swift RMF."""

    infile = make_data_path(RMFFILE)

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        rmf = io.read_rmf(infile)

    validate_replacement_warning(ws, 'RMF', infile)
    check_rmf(rmf)


@requires_data
@requires_fits
def test_roundtrip_rmf(make_data_path, tmp_path):
    """Can we write out a Swift RMF file and then read it back in?"""

    infile = make_data_path(RMFFILE)

    # Skip the warning check here.
    #
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")
        rmf1 = io.read_rmf(infile)

    outpath = tmp_path / "test.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, rmf1, clobber=True)

    # This should have no warnings as the ENERG_LO values were updated
    # when creating rmf1.
    #
    rmf2 = io.read_rmf(outfile)
    check_rmf(rmf2)


@requires_data
@requires_fits
def test_read_rmf_fails_pha(make_data_path):
    """Just check in we can't read in a PHA file as a RMF."""

    if backend_is("pyfits"):
        emsg = " does not appear to be an RMF"
        etype = IOErr
    elif backend_is("crates"):
        emsg = " does not contain a Response Matrix."
        etype = TypeError
    else:
        assert False, f"Internal error: unknown backend {io.backend}"

    infile = make_data_path(PHAFILE)
    with pytest.raises(etype) as excinfo:
        io.read_rmf(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_rmf_fails_arf(make_data_path):
    """Just check in we can't read in a ARF as a RMF."""

    if backend_is("pyfits"):
        emsg = " does not appear to be an RMF"
        etype = IOErr
    elif backend_is("crates"):
        emsg = " does not contain a Response Matrix."
        etype = TypeError
    else:
        assert False, f"Internal error: unknown backend {io.backend}"

    infile = make_data_path(ARFFILE)
    with pytest.raises(etype, match=emsg):
        io.read_rmf(infile)


@requires_data
@requires_fits
def test_can_use_swift_data(make_data_path, clean_astro_ui):
    """A basic check that we can read in and use the Swift data.

    Unlike the previous tests, that directly access the io module,
    this uses the ui interface.
    """

    # The Swift PHA file does not have the ANCRFILE/RESPFILE keywords
    # set up, so the responses have to be manually added.
    #
    ui.load_pha(make_data_path(PHAFILE))

    rmffile = make_data_path(RMFFILE)
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        ui.load_rmf(rmffile)

    validate_replacement_warning(ws, 'RMF', rmffile)

    arffile = make_data_path(ARFFILE)
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        ui.load_arf(arffile)

    validate_replacement_warning(ws, 'ARF', arffile)

    assert ui.get_analysis() == 'energy'

    arf = ui.get_arf()
    rmf = ui.get_rmf()
    assert arf.energ_lo[0] == EMIN
    assert rmf.energ_lo[0] == EMIN
    assert rmf.e_min[0] == 0.0

    ui.set_source(ui.powlaw1d.pl)
    ui.set_par('pl.ampl', 0.0003)

    stat = ui.calc_stat()

    # This check is purely a regression test, so the value has
    # not been externally validated.
    #
    assert stat == pytest.approx(58.2813692358182)

    # Pick an energy range which isn't affected by the first
    # bin.
    #
    # Unfortunately, using a range of 0.3-8.0 gives 771 bins
    # in XSPEC - channels 30 to 800 - but 770 bins in Sherpa,
    # channels 31 to 800.
    #
    # Note that the channel numbering starts at 0:
    # % dmlist target_sr.pha header,clean,raw | grep TLMIN
    # TLMIN1       = 0                    / Lowest legal channel number
    #
    # and so it's not clear when XSPEC says 30-800 what it
    # means. From https://github.com/sherpa/sherpa/issues/1211#issuecomment-881647128
    # we have that the first bin it is using is
    #     0.29-0.30
    # and the last bin is
    #     7.99-8.00
    # and I've checked with iplot that it has renumbered the
    # channels to 1-1024 from 0-1023
    #
    # % dmlist swxpc0to12s6_20130101v014.rmf.gz"[ebounds][channel=28:31]" data,clean
    #  CHANNEL    E_MIN                E_MAX
    #         28     0.28000000119209     0.28999999165535
    #         29     0.28999999165535     0.30000001192093
    #         30     0.30000001192093     0.31000000238419
    #         31     0.31000000238419     0.31999999284744
    # % dmlist swxpc0to12s6_20130101v014.rmf.gz"[ebounds][channel=798:801]" data,clean
    #  CHANNEL    E_MIN                E_MAX
    #        798         7.9800000191         7.9899997711
    #        799         7.9899997711                  8.0
    #        800                  8.0         8.0100002289
    #        801         8.0100002289         8.0200004578
    #
    # If I use ignore(None, 0.3); ignore(8.0, None) instead then the
    # result is 771 bins (channels 31 to 800). This is because the
    # e_min/max of the RMF has channel widths of 0.01 keV, starting at
    # 0, so both 0.3 and 8.0 fall on a bin boundary. So, it's either a
    # difference in < or <= (or > vs >=), or a rounding issue due to
    # floating-point conversion leading to one bin boundary being
    # slightly different in Sherpa vs XSPEC).
    #
    # When using ui.notice(0.3, 8.0); ui.get_indep(filter=True)
    # returns 770 channels, 31 to 800.
    #
    # Using ui.notice(0.3, 7.995) selects channels 31 to 800.
    # Using ui.notice(0.299, 8.0) selects channels 30 to 800.
    # Using ui.notice(0.299, 7.995) selects channels 30 to 800.
    #
    ui.notice(0.299, 8.0)

    # Check the selected range
    pha = ui.get_data()
    expected = np.zeros(1024, dtype=bool)
    expected[29:800] = True
    assert pha.mask == pytest.approx(expected)
    assert pha.get_mask() == pytest.approx(expected)

    # XSPEC 12.9.1b calculation of the statistic:
    #   chi sq = 203.88 from 771 bins with 769 dof
    #   cstat  = 568.52
    #
    # Histotically there have been issues with chi2xspecvar,
    # such as https://github.com/sherpa/sherpa/issues/356
    # so check with this setting.
    #
    ui.set_stat('chi2xspecvar')
    stat_xvar = ui.get_stat_info()

    assert len(stat_xvar) == 1
    stat_xvar = stat_xvar[0]
    assert stat_xvar.numpoints == 771
    assert stat_xvar.dof == 769
    assert stat_xvar.statval == pytest.approx(203.88,
                                              rel=0, abs=0.005)

    ui.set_stat('cstat')
    stat_cstat = ui.get_stat_info()

    assert len(stat_cstat) == 1
    stat_cstat = stat_cstat[0]
    assert stat_cstat.numpoints == 771
    assert stat_cstat.dof == 769
    assert stat_cstat.statval == pytest.approx(568.52,
                                               rel=0, abs=0.005)


@requires_fits
@requires_data
@pytest.mark.parametrize("mode", [["arf"], ["rmf"], ["arf", "rmf"]])
def test_1209_response(mode, make_data_path):
    """Do we pick up the header keywords from the response?

    This is related to issue #1209
    """

    # We could set up channels and counts, but let's not.
    #
    d = DataPHA("dummy", None, None)
    assert d.header["TELESCOP"] == "none"
    assert d.header["INSTRUME"] == "none"
    assert d.header["FILTER"] == "none"

    # We do not care about the warning messages here from
    # ENERG_LO replacement.
    #
    if "arf" in mode:
        infile = make_data_path(ARFFILE)
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("ignore")
            arf = io.read_arf(infile)

        d.set_arf(arf)

    if "rmf" in mode:
        infile = make_data_path(RMFFILE)
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("ignore")
            rmf = io.read_rmf(infile)

        d.set_rmf(rmf)

    # The PHA file contains a FILTER keyword but the responses do not.
    #
    assert d.header["TELESCOP"] == "SWIFT"
    assert d.header["INSTRUME"] == "XRT"
    assert d.header["FILTER"] == "none"


@requires_fits
@requires_data
def test_1209_background(make_data_path):
    """Do we pick up the header keywords from the background?

    This is related to issue #1209
    """

    # We need to set up the channels array to match the background.
    #
    d = DataPHA("dummy", np.arange(0, 1024, dtype=np.int16), None)
    assert d.header["TELESCOP"] == "none"
    assert d.header["INSTRUME"] == "none"
    assert d.header["FILTER"] == "none"

    infile = make_data_path(PHAFILE)
    bkg = io.read_pha(infile)
    d.set_background(bkg)

    # The PHA file contains a FILTER keyword but the responses do not.
    #
    assert d.header["TELESCOP"] == "SWIFT"
    assert d.header["INSTRUME"] == "XRT"
    assert d.header["FILTER"] == "NONE"


@requires_fits
@requires_data
def test_get_x_channel0(make_data_path):
    """Are we correctly matching up channel and energies?"""

    infile = make_data_path(PHAFILE)
    rmffile = make_data_path(RMFFILE)
    arffile = make_data_path(ARFFILE)

    pha = io.read_pha(infile)
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")
        rmf = io.read_rmf(rmffile)
        arf = io.read_arf(arffile)

    pha.set_rmf(rmf)
    pha.set_arf(arf)

    # get_x returns all data, no matter the filter. We can check this
    # by filtering out everything.
    pha.ignore()

    pha.set_analysis("channel")
    x_chan = pha.get_x()

    pha.set_analysis("energy")
    x_energy = pha.get_x()

    pha.set_analysis("wave")
    x_wave = pha.get_x()

    assert len(x_chan) == 1024
    assert len(x_energy) == 1024
    assert len(x_wave) == 1024

    # Check some expected values. The energy bins are (modulo rounding
    # issues):
    #
    #    0    -  0.01
    #    0.01 -  0.02
    #    0.02 -  0.03
    #    0.03 -  0.04
    #    0.04 -  0.05
    #    ...
    #   10.23 - 10.24
    #
    assert x_chan[0:5] == pytest.approx([0, 1, 2, 3, 4])
    assert x_energy[0:5] == pytest.approx([0.005, 0.015, 0.025, 0.035, 0.045])
    assert x_wave[0:5] == pytest.approx([2479.68380343, 826.56126781,
                                         495.93676069, 354.24054335,
                                         275.5204169])

    # Do not bother checking wavelength since it is calculated from
    # energy and we have shown above the wavelength mapping is
    # correct.
    #
    assert x_chan[1023] == pytest.approx(1023)
    assert x_energy[1023] == pytest.approx(10.235)

    # check ascending or descending
    assert (x_chan[1:] > x_chan[:-1]).all()
    assert (x_energy[1:] > x_energy[:-1]).all()
    assert (x_wave[1:] < x_wave[:-1]).all()
