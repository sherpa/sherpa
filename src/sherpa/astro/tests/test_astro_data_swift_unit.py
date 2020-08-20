#
#  Copyright (C) 2017, 2018, 2020  Smithsonian Astrophysical Observatory
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
since the high-level routins like sherpa.astro.ui.unpack_pha are a
very-thin wrapper around the io routines.

References
----------

.. [1] The Swift Gamma Ray-Burst Mission https://swift.gsfc.nasa.gov/

"""

import warnings

import numpy as np
from numpy.testing import assert_allclose

import pytest

from sherpa.testing import requires_data, requires_fits
from sherpa.utils.err import IOErr
from sherpa.astro.data import DataARF, DataPHA, DataRMF

# Should each test import io instead of this? Also, do we have a
# better way of determining what the backend is?
#
try:
    from sherpa.astro import io
    if io.backend.__name__ == "sherpa.astro.io.pyfits_backend":
        backend = "pyfits"
    elif io.backend.__name__ == "sherpa.astro.io.crates_backend":
        backend = "crates"
    else:
        # Should not happen, but do not want to error out here.
        # Leave io as whatever was loaded, which will likely cause
        # the tests to fail.
        backend = None

except ImportError:
    io = None
    backend = None

from sherpa.astro import ui


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


@requires_data
@requires_fits
def test_read_pha(make_data_path):
    """Can we read in a Swift PHA file."""

    infile = make_data_path(PHAFILE)
    pha = io.read_pha(infile)
    assert isinstance(pha, DataPHA)

    # assert_allclose defaults to a tolerance of 1e-7
    # which is close to the precision these values are stored
    # in the file.
    assert_allclose(pha.exposure, 4812.26699895)
    assert_allclose(pha.backscal, 1.148e-3)
    assert_allclose(pha.areascal, 1.0)

    # Although channel is an integer, it's treated as a
    # float. Note that the channels start at 0 in the file,
    # but they are promoted to 1 by the crates backend.
    # XSPEC 12.9.1b, on reading this file, reports channel numbers
    # between 1 and 1024 inclusive.
    #
    nchan = 1024
    assert_allclose(pha.channel, np.arange(1, nchan + 1))

    # Rather than check each element, use some simple summary
    # statistics.
    assert len(pha.counts) == nchan
    assert_allclose(pha.counts.min(), 0.0)
    assert_allclose(pha.counts.max(), 3.0)
    assert_allclose(pha.counts.sum(), 58.0)
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

    if backend == 'pyfits':
        emsg = " does not appear to be a PHA spectrum"
        etype = IOErr
    elif backend == 'crates':
        emsg = 'File must be a PHA file.'
        etype = TypeError
    else:
        assert False, "Internal error: unknown backend {}".format(backend)

    infile = make_data_path(ARFFILE)
    with pytest.raises(etype) as excinfo:
        io.read_pha(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_pha_fails_rmf(make_data_path):
    """Just check in we can't read in a RMF as a PHA file."""

    if backend == 'pyfits':
        emsg = " does not appear to be a PHA spectrum"
        etype = IOErr
    elif backend == 'crates':
        emsg = 'File must be a PHA file.'
        etype = TypeError
    else:
        assert False, "Internal error: unknown backend {}".format(backend)

    infile = make_data_path(RMFFILE)
    with pytest.raises(etype) as excinfo:
        io.read_pha(infile)

    assert emsg in str(excinfo.value)


def validate_replacement_warning(ws, rtype, label):

    assert len(ws) == 1
    w = ws[0]
    assert w.category == UserWarning

    emsg = "The minimum ENERG_LO in the {} '{}' ".format(rtype, label) + \
           "was 0 and has been replaced by {}".format(EMIN)
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

    assert isinstance(arf, DataARF)

    nbins = 2400

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert len(arf.energ_lo) == nbins
    assert_allclose(arf.energ_lo[1:], arf.energ_hi[:-1])
    assert_allclose(arf.energ_lo[0], EMIN)
    assert_allclose(arf.energ_hi[0], 0.005)
    assert_allclose(arf.energ_lo[-1], 11.995)
    assert_allclose(arf.energ_hi[-1], 12.0)

    # Rather than check each element, use some simple summary
    # statistics.
    assert len(arf.specresp) == nbins
    assert_allclose(arf.specresp.min(), 0.1506345272064209)
    assert_allclose(arf.specresp.max(), 145.38259887695312)
    assert_allclose(arf.specresp.sum(), 178696.1007297188)
    assert np.argmax(arf.specresp) == 309

    for field in ['bin_lo', 'bin_hi', 'exposure']:
        assert getattr(arf, field) is None


@requires_data
@requires_fits
def test_read_arf_fails_pha(make_data_path):
    """Just check in we can't read in a PHA file as an ARF."""

    if backend == 'pyfits':
        emsg = ' does not appear to be an ARF'
    elif backend == 'crates':
        emsg = "Required column 'ENERG_LO' not found in "
    else:
        assert False, "Internal error: unknown backend {}".format(backend)

    infile = make_data_path(PHAFILE)
    with pytest.raises(IOErr) as excinfo:
        io.read_arf(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_arf_fails_rmf(make_data_path):
    """Just check in we can't read in a RNF as an ARF."""

    if backend == 'pyfits':
        emsg = ' does not appear to be an ARF'
    elif backend == 'crates':
        emsg = "Required column 'SPECRESP' not found in "
    else:
        assert False, "Internal error: unknown backend {}".format(backend)

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
    assert_allclose(rmf.matrix.min(), 0.0)
    assert_allclose(rmf.matrix.max(), 0.05834391713142395)
    assert_allclose(rmf.matrix.sum(), 1133.8128197300257)
    assert np.argmax(rmf.matrix) == 269443

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert_allclose(rmf.energ_lo[1:], rmf.energ_hi[:-1])
    assert_allclose(rmf.e_min[1:], rmf.e_max[:-1])

    assert_allclose(rmf.energ_lo[0], EMIN)
    assert_allclose(rmf.energ_hi[0], 0.005)
    assert_allclose(rmf.energ_lo[-1], 11.995)
    assert_allclose(rmf.energ_hi[-1], 12.0)

    assert_allclose(rmf.e_min[0], 0.0)
    assert_allclose(rmf.e_max[0], 0.01)
    assert_allclose(rmf.e_min[-1], 10.23)
    assert_allclose(rmf.e_max[-1], 10.24)


@requires_data
@requires_fits
def test_read_rmf_fails_pha(make_data_path):
    """Just check in we can't read in a PHA file as a RMF."""

    if backend == 'pyfits':
        emsg = ' does not appear to be an RMF'
        etype = IOErr
    elif backend == 'crates':
        emsg = ' does not contain a Response Matrix.'
        etype = TypeError
    else:
        assert False, "Internal error: unknown backend {}".format(backend)

    infile = make_data_path(PHAFILE)
    with pytest.raises(etype) as excinfo:
        io.read_rmf(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_rmf_fails_arf(make_data_path):
    """Just check in we can't read in a ARF as a RMF."""

    if backend == 'pyfits':
        emsg = " does not have a 'DETCHANS' keyword"
        etype = IOErr
    elif backend == 'crates':
        emsg = ' does not contain a Response Matrix.'
        etype = TypeError
    else:
        assert False, "Internal error: unknown backend {}".format(backend)

    infile = make_data_path(ARFFILE)
    with pytest.raises(etype) as excinfo:
        io.read_rmf(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_can_use_swift_data(make_data_path):
    """A basic check that we can read in and use the Swift data.

    Unlike the previous tests, that directly access the io module,
    this uses the ui interface.
    """

    # QUS are there pytest fixtures that ensure the state is
    # clean on entry and exit?
    ui.clean()

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
    assert_allclose(stat, 58.2813692358182)

    # Pick an energy range which isn't affected by the first
    # bin.
    #
    # Unfortunately, using a range of 0.3-8.0 gives 771 bins
    # in XSPEC - channels 30 to 800 - but 772 bins in Sherpa.
    # If I use ignore(None, 0.3); ignore(8.0, None) instead
    # then the result is 771 bins. This is because the e_min/max
    # of the RMF has channel widths of 0.01 keV, starting at 0,
    # so both 0.3 and 8.0 fall on a bin boundary. So, it's either
    # a difference in < or <= (or > vs >=), or a rounding issue
    # due to floating-point conversion leading to one bin boundary
    # being slightly different in Sherpa vs XSPEC).
    #
    # When using ui.notice(0.3, 8.0); ui.get_indep(filter=True)
    # returns 772 channels, 30 to 801.
    #
    # Using ui.notice(0.3, 7.995) selects channels 30 to 800. So
    # this range is used. Alternatively, channel 801 could have been
    # excluded explicitly.
    #
    # ui.notice(0.3, 8.0)
    ui.notice(0.3, 7.995)

    # XSPEC 12.9.1b calculation of the statistic:
    #   chi sq = 203.88 from 771 bins with 769 dof
    #   cstat  = 568.52
    #
    # There are known differences between XSPEC and Sherpa
    # with chi2xspecvar. This only affects data sets where
    # there is background subtraction, which is not the case
    # here. See https://github.com/sherpa/sherpa/issues/356
    #
    ui.set_stat('chi2xspecvar')
    stat_xvar = ui.get_stat_info()

    assert len(stat_xvar) == 1
    stat_xvar = stat_xvar[0]
    assert stat_xvar.numpoints == 771
    assert stat_xvar.dof == 769
    assert_allclose(stat_xvar.statval, 203.88,
                    rtol=0, atol=0.005)

    ui.set_stat('cstat')
    stat_cstat = ui.get_stat_info()

    assert len(stat_cstat) == 1
    stat_cstat = stat_cstat[0]
    assert stat_cstat.numpoints == 771
    assert stat_cstat.dof == 769
    assert_allclose(stat_cstat.statval, 568.52,
                    rtol=0, atol=0.005)

    ui.clean()
