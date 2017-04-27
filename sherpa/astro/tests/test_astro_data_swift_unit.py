#
#  Copyright (C) 2017  Smithsonian Astrophysical Observatory
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

import numpy as np
from numpy.testing import assert_allclose

import pytest

from sherpa.utils import requires_data, requires_fits
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


# PHA and ARF are stored uncompressed while RMF is gzip-compressed.
# Use the same name (no .gz suffix) to access all files.
#
PHAFILE = 'target_sr.pha'
ARFFILE = 'target_sr.arf'
RMFFILE = 'swxpc0to12s6_20130101v014.rmf'


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
        etype == TypeError

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

    infile = make_data_path(RMFFILE)
    with pytest.raises(etype) as excinfo:
        io.read_pha(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_arf(make_data_path):
    """Can we read in a Swift ARF."""

    infile = make_data_path(ARFFILE)
    arf = io.read_arf(infile)
    assert isinstance(arf, DataARF)

    nbins = 2400

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert len(arf.energ_lo) == nbins
    assert_allclose(arf.energ_lo[1:], arf.energ_hi[:-1])
    assert_allclose(arf.energ_lo[0], 0.0)
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

    infile = make_data_path(RMFFILE)
    with pytest.raises(IOErr) as excinfo:
        io.read_arf(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_rmf(make_data_path):
    """Can we read in a Swift RMF."""

    infile = make_data_path(RMFFILE)
    rmf = io.read_rmf(infile)
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

    assert_allclose(rmf.energ_lo[0], 0.0)
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

    infile = make_data_path(ARFFILE)
    with pytest.raises(etype) as excinfo:
        io.read_rmf(infile)

    assert emsg in str(excinfo.value)
