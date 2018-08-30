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

"""Test handling of ROSAT [1]_ data.

At the moment the tests are of routines in the sherpa.astro.io module,
since the high-level routins like sherpa.astro.ui.unpack_pha are a
very-thin wrapper around the io routines.

References
----------

.. [1] The ROSAT Mission https://heasarc.gsfc.nasa.gov/docs/rosat/rosat3.html

"""

import numpy as np
from numpy.testing import assert_allclose

import pytest
import six

from sherpa.utils.testing import requires_data, requires_fits
from sherpa.utils.err import IOErr
from sherpa.astro.data import DataPHA, DataRosatRMF

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


@pytest.fixture(autouse=True)
def cleanup_astro_session(request):
    """Ensure sherpa.astro.ui is cleaned before and after the test."""
    ui.clean()

    def fin():
        ui.clean()

    request.addfinalizer(fin)


# The RMF includes the ARF (so is a "RSP" file).
#
PHAFILE = 'xrbg_xspec.pi'
RMFFILE = 'pspcc_gain1_256.rsp'

"""
Plotting the data - once the response has been read in - into
XSPEC 12.9.1p returns the following (this is the output of WDATA and
represents X value, X error, Y value, Y error).

First for X=channel, Y=normalized count/s/channel

    1 0.5 0 1
    2 0.5 2.80112115E-4 4.72556394E-6
    3 0.5 4.47507802E-4 5.81684435E-6
    4 0.5 4.62442949E-5 2.01607168E-6
    5 0.5 7.15584611E-5 2.57730176E-6
    6 0.5 8.68857605E-5 2.57510169E-6
    7 0.5 1.03978207E-4 2.74851345E-6
    8 0.5 5.3703181E-5 2.34389222E-6
    9 0.5 0 1
    ...
    63 0.5 0 1

Now for X=keV, Y=normalized count/s/keV

    4.96509969E-2 3.46949995E-2 0 14.4112997
    0.143822998 5.94770014E-2 2.35479348E-3 3.97259755E-5
    0.312340498 0.109040499 2.05202564E-3 2.66728621E-5
    0.470945001 4.95640039E-2 4.66510886E-4 2.03380623E-5
    0.609724522 8.92154872E-2 4.01042809E-4 1.44442511E-5
    0.803023994 0.104084015 4.17382835E-4 1.2370303E-5
    1.11032045 0.20321247 2.55836174E-4 6.76265927E-6
    1.66048098 0.346947968 7.73937136E-5 3.37787287E-6
    2.01238537 4.95660305E-3 0 100.875542
    ...
    2.54767656 4.95648384E-3 0 100.877968

"""


def validate_pha(pha, errors=False):
    """Has the PSPC PHA file been read in correctly?

    Parameters
    ----------
    pha : sherpa.astro.data.DataPHA
    errors : bool, optional
        Have the systematic errors been read in from the file?

    Notes
    ------
    XSPEC 12.9.1p reports this message when reading in this file:

    ***Warning: Unrecognized grouping for channel(s). It/they will be reset to 1.
    HOWEVER NOTE: Resets only occurred for channels of bad quality

    """

    assert isinstance(pha, DataPHA)

    # assert_allclose defaults to a relative tolerance of 1e-7.
    assert_allclose(pha.exposure, 1.0)
    assert_allclose(pha.backscal, 1.0)
    assert_allclose(pha.areascal, 1.0)

    # Although channel is an integer, it's treated as a
    # float.
    #
    nchan = 256
    assert_allclose(pha.channel, np.arange(1, nchan + 1))

    # The file has a RATE rather than COUNTS column; this gets
    # converted into a counts value, and the rate flag is set.
    # The exposure time here is 1.0, so the values are the
    # same (rate and counts) even if the (implied) units are not.
    #
    assert len(pha.counts) == nchan
    assert_allclose(pha.counts.min(), 0.0)
    assert_allclose(pha.counts.max(), 4.4750780216e-4)
    assert_allclose(pha.counts.sum(), 1.0899898225e-3)
    assert np.argmax(pha.counts) == 29

    idx = np.where(pha.counts > 0)[0]

    nonzerobins = [13, 29, 44, 59, 79, 110, 154]
    assert np.all(idx == nonzerobins)

    assert_allclose(pha.counts[nonzerobins],
                    [2.80112115e-4,
                     4.47507802e-4,
                     4.62442949e-5,
                     7.15584611e-5,
                     8.68857605e-5,
                     1.03978207e-4,
                     5.3703181e-5])

    # Simple statistics on the GROUPING and QUALITY channels
    #
    assert len(pha.quality) == nchan
    assert_allclose(pha.quality.min(), 0)
    assert_allclose(pha.quality.max(), 2)
    assert_allclose(pha.quality.sum(), 124)
    assert np.argmin(pha.quality) == 7
    assert np.argmax(pha.quality) == 0

    assert len(pha.grouping) == nchan
    assert_allclose(pha.grouping.min(), -1)
    assert_allclose(pha.grouping.max(), 1)
    assert_allclose(pha.grouping.sum(), -185)
    assert np.argmin(pha.grouping) == 1
    assert np.argmax(pha.grouping) == 0

    # TODO: check that grouping/quality is applied correctly

    emptyfields = ['syserror', 'bin_lo', 'bin_hi']
    if errors:
        assert len(pha.staterror) == nchan
        assert_allclose(pha.staterror.min(), 0.0)
        assert_allclose(pha.staterror.max(), 5.8168443502e-06)
        assert_allclose(pha.staterror.sum(), 2.2803289085e-05)
        assert np.argmin(pha.staterror) == 0
        assert np.argmax(pha.staterror) == 29

        idx = np.where(pha.counts > 0)[0]
        assert np.all(idx == nonzerobins)

        assert_allclose(pha.staterror[idx],
                        [4.72556394E-6,
                         5.81684435E-6,
                         2.01607168E-6,
                         2.57730176E-6,
                         2.57510169E-6,
                         2.74851345E-6,
                         2.34389222E-6])

    else:
        emptyfields.insert(0, 'staterror')

    for field in emptyfields:
        assert getattr(pha, field) is None

    assert pha.grouped is True
    assert pha.subtracted is False
    assert pha.units == 'channel'
    assert pha.rate is True
    assert pha.plot_fac == 0
    assert pha.response_ids == []
    assert pha.background_ids == []
    assert pha.get_background() is None


def validate_rmf(rmf):
    """Has the PSPC RMF file been read in correctly?

    Parameters
    ----------
    rmf : sherpa.astro.data.DataRMF

    Notes
    -----
    This is based on the values returned by Crates in Python 2.7,
    from CIAO 4.9, since it produces the same values as XSPEC 12.9.1
    when used to evaluate a power-law model.

    """

    assert isinstance(rmf, DataRosatRMF)

    # There are 10 bins which have N_GRP set to 0. These rows are filtered
    # out from the F_CHAN and N_CHAN arrays.
    #
    nchans = 256
    nbins_all = 729
    nbins_nonzero = nbins_all - 10

    assert rmf.detchans == nchans
    assert rmf.offset == 1

    for field in ['energ_lo', 'energ_hi', 'n_grp']:
        assert getattr(rmf, field).size == nbins_all

    for field in ['f_chan', 'n_chan']:
        assert getattr(rmf, field).size == nbins_nonzero

    nmatrix = 115269
    assert rmf.matrix.size == nmatrix

    for field in ['e_min', 'e_max']:
        assert getattr(rmf, field).size == nchans

    # Where are the N_GRP=0 bins?
    zeroidx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 728]
    assert np.all(np.where(rmf.n_grp == 0)[0] == zeroidx)

    # All n_grp values are 0 or 1
    assert rmf.n_grp.sum() == nbins_nonzero
    assert (rmf.n_grp[rmf.n_grp > 0] == 1).all()

    assert rmf.f_chan.min() == 8
    assert rmf.f_chan.max() == 194
    assert rmf.f_chan.sum() == 7732
    assert np.argmax(rmf.f_chan) == 718

    assert rmf.n_chan.min() == 2
    assert rmf.n_chan.max() == 245
    assert rmf.n_chan.sum() == nmatrix
    assert np.argmax(rmf.n_chan) == 410

    # Rather than check each element, use some simple summary
    # statistics.
    assert_allclose(rmf.matrix.min(), 1.0e-5)
    assert_allclose(rmf.matrix.max(), 11.146280288696289)
    assert_allclose(rmf.matrix.sum(), 42994.357679666558)
    assert np.argmin(rmf.matrix) == 0
    assert np.argmax(rmf.matrix) == 7561

    # Expect the upper edge of bin j to equal the lower
    # edge of bin j + 1.
    #
    assert_allclose(rmf.energ_lo[1:], rmf.energ_hi[:-1])
    assert_allclose(rmf.e_min[1:], rmf.e_max[:-1])

    assert_allclose(rmf.energ_lo[0], 0.054607998579740524)
    assert_allclose(rmf.energ_hi[0], 0.071539998054504395)
    assert_allclose(rmf.energ_lo[-1], 3.0)
    assert_allclose(rmf.energ_hi[-1], 3.0099999904632568)

    assert_allclose(rmf.e_min[0], 0.014956000261008739)
    assert_allclose(rmf.e_max[0], 0.024869000539183617)
    assert_allclose(rmf.e_min[-1], 2.5427200794219971)
    assert_allclose(rmf.e_max[-1], 2.5526330471038818)


@requires_data
@requires_fits
def test_read_pha(make_data_path):
    """Can we read in a ROSAT PHA file."""

    infile = make_data_path(PHAFILE)
    pha = io.read_pha(infile)
    validate_pha(pha, errors=False)


@requires_data
@requires_fits
def test_read_pha_errors(make_data_path):
    """Can we read in a ROSAT PHA file."""

    infile = make_data_path(PHAFILE)
    pha = io.read_pha(infile, use_errors=True)
    validate_pha(pha, errors=True)


# In the error checks (i.e. checking that invalid files fail),
# the error messages and types depend on the backend. An attempt
# is made to be specific when testing these, primarily to point
# out when any refactoring changes this (since in many ways it
# would be nice to have the same error messages between the
# backends when possible).
#
@requires_data
@requires_fits
def test_read_pha_fails_rmf(make_data_path):
    """Just check in we can't read in a RMF as a PHA file."""

    if backend == 'pyfits':
        emsg = " does not appear to be a PHA spectrum"
        etype = IOErr
    elif backend == 'crates':
        emsg = "dmKeyRead() could not find key. 'HDUCLAS1'"
        etype = ValueError
    else:
        assert False, "Internal error: unknown backend {}".format(backend)

    infile = make_data_path(RMFFILE)
    with pytest.raises(etype) as excinfo:
        io.read_pha(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_read_arf_fails_pha(make_data_path):
    """Just check in we can't read in a PHA file as an ARF."""

    if backend == 'pyfits':
        emsg = ' does not appear to be an ARF'
    elif backend == 'crates':
        emsg = " is not a filename or ARFCrate obj"
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
    infile = make_data_path(RMFFILE)
    rmf = io.read_rmf(infile)
    validate_rmf(rmf)


@requires_data
@requires_fits
def test_read_rmf_fails_pha(make_data_path):
    """Just check in we can't read in a PHA file as a RMF."""

    if backend == 'pyfits':
        emsg = ' does not appear to be an RMF'
    elif backend == 'crates':
        emsg = "Required column 'ENERG_LO' not found in "
    else:
        assert False, "Internal error: unknown backend {}".format(backend)

    infile = make_data_path(PHAFILE)
    with pytest.raises(IOErr) as excinfo:
        io.read_rmf(infile)

    assert emsg in str(excinfo.value)


@requires_data
@requires_fits
def test_can_use_pspc_data(make_data_path):
    """A basic check that we can read in and use the ROSAT PSPC data.

    Unlike the previous tests, that directly access the io module,
    this uses the ui interface.
    """

    # The PSPC PHA file does not have the ANCRFILE/RESPFILE keywords
    # set up, so the responses has to be manually added.
    #
    ui.load_pha(make_data_path(PHAFILE), use_errors=True)
    assert ui.get_analysis() == 'channel'

    ui.load_rmf(make_data_path(RMFFILE))
    assert ui.get_analysis() == 'energy'

    ui.set_source(ui.powlaw1d.pl)
    ui.set_par('pl.gamma', 1.7)
    ui.set_par('pl.ampl', 2e-6)

    s = ui.get_stat_info()[0]
    assert s.numpoints == 63
    assert s.dof == 61

    # Value obtained from XSPEC 12.9.1p; Sherpa returns
    # sexpected = 973.2270845920297
    sexpected = 973.23
    assert_allclose(s.statval, sexpected, rtol=0, atol=0.005)

    # apply an energy filter to remove the "bogus" points
    ui.ignore(None, 0.05)

    s = ui.get_stat_info()[0]
    assert s.numpoints == 62
    assert s.dof == 60
    assert_allclose(s.statval, sexpected, rtol=0, atol=0.005)

    ui.ignore(2.01, None)

    s = ui.get_stat_info()[0]
    assert s.numpoints == 7
    assert s.dof == 5

    assert_allclose(s.statval, sexpected, rtol=0, atol=0.005)
