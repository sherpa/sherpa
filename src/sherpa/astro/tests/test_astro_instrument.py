#
#  Copyright (C) 2017, 2020 - 2024
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

"""
Test the instrument-handling code without requiring external
data files.
"""

from dataclasses import dataclass
import logging
import warnings

import numpy as np
from numpy.testing import assert_allclose

import pytest  # pytest >= 3.0 is needed

from sherpa.astro import hc
from sherpa.astro.data import DataPHA, DataRMF, DataIMG
from sherpa.astro.instrument import ARF1D, ARFModelNoPHA, ARFModelPHA, \
    Response1D, RMF1D, RMFModelNoPHA, RMFModelPHA, \
    RSPModelNoPHA, RSPModelPHA, PSFModel, RMFMatrix, \
    create_arf, create_delta_rmf, create_non_delta_rmf, \
    rmf_to_matrix, rmf_to_image, has_pha_response
from sherpa.astro import io
from sherpa.data import Data1D
from sherpa.fit import Fit
from sherpa.models.basic import Box1D, Const1D, Gauss1D, Polynom1D, \
    PowLaw1D, TableModel
from sherpa.models.model import ArithmeticModel, \
    ArithmeticConstantModel, BinaryOpModel
from sherpa.utils.err import DataErr, PSFErr
from sherpa.utils.testing import requires_xspec, requires_data, requires_fits

try:
    from sherpa.astro.xspec import XSconstant
except ImportError:
    XSConstant = None


def validate_zero_replacement(ws, rtype, label, ethresh):
    """Ensure that there is one warning about replacing a 0 energy bin"""

    assert len(ws) == 1
    w = ws[0]
    assert w.category == UserWarning

    emsg = f"The minimum ENERG_LO in the {rtype} '{label}' " + \
           f"was 0 and has been replaced by {ethresh}"
    assert str(w.message) == emsg


def get_non_delta_matrix():
    """Return a 2D matrix representing the RMF in create_non_delta_rmf_local."""

    # X axis is channel number (e_min/max grid)
    # Y axis is energ_lo/hi grid
    #
    # Have nX=8 nY=20
    #
    matrix = np.zeros((20, 8), dtype=np.float32)
    matrix[0, 0] = 1.0
    matrix[1, 0] = 0.5
    matrix[1, 1] = 0.5
    matrix[2, 0] = 0.5
    matrix[2, 1] = 0.5
    matrix[3, 0] = 0.5
    matrix[3, 1] = 0.5
    matrix[4, 0] = 0.4
    matrix[4, 1] = 0.3
    matrix[4, 2] = 0.3
    matrix[5, 1] = 0.6
    matrix[5, 2] = 0.4
    matrix[6, 1] = 1.0
    matrix[7, 0] = 0.1
    matrix[7, 1] = 0.2
    matrix[7, 2] = 0.4
    matrix[7, 3] = 0.3
    matrix[8, 1] = 0.2
    matrix[8, 2] = 0.4
    matrix[8, 4] = 0.2
    matrix[8, 5] = 0.2
    matrix[9, 3] = 1.0
    matrix[10, 3] = 1.0
    matrix[11, 4] = 1.0
    matrix[12, 4] = 1.0
    matrix[13, 4] = 0.1
    matrix[13, 5] = 0.8
    matrix[13, 6] = 0.1
    matrix[14, 1] = 0.2
    matrix[14, 5] = 0.8
    matrix[15, 2] = 0.3
    matrix[15, 6] = 0.6
    matrix[15, 7] = 0.1
    matrix[16, 5] = 0.5
    matrix[16, 6] = 0.5
    matrix[17, 6] = 0.5
    matrix[17, 7] = 0.5
    matrix[18, 6] = 0.4
    matrix[18, 7] = 0.6
    matrix[19, 1] = 0.1
    matrix[19, 3] = 0.3
    matrix[19, 5] = 0.6

    # just check that the matrix is correct
    for row in matrix:
        assert row.sum() == pytest.approx(1.0)

    return matrix


# The _local suffix is to avoid a conflict with
# sherpa.astro.instrument.create_non_delta_rmf
#
def create_non_delta_rmf_local():
    """Create a RMF which does not have a delta-function response.

    This is hard-coded to have a range of behavior: some energies
    it is a delta function, some a small (width) response, and some multiple
    peaks.

    Returns
    -------
    rmf : DataRMF instance

    """

    startchan = 1

    # e_min/max define the "X" axis of the matrix (which is the
    # number of channels)
    # energ_lo/hi define the "Y" axis of the matrix
    #
    # The bins do not have to be constant-width
    echan = np.asarray([0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    e_min = echan[:-1]
    e_max = echan[1:]

    nchans = e_min.size

    energ = np.arange(0.05, 1.1, 0.05)
    energ_lo = energ[:-1]
    energ_hi = energ[1:]

    # Deconstruct the matrix to get the "condensed" representation
    # used by the RMF.
    #
    matrix = get_non_delta_matrix()
    n_grp = []
    n_chan = []
    f_chan = []
    for row in (matrix > 0):
        # simple run-length encoding, bsaed on code in
        # https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array
        #
        flag = np.hstack([[0], row, [0]])
        diffs = np.diff(flag, n=1)
        starts, = np.where(diffs > 0)
        ends, = np.where(diffs < 0)

        n_chan.extend(ends - starts)
        f_chan.extend(starts + 1)
        n_grp.append(len(starts))

    matrix = matrix.flatten()
    matrix = matrix[matrix > 0]

    n_grp = np.asarray(n_grp, dtype=np.int16)
    f_chan = np.asarray(f_chan, dtype=np.int16)
    n_chan = np.asarray(n_chan, dtype=np.int16)

    return DataRMF('non-delta-rmf', detchans=nchans,
                   energ_lo=energ_lo, energ_hi=energ_hi,
                   n_grp=n_grp, n_chan=n_chan,
                   f_chan=f_chan, matrix=matrix,
                   offset=startchan,
                   e_min=e_min, e_max=e_max)


def create_non_delta_specresp():
    """A SPECRESP column compatible with create_non_delta_rmf_local.

    Returns
    -------
    specresp : array
    """

    return np.asarray([20.0, 10.0, 0.0, 10.0, 15.0,
                       20.0, 20.0, 20.0, 20.0, 20.0,
                       18.0, 17.0, 0.0, 20.0, 20.0,
                       15.0, 0.0, 10.0, 20.0, 20.0],
                      dtype=np.float32)


def test_datarmf_get_y():
    """Just check the y accessor works.

    This is technically a test of sherpa.astro.data code
    but do it here as it's easier (since we've set up
    a "complex" RMF anyway).
    """

    rmf = create_non_delta_rmf_local()
    matrix = get_non_delta_matrix()
    y = matrix[matrix > 0]
    assert rmf.y == pytest.approx(y)


def test_arf1d_empty():

    arf = ARF1D(None)

    assert arf._arf is None
    assert arf._pha is None

    assert str(arf) == str(None)

    # Since there's no ARF to fall through, it should be an error
    # to access the name attribute.
    #
    with pytest.raises(AttributeError):
        arf.name


def test_rmf1d_empty():

    rmf = RMF1D(None)

    assert rmf._rmf is None
    assert rmf._pha is None

    assert str(rmf) == str(None)

    # Since there's no RMF to fall through, it should be an error
    # to access the name attribute.
    #
    with pytest.raises(AttributeError):
        rmf.name


def test_rsp1d_empty():

    # As there's currently no explicit check for the input arg
    # being set, the result of sending None in should be an
    # Attribute Error.
    #
    with pytest.raises(AttributeError):
        Response1D(None)


def test_rsp1d_pha_empty():

    # Create a PHA with no ARF or RMF
    channels = np.arange(1, 5, dtype=np.int16)
    counts = np.asarray([10, 5, 12, 7], dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=12.2)

    with pytest.raises(DataErr,
                       match="No instrument response found for dataset test-pha"):
        Response1D(pha)


def test_arf1d_no_pha_no_exposure_basic():
    "Can we create an ARF with no PHA?"

    egrid = np.arange(0.1, 0.6, 0.1)
    adata = create_arf(egrid[:-1], egrid[1:])
    arf = ARF1D(adata)

    assert arf._arf == adata
    assert arf._pha is None

    # Does the ARF1D pass through functionality to the DataARF?
    assert arf.name == 'user-arf'
    assert arf.exposure is None
    assert str(arf) == str(adata)

    assert dir(arf) == dir(adata)

    # Should be exact
    specresp = np.ones(egrid.size - 1, dtype=np.float32)
    assert (specresp == arf.specresp).all()


def test_rmf1d_no_pha_delta():
    "Can we create an RMF (delta function) with no PHA?"

    egrid = np.arange(0.1, 0.6, 0.1)
    rdata = create_delta_rmf(egrid[:-1], egrid[1:])
    rmf = RMF1D(rdata)

    assert rmf._rmf == rdata
    assert rmf._pha is None

    # Does the RMF1D pass through functionality to the DataRMF?
    assert rmf.name == 'delta-rmf'
    assert str(rmf) == str(rdata)

    assert dir(rmf) == dir(rdata)

    matrix = np.ones(egrid.size - 1, dtype=np.float32)
    assert (matrix == rmf.matrix).all()

    # Unlike the ARF, the RMF doesn't have an exposure field
    with pytest.raises(AttributeError):
        rmf.exposure


def test_rmf1d_no_pha_matrix():
    "Can we create an RMF (matrix) with no PHA?"

    rdata = create_non_delta_rmf_local()
    rmf = RMF1D(rdata)

    assert rmf._rmf == rdata
    assert rmf._pha is None

    # Does the RMF1D pass through functionality to the DataRMF?
    assert rmf.name == 'non-delta-rmf'
    assert str(rmf) == str(rdata)

    assert dir(rmf) == dir(rdata)

    # We do not do a full test, just some global checks and a few
    # spot checks
    rmatrix = rmf.matrix
    assert rmatrix.sum() == pytest.approx(20.0)

    expected = [1.0, 0.4, 0.2, 1.0]
    assert_allclose(rmatrix[[0, 7, 19, 23]], expected)

    elo = np.linspace(0.05, 1.0, 20)
    assert_allclose(rmf.energ_lo, elo)

    emin = [0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    assert_allclose(rmf.e_min, emin)

    assert (rmf.energ_hi > rmf.energ_lo).all()
    assert (rmf.e_max > rmf.e_min).all()

    # Unlike the ARF, the RMF doesn't have an exposure field
    with pytest.raises(AttributeError):
        rmf.exposure


@pytest.mark.parametrize("exposure", [None, 12.2])
def test_arf1d_no_pha_no_exposure_call(exposure):
    "Can we call an ARF (no PHA)"

    egrid = np.arange(0.1, 0.6, 0.1)
    specresp = np.asarray([1.1, 1.2, 1.3, 1.4])
    adata = create_arf(egrid[:-1], egrid[1:], specresp,
                       exposure=exposure)
    arf = ARF1D(adata)

    mdl = Const1D('flat')
    mdl.c0 = 2.3

    wrapped = arf(mdl)
    assert isinstance(wrapped, ARFModelNoPHA)

    wmdl = wrapped.model
    if exposure is None:
        assert wmdl == mdl
    else:
        # It looks like equality is not well defined for the model
        # classes, since the following assertion fails
        #     assert wmdl == (exposure * mdl)
        # so manually check
        #
        assert isinstance(wmdl, BinaryOpModel)
        assert wmdl.op == np.multiply
        assert wmdl.rhs == mdl
        assert isinstance(wmdl.lhs, ArithmeticConstantModel)
        assert wmdl.lhs.val == pytest.approx(exposure)


def test_rmf1d_simple_no_pha_call():
    "Can we call an RMF (delta function) with no PHA"

    egrid = np.arange(0.1, 0.6, 0.1)
    rdata = create_delta_rmf(egrid[:-1], egrid[1:])
    rmf = RMF1D(rdata)

    mdl = Const1D('flat')
    mdl.c0 = 2.3

    wrapped = rmf(mdl)
    assert isinstance(wrapped, RMFModelNoPHA)

    wmdl = wrapped.model
    assert wmdl == mdl


def test_rmf1d_matrix_no_pha_call():
    "Can we call an RMF (delta function) with no PHA"

    rdata = create_non_delta_rmf_local()
    rmf = RMF1D(rdata)

    mdl = Const1D('flat')
    mdl.c0 = 2.3

    wrapped = rmf(mdl)
    assert isinstance(wrapped, RMFModelNoPHA)

    wmdl = wrapped.model
    assert wmdl == mdl


def test_rmf1d_delta_arf_no_pha_call():
    "Can we call an RMF (delta function) with ARF (no PHA)"

    egrid = np.arange(0.1, 0.6, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]
    rdata = create_delta_rmf(elo, ehi)
    adata = create_arf(elo, ehi)
    rmf = RMF1D(rdata, arf=adata)

    mdl = Const1D('flat')
    mdl.c0 = 2.3

    wrapped = rmf(mdl)

    # It seems like the wrapper should match the following:
    #     assert isinstance(wrapped, RSPModelNoPHA)
    # but at the time the test was written (June 2017) it
    # does not.
    #
    assert isinstance(wrapped, RMFModelNoPHA)

    wmdl = wrapped.model
    assert wmdl == mdl


def test_rmf1d_matrix_arf_no_pha_call():
    "Can we call an RMF () with ARF (no PHA)"

    # NOTE: there is no check that the grids are compatible
    #       so this is probably not that useful a test
    #
    rdata = create_non_delta_rmf_local()
    elo = rdata.e_min
    ehi = rdata.e_max
    adata = create_arf(elo, ehi)
    rmf = RMF1D(rdata, arf=adata)

    mdl = Const1D('flat')
    mdl.c0 = 2.3

    wrapped = rmf(mdl)

    # It seems like the wrapper should match the following:
    #     assert isinstance(wrapped, RSPModelNoPHA)
    # but at the time the test was written (June 2017) it
    # does not.
    #
    assert isinstance(wrapped, RMFModelNoPHA)

    wmdl = wrapped.model
    assert wmdl == mdl


def test_arfmodelnopha_call():
    "What happens calling an arf with no pha?"

    # Note: the exposure is set in the ARF, but should not be
    #       used when evaluating the model; it's value has been
    #       set to a value that the test will fail it it is.
    #
    egrid = np.arange(0.01, 0.06, 0.01)
    svals = [1.1, 1.2, 1.3, 1.4]
    specresp = np.asarray(svals)
    adata = create_arf(egrid[:-1], egrid[1:], specresp,
                       exposure=200.1)

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    wrapped = ARFModelNoPHA(adata, mdl)

    # The model is evaluated on the ARF grid, not whatever
    # is sent in. It is also integrated across the bins,
    # which is why there is a multiplication by the
    # grid width (for this constant model).
    #
    de = egrid[1:] - egrid[:-1]
    expected = constant * np.asarray(svals) * de
    out = wrapped([4, 5])
    assert_allclose(out, expected)


def test_rmfmodelnopha_delta_call():
    "What happens calling an rmf (delta function) with no pha?"

    egrid = np.arange(0.01, 0.06, 0.01)
    rdata = create_delta_rmf(egrid[:-1], egrid[1:])

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    wrapped = RMFModelNoPHA(rdata, mdl)

    # The model is evaluated on the RMF grid, not whatever
    # is sent in. It is also integrated across the bins,
    # which is why there is a multiplication by the
    # grid width (for this constant model).
    #
    de = egrid[1:] - egrid[:-1]
    expected = constant * de
    out = wrapped([4, 5])
    assert_allclose(out, expected)


def test_rmfmodelnopha_matrix_call():
    "What happens calling an rmf (matrix) with no pha?"

    rdata = create_non_delta_rmf_local()

    # Do not use a flat model as it is not as useful a check
    # that the RMF is doing its job.
    #
    constant = 2.3
    slope = -0.1
    mdl = Polynom1D('mdl')
    mdl.c0 = constant
    mdl.c1 = slope

    wrapped = RMFModelNoPHA(rdata, mdl)

    # Calculate the model analytically.
    #
    modvals = mdl(rdata.energ_lo, rdata.energ_hi)
    matrix = get_non_delta_matrix()
    expected = np.matmul(modvals, matrix)

    out = wrapped([4, 5])
    assert_allclose(out, expected)

    # the RMF convolution shouldn't lose flux, although
    # given the previous check it's not clear if this really
    # adds any extra confidence.
    #
    assert out.sum() == pytest.approx(modvals.sum())


def test_rspmodelnopha_delta_call():
    "What happens calling an RMF (delta)+ARF with no pha?"

    exposure = 200.1
    egrid = np.arange(0.01, 0.06, 0.01)
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = np.asarray([1.2, 0.0, 0.5, 4.3])
    rdata = create_delta_rmf(elo, ehi)
    adata = create_arf(elo, ehi, specresp, exposure=exposure)

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    wrapped = RSPModelNoPHA(adata, rdata, mdl)

    # The model is evaluated on the RMF grid, not whatever
    # is sent in. It is also integrated across the bins,
    # which is why there is a multiplication by the
    # grid width (for this constant model).
    #
    de = egrid[1:] - egrid[:-1]
    expected = constant * specresp * de
    out = wrapped([4, 5])
    assert_allclose(out, expected)


def test_rspmodelnopha_matrix_call():
    "What happens calling an RMF (matrix)+ARF with no pha?"

    rdata = create_non_delta_rmf_local()
    exposure = 200.1
    specresp = np.asarray([200.0, 100.0, 0.0, 175.0, 300.0,
                           400.0, 350.0, 200.0, 250.0, 300.0,
                           200.0, 100.0, 100.0, 150.0, 175.0,
                           125.0, 100.0, 90.0, 80.0, 0.0])
    adata = create_arf(rdata.energ_lo, rdata.energ_hi,
                       specresp, exposure=exposure)

    constant = 2.3
    slope = -0.25
    mdl = Polynom1D('mdl')
    mdl.c0 = constant
    mdl.c1 = slope

    wrapped = RSPModelNoPHA(adata, rdata, mdl)

    # Calculate the model analytically. Note that the exposure
    # value is ignored.
    #
    modvals = specresp * mdl(rdata.energ_lo, rdata.energ_hi)
    matrix = get_non_delta_matrix()
    expected = np.matmul(modvals, matrix)

    out = wrapped([4, 5])
    assert_allclose(out, expected)


@pytest.mark.parametrize("ignore", [None, 0, 1, 2, 3])
def test_arfmodelpha_call(ignore):
    """What happens calling an arf with a pha?

    The ignore value indicates what channel to ignore (0 means
    nothing is ignored). The aim is to check edge effects,
    and as there are only a few channels, it was decided to
    test all channels.
    """

    # Note: the exposure is set in the PHA and ARF, but should not be
    #       used when evaluating the model; it's value has been
    #       set to a value that the test will fail it it is.
    #
    exposure = 200.1
    estep = 0.01
    egrid = np.arange(0.01, 0.06, estep)
    svals = [1.1, 1.2, 1.3, 1.4]
    specresp = np.asarray(svals)
    adata = create_arf(egrid[:-1], egrid[1:], specresp,
                       exposure=exposure)

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    channels = np.arange(1, 5, dtype=np.int16)
    counts = np.asarray([10, 5, 12, 7], dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_arf(adata)

    # force energy units (only needed if ignore is set)
    pha.set_analysis('energy')

    if ignore is not None:
        de = estep * 0.9
        e0 = egrid[ignore]
        pha.notice(lo=e0, hi=e0 + de, ignore=True)

        # The assert are intended to help people reading this
        # code rather than being a useful check that the code
        # is working.
        mask = [True, True, True, True]
        mask[ignore] = False
        assert (pha.mask == mask).all()

    wrapped = ARFModelPHA(adata, pha, mdl)

    # The model is evaluated on the ARF grid, not whatever
    # is sent in. It is also integrated across the bins,
    # which is why there is a multiplication by the
    # grid width (for this constant model).
    #
    # Note that the filter doesn't change the grid.
    #
    de = egrid[1:] - egrid[:-1]
    expected = constant * np.asarray(svals) * de
    out = wrapped([4, 5])
    assert_allclose(out, expected)


@pytest.mark.parametrize("analysis", ["energy", "wave"])
def test_rmfmodelpha_delta_no_ebounds(analysis, caplog):
    """What happens calling an rmf with a pha and no EBOUNDS is set

    create_delta_rmf, prior to 4.16.0, did not set the e_min/max
    values when called like this, which meant that we couldn't
    use an energy filter, so a test was created to test the
    behavior. It now sets the values so we just check it works
    as expected. See test_rmfmodelpha_delta_no_ebounds_manual.

    """

    egrid = np.arange(0.01, 0.07, 0.01)
    rdata = create_delta_rmf(egrid[:-1], egrid[1:])

    channels = np.arange(1, 6, dtype=np.int16)
    counts = np.asarray([10, 5, 12, 7, 3], dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts)
    pha.set_rmf(rdata)

    pha.set_analysis(analysis)
    with caplog.at_level(logging.INFO, logger='sherpa'):
        pha.notice(0.025, 0.045, ignore=False)

    assert len(caplog.records) == 0
    if analysis == "energy":
        assert pha.mask == pytest.approx([False, True, True, True, False])
    else:
        assert not pha.mask.any()


@pytest.mark.parametrize("analysis", ["energy", "wave"])
def test_rmfmodelpha_delta_no_ebounds_manual(analysis, caplog):
    """What happens calling an rmf with a pha and no EBOUNDS is set

    Ensure we can't filter on energy or wavelength since there's no
    EBOUNDS information. This behavior was seen when writing
    test_rmfmodelpha_call, so a test was written for it.

    The code used to raise a DataErr but now just displays a
    logged warning.
    """

    egrid = np.arange(0.01, 0.07, 0.01)
    nchans = egrid.size - 1
    matrix = np.ones(nchans, dtype=np.float32)
    dummy = np.ones(nchans, dtype=np.int16)
    f_chan = np.arange(1, nchans + 1, dtype=np.int16)

    rdata = DataRMF(name="f", detchans=nchans, energ_lo=egrid[:-1], energ_hi=egrid[1:],
                    n_grp=dummy, n_chan=dummy, f_chan=f_chan,
                    matrix=matrix, offset=1, e_min=None, e_max=None,
                    ethresh=None, header=None)

    channels = np.arange(1, 6, dtype=np.int16)
    counts = np.asarray([10, 5, 12, 7, 3], dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts)
    pha.set_rmf(rdata)

    pha.set_analysis(analysis)
    with caplog.at_level(logging.INFO, logger='sherpa'):
        pha.notice(0.025, 0.045, ignore=False)

    assert len(caplog.records) == 1
    log_name, log_level, message = caplog.record_tuples[0]
    assert log_name == 'sherpa.astro.data'
    assert log_level == logging.INFO
    assert message == 'Skipping dataset test-pha: RMF does not specify energy bins'


@pytest.mark.parametrize("analysis", ["energy", "wave"])
def test_rmfmodelpha_matrix_mismatch(analysis):
    """Check that an error is raised if there's a mismatch.

    """

    exposure = 200.1
    rdata = create_non_delta_rmf_local()

    # nchans should be rdata.e_min.size for the sizes to match
    nchans = rdata.energ_lo.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    with pytest.raises(DataErr,
                       match="RMF 'non-delta-rmf' is incompatible with PHA dataset 'test-pha'"):
        pha.set_analysis(analysis)


@pytest.mark.parametrize("analysis", ["energy", "wave"])
def test_rsp_normf_error(analysis):
    """Check that an error is raised on set_analysis

    """

    exposure = 200.1

    # rdata is only used to define the grids
    rdata = create_non_delta_rmf_local()
    specresp = create_non_delta_specresp()
    adata = create_arf(rdata.energ_lo,
                       rdata.energ_hi,
                       specresp,
                       exposure=exposure)

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_arf(adata)

    with pytest.raises(DataErr,
                       match="response incomplete for dataset test-pha, " +
                       "check the instrument model"):
        pha.set_analysis(analysis)


def test_rsp_norsp_error():
    """Check that an error is raised when creating a wrapped model

    """

    # rdata is only used to define the grids
    rdata = create_non_delta_rmf_local()

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts)

    with pytest.raises(DataErr,
                       match="No instrument response found for dataset test-pha"):
        Response1D(pha)


@pytest.mark.parametrize("ignore", [None, 0, 1, 17, 26])
def test_rmfmodelpha_delta_call(ignore):
    """What happens calling an rmf (delta) with a pha?

    The ignore value gives the channel to ignore (counting from 0).
    """

    exposure = 200.1
    estep = 0.025
    egrid = np.arange(0.1, 0.8, estep)
    elo = egrid[:-1]
    ehi = egrid[1:]
    rdata = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    nchans = elo.size

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    # force energy units (only needed if ignore is set)
    pha.set_analysis('energy')

    if ignore is not None:
        de = estep * 0.9
        e0 = egrid[ignore]
        pha.notice(lo=e0, hi=e0 + de, ignore=True)

        # The assert are intended to help people reading this
        # code rather than being a useful check that the code
        # is working.
        mask = [True] * nchans
        mask[ignore] = False
        assert (pha.mask == mask).all()

    wrapped = RMFModelPHA(rdata, pha, mdl)

    # The model is evaluated on the RMF grid, not whatever
    # is sent in. It is also integrated across the bins,
    # which is why there is a multiplication by the
    # grid width (for this constant model).
    #
    # Note that the filter doesn't change the grid.
    #
    de = egrid[1:] - egrid[:-1]
    expected = constant * de
    out = wrapped([4, 5])
    assert_allclose(out, expected)


@pytest.mark.parametrize("ignore", [None, 0, 1, 6, 7])
def test_rmfmodelpha_matrix_call(ignore):
    """What happens calling an rmf (matrix) with a pha?

    The ignore value gives the channel to ignore (counting from 0).
    """

    exposure = 200.1
    rdata = create_non_delta_rmf_local()
    elo = rdata.e_min
    ehi = rdata.e_max
    nchans = elo.size

    constant = 12.2
    slope = 0.01
    mdl = Polynom1D('not-flat')
    mdl.c0 = constant
    mdl.c1 = slope

    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    # force energy units (only needed if ignore is set)
    pha.set_analysis('energy')

    if ignore is not None:
        e0 = elo[ignore]
        e1 = ehi[ignore]
        de = 0.9 * (e1 - e0)
        pha.notice(lo=e0, hi=e0 + de, ignore=True)

        # The assert are intended to help people reading this
        # code rather than being a useful check that the code
        # is working.
        mask = [True] * nchans
        mask[ignore] = False
        assert (pha.mask == mask).all()

    wrapped = RMFModelPHA(rdata, pha, mdl)

    # Note that the evaluation ignores any filter we've applied.
    # and the exposure time is not used.
    #
    modvals = mdl(rdata.energ_lo, rdata.energ_hi)
    matrix = get_non_delta_matrix()
    expected = np.matmul(modvals, matrix)

    out = wrapped([4, 5])
    assert_allclose(out, expected)


@pytest.mark.parametrize("ignore", [None, 0, 1, 17, 26])
def test_rspmodelpha_delta_call(ignore):
    """What happens calling a rsp with a pha (RMF is a delta fn)?

    The ignore value gives the channel to ignore (counting from 0).
    """

    exposure = 200.1
    estep = 0.025
    egrid = np.arange(0.1, 0.8, estep)
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = 2.4 * np.ones(elo.size, dtype=np.float32)
    specresp[2:5] = 0.0
    specresp[16:19] = 3.2
    adata = create_arf(elo, ehi, specresp, exposure=exposure)
    rdata = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    nchans = elo.size

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    # force energy units (only needed if ignore is set)
    pha.set_analysis('energy')

    if ignore is not None:
        de = estep * 0.9
        e0 = egrid[ignore]
        pha.notice(lo=e0, hi=e0 + de, ignore=True)

        # The assert are intended to help people reading this
        # code rather than being a useful check that the code
        # is working.
        mask = [True] * nchans
        mask[ignore] = False
        assert (pha.mask == mask).all()

    wrapped = RSPModelPHA(adata, rdata, pha, mdl)

    # The model is evaluated on the RMF grid, not whatever
    # is sent in. It is also integrated across the bins,
    # which is why there is a multiplication by the
    # grid width (for this constant model).
    #
    # Note that the filter doesn't change the grid.
    #
    de = egrid[1:] - egrid[:-1]
    expected = constant * specresp * de
    out = wrapped([4, 5])
    assert_allclose(out, expected)


@pytest.mark.parametrize("ignore", [None, 0, 1, 5, 6, 7])
def test_rspmodelpha_matrix_call(ignore):
    """What happens calling a rsp with a pha (RMF is a matrix)?

    The ignore value gives the channel to ignore (counting from 0).
    """

    exposure = 200.1
    rdata = create_non_delta_rmf_local()
    specresp = create_non_delta_specresp()
    elo = rdata.energ_lo
    ehi = rdata.energ_hi

    adata = create_arf(elo, ehi, specresp, exposure=exposure)
    nchans = rdata.e_min.size

    constant = 22.3
    slope = -1.2
    mdl = Polynom1D('sloped')
    mdl.c0 = constant
    mdl.c1 = slope

    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    # force energy units (only needed if ignore is set)
    pha.set_analysis('energy')

    if ignore is not None:
        e0 = rdata.e_min[ignore]
        e1 = rdata.e_max[ignore]
        de = 0.9 * (e1 - e0)
        pha.notice(lo=e0, hi=e0 + de, ignore=True)

        # The assert are intended to help people reading this
        # code rather than being a useful check that the code
        # is working.
        mask = [True] * nchans
        mask[ignore] = False
        assert (pha.mask == mask).all()

    wrapped = RSPModelPHA(adata, rdata, pha, mdl)

    # The filter does not change the grid
    modvals = specresp * mdl(rdata.energ_lo, rdata.energ_hi)
    matrix = get_non_delta_matrix()
    expected = np.matmul(modvals, matrix)

    out = wrapped([4, 5])
    assert_allclose(out, expected)


def test_rspmodelpha_delta_call_wave():
    """What happens calling a rsp with a pha (RMF is a delta fn)? Wavelength.

    Unlike the energy case no bins are ignored, as this code path
    has already been tested.
    """

    exposure = 200.1
    estep = 0.025
    egrid = np.arange(0.1, 0.8, estep)
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = 2.4 * np.ones(elo.size, dtype=np.float32)
    specresp[2:5] = 0.0
    specresp[16:19] = 3.2
    adata = create_arf(elo, ehi, specresp, exposure=exposure)
    rdata = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    nchans = elo.size

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    pha.set_analysis('wave')

    wrapped = RSPModelPHA(adata, rdata, pha, mdl)

    # Note that this is a Sherpa model, so it's normalization is
    # per unit x axis, so when integrated here the bins are in
    # Angstroms, so the bin width to multiply by is
    # Angstroms, not keV.
    #
    dl = (hc / elo) - (hc / ehi)
    expected = constant * specresp * dl

    out = wrapped([4, 5])
    assert_allclose(out, expected)


def test_rspmodelpha_delta_call_channel():
    """What happens calling a rsp with a pha (RMF is a delta fn)? Channels.

    I am not convinced I understand the bin width calculation here,
    as it doesn't seem to match the wavelength case.
    """

    exposure = 200.1
    estep = 0.025
    egrid = np.arange(0.1, 0.8, estep)
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = 2.4 * np.ones(elo.size, dtype=np.float32)
    specresp[2:5] = 0.0
    specresp[16:19] = 3.2
    adata = create_arf(elo, ehi, specresp, exposure=exposure)
    rdata = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    nchans = elo.size

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    pha.set_analysis('channel')

    wrapped = RSPModelPHA(adata, rdata, pha, mdl)

    # Since this is channels you might expect the bin width to be 1,
    # but it is actually still dE.
    #
    de = ehi - elo
    expected = constant * specresp * de

    out = wrapped([4, 5])
    assert_allclose(out, expected)


@requires_xspec
def test_rspmodelpha_matrix_call_xspec():
    """Check XSPEC constant is invariant to wavelength/energy setting.

    As XSPEC models internally convert from Angstrom to keV,
    do a simple check here.
    """

    exposure = 200.1
    rdata = create_non_delta_rmf_local()
    specresp = create_non_delta_specresp()
    adata = create_arf(rdata.energ_lo,
                       rdata.energ_hi,
                       specresp,
                       exposure=exposure)

    constant = 2.3
    mdl = XSconstant('flat')
    mdl.factor = constant

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)

    # The set_arf call isn't necessary, but leave in
    pha.set_arf(adata)
    pha.set_rmf(rdata)

    # The XSPEC models are evaluated on an energy grid, even when
    # the analysis setting is wavelength. Also, unlike the Sherpa
    # Constant model, the XSPEC XSconstant model is defined
    # over the integrated bin, so no correction is needed for the
    # bin width.
    #
    modvals = constant * specresp
    matrix = get_non_delta_matrix()
    expected = np.matmul(modvals, matrix)

    wrapped = RSPModelPHA(adata, rdata, pha, mdl)

    pha.set_analysis('wave')
    out_wl = wrapped([4, 5])
    assert_allclose(out_wl, expected)

    pha.set_analysis('energy')
    out_en = wrapped([4, 5])
    assert_allclose(out_en, expected)


@pytest.mark.parametrize("analysis", ['channel', 'energy', 'wave'])
@pytest.mark.parametrize("arfexp", [True, False])
@pytest.mark.parametrize("phaexp", [True, False])
def test_rsp_matrix_call(analysis, arfexp, phaexp):
    """Check out Response1D with matrix.

    analysis is the analysis setting
    arfexp determines whether the arf has an exposure time
    phaexp determines whether the PHA has an exposure time
    """

    # Chose different exposure times for ARF and PHA to see which
    # gets picked up.
    #
    if arfexp:
        arf_exposure = 200.1
    else:
        arf_exposure = None

    if phaexp:
        pha_exposure = 220.9
    else:
        pha_exposure = None

    if phaexp:
        exposure = pha_exposure
        mdl_label = f'{exposure} * flat'
    elif arfexp:
        exposure = arf_exposure
        mdl_label = f'{exposure} * flat'
    else:
        exposure = 1.0
        mdl_label = 'flat'

    rdata = create_non_delta_rmf_local()
    specresp = create_non_delta_specresp()
    adata = create_arf(rdata.energ_lo,
                       rdata.energ_hi,
                       specresp,
                       exposure=arf_exposure)

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    # Turn off integration on this model, so that it is not integrated
    # across the bin width.
    #
    mdl.integrate = False

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=pha_exposure)

    pha.set_arf(adata)
    pha.set_rmf(rdata)

    rsp = Response1D(pha)
    wrapped = rsp(mdl)

    assert isinstance(wrapped, ArithmeticModel)

    expname = f'apply_rmf(apply_arf({mdl_label}))'
    assert wrapped.name == expname

    modvals = exposure * constant * specresp
    matrix = get_non_delta_matrix()
    expected = np.matmul(modvals, matrix)

    pha.set_analysis(analysis)
    out = wrapped([4, 5])
    assert_allclose(out, expected, rtol=2e-7)


@pytest.mark.parametrize("arfexp", [True, False])
@pytest.mark.parametrize("phaexp", [True, False])
def test_rsp_normf_call(arfexp, phaexp):
    """Check out Response1D with no RMF.

    analysis is the analysis setting
    arfexp determines whether the arf has an exposure time
    phaexp determines whether the PHA has an exposure time

    This only uses the channel setting
    """

    # Chose different exposure times for ARF and PHA to see which
    # gets picked up.
    #
    if arfexp:
        arf_exposure = 200.1
    else:
        arf_exposure = None

    if phaexp:
        pha_exposure = 220.9
    else:
        pha_exposure = None

    if phaexp:
        exposure = pha_exposure
        mdl_label = f'{exposure} * flat'
    elif arfexp:
        exposure = arf_exposure
        mdl_label = f'{exposure} * flat'
    else:
        exposure = 1.0
        mdl_label = 'flat'

    # rdata is only used to define the grids
    rdata = create_non_delta_rmf_local()
    specresp = create_non_delta_specresp()
    adata = create_arf(rdata.energ_lo,
                       rdata.energ_hi,
                       specresp,
                       exposure=arf_exposure)

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    # Turn off integration on this model, so that it is not integrated
    # across the bin width.
    #
    mdl.integrate = False

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=pha_exposure)

    pha.set_arf(adata)

    rsp = Response1D(pha)
    wrapped = rsp(mdl)

    assert isinstance(wrapped, ArithmeticModel)

    expname = f'apply_arf({mdl_label})'
    assert wrapped.name == expname

    expected = exposure * constant * specresp

    pha.set_analysis('channel')
    out = wrapped([4, 5])
    assert_allclose(out, expected)


@pytest.mark.parametrize("analysis", ['channel', 'energy', 'wave'])
@pytest.mark.parametrize("phaexp", [True, False])
def test_rsp_no_arf_matrix_call(analysis, phaexp):
    """Check out Response1D with matrix but no ARF

    analysis is the analysis setting
    arfexp determines whether the arf has an exposure time
    phaexp determines whether the PHA has an exposure time
    """

    if phaexp:
        pha_exposure = 220.9
    else:
        pha_exposure = None

    if phaexp:
        exposure = pha_exposure
        mdl_label = f'{exposure} * flat'
    else:
        exposure = 1.0
        mdl_label = 'flat'

    rdata = create_non_delta_rmf_local()

    constant = 2.3
    mdl = Const1D('flat')
    mdl.c0 = constant

    # Turn off integration on this model, so that it is not integrated
    # across the bin width.
    #
    mdl.integrate = False

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=pha_exposure)

    pha.set_rmf(rdata)

    rsp = Response1D(pha)
    wrapped = rsp(mdl)

    assert isinstance(wrapped, ArithmeticModel)

    expname = f'apply_rmf({mdl_label})'
    assert wrapped.name == expname

    modvals = exposure * constant * np.ones(rdata.energ_lo.size)
    matrix = get_non_delta_matrix()
    expected = np.matmul(modvals, matrix)

    pha.set_analysis(analysis)
    out = wrapped([4, 5])
    assert_allclose(out, expected)


class MyPowLaw1D(PowLaw1D):
    """A model that creates a NaN if the first bin starts <= 0"""

    def calc(self, pars, *args, **kwargs):
        out = PowLaw1D.calc(self, pars, *args, **kwargs)
        if args[0][0] <= 0.0:
            out[0] = np.nan
        return out


def test_arf1d_no_pha_zero_energy_bin():
    """What happens when the first bin starts at 0, no replacement

    This replicates a test in test_data.py; note that this test
    is left here, but other tests below only include the "with
    replacement" version, to avoid duplication.
    """

    exposure = 0.1
    egrid = np.asarray([0.0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = np.asarray([10.2, 9.8, 10.0, 12.0, 8.0, 10.0])
    with pytest.raises(DataErr,
                       match="The ARF 'user-arf' has an ENERG_LO value <= 0"):
        create_arf(elo, ehi, specresp, exposure=exposure)


def test_arf1d_no_pha_zero_energy_bin_replace():
    "What happens when the first bin starts at 0, with replacement"

    ethresh = 1e-5

    exposure = 0.1
    egrid = np.asarray([0.0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = np.asarray([10.2, 9.8, 10.0, 12.0, 8.0, 10.0])

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        adata = create_arf(elo, ehi, specresp, exposure=exposure,
                           ethresh=ethresh)

    validate_zero_replacement(ws, 'ARF', 'user-arf', ethresh)

    arf = ARF1D(adata)

    mdl = MyPowLaw1D()
    tmdl = PowLaw1D()

    wrapped = arf(mdl)

    out = wrapped([0.1, 0.2])

    elo[0] = ethresh
    expected = exposure * specresp * tmdl(elo, ehi)

    assert_allclose(out, expected)
    assert not np.isnan(out[0])


def test_arf1d_pha_zero_energy_bin():
    "What happens when the first bin starts at 0, with replacement"

    ethresh = 1.0e-10

    # Note: the two exposures are different to check which is
    #       used (the answer is neither, which seems surprising)
    #
    exposure1 = 0.1
    egrid = np.asarray([0.0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = np.asarray([10.2, 9.8, 10.0, 12.0, 8.0, 10.0])

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        adata = create_arf(elo, ehi, specresp, exposure=exposure1,
                           ethresh=ethresh)

    validate_zero_replacement(ws, 'ARF', 'user-arf', ethresh)

    arf = ARF1D(adata)

    exposure2 = 2.4
    channels = np.arange(1, 7, dtype=np.int16)
    counts = np.ones(6, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure2)
    pha.set_arf(adata)

    pha.set_analysis('energy')

    mdl = MyPowLaw1D()
    tmdl = PowLaw1D()

    wrapped = ARFModelPHA(arf, pha, mdl)

    out = wrapped([0.1, 0.2])
    elo[0] = ethresh
    expected = specresp * tmdl(elo, ehi)

    assert_allclose(out, expected)
    assert not np.isnan(out[0])


def test_rmf1d_delta_no_pha_zero_energy_bin():
    "What happens when the first bin starts at 0, no replacement"

    ethresh = None

    egrid = np.asarray([0.0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]

    with pytest.raises(DataErr,
                       match="The RMF 'delta-rmf' has an ENERG_LO value <= 0"):
        create_delta_rmf(elo, ehi, ethresh=ethresh)


def test_rmf1d_delta_no_pha_zero_energy_bin_replace():
    "What happens when the first bin starts at 0, with replacement"

    ethresh = 1e-8

    egrid = np.asarray([0.0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        rdata = create_delta_rmf(elo, ehi, ethresh=ethresh)

    validate_zero_replacement(ws, 'RMF', 'delta-rmf', ethresh)

    rmf = RMF1D(rdata)

    mdl = MyPowLaw1D()
    tmdl = PowLaw1D()

    wrapped = rmf(mdl)

    out = wrapped([0.1, 0.2])

    elo[0] = ethresh
    expected = tmdl(elo, ehi)

    assert_allclose(out, expected)
    assert not np.isnan(out[0])


def test_rmf1d_delta_pha_zero_energy_bin():
    "What happens when the first bin starts at 0, with replacement"

    ethresh = 2e-7

    egrid = np.asarray([0.0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        rdata = create_delta_rmf(elo, ehi, ethresh=ethresh)

    validate_zero_replacement(ws, 'RMF', 'delta-rmf', ethresh)

    exposure = 2.4
    channels = np.arange(1, 7, dtype=np.int16)
    counts = np.ones(6, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    pha.set_analysis('energy')

    mdl = MyPowLaw1D()
    tmdl = PowLaw1D()

    wrapped = RMFModelPHA(rdata, pha, mdl)

    out = wrapped([0.1, 0.2])

    elo[0] = ethresh
    expected = tmdl(elo, ehi)

    assert_allclose(out, expected)
    assert not np.isnan(out[0])


def test_rsp1d_delta_no_pha_zero_energy_bin():
    "What happens when the first bin starts at 0, with replacement"

    ethresh = 1.0e-9

    exposure = 0.1
    egrid = np.asarray([0.0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = np.asarray([10.2, 9.8, 10.0, 12.0, 8.0, 10.0])

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        adata = create_arf(elo, ehi, specresp, exposure=exposure,
                           ethresh=ethresh)

    validate_zero_replacement(ws, 'ARF', 'user-arf', ethresh)

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        rdata = create_delta_rmf(elo, ehi, ethresh=ethresh)

    validate_zero_replacement(ws, 'RMF', 'delta-rmf', ethresh)

    mdl = MyPowLaw1D()
    tmdl = PowLaw1D()

    wrapped = RSPModelNoPHA(adata, rdata, mdl)

    out = wrapped([0.1, 0.2])

    elo[0] = ethresh
    expected = specresp * tmdl(elo, ehi)

    assert_allclose(out, expected)
    assert not np.isnan(out[0])


def test_rsp1d_delta_pha_zero_energy_bin():
    "What happens when the first bin starts at 0, with replacement"

    ethresh = 2.0e-7

    # PHA and ARF have different exposure ties
    exposure1 = 0.1
    exposure2 = 2.4
    egrid = np.asarray([0.0, 0.1, 0.2, 0.4, 0.5, 0.7, 0.8])
    elo = egrid[:-1]
    ehi = egrid[1:]
    specresp = np.asarray([10.2, 9.8, 10.0, 12.0, 8.0, 10.0])

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        adata = create_arf(elo, ehi, specresp, exposure=exposure1,
                           ethresh=ethresh)

    validate_zero_replacement(ws, 'ARF', 'user-arf', ethresh)

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        rdata = create_delta_rmf(elo, ehi, ethresh=ethresh)

    validate_zero_replacement(ws, 'RMF', 'delta-rmf', ethresh)

    channels = np.arange(1, 7, dtype=np.int16)
    counts = np.ones(6, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure2)
    pha.set_rmf(rdata)
    pha.set_arf(adata)

    pha.set_analysis('energy')

    mdl = MyPowLaw1D()
    tmdl = PowLaw1D()

    wrapped = RSPModelPHA(adata, rdata, pha, mdl)

    out = wrapped([0.1, 0.2])

    elo[0] = ethresh
    expected = specresp * tmdl(elo, ehi)

    assert_allclose(out, expected)
    assert not np.isnan(out[0])


def test_rsp1d_matrix_pha_zero_energy_bin():
    """What happens when the first bin starts at 0, with replacement.

    Unlike test_rsp1d_delta_pha_zero_energy_bin this directly
    calls Response1D to create the model.
    """

    ethresh = 1.0e-5

    rdata = create_non_delta_rmf_local()

    # hack the first bin to have 0 energy
    rdata.energ_lo[0] = 0.0

    # PHA and ARF have different exposure ties
    exposure_arf = 0.1
    exposure_pha = 2.4

    specresp = create_non_delta_specresp()

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        adata = create_arf(rdata.energ_lo,
                           rdata.energ_hi,
                           specresp,
                           exposure=exposure_arf,
                           ethresh=ethresh)

    validate_zero_replacement(ws, 'ARF', 'user-arf', ethresh)

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure_pha)
    pha.set_rmf(rdata)
    pha.set_arf(adata)

    pha.set_analysis('energy')

    mdl = MyPowLaw1D()

    rsp = Response1D(pha)
    wrapped = rsp(mdl)

    # Evaluate the statistic / model. The value was calculated using
    # commit a65fb94004664eab219cc09652172ffe1dad80a6 on a linux
    # system (Ubuntu 17.04).
    #
    f = Fit(pha, wrapped)
    ans = f.calc_stat()
    assert ans == pytest.approx(37971.8716151947)


@requires_data
@requires_fits
def test_create_rmf(make_data_path):
    """It would have been nice to have the original RMF used to create
    the image, to compare against.

    Running 'dmlist test_rmfimg.fits rmfimg' returns

    rmfimg infile="acisf04938_000N002_r0043_rmf3.fits" outfile="test_rmfimg.fits" arf="" arfout="" product="no" verbose="0" clobber="no"

    From CSC 1.1 we get acisf04938_000N003_r0043_arf/pha/rmf3.fits.gz
    and running rmfimg on this gives the same pixel values as
    test_rmfmg.fits. Using this RMF to convolve a constant (c0=100)
    gives the checks against the out array below.

    """

    # The energy grid is 0.3-0.31, 0.31-0.32, .., 9.29-9.3,
    # which is 900 elements, but there are 1024 channels.
    #
    energ = np.arange(0.3, 9.301, 0.01)
    rmflo = energ[:-1]
    rmfhi = energ[1:]

    # We do not know the original grid, but could guess as it's from a
    # ACIS RMF. However, the exact values do not matter here.
    #
    ebounds = np.linspace(0.1, 12, 1025)
    elo = ebounds[:-1]
    ehi = ebounds[1:]

    # The header has NUMGRP=1039 and NUMELT=380384 so we can check these.
    #
    fname = make_data_path('test_rmfimg.fits')
    datarmf = create_non_delta_rmf(rmflo, rmfhi, fname, e_min=elo, e_max=ehi)
    assert len(datarmf.f_chan) == 1039
    assert len(datarmf.n_chan) == 1039
    assert len(datarmf.n_grp) == 900
    assert len(datarmf.matrix) == 380384
    # How much use is the minimum value check?
    assert datarmf.matrix.min() == pytest.approx(1.0000736e-06)
    assert datarmf.matrix.max() == pytest.approx(0.23166645)
    assert datarmf.matrix.sum() == pytest.approx(900.14188997)

    assert datarmf.detchans == 1024

    # Apply the model to a constant model.
    #
    wrapper = RMF1D(datarmf)
    orig = Const1D("orig")
    orig.c0 = 100
    conv = wrapper(orig)
    out = conv(np.arange(1, 1025, dtype=np.int16))

    # Compare to acisf04938_000N003_r0043_rmf3.fits.gz
    assert len(out) == 1024
    assert out.min() == 0.0
    assert out.max() == pytest.approx(1.6668367662320187)
    assert out.argmax() == 118
    assert out.sum() == pytest.approx(900.1419080498199)
    assert out[5] == pytest.approx(0.30729822633750686)
    assert out[600] == pytest.approx(1.4288386887490545)
    assert out[700] == pytest.approx(0.0005906268168658253)
    assert (out[765:] == 0).all()


@pytest.mark.parametrize("kwargs,msg",
                         [({"e_min": [1, 2], "e_max": None},
                           "e_min/max must both be set or empty"),
                          ({"e_min": None, "e_max": [1, 2]},
                           "e_min/max must both be set or empty"),
                          ({"e_min": [1, 2], "e_max": [2, 3, 4]},
                           "e_min/max mismatch in size: 2 vs 3"),
                          ({"e_min": [1, 2, 3], "e_max": [2, 3, 4]},
                           "detchans mismatch with e_min/max: 2 vs 3")
                          ])
def test_rmf_creation_fails_approximate_energy_bounds(kwargs, msg):
    """Check we error out when a validation step fails."""

    elo = np.asarray([1, 2])
    ehi = np.asarray([2, 3])
    n_grp = np.asarray([1, 1])
    f_chan = np.asarray([1, 2])
    n_chan = np.asarray([1, 1])
    matrix = np.asarray([1, 1])
    with pytest.raises(DataErr, match=f"^{msg}$"):
        DataRMF("test", 2, elo, ehi, n_grp, f_chan, n_chan, matrix, **kwargs)


@requires_data
@requires_fits
def test_recreate_img_from_rmf(make_data_path):
    """The inverse of test_create_rmf"""

    fname = make_data_path('test_rmfimg.fits')
    energ = np.arange(0.3, 9.301, 0.01)
    rmflo = energ[:-1]
    rmfhi = energ[1:]

    # fake up an energy range
    #
    # nchan = 1024  TODO it should be this, see #1889
    nchan = 900
    ebounds = np.linspace(0.1, 12, nchan + 1)
    elo = ebounds[:-1]
    ehi = ebounds[1:]

    rmf = create_non_delta_rmf(rmflo, rmfhi, fname,
                               e_min=elo, e_max=ehi)

    # What is the input image?
    #
    iblock, _ = io.backend.get_image_data(fname)
    expected = iblock.image
    assert expected.shape == (900, 1024)

    # Can we reconstruct the image in test_rmfimg.fits?
    #
    mat = rmf_to_matrix(rmf)
    assert isinstance(mat, RMFMatrix)

    matrix = mat.matrix
    xaxis = mat.channels
    yaxis = mat.energies
    assert not xaxis.is_integrated
    assert yaxis.is_integrated
    assert xaxis.x_axis.size == nchan
    assert yaxis.x_axis.size == 900
    assert xaxis.start == 1
    assert xaxis.end == nchan
    assert yaxis.start == pytest.approx(0.3)
    assert yaxis.end == pytest.approx(9.3)

    # Need to subset expected because of #1889
    assert matrix.shape == (900, nchan)
    # assert matrix == pytest.approx(expected)
    assert matrix == pytest.approx(expected[:, 0:900])


def check_rmf_as_image(infile, tmpfile,  nchan, nenergy, offset):
    """Read in RMF, convert to matrix and image. Check things.

    We round-trip the data: convert RMF to image, write out,
    read back in, and then compare the two RMFs by evaluating
    a model.
    """

    rmf = io.read_rmf(infile)

    mat = rmf_to_matrix(rmf)
    assert mat.matrix.shape == (nenergy, nchan)
    assert mat.channels.x_axis.size == nchan
    assert mat.energies.x_axis.size == nenergy
    assert mat.channels.start == offset
    assert mat.channels.end == (nchan - 1 + offset)

    img = rmf_to_image(rmf)
    assert isinstance(img, DataIMG)

    # Read the imge back to check the roundtrip behaviour.
    #
    io.write_image(tmpfile, img, ascii=False)

    elo, ehi = mat.energies.grid
    nrmf = create_non_delta_rmf(elo, ehi, tmpfile,
                                e_min=rmf.e_min, e_max=rmf.e_max)

    # Convolve a simple model with the two RMFs and check the result
    # is the same.
    #
    r1 = RMF1D(rmf)
    r2 = RMF1D(nrmf)
    mdl = PowLaw1D()
    mdl.gamma = 1.7
    mdl.ampl = 1000

    c1 = r1(mdl)
    c2 = r2(mdl)
    chans = mat.channels.grid
    assert c2(chans) == pytest.approx(c1(chans))


@requires_data
@requires_fits
@pytest.mark.parametrize("rmfname,nchan,nenergy,offset",
                         [("swxpc0to12s6_20130101v014.rmf.gz",
                           1024, 2400, 0),  # SWIFT
                          ("sis0.rmf", 1024, 1180, 1),  # ASCA
                          ("chandra_hrci/hrcf24564_000N030_r0001_rmf3.fits.gz",
                           1024, 16921, 1),  # CHANDRA/HRC-I
                          ("xmmrgs/P0112880201R1S004RSPMAT1003.FTZ",
                           3600, 4000, 1)  # XMM/RGS
                          ])

def test_rmf_to_matrix_mission(rmfname, nchan, nenergy, offset,
                               make_data_path, tmp_path, recwarn):
    """Check other missions for RMF

    Note that the SWIFT RMF is essentially stored as an image,
    i.e. each row is stored completely. It also starts at column 0
    and triggers a user warning when read in.

    The ASCA RMF has N_GRP=0 or 1.

    The Chandra HRC-I RMF is used because it has caused problems for
    the pyfits backend in the past.

    XMM/RGS has EBOUNDS in decreasing order.
    """

    infile = make_data_path(rmfname)
    outpath = tmp_path / "fake.img"
    outfile = str(outpath)
    check_rmf_as_image(infile, outfile, nchan, nenergy, offset)

    # We do not care about any warnings here, so clear the state so
    # that the capture_all_warnings fixture is not triggered.
    #
    recwarn.clear()


@requires_data
@requires_fits
def test_rmf_image_offset_0(make_data_path, tmp_path, recwarn):
    """The SWIFT RMF is stored almost as an image with offset=0.

    Use this to validate the "image" checks.
    """

    infile = make_data_path("swxpc0to12s6_20130101v014.rmf.gz")

    # Read the file in directly
    #
    rmf_direct = io.read_rmf(infile)

    # Messy way to get the MATRIX data.
    #
    _, blocks, _ = io.read_table_blocks(infile)
    matrix = blocks[2]["MATRIX"]
    assert matrix.shape == (2400, 1024)
    assert matrix.max() == pytest.approx(0.058343917)

    elo = blocks[2]["ENERG_LO"]
    ehi = blocks[2]["ENERG_HI"]
    e_min = blocks[3]["E_MIN"]
    e_max = blocks[3]["E_MAX"]

    # Create a DataIMG with this data and write it out, as
    # sherpa.astro.instrument.create_non_delta_rmf does not have a
    # "read the data from an array" option.
    #
    x1, x0 = np.mgrid[1:2401, 1:1025]
    img = DataIMG("rmfimg", x0.flatten(), x1.flatten(),
                  matrix.flatten(), shape=(2400, 1024))

    outpath = tmp_path / "rmf.img"
    outfile = str(outpath)
    io.write_image(outfile, img, ascii=False)

    # Create a RMF with this image.
    #
    rmf_image = create_non_delta_rmf(elo, ehi, outfile, offset=0,
                                     e_min=e_min, e_max=e_max,
                                     ethresh=1e-10)

    # A couple of simple checks. We can not directly compare
    # n_grp/n_chan/f_chan/matrix values. We can check some basic
    # matrix properties, such as the max value and the summation
    # (assuming the tolerance filter does not significantly change the
    # results).
    #
    assert rmf_image.detchans == rmf_direct.detchans
    assert rmf_image.offset == rmf_direct.offset

    matrix_direct = rmf_direct.matrix
    matrix_image = rmf_image.matrix
    assert matrix_image.max() == pytest.approx(matrix_direct.max())
    assert matrix_image.sum() == pytest.approx(matrix_direct.sum())

    # Pass a simple model through the full grid.
    #
    mvals = np.linspace(2, 10, 2400)

    y_direct = rmf_direct.apply_rmf(mvals)
    y_image = rmf_image.apply_rmf(mvals)
    assert y_image == pytest.approx(y_direct)

    # Now filter the RMF
    #
    chans = np.arange(200, 500, dtype=np.int16)
    selected_direct = rmf_direct.notice(chans)
    selected_image = rmf_image.notice(chans)

    # These filters are not the same, for some reason, but let's add
    # basic checks.
    #
    assert len(selected_direct) == 2400
    assert len(selected_image) == 2400

    assert selected_direct.all()

    # It is not clear what these values should be, so treat as a
    # regression test.
    where, = np.where(selected_image)
    expected = [146, 178, 184, 185, 186, 187, 188, 189, 190, 191]
    assert where[:10] == pytest.approx(expected)

    y_direct = rmf_direct.apply_rmf(mvals[selected_direct])
    y_image = rmf_image.apply_rmf(mvals[selected_image])

    assert len(y_direct) == 1024
    assert len(y_image) == 1024

    # Because the selected channel range is different the results will
    # be different, but we can check the position and value of the max
    # values as, in this case, they should be the same.
    #
    idx_direct = y_direct.argmax()
    idx_image = y_image.argmax()
    assert idx_image == idx_direct

    assert y_image[idx_image] == pytest.approx(y_direct[idx_direct])

    # We do not care about any warnings here, so clear the state so
    # that the capture_all_warnings fixture is not triggered.
    #
    recwarn.clear()


@requires_data
@requires_fits
def test_rmf_to_matrix_offset_0(make_data_path, recwarn):
    """The SWIFT RMF is stored almost as an image with offset=0."""

    infile = make_data_path("swxpc0to12s6_20130101v014.rmf.gz")
    rmf = io.read_rmf(infile)
    recwarn.clear()

    mat = rmf_to_matrix(rmf)
    assert mat.matrix.shape == (2400, 1024)
    assert mat.matrix.max() == pytest.approx(rmf.matrix.max())
    assert mat.matrix.sum() == pytest.approx(rmf.matrix.sum())

    assert len(mat.channels.grid) == 1
    assert mat.channels.grid[0] == pytest.approx(np.arange(0, 1024))

    assert len(mat.energies.grid) == 2
    assert mat.energies.grid[0] == pytest.approx(rmf.energ_lo)
    assert mat.energies.grid[1] == pytest.approx(rmf.energ_hi)


@requires_data
@requires_fits
def test_rmf_to_image_offset_0(make_data_path, recwarn):
    """The SWIFT RMF is stored almost as an image with offset=0."""

    infile = make_data_path("swxpc0to12s6_20130101v014.rmf.gz")
    rmf = io.read_rmf(infile)
    recwarn.clear()

    img = rmf_to_image(rmf)
    assert img.shape == (2400, 1024)

    assert img.y.max() == pytest.approx(rmf.matrix.max())
    assert img.y.sum() == pytest.approx(rmf.matrix.sum())

    # Note this starts at 1, not 0.
    #
    assert img.x0[0:10] == pytest.approx([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    assert img.x1[0:10] == pytest.approx([1] * 10)


# Several tests from sherpa/tests/test_instrument.py repeated to check out
# the astro-version of PSFModel
#
def test_psf1d_empty_pars():
    """What does .pars mean for an empty PSFModel?"""

    m = PSFModel()
    assert m.pars == ()


def test_psf1d_pars():
    """What does .pars mean for a PSFModel?"""

    b = Box1D()
    m = PSFModel(kernel=b)
    assert m.pars == ()


def test_psf1d_convolved_pars():
    """What does .pars mean for a PSFModel applied to a model?"""

    b1 = Box1D('b1')
    m = PSFModel(kernel=b1)
    b2 = Box1D('b2')
    c = m(b2)

    b1.xlow = 1
    b1.xhi = 10
    b1.ampl = 0.2

    b2.xlow = 4
    b2.xhi = 8
    b2.ampl = 0.4

    bpars = b1.pars + b2.pars
    assert len(bpars) == 6

    cpars = c.pars
    assert len(cpars) == 6
    for bpar, cpar in zip(bpars, cpars):
        assert cpar == bpar


def test_psfmodel_kernel_has_no_dimension():
    """It is expected that this will error out, but just where.

    It was intended to catch a different error condition but
    it is hard to trigger without writing very low-level test
    code, so just check what happens here.
    """

    x = np.arange(0.5, 10.5, 1)
    y = np.ones_like(x)
    data = Data1D("data-data", x, y)

    m = PSFModel(kernel=TableModel())
    with pytest.raises(PSFErr,
                       match="PSF model dimension must be <= 2"):
        m.get_kernel(data)


def test_psfmodel_fold_check_kernel_no_cdelt_warning(caplog):
    """Check behavior"""

    # The actual kernel does not matter, but do something that
    # is vaguely PSF-like.
    #
    kx1, kx0 = np.mgrid[1:4, 1:3]
    kshape = kx0.shape
    kx0 = kx0.flatten()
    kx1 = kx1.flatten()
    ky = 1 / np.abs((kx0 - 0.4) + (kx1 - 1.9))
    kdata = DataIMG("kernel-data", kx0, kx1, ky, kshape)

    x1, x0 = np.mgrid[-4:5, -5:-4]
    shape = x0.shape
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.ones_like(x0)
    data = DataIMG("data-data", x0, x1, y, shape)

    m = PSFModel(kernel=kdata)

    # The message is via a warning, not the logger.
    #
    assert len(caplog.records) == 0
    with pytest.warns(UserWarning,
                      match="^PSF Image does not have a pixel size. " +
                      "Sherpa will assume the pixel size is the same as the data$"):
        m.fold(data)

    assert len(caplog.records) == 0


@dataclass
class FakeWCS:
    cdelt: float = 1.0


def test_psfmodel_fold_check_data_no_cdelt_warning(caplog):
    """Check behavior"""

    # The actual kernel does not matter, but do something that
    # is vaguely PSF-like.
    #
    kx1, kx0 = np.mgrid[1:4, 1:3]
    kshape = kx0.shape
    kx0 = kx0.flatten()
    kx1 = kx1.flatten()
    ky = 1 / np.abs((kx0 - 0.4) + (kx1 - 1.9))
    kdata = DataIMG("kernel-data", kx0, kx1, ky, kshape)

    # For this we just need a kdata.sky.cdelt setting, we do not need
    # to bother with a full WCS object (at least at present; we could
    # require it in the future).
    #
    kdata.sky = FakeWCS()

    x1, x0 = np.mgrid[-4:5, -5:-4]
    shape = x0.shape
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.ones_like(x0)
    data = DataIMG("data-data", x0, x1, y, shape)

    m = PSFModel(kernel=kdata)

    # The message is via a warning, not the logger.
    #
    assert len(caplog.records) == 0
    with pytest.warns(UserWarning,
                      match="^Data Image does not have a pixel size. " +
                      "Sherpa will assume the pixel size is the same as the PSF$"):
        m.fold(data)

    assert len(caplog.records) == 0


def test_has_pha_response():
    """Check the examples from the docstring"""

    exposure = 200.1
    rdata = create_non_delta_rmf_local()
    specresp = create_non_delta_specresp()
    adata = create_arf(rdata.energ_lo,
                       rdata.energ_hi,
                       specresp,
                       exposure=exposure)

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)

    pha.set_arf(adata)
    pha.set_rmf(rdata)

    rsp = Response1D(pha)
    m1 = Gauss1D()
    m2 = PowLaw1D()

    assert not has_pha_response(m1)
    assert has_pha_response(rsp(m1))
    assert not has_pha_response(m1 + m2)
    assert has_pha_response(rsp(m1 + m2))
    assert has_pha_response(m1 + rsp(m2))

    # reflexivity check
    assert has_pha_response(rsp(m1) + m2)
    assert has_pha_response(rsp(m1) + rsp(m2))
