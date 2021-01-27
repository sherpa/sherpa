#
#  Copyright (C) 2017, 2020, 2021  Smithsonian Astrophysical Observatory
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

TODO:

  - add grouping to the PHA
  - ARF and RMF can have different grids; I believe this is supported
  - no testing of the startup/teardown methods; should this be done here?
  - no testing of the nested models (also not fully tested in the full suite?)
  - no testing of pileup
  - no testing of PSFModel (also not fully tested in the full suite?)

"""

import logging
import warnings

import numpy as np
from numpy.testing import assert_allclose

import pytest  # pytest >= 3.0 is needed

from sherpa.models.model import ArithmeticModel, \
    ArithmeticConstantModel, BinaryOpModel
from sherpa.astro.instrument import ARF1D, ARFModelNoPHA, ARFModelPHA, \
    Response1D, RMF1D, RMFModelNoPHA, RMFModelPHA, \
    RSPModelNoPHA, RSPModelPHA, create_arf, create_delta_rmf
from sherpa.fit import Fit
from sherpa.astro.data import DataPHA, DataRMF
from sherpa.models.basic import Const1D, Polynom1D, PowLaw1D
from sherpa.utils.err import DataErr
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

    emsg = "The minimum ENERG_LO in the {} '{}' ".format(rtype, label) + \
           "was 0 and has been replaced by {}".format(ethresh)
    assert str(w.message) == emsg


def get_non_delta_matrix():
    """Return a 2D matrix representing the RMF in create_non_delta_rmf."""

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


def create_non_delta_rmf():
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
    """A SPECRESP column compatible with create_non_delta_rmf.

    Returns
    -------
    specresp : array
    """

    return np.asarray([20.0, 10.0, 0.0, 10.0, 15.0,
                       20.0, 20.0, 20.0, 20.0, 20.0,
                       18.0, 17.0, 0.0, 20.0, 20.0,
                       15.0, 0.0, 10.0, 20.0, 20.0],
                      dtype=np.float32)


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

    with pytest.raises(DataErr) as exc:
        Response1D(pha)

    emsg = 'No instrument response found for dataset test-pha'
    assert str(exc.value) == emsg


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

    rdata = create_non_delta_rmf()
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

    rdata = create_non_delta_rmf()
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
    rdata = create_non_delta_rmf()
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

    rdata = create_non_delta_rmf()

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

    rdata = create_non_delta_rmf()
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

    Ensure we can't filter on energy or wavelength since there's no
    EBOUNDS information. This behavior was seen when writing
    test_rmfmodelpha_call, so a test was written for it.

    The code used to raise a DataErr but now just displays a
    logged warning.
    """

    estep = 0.01
    egrid = np.arange(0.01, 0.06, estep)
    rdata = create_delta_rmf(egrid[:-1], egrid[1:])

    channels = np.arange(1, 5, dtype=np.int16)
    counts = np.asarray([10, 5, 12, 7], dtype=np.int16)
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
    rdata = create_non_delta_rmf()

    # nchans should be rdata.e_min.size for the sizes to match
    nchans = rdata.energ_lo.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=exposure)
    pha.set_rmf(rdata)

    with pytest.raises(DataErr) as exc:
        pha.set_analysis(analysis)

    emsg = "RMF 'non-delta-rmf' is incompatible with PHA dataset 'test-pha'"
    assert str(exc.value) == emsg


@pytest.mark.parametrize("analysis", ["energy", "wave"])
def test_rsp_normf_error(analysis):
    """Check that an error is raised on set_analysis

    """

    exposure = 200.1

    # rdata is only used to define the grids
    rdata = create_non_delta_rmf()
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

    with pytest.raises(DataErr) as exc:
        pha.set_analysis(analysis)

    emsg = "response incomplete for dataset test-pha, " + \
           "check the instrument model"
    assert str(exc.value) == emsg


def test_rsp_norsp_error():
    """Check that an error is raised when creating a wrapped model

    """

    # rdata is only used to define the grids
    rdata = create_non_delta_rmf()

    nchans = rdata.e_min.size
    channels = np.arange(1, nchans + 1, dtype=np.int16)
    counts = np.ones(nchans, dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts)

    with pytest.raises(DataErr) as exc:
        Response1D(pha)

    emsg = "No instrument response found for dataset test-pha"
    assert str(exc.value) == emsg


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
    rdata = create_non_delta_rmf()
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
    rdata = create_non_delta_rmf()
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
    dl = (DataPHA._hc / elo) - (DataPHA._hc / ehi)
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
    rdata = create_non_delta_rmf()
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
        mdl_label = '({} * flat)'.format(exposure)
    elif arfexp:
        exposure = arf_exposure
        mdl_label = '({} * flat)'.format(exposure)
    else:
        exposure = 1.0
        mdl_label = 'flat'

    rdata = create_non_delta_rmf()
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

    expname = 'apply_rmf(apply_arf({}))'.format(mdl_label)
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
        mdl_label = '({} * flat)'.format(exposure)
    elif arfexp:
        exposure = arf_exposure
        mdl_label = '({} * flat)'.format(exposure)
    else:
        exposure = 1.0
        mdl_label = 'flat'

    # rdata is only used to define the grids
    rdata = create_non_delta_rmf()
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

    expname = 'apply_arf({})'.format(mdl_label)
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
        mdl_label = '({} * flat)'.format(exposure)
    else:
        exposure = 1.0
        mdl_label = 'flat'

    rdata = create_non_delta_rmf()

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

    expname = 'apply_rmf({})'.format(mdl_label)
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
    with pytest.raises(DataErr) as exc:
        create_arf(elo, ehi, specresp, exposure=exposure)

    emsg = "The ARF 'user-arf' has an ENERG_LO value <= 0"
    assert str(exc.value) == emsg


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

    with pytest.raises(DataErr) as exc:
        create_delta_rmf(elo, ehi, ethresh=ethresh)

    emsg = "The RMF 'delta-rmf' has an ENERG_LO value <= 0"
    assert str(exc.value) == emsg


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

    rdata = create_non_delta_rmf()

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
    from sherpa.astro.ui.utils import Session
    ui = Session()
    energ = np.arange(0.05, 1.1, 0.05)
    rmflo = energ[:-1]
    rmfhi = energ[1:]
    fname = make_data_path('test_rmfimg.fits')
    datarmf = ui.create_rmf(rmflo, rmfhi, fname=fname)
    assert len(datarmf._fch) == 1039
    assert len(datarmf._nch) == 1039
    assert len(datarmf.n_grp) == 900
    assert datarmf._rsp.shape[0] == 380384
