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

"""
Test handling of area/background scaling in PHA data sets.

A PHA file can have AREASCAL and BACKSCAL values (either a scalar or
column). The BACKSCAL column has tests in other parts of the system,
so this file concentrates on the area scaling column.

The handling of AREASCAL is scattered throughout the system: it's
not just the DataPHA object but also the ARF/RMF/RSP PHA instrument
objects.

The tests here focus on the object API, rather than the UI layer
(which, at present, doesn't do much extra with regard to area
scaling, so this should catch most problems).

This does *not* test any plots (e.g. to ensure that the correctly-scaled
data and/or model values are shown).

"""

import pytest

import numpy as np
from numpy.testing import assert_allclose

from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.astro.instrument import ARFModelPHA, RMFModelPHA, RSPModelPHA
from sherpa.models.basic import Const1D, StepHi1D
from sherpa.stats import Chi2DataVar, CStat
from sherpa.utils.err import ArgumentErr, DataErr

from sherpa.utils import requires_data, requires_fits

from sherpa.astro import ui


def expected_basic_areascal():
    """Return the expected areascal values."""

    areascal = np.ones(10)
    areascal[0] = 0.0
    areascal[6:] = 0.5
    return areascal.copy()


def expected_basic_areascal_bgnd():
    """Return the expected areascal values for the background."""

    return np.asarray([0.0, 2.0, 2.0, 2.0, 2.0,
                       4.0, 4.0, 4.0, 4.0, 4.0]).copy()


def expected_basic_counts(scale=False):
    """Return the expected count values.

    Parameters
    ----------
    scale : bool, optional
        If True then the counts are scaled by areascal, otherwise
        (when False) it is the value that should be stored in the
        DataPHA object.

    Notes
    -----
    The idea is that there's two constant regions, with values
    of 20 and 40: i.e.

       Const1D + StepHi1D

    where const1d.c0 = 20, stephi1d.xcut = 5.5,stephi1d.ampl = 20

    However, as areascal is 0.5 for the last 4 bins, the counts
    are all close to 20 (it is only after applying areascal that
    the step is obvious).
    """

    counts = np.asarray([0, 12, 21, 25, 18, 24, 23, 16, 20, 19],
                        dtype=np.int16)

    if scale:
        ascal = expected_basic_areascal()

        # the first channel has areascal=0, so leave (it is 0).
        counts[1:] = counts[1:] / ascal[1:]

        counts = counts.astype(np.int16)

    return counts.copy()


def expected_basic_counts_bgnd(scale=False):
    """Return the expected count values for the background.

    Parameters
    ----------
    scale : bool, optional
        If True then the result is scaled by the source AREASCAL
        value.

    Notes
    -----
    There's three regions: the first bin, which is bad so has
    an area scaling of 0, then the next four bins with a value
    of 2.0 and then the remaining 5 bins with a value of 4.0.
    The aim is to have a slightly-different structure to the
    "source" data set, and to use values larger than 1. The
    background level is constant, at 10, before any scaling;
    this means that the scaled values are not always integer
    values.
    """

    counts = np.zeros(10, dtype=np.int16) + 10
    counts[0] = 0
    if scale:
        ascal = expected_basic_areascal_bgnd()

        # the first channel has areascal=0, so leave (it is 0).
        counts = counts.astype(np.float64)
        counts[1:] = counts[1:] / ascal[1:]

        # counts = counts.astype(np.int16)

    return counts.copy()


def expected_basic_chisquare_errors(scale=False):
    """Return the expected error values (chi square).

    Parameters
    ----------
    scale : bool, optional
        If True then the results are scaled by areascal.

    Notes
    ------
    The calculation assumes Chi2DataVar statistic and
    ignore_bad() is assumed to have been called.
    """

    # Calculate the errors based on the counts. There are no
    # counts less than 1, so this is just the square root of
    # the observed data value, which then has to be scaled by
    # the area scaling.
    #
    # Since ignore_bad has been called we ignore the first bin.
    counts = expected_basic_counts(scale=False)[1:]
    expected = np.sqrt(counts)

    if scale:
        ascal = expected_basic_areascal()[1:]
        expected /= ascal

    return expected.copy()


def expected_basic_chisquare_errors_bgnd(scale=False):
    """Return the expected error values (chi square) after bg subtraction.

    Parameters
    ----------
    scale : bool, optional
        If True then the results are scaled by areascal.

    Notes
    -----
    The calculation assumes the Chi2DataVar statistic and
    ignore_bad() is assumed to have been called.

    It is not at all clear if this code is doing the right thing
    at the moment.
    """

    # Calculate the errors based on the counts. There are no
    # counts less than 1, so this is just the square root of
    # the observed data value, which then has to be scaled by
    # the area scaling.
    #
    # Since ignore_bad has been called we ignore the first bin.
    counts = expected_basic_counts(scale=False)[1:]
    ascal = expected_basic_areascal()[1:]

    bgcounts = expected_basic_counts_bgnd(scale=False)[1:]
    bgascal = expected_basic_areascal_bgnd()[1:]

    # It is unclear what to do when scale=False, since this
    # is saying "let's ignore the areascal", but we still
    # have to include the background areascal term when calculating
    # the background contribution to the area.
    #
    # This mess suggests the abstraction is not ideal.
    #
    if scale:
        expected = counts / (ascal * ascal)
    else:
        expected = counts

    if scale:
        expected += bgcounts / (bgascal * bgascal)
    else:
        expected = expected + bgcounts * ascal * ascal / (bgascal * bgascal)

    return np.sqrt(expected)


def expected_basic_chisquare_errors_scaling_bgnd(scale=False):
    """Return the expected error values (chi square) after bg subtraction.

    Parameters
    ----------
    scale : bool, optional
        If True then the results are scaled by areascal.
        It is unclear what this should be doing.

    Notes
    -----
    The errors are calculated using the Chi2DataVar algorithm.

    It assumes:
       source exposure = 20
       background exposure = 15

       source backscal = 0.01
       background backscal = 0.1

    ignore_bad() is assumed to have been called.

    It is not at all clear if this code is doing the right thing
    at the moment.
    """

    # Calculate the errors based on the counts. There are no
    # counts less than 1, so this is just the square root of
    # the observed data value, which then has to be scaled by
    # the exposure, backscal, and area scaling.
    #
    # Since ignore_bad has been called we ignore the first bin.
    counts = expected_basic_counts(scale=False)[1:]

    exposure_src = 20.0
    exposure_bg = 15.0

    backscal_src = 0.01
    backscal_bg = 0.1

    ascal = expected_basic_areascal()[1:]

    bgcounts = expected_basic_counts_bgnd(scale=False)[1:]
    bgascal = expected_basic_areascal_bgnd()[1:]

    # see comments in expected_basic_chisquare_errors_bgnd
    # about confusion of what to do when scale = False
    #
    if scale:
        expected = counts / (ascal * ascal)
    else:
        expected = 1.0 * counts

    # The background counts are re-scaled to match the source
    # settings.
    bgdenom = bgascal * backscal_bg * exposure_bg / \
        (backscal_src * exposure_src)

    if not scale:
        bgdenom /= ascal

    expected += bgcounts / (bgdenom * bgdenom)

    return np.sqrt(expected)


def setup_basic_dataset():
    """Create a basic PHA data set with an AREASCAL value.

    Returns
    -------
    dataset : sherpa.astro.data.DataPHA instance
        The first channel has non-zero quality, but is not
        masked out (so the caller needs to call ignore_bad
        to ignore it).
    """

    channels = np.arange(1, 11)
    counts = expected_basic_counts(scale=False)

    quality = np.zeros(10, dtype=np.int16)
    quality[0] = 1

    areascal = expected_basic_areascal()

    return DataPHA('test', channel=channels, counts=counts,
                   quality=quality, areascal=areascal)


def setup_basic_dataset_bgnd():
    """Create a basic PHA data set with an AREASCAL value and a background.

    Returns
    -------
    dataset : sherpa.astro.data.DataPHA instance
        The first channel has non-zero quality, but is not
        masked out (so the caller needs to call ignore_bad
        to ignore it).
    """

    dset = setup_basic_dataset()

    channels = np.arange(1, 11)
    counts = expected_basic_counts_bgnd(scale=False)

    quality = np.zeros(10, dtype=np.int16)
    quality[0] = 1

    areascal = expected_basic_areascal_bgnd()

    bset = DataPHA('testbg', channel=channels, counts=counts,
                   quality=quality, areascal=areascal)

    dset.set_background(bset)
    return dset


# The first few tests are really of the DataPHA class, and so
# should not be necessary here, but they are included just in
# case.
#
def test_analysis_is_channel():
    """There's no response, so we have to be using channels."""

    dset = setup_basic_dataset()
    assert dset.get_analysis() == 'channel'


def test_counts_is_set():
    """Is the counts column set correctly?"""

    dset = setup_basic_dataset()
    expected = expected_basic_counts(scale=False)
    assert_allclose(dset.counts, expected)


def test_areascal_is_set():
    """Is the areascal column set correctly?"""

    dset = setup_basic_dataset()
    expected = expected_basic_areascal()
    assert_allclose(dset.areascal, expected)


def test_staterror_is_not_set():
    """Is the staterror column not set?"""

    dset = setup_basic_dataset()
    assert dset.staterror is None


def test_get_dep():
    """What does get_dep return"""

    dset = setup_basic_dataset()
    expected = expected_basic_counts(scale=False)
    expected = expected.astype(np.float64)

    assert_allclose(dset.get_dep(), expected)


def test_get_y():
    """What does get_y return"""

    dset = setup_basic_dataset()
    expected = expected_basic_counts(scale=True)
    expected = expected.astype(np.float64)

    assert_allclose(dset.get_y(), expected)


def test_get_staterror():
    """What does get_staterror return?

    This uses the data-variance calculation.
    """

    dset = setup_basic_dataset()
    dset.ignore_bad()

    stat = Chi2DataVar()
    errors = dset.get_staterror(filter=True,
                                staterrfunc=stat.calc_staterror)

    expected = expected_basic_chisquare_errors(scale=False)
    assert_allclose(errors, expected)


def test_get_yerr():
    """What does get_yerr return?

    This uses the data-variance calculation.
    """

    dset = setup_basic_dataset()
    dset.ignore_bad()

    stat = Chi2DataVar()
    errors = dset.get_yerr(filter=True,
                           staterrfunc=stat.calc_staterror)

    expected = expected_basic_chisquare_errors(scale=True)
    assert_allclose(errors, expected)


# Ensure that the interfaces used by the statistic object
# are behaving correctly.
#
def test_chisquare():
    """Is the chi square correct?

    This uses the data-variance calculation.
    """

    dset = setup_basic_dataset()
    dset.ignore_bad()

    cpt1 = Const1D()
    cpt2 = StepHi1D()
    cpt1.c0 = 20
    cpt2.ampl = 20
    cpt2.xcut = 6.5
    mdl = cpt1 + cpt2

    # Since the model does not contain a *PHA instrument model
    # it will not include the area-scaling, so the data and
    # errors should not be scaled.
    #
    scale = False
    counts = expected_basic_counts(scale=scale)[1:]
    errors = expected_basic_chisquare_errors(scale=scale)
    mvals = mdl(dset.channel[1:])

    expected = (counts - mvals)**2 / (errors**2)
    expected = expected.sum()

    stat = Chi2DataVar()
    sval = stat.calc_stat(dset, mdl)

    assert_allclose(sval[0], expected)


def make_arf(energ_lo, energ_hi, specresp=None, exposure=1.0,
             name='arf'):
    """A simple in-memory representation of an ARF.

    Parameters
    ----------
    energ_lo, energ_hi : array
        The energy grid over which the ARF is defined. The units are
        keV and each bin has energ_hi > energ_lo. The arrays are
        assumed to be ordered, but it is not clear yet whether they
        have to be in ascending order.
    specresp : array or None, optional
        The spectral response (effective area) for each bin, in cm^2.
        If not given then a value of 1.0 per bin is used.
    exposure : number, optional
        The exposure time, in seconds. It must be positive.
    name : str, optional
        The name to give to the ARF instance.

    Returns
    -------
    arf : sherpa.astro.data.DataARF
        The ARF.
    """

    elo = np.asarray(energ_lo)
    ehi = np.asarray(energ_hi)
    if elo.size != ehi.size:
        raise DataErr('mismatch', 'energ_lo', 'energ_hi')

    if specresp is None:
            specresp = np.ones(elo.size, dtype=np.float32)
    else:
        specresp = np.asarray(specresp)
        if specresp.size != elo.size:
            raise DataErr('mismatch', 'energy grid', 'effarea')

    if exposure <= 0.0:
        raise ArgumentErr('bad', 'exposure', 'value must be positive')

    return DataARF(name=name, energ_lo=elo, energ_hi=ehi,
                   specresp=specresp, exposure=exposure)


def make_ideal_rmf(e_min, e_max, offset=1, name='rmf'):
    """A simple in-memory representation of an ideal RMF.

    This RMF represents a 1-to-1 mapping from channel to energy
    bin (i.e. there's no blurring or secondary channels).

    Parameters
    ----------
    e_min, e_max : array
        The energy ranges corresponding to the channels. The units are
        in keV and each bin has energ_hi > energ_lo. The arrays are
        assumed to be ordered, but it is not clear yet whether they
        have to be in ascending order. The sizes must match each
        other. This corresponds to the E_MIN and E_MAX columns of
        the EBOUNDS extension of the RMF file format.
    offset : int, optional
        The value of the first channel (corresponding to the TLMIN
        value of the F_CHAN column of the 'MATRIX' or 'SPECRESP
        MATRIX' block. It is expected to be 0 or 1, but the only
        restriction is that it is 0 or greater.
    name : str, optional
        The name to give to the RMF instance.

    Returns
    -------
    rmf : sherpa.astro.data.DatAMRF
        The RMF.

    """

    elo = np.asarray(e_min)
    ehi = np.asarray(e_max)
    if elo.size != ehi.size:
        raise DataErr('mismatch', 'e_min', 'e_max')
    detchans = elo.size

    if offset < 0:
        raise ArgumentErr('bad', 'offset', 'value can not be negative')

    # The "ideal" matrix is the identity matrix, which, in compressed
    # form, is an array of 1.0's (matrix) and an array of locations
    # giving the column where the element is 1 (fchan). It appears
    # that this uses 1 indexing.
    #
    dummy = np.ones(detchans, dtype=np.int16)
    matrix = np.ones(detchans, dtype=np.float32)
    fchan = np.arange(1, detchans + 1, dtype=np.int16)

    return DataRMF(name=name, detchans=detchans,
                   energ_lo=elo, energ_hi=ehi,
                   n_grp=dummy, n_chan=dummy,
                   f_chan=fchan, matrix=matrix,
                   offset=offset)


def setup_likelihood(scale=False):
    """Create data and models for likelihood-based tests.

    Parameters
    ----------
    scale : bool, optional
        If True the expected statistic value includes the
        AREASCAL values, otherwise it doesn't.
    """

    dset = setup_basic_dataset()
    dset.ignore_bad()

    cpt1 = Const1D()
    cpt2 = StepHi1D()
    cpt1.c0 = 20
    cpt2.ampl = 20
    cpt2.xcut = 6.5
    mdl = cpt1 + cpt2

    # The asserts are used to ensure that any changes to parameter
    # values (or the models) don't lead to invalid data values
    # for the CSTAT calculation (i.e. this code does not have
    # to deal with the truncation value). As this routine is
    # called multiple times it's a bit wasteful.
    #
    counts = expected_basic_counts(scale=False)[1:]
    assert (counts > 0).all()

    # The model calculation uses the "integrated" version
    # when scale=True, since this compares the results to
    # the model evaluated through a response. For the
    # scale=False case, when the model is evaluated directly,
    # then the model is evaluated at the bin value, not
    # integrated over the bin (this is because the instrument
    # models actually ignore the channel grid sent to them
    # and do the calculation over the bins, whereas when there
    # is no instrument model the input arguments are used,
    # i.e. mdl([1,2,3]) rather than mdl([1,2,3], [2,3,4])
    #
    chanlo = dset.channel[1:]
    chanhi = 1 + chanlo
    mvals_noascal = mdl(chanlo)
    assert (mvals_noascal > 0).all()

    def calc_stat(m):
        statval = m - counts + counts * np.log(counts / m)
        return 2.0 * statval.sum()

    # Calculate the expected statistics; this should really be
    # a hard-coded value but is left as a calculation while the
    # tests are being worked on.
    #
    expected_noascal = calc_stat(mvals_noascal)

    ascal = dset.get_areascal(filter=True)
    mvals_ascal = mdl(chanlo, chanhi) * ascal
    assert (mvals_ascal > 0).all()

    # Calculate the expected statistics
    #
    expected_ascal = calc_stat(mvals_ascal)

    # TODO: remove once the values are hard-coded.
    #
    # Just to check that the area scaling is working in the
    # expected direction.
    assert (expected_noascal > expected_ascal)

    if scale:
        expected = expected_ascal
    else:
        expected = expected_noascal

    return dset, mdl, expected


def test_cstat_nophamodel():
    """What does CSTAT calculate when there is no PHA instrument model.

    The value here is technically wrong, in that the AREASCAL value
    is not being included in the calculation, but is included as a
    test to validate the current approach.

    See Also
    --------
    test_cstat_arfpha, test_cstat_rmfpha, test_cstat_rsppha
    """

    dset, mdl, expected = setup_likelihood(scale=False)

    stat = CStat()
    sval_noascal = stat.calc_stat(dset, mdl)

    assert_allclose(sval_noascal[0], expected)


def test_cstat_arfpha():
    """What does CSTAT calculate when there is an ARF+PHA instrument model.

    The value here is technically wrong, in that the AREASCAL value
    is not being included in the calculation, but is included as a
    test to validate the current approach.

    See Also
    --------
    test_cstat_nophamodel, test_cstat_rmfpha, test_cstat_rsppha
    """

    dset, mdl, expected = setup_likelihood(scale=True)

    # Use the channel grid as the "energy axis".
    #
    arf = make_arf(energ_lo=dset.channel,
                   energ_hi=dset.channel + 1)
    mdl_ascal = ARFModelPHA(arf, dset, mdl)

    stat = CStat()
    sval_ascal = stat.calc_stat(dset, mdl_ascal)

    assert_allclose(sval_ascal[0], expected)


def test_cstat_rmfpha():
    """What does CSTAT calculate when there is an RMF+PHA instrument model.

    This includes the AREASCAL when evaluating the model.

    See Also
    --------
    test_cstat_nophamodel, test_cstat_arfpha, test_cstat_rsppha
    """

    dset, mdl, expected = setup_likelihood(scale=True)

    # use the full channel grid; the energy grid has to be
    # "the same" as the channel values since the model
    # has a dependency on the independent axis
    #
    egrid = 1.0 * np.concatenate((dset.channel,
                                  [dset.channel.max() + 1]))
    rmf = make_ideal_rmf(egrid[:-1], egrid[1:])

    mdl_ascal = RMFModelPHA(rmf, dset, mdl)

    stat = CStat()
    sval_ascal = stat.calc_stat(dset, mdl_ascal)

    assert_allclose(sval_ascal[0], expected)


def test_cstat_rsppha():
    """What does CSTAT calculate when there is an RSP+PHA instrument model.

    This includes the AREASCAL when evaluating the model.

    See Also
    --------
    test_cstat_nophamodel, test_cstat_arfpha, test_cstat_rmfpha
    """

    dset, mdl, expected = setup_likelihood(scale=True)

    # use the full channel grid; the energy grid has to be
    # "the same" as the channel values since the model
    # has a dependency on the independent axis
    #
    egrid = 1.0 * np.concatenate((dset.channel,
                                  [dset.channel.max() + 1]))
    arf = make_arf(energ_lo=egrid[:-1],
                   energ_hi=egrid[1:])
    rmf = make_ideal_rmf(e_min=egrid[:-1], e_max=egrid[1:])

    mdl_ascal = RSPModelPHA(arf, rmf, dset, mdl)

    stat = CStat()
    sval_ascal = stat.calc_stat(dset, mdl_ascal)

    assert_allclose(sval_ascal[0], expected)


def test_get_dep_no_bgnd():
    """What does get_dep return: background but not subtracted"""

    # As no background subtracted, the same as test_get_dep,
    # except that the bad channel is now ignored.
    #
    dset = setup_basic_dataset_bgnd()
    dset.ignore_bad()

    expected = expected_basic_counts(scale=False)
    expected = expected.astype(np.float64)

    assert_allclose(dset.get_dep(), expected)


def test_get_y_no_bgnd():
    """What does get_y return: background but not subtracted"""

    # As no background subtracted, the same as test_get_y,
    # except that the bad channel is now ignored.
    #
    dset = setup_basic_dataset_bgnd()
    dset.ignore_bad()

    expected = expected_basic_counts(scale=True)
    expected = expected.astype(np.float64)

    assert_allclose(dset.get_y(), expected)


def test_get_dep_bgnd():
    """What does get_dep return: background subtracted"""

    dset = setup_basic_dataset_bgnd()
    dset.subtract()

    src = expected_basic_counts(scale=False)
    bg = expected_basic_counts_bgnd(scale=True)
    expected = src - bg

    assert_allclose(dset.get_dep(), expected)


def test_get_y_bgnd():
    """What does get_y return: background but not subtracted"""

    # As no background subtracted, the same as test_get_dep.
    dset = setup_basic_dataset_bgnd()
    dset.subtract()

    src = expected_basic_counts(scale=False)
    bg = expected_basic_counts_bgnd(scale=True)
    expected = src - bg
    ascal = expected_basic_areascal()
    ascal[ascal <= 0] = 1.0
    expected /= ascal

    assert_allclose(dset.get_y(), expected)


def test_get_staterror_no_bgnd():
    """What does get_staterror return when bgnd is not subtracted.

    This is the same as test_get_staterror, as the background is
    ignored in this case.
    """

    dset = setup_basic_dataset_bgnd()
    dset.ignore_bad()

    stat = Chi2DataVar()
    errors = dset.get_staterror(filter=True,
                                staterrfunc=stat.calc_staterror)

    expected = expected_basic_chisquare_errors(scale=False)
    assert_allclose(errors, expected)


def test_get_yerr_no_bgnd():
    """What does get_yerr return when bgnd is not subtracted.

    This is the same as test_get_yerr, as the background is
    ignored in this case.
    """

    dset = setup_basic_dataset_bgnd()
    dset.ignore_bad()

    stat = Chi2DataVar()
    errors = dset.get_yerr(filter=True,
                           staterrfunc=stat.calc_staterror)

    expected = expected_basic_chisquare_errors(scale=True)
    assert_allclose(errors, expected)


def test_get_staterror_bgnd():
    """What does get_staterror return when bgnd is subtracted.

    This is the easy case, since BACKSCAL and EXPOSURE are 1.
    """

    dset = setup_basic_dataset_bgnd()
    dset.ignore_bad()
    dset.subtract()

    stat = Chi2DataVar()
    errors = dset.get_staterror(filter=True,
                                staterrfunc=stat.calc_staterror)

    expected = expected_basic_chisquare_errors_bgnd(scale=False)
    assert_allclose(errors, expected)


def test_get_yerr_bgnd():
    """What does get_yerr return when bgnd is subtracted.

    This is the easy case, since BACKSCAL and EXPOSURE are 1.
    """

    dset = setup_basic_dataset_bgnd()
    dset.ignore_bad()
    dset.subtract()

    stat = Chi2DataVar()
    errors = dset.get_yerr(filter=True,
                           staterrfunc=stat.calc_staterror)

    expected = expected_basic_chisquare_errors_bgnd(scale=True)
    assert_allclose(errors, expected)


def test_get_staterror_scaling_bgnd():
    """What does get_staterror return when bgnd is subtracted.

    This expands on test_get_staterror_bgnd by having different
    BACKSCAL and EXPOSURE values in the source and background
    data sets.
    """

    dset = setup_basic_dataset_bgnd()
    dset.exposure = 20.0
    dset.backscal = 0.01
    dset.ignore_bad()
    dset.subtract()

    bset = dset.get_background()
    bset.exposure = 15.0
    bset.backscal = 0.1

    stat = Chi2DataVar()
    errors = dset.get_staterror(filter=True,
                                staterrfunc=stat.calc_staterror)

    expected = expected_basic_chisquare_errors_scaling_bgnd(scale=False)
    assert_allclose(errors, expected)


# Marked as xfail since I don't have the energy to resolve this just
# now, and I want to get a code commit in before something goes too
# wrong.
#
@pytest.mark.xfail
def test_get_yerr_scaling_bgnd():
    """What does get_yerr return when bgnd is subtracted.

    This expands on test_get_yerr_bgnd by having different
    BACKSCAL and EXPOSURE values in the source and background
    data sets.
    """

    dset = setup_basic_dataset_bgnd()
    dset.exposure = 20.0
    dset.backscal = 0.01
    dset.ignore_bad()
    dset.subtract()

    bset = dset.get_background()
    bset.exposure = 15.0
    bset.backscal = 0.1

    stat = Chi2DataVar()
    errors = dset.get_yerr(filter=True,
                           staterrfunc=stat.calc_staterror)

    expected = expected_basic_chisquare_errors_scaling_bgnd(scale=True)
    assert_allclose(errors, expected)


def create_xspec_comparison_dataset(make_data_path,
                                    keep_background=False):
    """Hack up the data set used in the XSPEC comparison tests."""

    infile = make_data_path('3c273.pi')
    dset = ui.unpack_pha(infile)

    # Remove background (if asked), ignore grouping, add bad channels and
    # an AREASCAL column. Sherpa ignores the error arrays in
    # the PHA file by default, so should have been read in,
    # hence the asserts.
    #
    assert dset.staterror is None
    assert dset.syserror is None

    dset.grouped = False

    chans = dset.channel
    quality = dset.quality

    ascal = np.ones(chans.size)
    ascal[(chans > 200) & (chans < 300)] = 0.1

    for cnum in [150, 151, 250, 251, 252]:
        quality[chans == cnum] = 1
        ascal[chans == cnum] = 0

    dset.areascal = ascal

    if not keep_background:
        dset.delete_background()
        return dset

    # Adjust the background
    #
    bkg = dset.get_background()

    ascal = np.ones(chans.size) * 1.2
    ascal[(chans > 150) & (chans < 250)] = 1.4

    bkg.areascal = ascal

    return dset


def validate_xspec_result(l, h, npts, ndof, statval):
    """Check that the stat results match those from XSPEC.

    This assumes the first data set returned by get_stat_info
    should be used.
    """

    ui.notice(None, None)
    ui.ignore(None, l)
    ui.ignore(h, None)
    ui.ignore_bad()
    sinfo = ui.get_stat_info()[0]
    assert sinfo.numpoints == npts
    assert sinfo.dof == ndof

    # XSPEC displays results to ~ 2dp, so the tolerance
    # is quite forgiving here. Or I could use pyxspec to
    # calculate and display this.
    #
    assert_allclose(sinfo.statval, statval, rtol=0, atol=0.005)


@requires_data
@requires_fits
def test_cstat_comparison_xspec(make_data_path):
    """Compare CSTAT values for a data set to XSPEC.

    This checks that the "UI layer" works, although ideally there
    should be a file that can be read in rather than having to
    manipulate it (the advantage here is that it means there is
    no messing around with adding a file to the test data set).

    The XSPEC version used was 12.9.0o.
    """

    dset = create_xspec_comparison_dataset(make_data_path,
                                           keep_background=False)

    ui.clean()
    ui.set_data(dset)
    # use powlaw1d rather than xspowerlaw so do not need XSPEC
    ui.set_source(ui.powlaw1d.pl)
    ui.set_par('pl.ampl', 1e-4)

    ui.set_stat('cstat')
    ui.set_analysis('channel')

    # Test three ranges. These cover
    #  a) areascal = 1, with bad channels
    #  b) areascal = 0.1, with bad channels
    #  c) combination of both
    #
    # Note that XSPEC and Sherpa deal differently with bins where the model
    # value is small. Sherpa replaces model <= 0 by a set value, whereas
    # XSPEC uses a different value, and a slightly-different replacement
    # scheme. That means that the ranges used here are chosen to
    # avoid areas where this replacement happens.
    #
    # Since XSPEC ignores bad channels, the equivalent ignore
    # lines for these ranges are
    #   ignore **-99, 151-**
    #   ignore **-209,291-**
    #   ignore **-7,768-**
    #
    for args in [(99, 153, 51, 49, 62.31),
                 (211, 296, 81, 79, 244.80),
                 (7, 773, 760, 758, 1063.09)]:
        validate_xspec_result(*args)

    ui.clean()


@requires_data
@requires_fits
def test_wstat_comparison_xspec(make_data_path):
    """Compare WSTAT values for a data set to XSPEC.

    See test_cstat_comparison_xspec.

    The XSPEC version used was 12.9.0o.
    """

    dset = create_xspec_comparison_dataset(make_data_path,
                                           keep_background=True)

    ui.clean()
    ui.set_data(dset)
    ui.set_source(ui.powlaw1d.pl)
    ui.set_par('pl.ampl', 1e-4)

    ui.set_stat('wstat')
    ui.set_analysis('channel')

    # Since XSPEC ignores bad channels, the equivalent ignore
    # lines for these ranges are
    #   ignore **-99, 151-**
    #   ignore **-209,291-**
    #   ignore **-7,768-**
    #
    for args in [(99, 153, 51, 49, 61.82),
                 (211, 296, 81, 79, 242.89),
                 (7, 773, 760, 758, 1043.85)]:
        validate_xspec_result(*args)

    ui.clean()
