#
#  Copyright (C) 2017, 2018, 2020, 2021
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
data and/or model values are shown), but does test some of the accessor
methods the plotting code uses (e.g. get_y).

The unit-level tests here (mostly) focus on just the area scaling; the
more integration-level tests (i.e. those that read in data files)
are used to check that the scaling is being applied to cases
where the exposure time or background scaling may not be unity.

"""

import numpy as np
from numpy.testing import assert_allclose

import pytest

from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.astro.instrument import ARFModelPHA, RMFModelPHA, RSPModelPHA, \
    create_arf, create_delta_rmf
from sherpa.models.basic import Const1D, StepHi1D
from sherpa.stats import Chi2DataVar, CStat
from sherpa.utils.err import ArgumentErr, DataErr

from sherpa.utils.testing import requires_data, requires_fits, requires_group

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
        If True then the result is scaled by the ratio of the
        source and background AREASCAL values.

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
        # the first channel has areascal=0, so leave (it is 0).
        #
        src_ascal = expected_basic_areascal()[1:]
        bkg_ascal = expected_basic_areascal_bgnd()[1:]

        counts = counts.astype(np.float64)
        counts[1:] = counts[1:] * src_ascal / bkg_ascal

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

    # This could be simplified, but left as is to keep the
    # intent clear.
    #
    var = counts + bgcounts * ascal * ascal / (bgascal * bgascal)
    ans = np.sqrt(var)
    if scale:
        ans = ans / ascal

    return ans


def expected_basic_chisquare_errors_scaling_bgnd(scale=False):
    """Return the expected error values (chi square) after bg subtraction.

    Parameters
    ----------
    scale : bool, optional
        If True then the results are scaled by areascal.

    Notes
    -----
    The errors are calculated using the Chi2DataVar algorithm.

    It assumes:
       source exposure = 20
       background exposure = 15

       source backscal = 0.01
       background backscal = 0.1

    ignore_bad() is assumed to have been called.

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

    ratio = backscal_src * exposure_src * ascal / \
        (backscal_bg * exposure_bg * bgascal)
    var = counts + bgcounts * ratio * ratio
    ans = np.sqrt(var)
    if scale:
        ans = ans / ascal

    return ans


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
    test to validate this use case.

    See Also
    --------
    test_cstat_nophamodel, test_cstat_rmfpha, test_cstat_rsppha
    """

    dset, mdl, expected = setup_likelihood(scale=True)

    # Use the channel grid as the "energy axis".
    #
    arf = create_arf(dset.channel, dset.channel + 1)
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
    rmf = create_delta_rmf(egrid[:-1], egrid[1:])

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
    arf = create_arf(egrid[:-1], egrid[1:])
    rmf = create_delta_rmf(egrid[:-1], egrid[1:])

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
    expected = expected / dset.exposure
    assert_allclose(errors, expected)


def create_xspec_comparison_dataset(make_data_path,
                                    keep_background=False):
    """Hack up the data set used in the XSPEC comparison tests.

    This is to avoid adding a new data set to the sherpa-test-data
    repository.
    """

    infile = make_data_path('3c273.pi')
    dset = ui.unpack_pha(infile)

    # Remove background (if asked), ignore grouping, add bad channels and
    # an AREASCAL column. Sherpa ignores the error arrays in
    # the PHA file by default; the asserts are to check this.
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
@requires_data
@requires_fits
@pytest.mark.parametrize("l,h,ndp,ndof,statval",
                         [(99, 153, 51, 49, 62.31),
                          (211, 296, 81, 79, 244.80),
                          (7, 773, 760, 758, 1063.09)])
def test_cstat_comparison_xspec(make_data_path, l, h, ndp, ndof, statval):
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

    validate_xspec_result(l, h, ndp, ndof, statval)
    ui.clean()


@requires_data
@requires_fits
@pytest.mark.parametrize("l,h,ndp,ndof,statval",
                         [(99, 153, 51, 49, 61.82),
                          (211, 296, 81, 79, 242.89),
                          (7, 773, 760, 758, 1043.85)])
def test_wstat_comparison_xspec(make_data_path, l, h, ndp, ndof, statval):
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

    validate_xspec_result(l, h, ndp, ndof, statval)
    ui.clean()


# This uses energy ranges to filter the data, which makes XSPEC
# and Sherpa comparison easier (the ranges were chosen so that
# there is no difference in the start/end bin as it is not 100%
# clear what the respective filters in Sherpa and XSPEC do
# for range comparisons [inclusive or exclusive], or if there's
# some numerical differences in play at times).
#
@requires_data
@requires_fits
@pytest.mark.parametrize("l,h,ndp,ndof,statval",
                         [(0.5, 2.0, 101, 99, 201.21),
                          (3.2, 4.1, 57, 55, 255.23),
                          (0.5, 7.0, 439, 437, 904.91)])
def test_xspecvar_no_grouping_no_bg_comparison_xspec(make_data_path,
                                                     l, h, ndp, ndof, statval):
    """Compare chi2xspecvar values for a data set to XSPEC.

    The data set has no background.

    See test_cstat_comparison_xspec. Note that at present
    Sherpa and XSPEC treat bins with 0 values in them differently:
    see https://github.com/sherpa/sherpa/issues/356
    so for this test all bins are forced to have at least one
    count in them (source -> 5 is added per channel,background ->
    3 is added per channel).

    The XSPEC version used was 12.9.0o.
    """

    dset = create_xspec_comparison_dataset(make_data_path,
                                           keep_background=False)

    # Lazy, so add it to "bad" channels too
    dset.counts += 5

    ui.clean()
    ui.set_data(dset)

    ui.set_source(ui.powlaw1d.pl)
    ui.set_par('pl.ampl', 5e-4)

    ui.set_stat('chi2xspecvar')
    ui.set_analysis('energy')

    validate_xspec_result(l, h, ndp, ndof, statval)
    ui.clean()


@requires_data
@requires_fits
@pytest.mark.parametrize("l,h,ndp,ndof,statval",
                         [(0.5, 2.0, 101, 99, 228.21),
                          (3.2, 4.1, 57, 55, 251.95),
                          (0.5, 7.0, 439, 437, 960.56)])
def test_xspecvar_no_grouping_comparison_xspec(make_data_path,
                                               l, h, ndp, ndof, statval):
    """Compare chi2xspecvar values for a data set to XSPEC.

    The data set has a background. See
    test_xspecvar_no_grouping_no_bg_comparison_xspec

    The XSPEC version used was 12.9.0o.
    """

    dset = create_xspec_comparison_dataset(make_data_path,
                                           keep_background=True)

    # Lazy, so add it to "bad" channels too
    dset.counts += 5
    dset.get_background().counts += 3

    ui.clean()
    ui.set_data(dset)
    ui.subtract()

    ui.set_source(ui.powlaw1d.pl)
    ui.set_par('pl.ampl', 5e-4)

    ui.set_stat('chi2xspecvar')
    ui.set_analysis('energy')

    validate_xspec_result(l, h, ndp, ndof, statval)
    ui.clean()


# This is a regression test, not calculated from first principles
STATERROR_3C273_SRC = np.array([4.12310563, 3.87298335, 4.        , 3.87298335, 4.        ,
                                3.87298335, 4.24264069, 4.24264069, 3.87298335, 4.24264069,
                                3.87298335, 3.87298335, 4.35889894, 3.87298335, 3.87298335,
                                4.12310563, 4.        , 4.        , 4.12310563, 3.87298335,
                                4.35889894, 3.87298335, 4.        , 3.87298335, 4.        ,
                                4.12310563, 3.87298335, 4.24264069, 4.        , 3.87298335,
                                3.87298335, 4.        , 3.87298335, 3.87298335, 3.87298335,
                                4.        , 4.        , 3.87298335, 3.87298335, 4.        ,
                                4.        , 3.87298335, 4.        , 3.87298335, 3.87298335,
                                4.47213595])

STATERROR_3C273_BKG = np.array([ 1.41421356,  1.        ,  0.        ,  1.        ,  1.        ,
                                 0.        ,  0.        ,  0.        ,  1.        ,  1.        ,
                                 0.        ,  1.73205081,  1.        ,  1.        ,  1.        ,
                                 0.        ,  0.        ,  1.41421356,  1.41421356,  0.        ,
                                 1.73205081,  1.        ,  1.73205081,  1.        ,  0.        ,
                                 1.73205081,  0.        ,  1.73205081,  0.        ,  0.        ,
                                 2.        ,  0.        ,  1.73205081,  1.73205081,  2.        ,
                                 1.        ,  1.73205081,  2.23606798,  2.23606798,  2.23606798,
                                 2.44948974,  2.23606798,  2.44948974,  2.64575131,  4.35889894,
                                 10.44030651])

STATERROR_3C273_BSCALE = 0.134920643888096

STATERROR_3C273_BKG_REGROUP = np.array([0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        1.        , 1.        , 0.        , 0.        , 0.        ,
                                        2.23606798, 3.16227766, 3.31662479, 2.64575131, 2.44948974,
                                        2.82842712, 3.16227766, 3.31662479, 2.82842712, 2.44948974,
                                        2.64575131, 2.        , 0.        , 1.73205081, 1.73205081,
                                        2.        , 3.16227766, 3.16227766, 2.44948974, 1.        ,
                                        1.        , 0.        , 0.        , 1.        , 0.        ,
                                        0.        , 0.        , 0.        , 1.        , 0.        ,
                                        0.        , 1.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        1.        , 1.        , 0.        , 0.        , 0.        ,
                                        1.        , 0.        , 0.        , 1.        , 0.        ,
                                        1.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 1.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 1.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        1.        , 0.        , 0.        , 0.        , 1.41421356,
                                        0.        , 0.        , 1.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 1.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        1.41421356, 1.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 1.        , 0.        ,
                                        0.        , 1.        , 0.        , 0.        , 0.        ,
                                        1.        , 0.        , 0.        , 0.        , 0.        ,
                                        1.41421356, 1.        , 0.        , 1.        , 0.        ,
                                        0.        , 1.        , 0.        , 0.        , 1.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 1.        , 0.        ,
                                        0.        , 1.        , 1.        , 0.        , 0.        ,
                                        0.        , 1.        , 1.        , 1.        , 1.        ,
                                        0.        , 0.        , 0.        , 1.        , 0.        ,
                                        1.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 1.        , 0.        , 1.        ,
                                        1.        , 0.        , 1.        , 0.        , 0.        ,
                                        0.        , 0.        , 1.        , 0.        , 0.        ,
                                        0.        , 1.        , 0.        , 1.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        1.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 0.        , 0.        , 0.        , 0.        ,
                                        0.        , 1.        , 0.        , 6.164414  ])


# This is the mask to match
#    notice()
#    notice(0.5, 2)
#    notice(5, 7)
#
STATERROR_3C273_MASK = np.array([False, False, False,  True,  True,  True,  True,  True,  True,
                                 True,  True,  True,  True,  True,  True,  True,  True,  True,
                                 True,  True,  True,  True,  True,  True,  True,  True, False,
                                 False, False, False, False, False, False, False, False, False,
                                 False, False, False, False,  True,  True,  True,  True,  True,
                                 False])


# The background dataset now has a different mapping to groups
# than the source region but still using the "same" energy
# filter (modulo different boundaries). The tabStops means
# that the first 20 and last ~200 channels map to individual
# groups, hence the large number of elements.
#
STATERROR_3C273_MASK_REGROUP = np.zeros(264, bool)
STATERROR_3C273_MASK_REGROUP[20:23] = True
STATERROR_3C273_MASK_REGROUP[27:37] = True


def loadup_3c273(use_errors, subt, noticed, make_data_path):
    """Load up the 3C273 data.

    Requires @requires_fits and @requires_data on the test.
    """

    # We could create a PHA object but it's easiest to use
    # an on-disk version to setup everything.
    #
    import sherpa.astro.io

    path = make_data_path('3c273.pi')
    pha = sherpa.astro.io.read_pha(path, use_errors=use_errors)

    if subt:
        pha.subtract()
    else:
        pha.unsubtract()

    if noticed:
        # After this, we have both
        #   pha.get_filter(format='%.3f')
        #   pha.get_background().get_filter(format='%.3f')
        # returning
        #   '0.518:1.986,4.869:8.220'
        #
        # In channels this is
        #   '36:136,334:563'
        #
        pha.notice(0.5, 2)
        pha.notice(5, 7)

    return pha


@requires_data
@requires_fits
@pytest.mark.parametrize("filt", [False, True])
@pytest.mark.parametrize("subt", [False, True])
@pytest.mark.parametrize("noticed", [False, True])
def test_get_staterror_file_no_errors(filt, subt, noticed, make_data_path):
    """Simple test of staterror to improve test coverage.

    Checks of background handling.

    Done during https://github.com/sherpa/sherpa/pull/1151
    """

    # With use_errors=False there's no statistical error
    pha = loadup_3c273(False, subt, noticed, make_data_path)
    ans = pha.get_staterror(filter=filt)
    assert ans is None


@requires_data
@requires_fits
@pytest.mark.parametrize("filt", [False, True])
@pytest.mark.parametrize("noticed", [False, True])
def test_get_staterror_file_errors_unsubtracted(filt, noticed, make_data_path):
    """Simple test of staterror to improve test coverage.

    Checks of background handling.

    Done during https://github.com/sherpa/sherpa/pull/1151
    """

    pha = loadup_3c273(True, False, noticed, make_data_path)

    ans = pha.get_staterror(filter=filt)
    expected = STATERROR_3C273_SRC
    if filt and noticed:
        expected = expected[STATERROR_3C273_MASK]

    assert ans == pytest.approx(expected)


@requires_data
@requires_fits
@pytest.mark.parametrize("filt", [False, True])
@pytest.mark.parametrize("noticed", [False, True])
def test_get_staterror_file_errors_subtracted(filt, noticed, make_data_path):
    """Simple test of staterror to improve test coverage.

    Checks of background handling.

    Done during https://github.com/sherpa/sherpa/pull/1151
    """

    pha = loadup_3c273(True, True, noticed, make_data_path)

    ans = pha.get_staterror(filter=filt)
    expected = STATERROR_3C273_SRC**2 + (STATERROR_3C273_BSCALE * STATERROR_3C273_BKG)**2
    expected = np.sqrt(expected)
    if filt and noticed:
        expected = expected[STATERROR_3C273_MASK]

    assert ans == pytest.approx(expected)


@requires_data
@requires_fits
@pytest.mark.parametrize("filt", [False, True])
@pytest.mark.parametrize("noticed", [False, True])
def test_get_staterror_file_errors_bg(filt, noticed, make_data_path):
    """Simple test of staterror to improve test coverage.

    Checks of background handling.

    Done during https://github.com/sherpa/sherpa/pull/1151
    """

    pha = loadup_3c273(True, False, noticed, make_data_path)

    # We use the "source" grouping by default
    ans = pha.get_background().get_staterror(filter=filt)
    expected = STATERROR_3C273_BKG
    if filt and noticed:
        expected = expected[STATERROR_3C273_MASK]

    assert ans == pytest.approx(expected)


@requires_data
@requires_fits
@requires_group
@pytest.mark.parametrize("filt", [False, True])
@pytest.mark.parametrize("noticed", [False, True])
def test_get_staterror_file_errors_bg_regrouped(filt, noticed, make_data_path):
    """Simple test of staterror to improve test coverage.

    Checks of background handling. This time with a separate
    grouping to the source.

    Done during https://github.com/sherpa/sherpa/pull/1151
    """

    pha = loadup_3c273(True, False, noticed, make_data_path)

    bpha = pha.get_background()

    # Manually create some tab stops
    tabStops = np.zeros(1024, dtype=bool)
    tabStops[0:20] = True
    tabStops[800:] = True
    bpha.group_width(40, tabStops=tabStops)

    # After this, we have
    #  pha.get_filter(format='%.3f')
    #    '0.518:1.986,4.869:8.220'
    # and
    #  bpha.get_filter(format='%.3f')
    #    '0.584:1.752,4.672:9.928'
    #  which in channels os
    #    '40:120,320:680'
    #
    ans = pha.get_background().get_staterror(filter=filt)
    expected = STATERROR_3C273_BKG_REGROUP
    if filt and noticed:
        expected = expected[STATERROR_3C273_MASK_REGROUP]

    assert ans == pytest.approx(expected)
