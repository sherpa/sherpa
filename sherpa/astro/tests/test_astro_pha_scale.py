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

The area scaling is applied to the counts in that channel, rather than
being factored into the ARF. This complicates data access.

The tests here focus on the object API, rather than the UI layer
(which, at present, doesn't do much extra with regard to area
scaling, so this should catch most problems).

This does *not* test any plots (e.g. to ensure that the correctly-scaled
data and/or model values are shown).

"""

# import pytest
import six
import warnings

import numpy as np
from numpy.testing import assert_allclose

from sherpa.astro.data import DataPHA
from sherpa.models.basic import Const1D, StepHi1D
from sherpa.stats import Chi2DataVar


def expected_basic_areascal():
    """Return the expected areascal values."""

    areascal = np.ones(10)
    areascal[0] = 0.0
    areascal[6:] = 0.5
    return areascal.copy()


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


def expected_basic_chisquare_errors():
    """Return the expected error values (chi square).

    The calculation is sqrt(observed) / areascaling, so should
    be compared to the Chi2DataVar statistic.
    """

    # Calculate the errors based on the counts. There are no
    # counts less than 1, so this is just the square root of
    # the observed data value, which then has to be scaled by
    # the area scaling.
    #
    # Since ignore_bad has been called we ignore the first bin.
    counts = expected_basic_counts(scale=False)[1:]
    ascal = expected_basic_areascal()[1:]

    expected = np.sqrt(counts) / ascal
    return expected.copy()


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
    expected = expected_basic_counts(scale=True)
    expected = expected.astype(np.float64)
    expected[0] = np.nan

    assert_allclose(dset.get_dep(), expected,
                    equal_nan=True)


def test_get_y():
    """What does get_y return"""

    dset = setup_basic_dataset()
    expected = expected_basic_counts(scale=True)
    expected = expected.astype(np.float64)
    expected[0] = np.nan

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

    expected = expected_basic_chisquare_errors()
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

    counts = expected_basic_counts(scale=True)[1:]
    errors = expected_basic_chisquare_errors()
    mvals = mdl(dset.channel[1:])

    expected = (counts - mvals)**2 / (errors**2)
    expected = expected.sum()

    stat = Chi2DataVar()
    sval = stat.calc_stat(dset, mdl)

    assert_allclose(sval[0], expected)
