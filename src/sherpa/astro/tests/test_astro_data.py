#
#  Copyright (C) 2007, 2015, 2017, 2018, 2020 - 2024
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

import logging
import warnings

import numpy as np

import pytest

from sherpa.astro.ui.utils import Session
from sherpa.astro.data import Data1D, DataARF, DataPHA, DataRMF, DataIMG, DataIMGInt
from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa.models.basic import Gauss1D, Gauss2D
from sherpa.utils import parse_expr
from sherpa.utils.err import ArgumentTypeErr, DataErr
from sherpa.utils.numeric_types import SherpaFloat
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_group, requires_region


EMPTY_DATA_OBJECTS_1D = [(DataPHA, [None] * 2)]
EMPTY_DATA_OBJECTS_2D = [(DataIMG, [None] * 3),
                         (DataIMGInt, [None] * 5)]
EMPTY_DATA_OBJECTS = EMPTY_DATA_OBJECTS_1D + EMPTY_DATA_OBJECTS_2D

def _monotonic_warning(response_type, filename):
    return UserWarning("The {} '{}' has a non-monotonic ENERG_LO array".format(response_type, filename))


def _bin_warning(response_type, filename):
    return UserWarning("The {} '{}' has at least one bin with ENERG_HI < ENERG_LO".format(response_type, filename))


def _assert_userwarning(expected_warnings, observed_warnings):

    expected_warnings_set = {w.args for w in expected_warnings}
    observed_warnings_set = {w.message.args for w in observed_warnings}

    assert observed_warnings_set == expected_warnings_set


@pytest.fixture
def setUp1():

    emin = np.array([
        1.46000006e-03, 2.48199999e-01, 3.06600004e-01, 4.67200011e-01,
        5.69400012e-01, 6.42400026e-01, 7.00800002e-01, 7.44599998e-01,
        7.88399994e-01, 8.17600012e-01, 8.61400008e-01, 8.90600026e-01,
        9.49000001e-01, 9.92799997e-01, 1.03659999e+00, 1.09500003e+00,
        1.13880002e+00, 1.19719994e+00, 1.28480005e+00, 1.40160000e+00,
        1.47459996e+00, 1.60599995e+00, 1.69360006e+00, 1.81040001e+00,
        1.89800000e+00, 1.94180000e+00, 2.02940011e+00, 2.08780003e+00,
        2.19000006e+00, 2.27760005e+00, 2.39439988e+00, 2.58419991e+00,
        2.71560001e+00, 2.86159992e+00, 3.08060002e+00, 3.38720012e+00,
        3.56240010e+00, 3.79600000e+00, 4.02960014e+00, 4.24860001e+00,
        4.71579981e+00, 5.02239990e+00, 5.37279987e+00, 5.89839983e+00,
        6.57000017e+00, 9.86960030e+00])

    emax = np.array([
        0.2482, 0.3066, 0.46720001, 0.56940001, 0.64240003,
        0.7008, 0.7446, 0.78839999, 0.81760001, 0.86140001,
        0.89060003, 0.949, 0.9928, 1.03659999, 1.09500003,
        1.13880002, 1.19719994, 1.28480005, 1.4016, 1.47459996,
        1.60599995, 1.69360006, 1.81040001, 1.898, 1.9418,
        2.02940011, 2.08780003, 2.19000006, 2.27760005, 2.39439988,
        2.58419991, 2.71560001, 2.86159992, 3.08060002, 3.38720012,
        3.5624001, 3.796, 4.02960014, 4.24860001, 4.71579981,
        5.0223999, 5.37279987, 5.89839983, 6.57000017, 9.8696003,
        14.95040035])

    pha = DataPHA('', np.arange(46, dtype=float) + 1.,
                  np.zeros(46),
                  bin_lo=emin,
                  bin_hi=emax)
    pha.units = "energy"
    return pha


def test_filter_energy_grid_notice(setUp1):
    pha = setUp1
    pha.notice()
    pha.notice(0.0, 6.0)

    # Use approx to make it easy to check an array
    expected = np.ones(46, dtype=bool)
    expected[44:46] = False
    assert pha.mask == pytest.approx(expected)


def test_filter_energy_grid_ignore(setUp1):
    pha = setUp1
    pha.notice()
    pha.ignore(0.0, 1.0)
    pha.ignore(3.0, 15.0)

    expected = np.zeros(46, dtype=bool)
    expected[14:33] = True
    assert pha.mask == pytest.approx(expected)


@pytest.fixture
def setUp2():

    emin = np.array([
        2.39196181, 2.35973215, 2.34076023, 2.30973101, 2.2884388,
        2.25861454, 2.22371697, 2.20662117, 2.18140674, 2.14317489,
        2.12185216, 2.09055495, 2.06256914, 2.04509854, 2.02788448,
        2.00133967, 1.97772908, 1.96379483, 1.93868744, 1.91855776,
        1.89444292, 1.87936974, 1.85819471, 1.84568763, 1.82923627,
        1.78920078, 1.77360916, 1.76206875, 1.74499893, 1.73006463,
        1.70084822, 1.6883322, 1.67772949, 1.65171933, 1.63476169,
        1.59687376, 1.5745424, 1.55736887, 1.54051399, 1.52546024,
        1.50043869, 1.48890531, 1.47329199, 1.46072423, 1.44289041,
        1.43344045, 1.41616774, 1.40441585, 1.3979584, 1.38773119,
        1.37138033, 1.35170007, 1.33725214, 1.33249414, 1.31839108,
        1.30797839, 1.29657102, 1.28310275, 1.26550889, 1.25471842,
        1.24513853, 1.23672664, 1.22944438, 1.21509433, 1.21003771,
        1.20401597, 1.19705439, 1.18722582, 0.90194935, 0.89519638,
        0.88912934, 0.88492262, 0.87837797, 0.87366825, 0.8689999,
        0.86437255, 0.85693878, 0.84793305, 0.84404182, 0.83580172,
        0.82876647, 0.82395256, 0.81865752, 0.81185687, 0.80004948,
        0.79450154, 0.78852075, 0.77920061, 0.77340651, 0.76626247,
        0.76202762, 0.75783074, 0.75413191, 0.74727529, 0.74321008,
        0.73474538, 0.73166627, 0.72687, 0.71785438, 0.71488959,
        0.71068853, 0.70199603, 0.69832331, 0.69387686, 0.68788701,
        0.68354762, 0.67847627, 0.67117327, 0.66512167, 0.66175646,
        0.65620857, 0.6518243, 0.64605182, 0.64142239, 0.63754696,
        0.63128632, 0.62478495, 0.62006336, 0.61440694, 0.60915887,
        0.60591549, 0.60078359, 0.5938406, 0.59103745, 0.58488411,
        0.58124125, 0.57883304, 0.57406437, 0.57023615, 0.56442606,
        0.56041539, 0.55701393, 0.55392498, 0.55030966, 0.54346251,
        0.53728294, 0.53515989, 0.5291304, 0.52448714, 0.51990861,
        0.51589233, 0.50996011, 0.50509953, 0.49889025, 0.49512967,
        0.49003205, 0.48888513, 0.48524383, 0.48164544, 0.47720695,
        0.47283325, 0.46916556, 0.46660379, 0.46280268, 0.45925769,
        0.45514211, 0.45290345, 0.44987884, 0.44589564, 0.44333643,
        0.44099477, 0.43790293, 0.43446559, 0.43088335, 0.42605683,
        0.42131537, 0.41826019, 0.41506338, 0.41155648, 0.40895697,
        0.40502119, 0.40400422, 0.40164718, 0.39864835, 0.39584854,
        0.39389083, 0.39130434, 0.38890362, 0.38526753, 0.38292497,
        0.38075879, 0.37891743, 0.37648395, 0.37557775, 0.37347662,
        0.37154216, 0.36742872, 0.3641032, 0.36167556, 0.35983625,
        0.35634032, 0.35248783, 0.35085678, 0.34843227, 0.34669766,
        0.34418666, 0.33912122, 0.33720407, 0.33505177, 0.33279634,
        0.33081138, 0.32847831, 0.32592943, 0.3111549], float)

    emax = np.array([
        3.06803656, 2.39196181, 2.35973215, 2.34076023, 2.30973101,
        2.2884388, 2.25861454, 2.22371697, 2.20662117, 2.18140674,
        2.14317489, 2.12185216, 2.09055495, 2.06256914, 2.04509854,
        2.02788448, 2.00133967, 1.97772908, 1.96379483, 1.93868744,
        1.91855776, 1.89444292, 1.87936974, 1.85819471, 1.84568763,
        1.82923627, 1.78920078, 1.77360916, 1.76206875, 1.74499893,
        1.73006463, 1.70084822, 1.6883322, 1.67772949, 1.65171933,
        1.63476169, 1.59687376, 1.5745424, 1.55736887, 1.54051399,
        1.52546024, 1.50043869, 1.48890531, 1.47329199, 1.46072423,
        1.44289041, 1.43344045, 1.41616774, 1.40441585, 1.3979584,
        1.38773119, 1.37138033, 1.35170007, 1.33725214, 1.33249414,
        1.31839108, 1.30797839, 1.29657102, 1.28310275, 1.26550889,
        1.25471842, 1.24513853, 1.23672664, 1.22944438, 1.21509433,
        1.21003771, 1.20401597, 1.19705439, 1.18722582, 0.90194935,
        0.89519638, 0.88912934, 0.88492262, 0.87837797, 0.87366825,
        0.8689999, 0.86437255, 0.85693878, 0.84793305, 0.84404182,
        0.83580172, 0.82876647, 0.82395256, 0.81865752, 0.81185687,
        0.80004948, 0.79450154, 0.78852075, 0.77920061, 0.77340651,
        0.76626247, 0.76202762, 0.75783074, 0.75413191, 0.74727529,
        0.74321008, 0.73474538, 0.73166627, 0.72687, 0.71785438,
        0.71488959, 0.71068853, 0.70199603, 0.69832331, 0.69387686,
        0.68788701, 0.68354762, 0.67847627, 0.67117327, 0.66512167,
        0.66175646, 0.65620857, 0.6518243, 0.64605182, 0.64142239,
        0.63754696, 0.63128632, 0.62478495, 0.62006336, 0.61440694,
        0.60915887, 0.60591549, 0.60078359, 0.5938406, 0.59103745,
        0.58488411, 0.58124125, 0.57883304, 0.57406437, 0.57023615,
        0.56442606, 0.56041539, 0.55701393, 0.55392498, 0.55030966,
        0.54346251, 0.53728294, 0.53515989, 0.5291304, 0.52448714,
        0.51990861, 0.51589233, 0.50996011, 0.50509953, 0.49889025,
        0.49512967, 0.49003205, 0.48888513, 0.48524383, 0.48164544,
        0.47720695, 0.47283325, 0.46916556, 0.46660379, 0.46280268,
        0.45925769, 0.45514211, 0.45290345, 0.44987884, 0.44589564,
        0.44333643, 0.44099477, 0.43790293, 0.43446559, 0.43088335,
        0.42605683, 0.42131537, 0.41826019, 0.41506338, 0.41155648,
        0.40895697, 0.40502119, 0.40400422, 0.40164718, 0.39864835,
        0.39584854, 0.39389083, 0.39130434, 0.38890362, 0.38526753,
        0.38292497, 0.38075879, 0.37891743, 0.37648395, 0.37557775,
        0.37347662, 0.37154216, 0.36742872, 0.3641032, 0.36167556,
        0.35983625, 0.35634032, 0.35248783, 0.35085678, 0.34843227,
        0.34669766, 0.34418666, 0.33912122, 0.33720407, 0.33505177,
        0.33279634, 0.33081138, 0.32847831, 0.32592943], float)

    pha = DataPHA('', np.arange(204, dtype=float) + 1.,
                  np.zeros(204),
                  bin_lo=emin,
                  bin_hi=emax)
    pha.units = "energy"
    return pha


def test_test_energy_grid_reversed_notice(setUp2):
    pha = setUp2

    pha.notice()
    pha.notice(4., 8.3)

    expected = np.zeros(204, dtype=bool)
    expected[0:42] = True
    assert pha.mask == pytest.approx(expected)

def test_test_energy_grid_reversed_ignore(setUp2):
    pha = setUp2

    pha.notice()
    pha.ignore(10.3, 13.8)
    pha.ignore(4.6, 6.2)

    expected = np.ones(204, dtype=bool)
    expected[66:70] = False
    expected[0:17] = False
    assert pha.mask == pytest.approx(expected)


@pytest.fixture
def setUp3():

    emin = np.arange(205.7875, 0.9875, -0.0125)
    emax = emin + 0.0125

    pha = DataPHA('', np.arange(16384, dtype=float) + 1,
                  np.zeros(16384),
                  bin_lo=emin,
                  bin_hi=emax)
    pha.units = 'wavelength'
    return pha


def test_filter_wave_grid_notice(setUp3):
    pha = setUp3

    pha.notice()
    pha.notice(100.0, 225.0)

    expected = np.ones(16384, dtype=bool)
    expected[8464:16384] = False
    assert pha.mask == pytest.approx(expected)


def test_filter_wave_grid_ignore(setUp3):
    pha = setUp3

    pha.notice()
    pha.ignore(30.01, 225.0)
    pha.ignore(0.1, 6.0)

    expected = np.zeros(16384, dtype=bool)
    expected[14064:15984] = True
    assert pha.mask == pytest.approx(expected)


# It would be nice to add some unit testing here, but it's not trivial
# and time doesn't allow.
@requires_data
@requires_fits
def test_bug_275(make_data_path):
    session = Session()
    session.load_data(make_data_path('3c273.pi'))
    str(session.get_data())
    str(session.get_rmf())
    str(session.get_arf())

    session.load_data(make_data_path('img.fits'))
    str(session.get_data())


# Test some simple "invalid input" cases. Unfortunately some of them
# are seen with released data products, so it is not sensible to
# error out for all errors.
#
@pytest.mark.parametrize("ethresh", [0.0, -1e-10, -100])
def test_arf_with_non_positive_thresh(ethresh):
    """Check the error-handling works when ethresh <= 0"""

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]
    specresp = energ_lo * 0 + 1.0

    with pytest.raises(ValueError,
                       match="ethresh is None or > 0"):
        create_arf(energ_lo, energ_hi, specresp, ethresh=ethresh)


@pytest.mark.parametrize("idx", [0, 1, 5, -2, -1])
def test_arf_with_swapped_energy_bounds(idx):
    """What happens if elo >= ehi?

    The bin edges are swapped at position idx.
    """

    # Ensure energy grid starts > 0
    energy = 0.001 + np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]
    specresp = energ_lo * 0 + 1.0

    # test energ_hi < energ_lo
    energ_lo[idx], energ_hi[idx] = energ_hi[idx], energ_lo[idx]

    if idx != -1:
        expected_warnings = [_bin_warning('ARF', 'user-arf'), _monotonic_warning('ARF', 'user-arf')]
    else:
        expected_warnings = [_bin_warning('ARF', 'user-arf')]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_arf(energ_lo, energ_hi, specresp)

    _assert_userwarning(expected_warnings, ws)

    # test energ_hi == energ_lo
    energ_lo[idx] = energ_hi[idx]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_arf(energ_lo, energ_hi, specresp)

    _assert_userwarning(expected_warnings, ws)


@pytest.mark.parametrize("idx", [0, 1, 5, -3, -2])
def test_arf_with_non_monotonic_grid(idx):
    """What happens if the grid is not monotonic?"""

    # For this test we want the ehi values to be larger than the
    # corresponding elo values (otherwise a different condition)
    # is triggered, but for energ_lo or energ_hi itself to
    # be non-monotonic. A non-consecutive array is picked as
    # this is a form not used much in these tests.
    #
    energ_lo = np.asarray([0.1, 0.4, 0.8, 1.0, 1.2, 2.0, 3.0, 4.1, 4.8, 5.2])
    energ_hi = np.asarray([0.3, 0.7, 0.9, 1.1, 1.9, 2.1, 3.5, 4.6, 5.1, 5.4])
    specresp = energ_lo * 0 + 1.0

    idx1 = idx + 1
    energ_lo[idx], energ_lo[idx1] = energ_lo[idx1], energ_lo[idx]
    energ_hi[idx], energ_hi[idx1] = energ_hi[idx1], energ_hi[idx]

    expected_warnings = [_monotonic_warning('ARF', 'user-arf')]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_arf(energ_lo, energ_hi, specresp)

    _assert_userwarning(expected_warnings, ws)

    # now make the two consecutive bin edges be the same
    #
    energ_lo[idx] = energ_lo[idx1]
    energ_hi[idx] = energ_hi[idx1]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_arf(energ_lo, energ_hi, specresp)

    _assert_userwarning(expected_warnings, ws)


def test_arf_with_zero_energy_elem():
    """What happens creating an ARf with a zero-energy element.

    This is for the case where the first bin starts at E=0 keV.
    """

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]
    specresp = energ_lo * 0 + 1.0

    with pytest.raises(DataErr,
                       match="The ARF 'user-arf' has an ENERG_LO value <= 0"):
        create_arf(energ_lo, energ_hi, specresp)


def test_arf_with_zero_energy_elem_replace():
    """What happens creating an ARf with a zero-energy element?

    This is for the case where the first bin starts at E=0 keV.
    In this case the ARF is allowed to replace the 0 value.
    """

    ethresh = 1.0e-5

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]
    specresp = energ_lo * 0 + 1.0

    expected_warnings = [UserWarning("The minimum ENERG_LO in the ARF 'user-arf' was 0 " +
                                     "and has been replaced by {}".format(ethresh))]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        adata = create_arf(energ_lo, energ_hi, specresp, ethresh=ethresh)

    _assert_userwarning(expected_warnings, ws)

    assert isinstance(adata, DataARF)
    assert adata.energ_lo[0] == pytest.approx(ethresh)


def test_arf_with_grid_below_thresh():
    """The first bin starts above 0 but ends below the threshold.

    This is a valid grid (other than the fact it is not
    consecutive), so the ARF can be created. See
    test_arf_with_grid_below_thresh_zero() for the
    case when the bin starts at 0.
    """

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    energ_lo[0] = 1e-7
    energ_hi[0] = 2e-7

    # The test is to make sure that the call does not
    # error out
    adata = create_arf(energ_lo, energ_hi, ethresh=1e-5)
    assert isinstance(adata, DataARF)


def test_arf_with_grid_below_thresh_zero():
    """The first bin starts at 0 but ends below the threshold."""

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    energ_lo[0] = 0.0
    energ_hi[0] = 1e-7

    with pytest.raises(DataErr,
                       match="The ARF 'user-arf' has an ENERG_HI value <= " +
                       "the replacement value of 1e-05"):
        create_arf(energ_lo, energ_hi, ethresh=1e-5)


def test_arf_with_decreasing_energies():
    """ENERG_LO < ENERG_HI for each row, but in decreasing order.

    This does not appear to be a valid OGIP file: are there
    examples in the real world that do this?
    """

    energy = np.arange(0.1, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    # Reverse the arrays
    energ_lo = energ_lo[::-1]
    energ_hi = energ_hi[::-1]

    # Programmer sanity check...
    assert energ_lo[1] < energ_lo[0]
    assert energ_hi[1] < energ_hi[0]

    adata = create_arf(energ_lo, energ_hi)
    assert isinstance(adata, DataARF)
    assert np.all(adata.energ_lo == energ_lo)
    assert np.all(adata.energ_hi == energ_hi)


@pytest.mark.parametrize("ethresh", [0.0, -1e-10, -100])
def test_rmf_with_non_positive_thresh(ethresh):
    """Check the error-handling works when ethresh <= 0"""

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    with pytest.raises(ValueError,
                       match="ethresh is None or > 0"):
        create_delta_rmf(energ_lo, energ_hi, ethresh=ethresh)


@pytest.mark.parametrize("idx", [0, 1, 5, -2, -1])
def test_rmf_with_swapped_energy_bounds(idx):
    """What happens if elo >= ehi?

    The bin edges are swapped at position idx.
    """

    # Ensure energy grid starts > 0
    energy = 0.001 + np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    # test energ_hi < energ_lo
    if idx != -1:
        expected_warnings = [_monotonic_warning('RMF', 'delta-rmf'), _bin_warning('RMF', 'delta-rmf')]
    else:
        expected_warnings = [_bin_warning('RMF', 'delta-rmf')]

    energ_lo[idx], energ_hi[idx] = energ_hi[idx], energ_lo[idx]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_delta_rmf(energ_lo, energ_hi)

    _assert_userwarning(expected_warnings, ws)

    # test energ_hi == energ_lo
    energ_lo[idx] = energ_hi[idx]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_delta_rmf(energ_lo, energ_hi)

    _assert_userwarning(expected_warnings, ws)


@pytest.mark.parametrize("idx", [0, 1, 5, -3, -2])
def test_rmf_with_non_monotonic_grid(idx):
    """What happens if the grid is not monotonic?"""

    # For this test we want the ehi values to be larger than the
    # corresponding elo values (otherwise a different condition)
    # is triggered, but for energ_lo or energ_hi itself to
    # be non-monotonic. A non-consecutive array is picked as
    # this is a form not used much in these tests.
    #
    energ_lo = np.asarray([0.1, 0.4, 0.8, 1.0, 1.2, 2.0, 3.0, 4.1, 4.8, 5.2])
    energ_hi = np.asarray([0.3, 0.7, 0.9, 1.1, 1.9, 2.1, 3.5, 4.6, 5.1, 5.4])

    idx1 = idx + 1
    energ_lo[idx], energ_lo[idx1] = energ_lo[idx1], energ_lo[idx]
    energ_hi[idx], energ_hi[idx1] = energ_hi[idx1], energ_hi[idx]

    expected_warnings = [_monotonic_warning('RMF', 'delta-rmf')]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_delta_rmf(energ_lo, energ_hi)

    _assert_userwarning(expected_warnings, ws)

    # now make the two consecutive bin edges be the same
    #
    energ_lo[idx] = energ_lo[idx1]
    energ_hi[idx] = energ_hi[idx1]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_delta_rmf(energ_lo, energ_hi)

    _assert_userwarning(expected_warnings, ws)


def test_rmf_with_zero_energy_elem():
    """What happens creating a RMf with a zero-energy element.

    This is for the case where the first bin starts at E=0 keV.
    """

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    with pytest.raises(DataErr,
                       match="The RMF 'delta-rmf' has an ENERG_LO value <= 0"):
        create_delta_rmf(energ_lo, energ_hi)


def test_rmf_with_zero_energy_elem_replace():
    """What happens creating a RMf with a zero-energy element.

    This is for the case where the first bin starts at E=0 keV.
    In this case the RMF is allowed to replace the 0 value.
    """

    ethresh = 1.0e-4

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    expected_warnings = [UserWarning("The minimum ENERG_LO in the RMF 'delta-rmf' was 0 " +
                                     "and has been replaced by {}".format(ethresh))]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        rdata = create_delta_rmf(energ_lo, energ_hi, ethresh=ethresh)

    _assert_userwarning(expected_warnings, ws)

    assert isinstance(rdata, DataRMF)
    assert rdata.energ_lo[0] == pytest.approx(ethresh)


def test_arf_with_negative_energy_elem():
    """What happens creating an ARf with negative energies.

    Hopefully we do not have files like this in use.
    """

    # Special case it so that the first in ends at 0 keV
    energy = np.arange(0, 1.0, 0.1, dtype=np.float32)
    energy = energy - energy[1]

    energ_lo = energy[:-1]
    energ_hi = energy[1:]
    specresp = energ_lo * 0 + 1.0

    with pytest.raises(DataErr,
                       match="The ARF 'user-arf' has an ENERG_LO value <= 0"):
        create_arf(energ_lo, energ_hi, specresp)


def test_arf_with_negative_energy_elem_replace():
    """What happens creating an ARf with negative energies.

    Hopefully we do not have files like this in use. Note that
    this errors out even with the replacement value set.
    """

    ethresh = 1.0e-5

    # Special case it so that the first in ends at 0 keV
    energy = np.arange(0, 1.0, 0.1, dtype=np.float32)
    energy = energy - energy[1]

    energ_lo = energy[:-1]
    energ_hi = energy[1:]
    specresp = energ_lo * 0 + 1.0

    expected_warnings = [UserWarning("The ARF 'user-arf' has an ENERG_LO value < 0")]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_arf(energ_lo, energ_hi, specresp, ethresh=ethresh)

    _assert_userwarning(expected_warnings, ws)


def test_rmf_with_negative_energy_elem():
    """What happens creating an ARf with negative energies.

    Hopefully we do not have files like this in use.
    """

    # Special case it so that the first in ends at 0 keV
    energy = np.arange(0, 1.0, 0.1, dtype=np.float32)
    energy = energy - energy[1]

    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    with pytest.raises(DataErr,
                       match="The RMF 'delta-rmf' has an ENERG_LO value <= 0"):
        create_delta_rmf(energ_lo, energ_hi)


def test_rmf_with_negative_energy_elem_replace():
    """What happens creating an ARf with negative energies.

    Hopefully we do not have files like this in use.
    """

    ethresh = 0.001

    # Special case it so that the first in ends at 0 keV
    energy = np.arange(0, 1.0, 0.1, dtype=np.float32)
    energy = energy - energy[1]

    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    expected_warnings = [UserWarning("The RMF 'delta-rmf' has an ENERG_LO value < 0")]

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        create_delta_rmf(energ_lo, energ_hi, ethresh=ethresh)

    _assert_userwarning(expected_warnings, ws)


def test_rmf_with_grid_below_thresh():
    """The first bin starts above 0 but ends below the threshold.

    This is a valid grid (other than the fact it is not
    consecutive), so the RMF can be created. See
    test_rmf_with_grid_below_thresh_zero() for the
    case when the bin starts at 0.
    """

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    energ_lo[0] = 1e-7
    energ_hi[0] = 2e-7

    # The test is to make sure that the call does not
    # error out
    rdata = create_delta_rmf(energ_lo, energ_hi, ethresh=1e-5)
    assert isinstance(rdata, DataRMF)


def test_rmf_with_grid_below_thresh_zero():
    """The first bin starts at 0 but ends below the threshold."""

    energy = np.arange(0.0, 1.0, 0.1, dtype=np.float32)
    energ_lo = energy[:-1]
    energ_hi = energy[1:]

    energ_lo[0] = 0.0
    energ_hi[0] = 1e-7

    with pytest.raises(DataErr,
                       match="The RMF 'delta-rmf' has an ENERG_HI value <= " +
                       "the replacement value of 1e-05"):
        create_delta_rmf(energ_lo, energ_hi, ethresh=1e-5)


# Bug https://github.com/sherpa/sherpa/issues/572
@requires_data
@requires_fits
def test_arf_rmf_get_x(make_data_path):
    arf_name = make_data_path('3c120_heg_-1.arf')
    rmf_name = make_data_path('3c120_heg_-1.rmf')

    session = Session()
    arf = session.unpack_arf(arf_name)
    rmf = session.unpack_rmf(rmf_name)

    expected_array_10 = [0.57724115, 0.57730836, 0.57737556, 0.5774428, 0.57751006,
                         0.57757729, 0.57764456, 0.57771185, 0.57777914, 0.57784647]
    actual_arf_10 = arf.get_x()[0:10]
    actual_rmf_10 = rmf.get_x()[0:10]

    np.testing.assert_array_almost_equal(expected_array_10, actual_arf_10)
    np.testing.assert_array_almost_equal(expected_array_10, actual_rmf_10)


def test_arf_get_x_unit():
    session = Session()
    arf_x_lo, arf_x_hi = np.array([12.0, 12.1, 12.2]), np.array([12.1, 12.2, 12.3])
    arf = session.create_arf(arf_x_lo, arf_x_hi)
    expected_arf_x = (arf_x_hi + arf_x_lo)/2
    actual_arf_x = arf.get_x()
    np.testing.assert_array_almost_equal(expected_arf_x, actual_arf_x)


def test_rmf_get_x_unit():
    session = Session()
    rmf_x_lo, rmf_x_hi = np.array([21.0, 21.1, 21.2]), np.array([21.1, 21.2, 21.3])
    rmf = session.create_rmf(rmf_x_lo, rmf_x_hi)
    expected_rmf_x = (rmf_x_hi + rmf_x_lo)/2
    actual_rmf_x = rmf.get_x()
    np.testing.assert_array_almost_equal(expected_rmf_x, actual_rmf_x)


# https://github.com/sherpa/sherpa/pull/766
@requires_group
def test_ungroup():
    '''Make sure that ungrouped data can be ungrouped.

    This test just groups and ungroups a few times.
    '''
    session = Session()
    testdata = DataPHA('testdata', np.arange(50, dtype=float) + 1.,
                       np.zeros(50))
    session.set_data(1, testdata)
    session.ungroup(1)
    session.group_bins(1, 5)
    assert np.all(session.get_grouping(1)[::10] == 1)
    assert testdata.grouped
    # test it can be ungrouped
    session.ungroup(1)
    assert not testdata.grouped
    # test ungrouped data can be ungrouped without altering
    # the grouping
    session.ungroup(1)
    assert not testdata.grouped


# https://github.com/sherpa/sherpa/pull/766
def test_unsubtract():
    '''Make sure that unsubtracted data can be unsubtracted.

    This test just subtracts and unsubtracts a few times.
    '''
    session = Session()
    testdata = DataPHA('testdata', np.arange(1, 51, dtype=np.int16),
                       np.zeros(50))
    testbkg = DataPHA('testbkg', np.arange(1, 51, dtype=np.int16),
                      np.zeros(50))
    session.set_data(1, testdata)
    session.set_bkg(1, testbkg)
    session.unsubtract(1)
    session.subtract(1)
    assert testdata.subtracted
    # test it can be ungrouped
    session.unsubtract(1)
    assert not testdata.subtracted
    # test ungrouped data can be ungrouped without altering
    # the grouping
    session.unsubtract(1)
    assert not testdata.subtracted


@requires_data
@requires_fits
@pytest.mark.parametrize("infile", ["9774.pi", "3c273.pi"])
def test_pha_get_filter_none(infile, make_data_path):
    """Check get_filter with no filter

    It would be nice to do this with faked data, but easier this
    way.
    """

    import sherpa.astro.io

    pha = sherpa.astro.io.read_pha(make_data_path(infile))

    assert pha.get_filter(format='%.4f') == '0.0015:14.9504'


@requires_data
@requires_fits
@pytest.mark.parametrize("infile,expected",
                         [("9774.pi", '0.4964:1.9856,3.0076:7.0080'),
                          ("3c273.pi", '0.4672:1.9418,3.0806:9.8696')])
def test_pha_get_filter_filter(infile, expected, make_data_path):
    """Check get_filter with simple-ish filter

    It would be nice to do this with faked data, but easier this
    way.
    """

    import sherpa.astro.io

    pha = sherpa.astro.io.read_pha(make_data_path(infile))
    pha.notice(0.5, 7)
    pha.ignore(2, 3)

    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("infile,expected",
                         [("9774.pi", '0.4964:0.6132'),
                          ("3c273.pi", '0.4672:0.6424')])
def test_pha_get_filter_edgecase(infile, expected, make_data_path):
    """Check get_filter with an edge case

    Pick something that has caused problems (related to #917).
    """

    import sherpa.astro.io

    pha = sherpa.astro.io.read_pha(make_data_path(infile))
    pha.notice(0.501, 0.6)

    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("infile,expected",
                         [("9774.pi", '0.4964:7.0080'),
                          ("3c273.pi", '0.4672:9.8696')])
def test_pha_get_filter_false(infile, expected, make_data_path):
    """get_filter(group=False) with no grouping."""

    import sherpa.astro.io

    pha = sherpa.astro.io.read_pha(make_data_path(infile))
    pha.notice(0.5, 7)

    assert pha.get_filter(group=False, format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("infile", ["9774.pi", "3c273.pi"])
def test_pha_mask_default(infile, make_data_path):
    """Make sure we have some tests (these may exist elsewhere)"""

    import sherpa.astro.io

    pha = sherpa.astro.io.read_pha(make_data_path(infile))
    assert pha.mask


@requires_data
@requires_fits
@pytest.mark.parametrize("infile,size,nset",
                         [("9774.pi", 1024, 376),
                          ("3c273.pi", 46, 33)])
def test_pha_mask_filtered(infile, size, nset, make_data_path):
    """Make sure we have some tests (these may exist elsewhere)"""

    import sherpa.astro.io

    pha = sherpa.astro.io.read_pha(make_data_path(infile))
    pha.notice(0.5, 7)
    pha.ignore(2, 3)

    assert pha.mask.dtype == bool
    assert pha.mask.size == size
    assert pha.mask.sum() == nset


def test_sum_background_data_missing():
    """Check we error out if there's no background data"""

    d = DataPHA('tmp', np.arange(3), np.arange(3))
    with pytest.raises(DataErr,
                       match="data set 'tmp' does not have any associated backgrounds"):
        d.sum_background_data()


@requires_data
@requires_fits
def test_get_filter_channel_ungrouped(make_data_path):
    """What does get_filter return for ungrouped channel data.

    This should create the PHA but easier to use a file.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('9774.pi'))
    assert not pha.grouped
    pha.set_analysis('channel')

    assert pha.get_filter() == '1:1024'

    pha.ignore(400, 500)
    assert pha.get_filter() == '1:399,501:1024'


@requires_data
@requires_fits
@requires_group
def test_get_filter_channel_grouped(make_data_path):
    """What does get_filter return for grouped channel data.

    This is related to bug #920

    This should create the PHA but easier to use a file.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    assert pha.grouped
    pha.set_analysis('channel')

    # This returns channels (now)
    assert pha.get_filter() == '1:1024'

    # Reset the grouping to use an easier-to-check scheme: groups
    # have a fixed number of channels, in this case 50.
    #
    pha.group_width(50)
    assert pha.get_filter() == '1:1024'

    # What units does ignore use? It appears to be channels.
    pha.ignore(151, 300)
    assert pha.get_filter() == '1:150,301:1024'


@requires_data
@requires_fits
@requires_group
def test_get_filter_channel_grouped_prefiltered(make_data_path):
    """Add an energy filter before switching to channel space

    This is related to bug #920
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.group_width(50)

    # notice the 1.0-7 keV range, which is (ungrouped)
    # channels 69 - 480, so with groups of width 50
    # this is groups 2 - 10.
    pha.notice(1.0, 7.0)

    pha.set_analysis('channel')

    assert pha.get_filter() == '51:500'

    # What units does ignore use? It appears to be channels.
    pha.ignore(150, 300)
    assert pha.get_filter() == '51:100,301:500'


@requires_data
@requires_fits
@requires_group
@pytest.mark.parametrize("analysis", ["energy", "wavelength", "channel"])
def test_grouping_nofilter(analysis, make_data_path):
    """Can we change grouping (no filter).

    This is related to bug #920

    This should create the PHA but easier to use a file.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    # The expected dependent axis
    dep = np.array([17, 15, 16, 15, 16, 15, 18, 18, 15, 18,
                    15, 15, 19, 15, 15, 17, 16, 16, 17, 15,
                    19, 15, 16, 15, 16, 17, 15, 18, 16, 15,
                    15, 16, 15, 15, 15, 16, 16, 15, 15, 16,
                    16, 15, 16, 15, 15, 20])
    assert pha.get_dep(filter=True) == pytest.approx(dep)

    pha.set_analysis(analysis)
    pha.group_width(50)
    dep = np.array([105, 213, 136,  79,  47,  47,  29,  27,
                    18, 12, 0, 2, 0, 1, 3, 3, 1, 2, 1, 2, 8])
    assert pha.get_dep(filter=True) == pytest.approx(dep)


@requires_data
@requires_fits
def test_get_filter_group_bug(make_data_path):
    """This should be the same problem as seen in
    test_grouping_filter with analysis=channel.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.ignore(None, 1)
    pha.ignore(7, None)

    en1 = 1.036599993706
    en2 = 6.570000171661

    n1 = 71
    n2 = 449

    # this is to check we get the expected results
    elo, ehi = pha._get_ebins(group=False)
    assert elo[n1] == pytest.approx(en1)
    assert ehi[n2] == pytest.approx(en2)

    filters1 = parse_expr(pha.get_filter(group=False))
    assert len(filters1) == 1
    assert len(filters1[0]) == 2
    assert filters1[0][0] == pytest.approx(en1)
    assert filters1[0][1] == pytest.approx(en2)

    pha.set_analysis('channel')

    # The filter is in channel units, which is 1 + n
    # when n is used to access elements of elo/ehi).
    #
    filters2 = parse_expr(pha.get_filter(group=False))
    assert len(filters2) == 1
    assert len(filters2[0]) == 2
    assert filters2[0][0] == n1 + 1
    assert filters2[0][1] == n2 + 1


@requires_data
@requires_fits
@pytest.mark.parametrize("analysis", ["energy", "wavelength", "channel"])
def test_get_noticed_channels(analysis, make_data_path):
    """Check get_noticed_channels when analysis=channel.

    This was used when tracking down bug #920
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.ignore(None, 1)
    pha.ignore(2, 3)
    pha.ignore(7, None)

    pha.set_analysis(analysis)

    expected = np.concatenate((np.arange(72, 134), np.arange(212, 451)))
    assert pha.get_noticed_channels() == pytest.approx(expected)


@requires_data
@requires_fits
@requires_group
@pytest.mark.parametrize("analysis", ["energy", "wavelength", "channel"])
def test_grouping_filter(analysis, make_data_path):
    """Can we change grouping with energy units.

    This is related to bug #920
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    assert pha.get_analysis() == 'energy'
    pha.ignore(hi=1)
    pha.ignore(lo=7)

    pha.set_analysis(analysis)

    dep = np.array([15, 17, 16, 16, 17, 15, 19, 15, 16, 15,
                    16, 17, 15, 18, 16, 15, 15, 16, 15, 15,
                    15, 16, 16, 15, 15, 16, 16, 15, 16, 15])
    assert pha.get_dep(filter=True) == pytest.approx(dep)

    # The group mapping for group_width of 50 is listed in
    # test_grouping_filtering_binning. However, as of 4.16.0 the
    # group_width call will now use the inverse of the current mask as
    # the tabStops array. As of 4.16, this only groups over the range
    # pha.counts[71:450].
    #
    pha.group_width(50)
    dep = np.array([136, 129, 61, 48, 35, 35, 19, 11])
    assert pha.get_dep(filter=True) == pytest.approx(dep)

    qual = np.zeros(1024, dtype=int)
    qual[421:450] = 2
    assert pha.quality == pytest.approx(qual)

    # Prior to 4.16.0 the call would have acted like follows
    #
    pha.group_width(50, tabStops=[0] * 1024)
    dep = np.array([213, 136,  79,  47,  47,  29,  27, 18])
    assert pha.get_dep(filter=True) == pytest.approx(dep)

    qual = np.zeros(1024, dtype=int)
    qual[1000:1024] = 2
    assert pha.quality == pytest.approx(qual)


@requires_data
@requires_fits
@requires_group
@pytest.mark.parametrize("analysis", ["energy", "wavelength", "channel"])
def test_grouping_filtering_binning(analysis, make_data_path):
    """Low-level testing of test_grouping_filtering.

    Check that the grouping has created the results we
    expect.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.ignore(hi=1)
    pha.ignore(lo=7)

    # As of 4.16.0, the tabStops defaults to the inverse of
    # pha.get_mask, so to avoid changing this test from the original
    # results we need to override the setting.
    #
    pha.set_analysis(analysis)
    pha.group_width(50, tabStops=np.zeros(1024, dtype=np.int8))

    # We expect 1, 49 * -1, repeated and then the last bin.
    #
    gbin = [1] + [-1] * 49
    gend = [1] + [-1] * 23
    expected = np.concatenate((np.tile(gbin, 20), gend))
    assert (pha.grouping == expected).all()

    # This is based on the energy results. The grouping has
    # energy ranges:
    #
    #   i= 0   0.0015 -  0.7300
    #   i= 1   0.7300 -  1.4600
    #   i= 2   1.4600 -  2.1900
    #   i= 3   2.1900 -  2.9200
    #   i= 4   2.9200 -  3.6500
    #   i= 5   3.6500 -  4.3800
    #   i= 6   4.3800 -  5.1100
    #   i= 7   5.1100 -  5.8400
    #   i= 8   5.8400 -  6.5700
    #   i= 9   6.5700 -  7.3000
    #   i=10   7.3000 -  8.0300
    #   i=11   8.0300 -  8.7600
    #   i=12   8.7600 -  9.4900
    #   i=13   9.4900 - 10.2200
    #   i=14  10.2200 - 10.9500
    #   i=15  10.9500 - 11.6800
    #   i=16  11.6800 - 12.4100
    #   i=17  12.4100 - 13.1400
    #   i=18  13.1400 - 13.8700
    #   i=19  13.8700 - 14.6000
    #   i=20  14.6000 - 14.9504
    #
    expected = np.zeros(21, dtype=bool)
    expected[1:9] = True
    assert (pha.mask == expected).all()

    # For the ungrouped-data we have, selecting
    # a few ranges related to how the code could
    # work:
    #
    #   i=49   0.7154 -  0.7300
    #   i=50   0.7300 -  0.7446
    #   i=51   0.7446 -  0.7592
    #   ...
    #   i=68   0.9928 -  1.0074
    #   i=69   1.0074 -  1.0220
    #   i=70   1.0220 -  1.0366
    #   ...
    #   i=99   1.4454 -  1.4600
    #   i=100   1.4600 -  1.4746
    #   i=101   1.4746 -  1.4892
    #   ...
    #   i=448   6.5408 -  6.5554
    #   i=449   6.5554 -  6.5700
    #   i=450   6.5700 -  6.5846
    #   ...
    #   i=478   6.9788 -  6.9934
    #   i=479   6.9934 -  7.0080
    #   i=480   7.0080 -  7.0226
    #
    expected = np.zeros(1024, dtype=bool)
    expected[50:450] = True
    assert (pha.get_mask() == expected).all()


@requires_data
@requires_fits
@requires_group
@pytest.mark.parametrize("analysis", ["energy", "wavelength", "channel"])
def test_grouping_filtering_binning_416(analysis, make_data_path):
    """Low-level testing of test_grouping_filtering without tabStops override

    In 4.16 the behavior of the grpoup_xxx call changed to include
    the existing mask as a tabStops array (after  inversion) which
    changes the results of this test, so we now have two versions

      test_grouping_filtering_binning
      test_grouping_filtering_binning_416

    where the former keeps the pre-4.16 behavior (by setting tabStops
    to [0] * 1024) and this one, which checks the new behavior.

    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.ignore(hi=1)
    pha.ignore(lo=7)

    # As of 4.16.0, the tabStops defaults to the inverse of
    # pha.get_mask, so to avoid changing this test from the original
    # results we need to override the setting.
    #
    pha.set_analysis(analysis)
    pha.group_width(50)

    # The valid channel range is (starting at 0), 71 to 449
    # (inclusive), which is 379 channels, which means 8 50-channel
    # groups (the last partly filled).
    #
    gbin = [1] + [-1] * 49
    gend = [1] + [-1] * 28
    expected = np.concatenate(([0] * 71, np.tile(gbin, 7),
                               gend, [0] * 574))
    assert (pha.grouping == expected).all()

    expected = np.zeros(653, dtype=bool)
    expected[71:79] = True
    assert (pha.mask == expected).all()

    expected = np.zeros(1024, dtype=bool)
    expected[71:450] = True
    assert (pha.get_mask() == expected).all()


# The channel ranges associated with the groups can be
# found with:
#   pha.apply_grouping(pha.channel, pha._min)
#   pha.apply_grouping(pha.channel, pha._max)
#
# For the 3c273.pi file we have:
#
# group  1  2  3  4  5  6  7  8  9 10 11 12 13 14
# low    1 18 22 33 40 45 49 52 55 57 60 62 66 69
# high  17 21 32 39 44 48 51 54 56 59 61 65 68 71
#
# group 15 16 17 18 19  20  21  22  23  24  25  26
# low   72 76 79 83 89  97 102 111 117 125 131 134
# high  75 78 82 88 96 101 110 116 124 130 133 139
#
# group  27  28  29  30  31  32  33  34  35  36
# low   140 144 151 157 165 178 187 197 212 233
# high  143 150 156 164 177 186 196 211 232 244
#
# group  37  38  39  40  41  42  43  44  45  46
# low   245 261 277 292 324 345 369 405 451 677
# high  260 276 291 323 344 368 404 450 676 1024
#
# So the full dataspace for groups 1:46 is (mid-points)
#      (1 + 17) / 2  -  (677 + 1024) / 2
#                 9  -  850.5
#                 9  -  850  (rounding down)
#
# Channel 70  lies in group 14 which covers  69 -  71
#         390               43              369 - 404
#
# For reference:
#   channel  70: 1.0147 keV (1.0074 - 1.0220)
#           390: 5.6867      5.6794 - 5.6940
#
# and the actual edges of the groups
#            69: 1.0001 keV
#           404: 5.8911
#
# group 14: 1.0147 keV = (elo[68] + ehi[70]) / 2
# group 15: 1.0658     = (elo[71] + ehi[74]) / 2
# group 42: 5.1976     = (elo[364] + ehi[367]) / 2
# group 43: 5.6356     = (elo[368] + ehi[403]) / 2
#
# where
#    pha.set_analysis('energy')
#    elo, ehi = pha._get_ebins(group=False)
#

@requires_data
@requires_fits
def test_notice_energy_grouping(make_data_path):
    """Check that filtering selects the right bins with channels"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    # This is groups 14 - 43.
    #
    # I have chosen the mid-points of the bin so it
    # should match the get_filter call (prior to Sherpa 4.14.0).
    #
    pha.notice(1.0147, 5.6356)

    # So, does the mask include the start/end bins
    # which we overlap, or not? It does - we have
    # channels 69 - 404, which maps to indices 68,..,403
    # or 68:404.
    #
    mask = np.zeros(1024, dtype=bool)
    mask[68:404] = True

    assert pha.get_mask() == pytest.approx(mask)

    expected = '0.9928:5.8984'
    assert pha.get_filter(format='%.4f') == expected

    # This gives the mid-points of the first and last channels
    # covered by the groups, so channels 69 and 404.
    #
    assert pha.get_filter(format='%.4f', group=False) == expected


@requires_data
@requires_fits
def test_notice_channel_grouping(make_data_path):
    """Check that filtering selects the right bins with channels"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.set_analysis('channel')

    # this covers the same groups as the 1.0147 - 5.6356 keV
    # call in test_notice_energy_grouping (14-43).
    #
    pha.notice(70, 390)

    # This now creates the expected range.
    #
    mask = np.zeros(1024, dtype=bool)
    mask[68:404] = True

    assert pha.get_mask() == pytest.approx(mask)

    # Returns the channel numbers (groups 14 to 43
    # is 69 - 404).
    assert pha.get_filter(format='%.4f') == '69:404'
    assert pha.get_filter(format='%.4f', group=False) == '69:404'


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(-5, 2000, '1:1024'),
                          (30, 2000, '22:1024'),
                          (-5, 350, '1:368'),
                          (-20, -5, ''),
                          (2000, 3000, '')])
def test_notice_channel_grouping_outofbounds(lo, hi, expected, make_data_path):
    """Check what happens with silly results

    Note that filters outside the range end up including the first or
    last bin (because it gets reset to the limit), which we
    probably do not want (for the low edge, for the upper
    edge we are probably lucky due to < rather than >=).

    The groups are such that the first group has channels
    1-17 and the last group 677-1024, which have mid-points
    9 and 850.5, hence the 9:850 as the default filter.

    >>> pha.apply_grouping(pha.channel, pha._min)
    array([  1.,  18.,  22.,  33.,  40.,  45.,  49.,  52.,  55.,  57.,  60.,
            62.,  66.,  69.,  72.,  76.,  79.,  83.,  89.,  97., 102., 111.,
           117., 125., 131., 134., 140., 144., 151., 157., 165., 178., 187.,
           197., 212., 233., 245., 261., 277., 292., 324., 345., 369., 405.,
           451., 677.])

    >>> pha.apply_grouping(pha.channel, pha._max)
    array([  17.,   21.,   32.,   39.,   44.,   48.,   51.,   54.,   56.,
             59.,   61.,   65.,   68.,   71.,   75.,   78.,   82.,   88.,
             96.,  101.,  110.,  116.,  124.,  130.,  133.,  139.,  143.,
            150.,  156.,  164.,  177.,  186.,  196.,  211.,  232.,  244.,
            260.,  276.,  291.,  323.,  344.,  368.,  404.,  450.,  676.,
           1024.])

    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.set_analysis('channel')

    pha.notice(lo, hi)
    assert pha.get_filter() == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(-5, 2000, '1:1024'),
                          (30, 2000, '30:1024'),
                          (-5, 350, '1:350'),
                          (-20, -5, ''),
                          (2000, 3000, '')])
def test_notice_channel_grouping_outofbounds_ungrouped(lo, hi, expected, make_data_path):
    """Check what happens with silly results
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.ungroup()
    pha.set_analysis('channel')

    pha.notice(lo, hi)
    assert pha.get_filter() == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(0, 2000, '0.0015:14.9504'),
                          (0.7, 2000, '0.6424:14.9504'),
                          (0, 4.2, '0.0015:4.2486'),
                          (0, 1e-10, ''),
                          (2000, 3000, '')])
def test_notice_energy_grouping_outofbounds(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.set_analysis('energy')

    pha.notice(lo, hi)
    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(0, 2000, '0.0015:14.9504'),
                          (0.7, 2000, '0.6862:14.9504'),
                          (0, 4.2, '0.0015:4.2048'),
                          (0, 1e-10, ''),
                          (2000, 3000, '')])
def test_notice_energy_grouping_outofbounds_ungrouped(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.ungroup()
    pha.set_analysis('energy')

    pha.notice(lo, hi)
    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(0, 8000, '0.8293:8492.0673'),
                          (20, 8000, '19.3002:8492.0673'),
                          (0, 15, '0.8293:15.1644'),
                          (0, 1e-10, '0.8293:8492.0673'),
                          (8000, 9000, '49.9533:8492.0673')])
def test_notice_wave_grouping_outofbounds(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.set_analysis('wave')

    pha.notice(lo, hi)
    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(0, 8000, '0.8293:8492.0673'),
                          (20, 8000, '19.7490:8492.0673'),
                          (0, 15, '0.8293:15.1644'),
                          (0, 1e-10, '0.8293:8492.0673'),
                          (8000, 9000, '849.2067:8492.0673')])
def test_notice_wave_grouping_outofbounds_ungrouped(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.ungroup()
    pha.set_analysis('wave')

    pha.notice(lo, hi)
    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(-5, 2000, ''),
                          (30, 2000, '1:21'),
                          (-5, 350, '369:1024'),
                          (-20, -5, '1:1024'),
                          (2000, 3000, '1:1024')])
def test_ignore_channel_grouping_outofbounds(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.set_analysis('channel')

    pha.ignore(lo, hi)
    assert pha.get_filter() == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(-5, 2000, ''),
                          (30, 2000, '1:29'),
                          (-5, 350, '351:1024'),
                          (-20, -5, '1:1024'),
                          (2000, 3000, '1:1024')])
def test_ignore_channel_grouping_outofbounds_ungrouped(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.ungroup()
    pha.set_analysis('channel')

    pha.ignore(lo, hi)
    assert pha.get_filter() == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(0, 2000, ''),
                          (0.8, 2000, '0.0015:0.7884'),
                          (0, 3.5, '3.5624:14.9504'),
                          (0, 1e-10, '0.0015:14.9504'),
                          (2000, 3000, '0.0015:14.9504')])
def test_ignore_energy_grouping_outofbounds(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.set_analysis('energy')

    pha.ignore(lo, hi)
    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(0, 2000, ''),
                          (0.8, 2000, '0.0015:0.7884'),
                          (0, 3.5, '3.5040:14.9504'),
                          (0, 1e-10, '0.0015:14.9504'),
                          (2000, 3000, '0.0015:14.9504')])
def test_ignore_energy_grouping_outofbounds_ungrouped(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.ungroup()
    pha.set_analysis('energy')

    pha.ignore(lo, hi)
    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(0, 2000, ''),
                          (20, 2000, '0.8293:19.3002'),
                          (0, 15, '15.1644:8492.0673'),
                          (0, 1e-10, '0.8293:8492.0673'),
                          (2000, 3000, '0.8293:49.9533')])
def test_ignore_wave_grouping_outofbounds(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.set_analysis('wave')

    pha.ignore(lo, hi)
    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("lo,hi,expected",
                         [(0, 2000, ''),
                          (20, 2000, '0.8293:19.7490'),
                          (0, 15, '15.1644:8492.0673'),
                          (0, 1e-10, '0.8293:8492.0673'),
                          (2000, 3000, '0.8293:849.2067')])
def test_ignore_wave_grouping_outofbounds_ungrouped(lo, hi, expected, make_data_path):
    """Check what happens with silly results"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.ungroup()
    pha.set_analysis('wave')

    pha.ignore(lo, hi)
    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("units", ["energy", "wave"])
@pytest.mark.parametrize("notice", [True, False])
@pytest.mark.parametrize("lo,hi",
                         [(-5, -2), (-5, 2), (-5, None), (0, -2), (None, -2)])
def test_pha_validates_limits(units, notice, lo, hi, make_data_path):
    """Check the limits are validated."""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.units = units
    func = pha.notice if notice else pha.ignore

    with pytest.raises(DataErr,
                       match="^unknown "):
        func(lo, hi)


@requires_data
@requires_fits
@pytest.mark.parametrize("notice", [True, False])
def test_pha_validates_range(notice, make_data_path):
    """Ensure lo <= hi check is made"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    func = pha.notice if notice else pha.ignore

    with pytest.raises(DataErr,
                       match="unknown hi argument: 'must be >= lo'"):
        func(3, 2)


@requires_data
@requires_fits
def test_channel_changing_limits(make_data_path):
    """Test behavior seen while fixing #920

    Check we can ungroup and regroup and get back
    to where we were.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    pha.set_analysis('channel')

    assert pha.mask is True

    # selects
    #    group 11 (channels 60-61, mid=60)
    # to
    #    group 42 (channsls 345-368, mid=356)
    #
    pha.notice(60, 350)

    mexpected1 = np.zeros(46, dtype=bool)
    mexpected1[10:42] = True
    mexpected2 = np.zeros(1024, dtype=bool)
    mexpected2[59:368] = True
    assert pha.mask == pytest.approx(mexpected1)
    assert pha.get_mask() == pytest.approx(mexpected2)

    expected = '60:368'
    assert pha.get_filter() == expected
    assert pha.get_filter(group=False) == expected

    # We now return the full filter range of the
    # group, even when group=False.
    #
    pha.ungroup()
    assert pha.mask == pytest.approx(mexpected2)
    assert pha.get_mask() == pytest.approx(mexpected2)
    assert pha.get_filter() == expected
    assert pha.get_filter(group=False) == expected

    # We go back to the original filter
    pha.group()
    assert pha.mask == pytest.approx(mexpected1)
    assert pha.get_mask() == pytest.approx(mexpected2)
    assert pha.get_filter() == expected
    assert pha.get_filter(group=False) == expected


@requires_data
@requires_fits
def test_energy_filter_notice_ignore(make_data_path):
    """Add a simple check of notice/ignore handling.

    This may be tested elsewhere  notice then ignore - but make
    sure we have it covered.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.set_analysis('energy')

    pha.notice(0.5, 7)
    pha.ignore(1, 2)

    # Use pytest.approx to make it an easy check, and as values are True/False
    # it should be sufficient.
    #
    expected = np.zeros(46, dtype=bool)
    expected[3:13] = True
    expected[26:45] = True
    assert pha.mask == pytest.approx(expected)

    # This is now the ungrouped mask
    expected = np.zeros(1024, dtype=bool)
    expected[32:68] = True
    expected[139:676] = True
    assert pha.get_mask() == pytest.approx(expected)


@requires_data
@requires_fits
def test_energy_filter_notice_ignore_ungrouped(make_data_path):
    """Add a simple check of notice/ignore handling.

    This may be tested elsewhere  notice then ignore - but make
    sure we have it covered.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.set_analysis('energy')

    pha.notice(0.5, 7)
    pha.ignore(1, 2)

    pha.ungroup()

    # Check the two masks are the same
    assert pha.mask == pytest.approx(pha.get_mask())

    expected = np.zeros(1024, dtype=bool)
    expected[32:68] = True
    expected[139:676] = True
    assert pha.get_mask() == pytest.approx(expected)


@requires_data
@requires_fits
def test_energy_filter_ungrouped_notice_ignore_ungrouped(make_data_path):
    """Add a simple check of notice/ignore handling.

    This is similar to test_energy_filter_notice_ignore_ungrouped
    but calls ungroup first.
    """

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    pha.set_analysis('energy')

    pha.ungroup()

    pha.notice(0.5, 7)
    pha.ignore(1, 2)

    # Check the two masks are the same
    assert pha.mask == pytest.approx(pha.get_mask())

    expected = np.zeros(1024, dtype=bool)
    expected[34:68] = True
    expected[137:480] = True
    assert pha.get_mask() == pytest.approx(expected)


@requires_data
@requires_fits
def test_energy_filter_roundtrip(make_data_path):
    """If you ungroup and then group you should get back the original filter."""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))

    fall = pha.get_filter(format='%.5f')
    assert fall == '0.00146:14.95040'

    pha.notice(0.5, 7)
    pha.ignore(1, 2)

    expected = '0.46720:0.99280,2.02940:9.86960'
    f1 = pha.get_filter(format='%.5f')
    assert f1 == expected

    pha.ungroup()
    f2 = pha.get_filter(format='%.5f')
    assert f2 == expected

    pha.group()
    f3 = pha.get_filter(format='%.5f')
    assert f3 == expected


@pytest.mark.xfail
@requires_data
@requires_fits
def test_energy_filter_ordering(make_data_path):
    """If you ungroup then filter or filter than ungroup you should get back the same filter."""

    from sherpa.astro.io import read_pha

    infile = make_data_path('3c273.pi')
    pha1 = read_pha(infile)
    pha2 = read_pha(infile)

    pha2.ungroup()
    for pha in [pha1, pha2]:
        pha.notice(0.5, 6)
        pha.ignore(1.1, 1.2)
        pha.ignore(4, 4.1)

    pha1.ungroup()

    assert pha1.get_filter() == pha2.get_filter()

    mask1 = pha1.mask
    mask2 = pha2.mask
    assert mask1 == pytest.approx(mask2)

    mask1 = pha1.get_mask()
    mask2 = pha2.get_mask()
    assert mask1 == pytest.approx(mask2)


@pytest.mark.parametrize('units', ['bin', 'chan', 'energy', 'wave'])
def test_pha_get_ebins_internal_no_response(units):
    """Check that _get_ebins has an unlikely-used path checked.

    It's  not clear what we are meant to return here - i.e. no
    response but a units value has been set - and maybe there
    should be an error instead (for non-channel settings).
    """

    chans = np.arange(1, 10, dtype=np.int16)
    counts = np.ones(9, dtype=np.int16)
    pha = DataPHA('tst', chans, counts)
    pha.units = units
    lo, hi = pha._get_ebins()
    assert lo == pytest.approx(chans)
    assert hi == pytest.approx(chans + 1)


def test_get_background_scale_is_none():
    """We get None when there's no background"""

    d = DataPHA('tmp', np.arange(3), np.arange(3))
    assert d.get_background_scale() is None


@pytest.mark.parametrize("is_bkg,option,expected",
                         [(True, 'exposure', 80.0),
                          (False, 'exposure', 0.002),
                          (True, 'backscal', 0.04),
                          (False, 'backscal', 2.0),
                          (True, 'areascal', 0.1),
                          (False, 'areascal', 0.25)
                          ])
def test_get_background_scale_missing_option(is_bkg, option, expected):
    """Check we can calculate the scaling when an option isn't present.

    We calculate ratio(EXPOSURE) * ratio(BACKSCAL) * ratio(AREASCAL)
    and allow values to not be present (the default to 1).

    """

    d = DataPHA('tmp', np.arange(3), np.arange(3))
    b = DataPHA('tmp', np.arange(3), np.arange(3))
    d.set_background(b)

    d.exposure = 100.0
    d.backscal = 0.1
    d.areascal = 0.8

    b.exposure = 400.0
    b.backscal = 0.2
    b.areascal = 0.5

    obj = b if is_bkg else d
    setattr(obj, option, None)

    scale = d.get_background_scale()
    assert scale == pytest.approx(expected)


@pytest.mark.parametrize("lo,hi,emsg", [("1:20", None, 'lower'), (None, "2", 'upper'), ("0.5", "7", 'lower')])
@pytest.mark.parametrize("ignore", [False, True])
def test_pha_notice_errors_out_on_string_range(lo, hi, emsg, ignore):
    """Check we get an error if lo or hi are strings."""

    d = DataPHA('tmp', np.arange(3), np.arange(3))
    with pytest.raises(DataErr,
                       match=f'strings not allowed in {emsg} bound list'):
        d.notice(lo, hi, ignore=ignore)


def test_pha_creation_warns_about_non_numpy_channel():
    """What happens if the channel array is not a NumPy array?

    There used to be a UserWarning but this has now been removed.
    """

    chans = [1, 2, 3]
    counts = np.ones(3)
    d = DataPHA('tmp', chans, counts)

    assert isinstance(d.x, np.ndarray)
    assert isinstance(d.channel, np.ndarray)

    assert d.channel == pytest.approx(chans)


def test_pha_creation_warns_about_non_numpy_counts():
    """What happens if the counts array is not a NumPy array?

    There used to be a UserWarning but this has now been removed.
    """

    chans = np.asarray([1, 2, 3])
    counts = [1, 0, 1]
    d = DataPHA('tmp', chans, counts)

    assert isinstance(d.y, np.ndarray)
    assert isinstance(d.counts, np.ndarray)

    assert d.counts == pytest.approx(counts)


@requires_fits
@requires_data
def test_pha_check_filter(make_data_path):
    """Added a test found useful when changing get_filter."""

    import sherpa.astro.io

    infile = make_data_path('3c273.pi')
    pha = sherpa.astro.io.read_pha(infile)

    pha.notice(0.5, 7)
    assert pha.get_filter(format='%.4f') == '0.4672:9.8696'

    pha.ignore(None, 1)
    assert pha.get_filter(format='%.4f') == '1.0366:9.8696'

    pha.ignore(5, None)
    assert pha.get_filter(format='%.4f') == '1.0366:4.7158'

    plot = pha.to_plot()
    assert plot[0].size == 26

    pha.ungroup()

    assert pha.get_filter(format='%.4f') == '1.0366:4.7158'

    plot = pha.to_plot()
    assert plot[0].size == 252


@requires_fits
@requires_data
def test_pha_check_filter_channel(make_data_path):
    """test_pha_check_filter with channel units"""

    import sherpa.astro.io

    infile = make_data_path('3c273.pi')
    pha = sherpa.astro.io.read_pha(infile)
    pha.units = 'channel'

    # The data is grouped so it doesn't quite match.
    # See the comments above test_notice_energy_grouping
    # for the channel/group/energy range mapping for
    # this file.
    #
    pha.notice(35, 480)
    assert pha.get_filter(format='%i') == '33:676'
    pha.units = 'energy'
    assert pha.get_filter(format='%.4f') == '0.4672:9.8696'
    pha.units = 'channel'

    pha.ignore(None, 69)
    assert pha.get_filter(format='%i') == '72:676'
    pha.units = 'energy'
    assert pha.get_filter(format='%.4f') == '1.0366:9.8696'
    pha.units = 'channel'

    pha.ignore(343, None)
    assert pha.get_filter(format='%i') == '72:323'
    pha.units = 'energy'
    assert pha.get_filter(format='%.4f') == '1.0366:4.7158'
    pha.units = 'channel'

    plot = pha.to_plot()
    assert plot[0].size == 26

    pha.ungroup()

    assert pha.get_filter(format='%i') == '72:323'
    pha.units = 'energy'
    assert pha.get_filter(format='%.4f') == '1.0366:4.7158'
    pha.units = 'channel'

    plot = pha.to_plot()
    assert plot[0].size == 252


@pytest.mark.parametrize('ignore,lo,hi,expected',
                         [(False, 0, 0, None),
                          (True, 0, 0, True),    # why is this not 'not ignore'?
                          (False, 0, None, None),
                          (True, 0, None, None),
                          (False, None, 0, None),
                          (True, None, 0, None)])
def test_pha_filter_wave_0limit(ignore, lo, hi, expected):
    """Edge case checks for wavelength handling: 0 values

    Explicit handling of the wavelength handling when a filter limit
    is 0. It's not really obvious what we want these limits to mean
    so I am treating this as a regression test.

    """

    chans = np.arange(1, 7, dtype=int)
    pha = DataPHA('d', chans, np.zeros_like(chans))

    rmf = create_delta_rmf(chans, chans + 1,
                           e_min=chans, e_max=chans + 1)
    pha.set_rmf(rmf)
    pha.units = 'wave'

    assert pha.mask is True

    func = pha.ignore if ignore else pha.notice
    func(lo, hi)

    if expected is None:
        assert pha.mask is not ignore
    else:
        assert pha.mask == pytest.approx(expected)


def test_pha_filter_simple_channel1():
    """Simple tests of get_filter

    See also test_pha_filter_simple_energy
    """

    chans = np.arange(1, 7, dtype=int)
    pha = DataPHA('d', chans, np.zeros_like(chans))

    assert pha.get_filter() == '1:6'

    # Fake filter/notice calls
    pha.mask = np.ones(6, dtype=bool)
    assert pha.get_filter() == '1:6'

    pha.mask[1] = False
    pha.mask[4] = False
    assert pha.get_filter() == '1,3:4,6'

    pha.mask = ~pha.mask
    assert pha.get_filter() == '2,5'

    pha.mask = np.zeros(6, dtype=bool)
    assert pha.get_filter() == ''

    pha.mask[0] = True
    assert pha.get_filter() == '1'

    pha.mask[-1] = True
    assert pha.get_filter() == '1,6'

    pha.mask[0] = False
    assert pha.get_filter() == '6'


def test_pha_filter_simple_channel0():
    """Simple tests of get_filter

    See also test_pha_filter_simple_energy
    """

    chans = np.arange(0, 6, dtype=int)
    pha = DataPHA('d', chans, np.zeros_like(chans))

    assert pha.get_filter() == '0:5'

    # Fake filter/notice calls
    pha.mask = np.ones(6, dtype=bool)
    assert pha.get_filter() == '0:5'

    pha.mask[1] = False
    pha.mask[4] = False
    assert pha.get_filter() == '0,2:3,5'

    pha.mask = ~pha.mask
    assert pha.get_filter() == '1,4'

    pha.mask = np.zeros(6, dtype=bool)
    assert pha.get_filter() == ''

    pha.mask[0] = True
    assert pha.get_filter() == '0'

    pha.mask[-1] = True
    assert pha.get_filter() == '0,5'

    pha.mask[0] = False
    assert pha.get_filter() == '5'


def test_pha_filter_simple_energy1():
    """Simple tests of get_filter

    See also test_pha_filter_simple_channel1
    """

    chans = np.arange(1, 7, dtype=int)
    pha = DataPHA('d', chans, np.zeros_like(chans))

    rmf = create_delta_rmf(chans, chans + 1,
                           e_min=chans, e_max=chans + 1)
    pha.set_rmf(rmf)
    pha.units = 'energy'

    assert pha.get_filter(format='%.1f') == '1.0:7.0'

    # Fake filter/notice calls
    pha.mask = np.ones(6, dtype=bool)
    assert pha.get_filter(format='%.1f') == '1.0:7.0'

    pha.mask[1] = False
    pha.mask[4] = False
    assert pha.get_filter(format='%.1f') == '1.0:2.0,3.0:5.0,6.0:7.0'

    pha.mask = ~pha.mask
    assert pha.get_filter(format='%.1f') == '2.0:3.0,5.0:6.0'

    pha.mask = np.zeros(6, dtype=bool)
    assert pha.get_filter(format='%.1f') == ''

    pha.mask[0] = True
    assert pha.get_filter(format='%.1f') == '1.0:2.0'

    pha.mask[-1] = True
    assert pha.get_filter(format='%.1f') == '1.0:2.0,6.0:7.0'

    pha.mask[0] = False
    assert pha.get_filter(format='%.1f') == '6.0:7.0'


def test_pha_filter_simple_energy0():
    """Simple tests of get_filter

    See also test_pha_filter_simple_channel0
    """

    chans = np.arange(0, 6, dtype=int)
    pha = DataPHA('d', chans, np.zeros_like(chans))

    # use integer bins as easy to check but ensure
    # the first bin is not 0
    rmf = create_delta_rmf(chans + 10, chans + 11, offset=0,
                           e_min=chans + 10, e_max=chans + 11)
    pha.set_rmf(rmf)
    pha.units = 'energy'

    assert pha.get_filter(format='%.1f') == '10.0:16.0'

    # Fake filter/notice calls
    pha.mask = np.ones(6, dtype=bool)
    assert pha.get_filter(format='%.1f') == '10.0:16.0'

    pha.mask[1] = False
    pha.mask[4] = False
    assert pha.get_filter(format='%.1f') == '10.0:11.0,12.0:14.0,15.0:16.0'

    pha.mask = ~pha.mask
    assert pha.get_filter(format='%.1f') == '11.0:12.0,14.0:15.0'

    pha.mask = np.zeros(6, dtype=bool)
    assert pha.get_filter(format='%.1f') == ''

    pha.mask[0] = True
    assert pha.get_filter(format='%.1f') == '10.0:11.0'

    pha.mask[-1] = True
    assert pha.get_filter(format='%.1f') == '10.0:11.0,15.0:16.0'

    pha.mask[0] = False
    assert pha.get_filter(format='%.1f') == '15.0:16.0'


@pytest.mark.parametrize('ignore', [False, True])
@pytest.mark.parametrize('lo,hi,evals',
                         [(0.5, 2.3, (0, 10, 0)),
                          (0.7, 2.1, (1, 8, 1)),
                          (0.5, 0.7, (0, 2, 8)),
                          (1.1, 1.3, (3, 2, 5)),
                          (2.1, 2.3, (8, 2, 0)),
                          # special case filters that are within a single bin
                          (0.45, 0.55, (0, 1, 9)),
                          (0.65, 0.75, (1, 1, 8)),
                          (1.05, 1.15, (3, 1, 6)),
                          (2.25, 2.35, (9, 1, 0)),
                          # outside the limits
                          (0.1, 0.3, (10, 0, 0)),
                          (0.1, 0.4, (10, 0, 0)),
                          (0.1, 0.5, (0, 1, 9)),
                          (None, 0.5, (0, 1, 9)),
                          (2.3, 3.0, (9, 1, 0)),
                          (2.3, None, (9, 1, 0)),
                          (2.41, 2.8, (10, 0, 0)),
                          # Now queries on the edge of each bin; these would ideally
                          # only match 1 bin
                          (0.4, 0.6, (0, 1, 9)),
                          (0.6, 0.8, (1, 1, 8)),
                          (0.8, 1.0, (2, 1, 7)),
                          (1.0, 1.2, (3, 1, 6)),
                          (1.2, 1.4, (4, 1, 5)),
                          (1.4, 1.6, (5, 1, 4)),
                          (1.6, 1.8, (6, 1, 3)),
                          (1.8, 2.0, (7, 1, 2)),
                          (2.0, 2.2, (8, 1, 1)),
                          (2.2, 2.4, (9, 1, 0)),
                          # check last upper limit
                          (2.4, 2.6, (10, 0, 0)),
                          # exact limit
                          (0.4, 2.4, (0, 10, 0)),
                          # multiple bins
                          (1.0, 1.4, (3, 2, 5)),
                          (0.6, 2.0, (1, 7, 2)),
                          (1.4, 2.2, (5, 4, 1)),
                          # 0 values
                          (0, 0, (10, 0, 0)),
                          (0, 0.01, (10, 0, 0)),
                          (0, 0.5, (0, 1, 9)),
                          (0, 1, (0, 3, 7)),
                          (0, 2.4, (0, 10, 0)),
                          (0, 3.0, (0, 10, 0))
                         ])
def test_pha_check_limit(ignore, lo, hi, evals):
    """What happens when we hit values at bin edges [energy]?

    This includes some non-bin-edge tests for fun. Fortunately
    ignore and notice are inverses here so we can use the
    same pattern to generate the expected mask signal.
    """

    chans = np.arange(1, 11, dtype=int)
    counts = chans * 2
    pha = DataPHA('example', chans, counts)

    egrids = 0.2 + 0.2 * np.arange(1, 12)
    arf = DataARF('arf', egrids[:-1], egrids[1:],
                  np.ones(10))
    pha.set_arf(arf)
    pha.units = 'energy'

    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 10)

    func = pha.ignore if ignore else pha.notice
    func(lo, hi)

    if ignore:
        vout = True
        vin = False
    else:
        vout = False
        vin = True

    c1, c2, c3 = evals
    expected = [vout] * c1 + [vin] * c2 + [vout] * c3
    assert pha.mask == pytest.approx(pha.get_mask())
    assert pha.mask == pytest.approx(expected)


@pytest.mark.parametrize('ignore', [False, True])
@pytest.mark.parametrize('lo,hi,evals',
                         [(1, 10, (0, 10, 0)),
                          (2, 9, (1, 8, 1)),
                          (1, 2, (0, 2, 8)),
                          (4, 5, (3, 2, 5)),
                          (9, 10, (8, 2, 0)),
                          # outside the limits
                          (0, 0, (10, 0, 0)),
                          (0, 1, (0, 1, 9)),
                          (None, 1, (0, 1, 9)),
                          (10, 12, (9, 1, 0)),
                          (10, None, (9, 1, 0)),
                          (11, 12, (10, 0, 0)),
                          # Now queries on the edge of each bin; these would ideally
                          # only match 1 bin
                          (1, 1, (0, 1, 9)),
                          (2, 2, (1, 1, 8)),
                          (3, 3, (2, 1, 7)),
                          (4, 4, (3, 1, 6)),
                          (5, 5, (4, 1, 5)),
                          (6, 6, (5, 1, 4)),
                          (7, 7, (6, 1, 3)),
                          (8, 8, (7, 1, 2)),
                          (9, 9, (8, 1, 1)),
                          (10, 10, (9, 1, 0)),
                          # check last upper limit
                          (10, 11, (9, 1, 0)),
                          # exact limit
                          (1, 10, (0, 10, 0)),
                          # multiple bins
                          (4, 5, (3, 2, 5)),
                          (2, 8, (1, 7, 2)),
                          (6, 9, (5, 4, 1)),
                          # 0 values
                          (0, 0, (10, 0, 0)),
                          (0, 3, (0, 3, 7)),
                          (0, 10, (0, 10, 0)),
                          (0, 12, (0, 10, 0))
                         ])
def test_pha_check_limit_channel(ignore, lo, hi, evals):
    """What happens when we hit values at bin edges [channel]?

    Channel filtering behaves differently to energy/wavelength
    filtering so check.  It's not quite the same since I am
    restricting the channel ranges to integer values so some of the
    energy checks don't make sense here.

    The ARF is not necessary but keep in to better match
    test_pha_check_limit.

    """

    chans = np.arange(1, 11, dtype=int)
    counts = chans * 2
    pha = DataPHA('example', chans, counts)

    egrids = 0.2 + 0.2 * np.arange(1, 12)
    arf = DataARF('arf', egrids[:-1], egrids[1:],
                  np.ones(10))
    pha.set_arf(arf)
    pha.units = 'channel'

    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 10)

    func = pha.ignore if ignore else pha.notice
    func(lo, hi)

    if ignore:
        vout = True
        vin = False
    else:
        vout = False
        vin = True

    c1, c2, c3 = evals
    expected = [vout] * c1 + [vin] * c2 + [vout] * c3
    assert pha.mask == pytest.approx(pha.get_mask())
    assert pha.mask == pytest.approx(expected)


def test_pha_channel0_filtering():
    """If channel starts at 0 does filtering still work?

    This would ideally be done with a response to check the
    mapping but let's use a simple check for now.
    """

    chans = np.arange(1, 10, dtype=np.int16)
    counts = np.asarray([3, 4, 1, 0, 2, 5, 2, 2, 2], dtype=np.int16)
    p1 = DataPHA('p1', chans, counts)
    p0 = DataPHA('p0', chans - 1, counts)

    assert p1.channel[0] == 1
    assert p0.channel[0] == 0

    assert not p1.subtracted
    assert not p0.subtracted

    assert not p1.grouped
    assert not p0.grouped

    assert p1.mask is True
    assert p0.mask is True

    assert p1.get_dep(filter=True) == pytest.approx(counts)
    assert p0.get_dep(filter=True) == pytest.approx(counts)

    # The filter has to be different since the channel values are
    # different.
    #
    p1.notice(2, 6)
    p0.notice(1, 5)

    expected = np.zeros(9, dtype=bool)
    expected[1:6] = True
    assert p1.mask == pytest.approx(expected)
    assert p0.mask == pytest.approx(expected)

    assert p1.get_dep(filter=True) == pytest.approx(counts[1:6])
    assert p0.get_dep(filter=True) == pytest.approx(counts[1:6])


@requires_group
def test_pha_channel0_grouping():
    """If channel starts at 0 does grouping still work?"""

    chans = np.arange(1, 10, dtype=np.int16)
    counts = np.asarray([3, 4, 1, 0, 2, 5, 2, 2, 2], dtype=np.int16)
    p1 = DataPHA('p1', chans, counts)
    p0 = DataPHA('p0', chans - 1, counts)

    assert p1.channel[0] == 1
    assert p0.channel[0] == 0

    for p in [p0, p1]:
        p.group_counts(3)

    assert p1.grouped
    assert p0.grouped

    expected_grp = np.array([1, 1, 1, -1, -1, 1, 1, -1, 1], dtype=np.int16)
    assert p1.grouping == pytest.approx(expected_grp)
    assert p0.grouping == pytest.approx(expected_grp)

    expected_qual = np.zeros(9)
    expected_qual[-1] = 2
    assert p1.quality == pytest.approx(expected_qual)
    assert p0.quality == pytest.approx(expected_qual)

    # This is the grouped data.
    #
    expected = np.array([3, 4, 3, 5, 4, 2], dtype=np.int16)
    assert p1.get_dep(filter=True) == pytest.approx(expected)
    assert p0.get_dep(filter=True) == pytest.approx(expected)

    # The filter has to be different since the channel values are
    # different.
    #
    p1.notice(2, 6)
    p0.notice(1, 5)

    expected_mask = np.zeros(9, dtype=bool)
    expected_mask[1:6] = True
    assert p1.get_mask() == pytest.approx(expected_mask)
    assert p0.get_mask() == pytest.approx(expected_mask)

    assert p1.get_dep(filter=True) == pytest.approx(expected[1:4])
    assert p0.get_dep(filter=True) == pytest.approx(expected[1:4])


@requires_group
def test_pha_channel0_subtract():
    """If channel starts at 0 can we subtract the background?"""

    chans = np.arange(1, 10, dtype=np.int16)
    counts = np.asarray([3, 4, 1, 0, 2, 5, 2, 2, 2], dtype=np.int16)
    bcounts = np.asarray([0, 1, 1, 1, 0, 0, 2, 1, 1], dtype=np.int16)
    p1 = DataPHA('p1', chans, counts)
    p0 = DataPHA('p0', chans - 1, counts)

    b1 = DataPHA('b1', chans, bcounts)
    b0 = DataPHA('b0', chans - 1, bcounts)

    p1.set_background(b1)
    p0.set_background(b0)

    for p in [p0, p1]:
        p.subtract()
        p.group_counts(3)

    p1.notice(2, 6)
    p0.notice(1, 5)

    assert p1.channel[0] == 1
    assert b1.channel[0] == 1
    assert p0.channel[0] == 0
    assert b0.channel[0] == 0

    assert p1.subtracted
    assert p0.subtracted

    expected = counts - bcounts
    assert p1.get_dep() == pytest.approx(expected)
    assert p0.get_dep() == pytest.approx(expected)

    expected = np.array([3, 4, 3, 5, 4, 2], dtype=np.int16) - \
        np.array([0, 1, 2, 0, 3, 1], dtype=np.int16)
    assert p1.get_dep(filter=True) == pytest.approx(expected[1:4])
    assert p0.get_dep(filter=True) == pytest.approx(expected[1:4])


def test_set_channel_sets_independent_axis():
    """What happens if the channel attribute is set?"""

    chans = np.arange(1, 5)
    counts = chans + 1
    d = DataPHA("x", chans, counts)

    chans2 = np.arange(10, 20)
    with pytest.raises(DataErr,
                       match="independent axis can not change size: 4 to 10"):
        d.channel = chans2


def test_set_counts_sets_y_axis():
    """What happens if the counts attribute is set?

    This is meant to check if the dependent axis is also
    changed to match the new counts setting.

    """
    chans = np.arange(1, 5)
    counts = chans + 1
    d = DataPHA("x", chans, counts)

    counts2 = chans + 5
    d.counts = counts2
    assert len(d.indep) == 1
    assert d.indep[0] == pytest.approx(chans)
    assert d.y == pytest.approx(counts2)


@pytest.mark.parametrize("vals", [10, [10, 10, 10], (10, 10, 10)])
def test_pha_can_set_dep(vals):
    """Ensure we can call set_dep with a scalar or sequence for DataPHA"""

    chans = np.arange(1, 4)
    counts = chans[::-1]
    data = DataPHA("simple", chans, counts)

    assert data.get_dep() == pytest.approx(counts)

    data.set_dep(vals)

    expected = 10 * np.ones(3)
    assert data.get_dep() == pytest.approx(expected)
    assert data.counts == pytest.approx(expected)
    assert data.y == pytest.approx(expected)


# We don't really care if the arguments don't make much sense, at least
# not until we add validation code which will mean these need fixing up.
#
ELO = np.array([0.1, 0.2, 0.3])
EHI = np.array([0.2, 0.3, 0.4])
ONES = np.ones(3)
CHANS = np.arange(1, 4)
X0 = np.array([1, 2, 3, 4] * 3)
X1 = np.array([1] * 4 + [2] * 4 + [3] * 4)
IMG_ONES = np.ones(12)
ARF_ARGS = DataARF, ("arf", ELO, EHI, ONES)
RMF_ARGS = DataRMF, ("emf", 3, ELO, EHI, ONES, ONES, ONES, ONES)
PHA_ARGS = DataPHA, ("pha", CHANS, ONES)
IMG_ARGS = DataIMG, ("img", X0, X1, IMG_ONES)
IMGINT_ARGS = DataIMGInt, ("imgint", X0, X1, X0 + 1, X1 + 1, IMG_ONES)


def test_is_mask_reset_pha(caplog):
    """What happens to the mask attribute after the independent axis is changed? PHA"""

    data = PHA_ARGS[0](*PHA_ARGS[1])

    # Pick a value somewhere within the independent axis
    assert data.mask is True
    data.ignore(None, 2)
    assert data.mask == pytest.approx([False, False, True])

    # Change the independent axis, but to something of the same
    # length.
    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        data.indep = (data.channel + 20, )

    assert len(caplog.records) == 0

    # The mask has *not* been cleared
    assert data.mask == pytest.approx([False, False, True])


def test_is_mask_reset_pha_channel(caplog):
    """What happens to the mask attribute after the channel is changed? PHA

    Extends test_is_mask_reset_pha
    """

    data = PHA_ARGS[0](*PHA_ARGS[1])
    data.ignore(None, 2)

    # Change the independent axis via the channel field
    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        data.channel = data.channel + 40

    assert len(caplog.records) == 0

    # The mask has not been cleared
    assert data.mask == pytest.approx([False, False, True])


@requires_region
@pytest.mark.parametrize("data_args",
                         [IMG_ARGS, IMGINT_ARGS])
def test_is_mask_reset_img(data_args, caplog):
    """What happens to the mask attribute after the independent axis is changed? img/imgint"""

    data_class, args = data_args
    data = data_class(*args)

    # Pick a value somewhere within the independent axis
    assert data.mask is True
    data.notice2d("rect(0,0,3,3)", ignore=True)

    # The filter depends on the interpretation of the bin edges so just
    # check something has happened (the values are for the img and
    # imgint variants).
    #
    omask = data.mask.copy()
    assert omask.sum() in [3, 8]

    # Change the independent axis, but to something of the same
    # length.
    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        indep = [x + 40 for x in data.indep]
        data.indep = tuple(indep)

    assert len(caplog.records) == 1

    # The mask has been cleared
    assert data.mask

    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.astro.data"
    assert r[1] == logging.WARN
    assert r[2].startswith("Region filter has been removed from 'img")


@pytest.mark.parametrize("data_args",
                         [ARF_ARGS, RMF_ARGS, PHA_ARGS, IMG_ARGS, IMGINT_ARGS])
def test_invalid_independent_axis(data_args):
    """What happens if we use the wrong number of independent axes?

    We just duplicate the current axes.
    """

    data_class, args = data_args
    data = data_class(*args)
    indep = data.indep
    with pytest.raises(DataErr,
                       match="^data set '.*' sent wrong tuple size for the independent axis: [124] not [248]$"):
        data.indep = tuple(list(indep) * 2)


@pytest.mark.parametrize("data_args",
                         [ARF_ARGS, RMF_ARGS, IMG_ARGS, IMGINT_ARGS])
def test_invalid_independent_axis_component(data_args):
    """What happens if we use mis-matched sizes?

    We remove one entry from the second component,
    """

    data_class, args = data_args
    data = data_class(*args)
    indep = list(data.indep)
    indep[1] = indep[1][:-1]

    with pytest.raises(DataErr,
                       match=r"^size mismatch between (lo|x0) and (hi|x1): (3|11|12) vs (2|11|12)$"):
        data.indep = tuple(indep)


@pytest.mark.parametrize("data_args",
                         [ARF_ARGS, PHA_ARGS, IMG_ARGS, IMGINT_ARGS])
def test_invalid_dependent_axis(data_args):
    """What happens if the dependent axis does not match the independent axis?

    We do not include DataRMF as it's not clear what the dependent axis is
    here.
    """

    data_class, args = data_args
    data = data_class(*args)
    with pytest.raises(DataErr,
                       match=r"^size mismatch between independent axis and y: (3|12) vs (10?)$"):
        data.y = data.y[:-2]


@pytest.mark.parametrize("data_class,args",
                         [(DataARF, (ELO, EHI, np.ones(10))),
                          (DataPHA, (CHANS, np.ones(10))),
                          (DataIMG, (X0, X1, np.ones(34))),
                          (DataIMGInt, (X0, X1, X0 + 1, X1 + 1, np.ones(2))),
                          ])
def test_make_invalid_dependent_axis(data_class, args):
    """What happens if call constructor with invalid independent axis?

    We do not include DataRMF as it's not clear what the dependent axis is
    here.
    """

    with pytest.raises(DataErr,
                       match=r"^size mismatch between independent axis and y: (3|12) vs (2|10|34)$"):
        data_class("wrong", *args)


def test_pha_fails_when_bin_lo_is_invalid():
    """What happens when bin_lo has a different size to the data?

    This is to test the order of calls in the __init__ method.
    """

    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and bin_lo: 3 vs 2"):
        DataPHA("something", CHANS, ONES, bin_lo=[1, 3])


def test_pha_fails_when_bin_hi_is_invalid():
    """What happens when bin_hi has a different size to the data?

    This is to test the order of calls in the __init__ method.
    """

    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and bin_hi: 3 vs 2"):
        DataPHA("something", CHANS, ONES, bin_hi=[1, 3])


@pytest.mark.parametrize("data_args",
                         [ARF_ARGS, PHA_ARGS, IMG_ARGS, IMGINT_ARGS])
def test_set_independent_axis_to_none(data_args):
    """What happens if we clear the independent axis?"""

    data_class, args = data_args
    data = data_class(*args)

    assert all(d is not None for d in data.indep)

    indep = [None for d in data.indep]
    with pytest.raises(DataErr,
                       match="independent axis can not be cleared"):
        data.set_indep(tuple(indep))


def test_set_independent_axis_to_none_pha_channel():
    """What happens if we clear the independent axis via channel?"""

    data = PHA_ARGS[0](*PHA_ARGS[1])
    assert data.channel is not None
    assert data.indep[0] is not None

    with pytest.raises(DataErr,
                       match="independent axis can not be cleared"):
        data.channel = None


# Should we remove support for the error columns for DataARF/DataRMF?
#
@pytest.mark.parametrize("data_args",
                         [ARF_ARGS, RMF_ARGS, PHA_ARGS, IMG_ARGS, IMGINT_ARGS])
@pytest.mark.parametrize("column", ["staterror", "syserror"])
def test_set_error_axis_wrong_length(data_args, column):
    """What happens if the column is set to the wrong length?"""

    data_class, args = data_args
    data = data_class(*args)

    col = getattr(data, column)
    assert col is None

    with pytest.raises(DataErr,
                       match=rf"^size mismatch between independent axis and {column}: \d+ vs 2$"):
        setattr(data, column, [1, 2])


@pytest.mark.parametrize("related", ["y", "counts"])
def test_pha_dependent_field_can_not_be_a_scalar(related):
    """This is to contrast with test_pha_related_fields_can_not_be_a_scalar.

    This is tested elsewhere but leave here to point out that the related
    fields are not all handled the same.
    """

    data_class, args = PHA_ARGS
    data = data_class(*args)
    with pytest.raises(DataErr,
                       match="Array must be a sequence or None"):
        setattr(data, related, 2)


@pytest.mark.parametrize("vals", [[4], (2, 3, 4, 5)])
def test_pha_bin_field_must_match_initialization(vals):
    """The bin_lo/hi must match the data"""

    data_class, args = PHA_ARGS
    with pytest.raises(DataErr,
                       match=r"^size mismatch between independent axis and bin_hi: 3 vs [14]$"):
        data_class(*args, bin_lo=[1, 2, 3], bin_hi=vals)


@pytest.mark.parametrize("vals", [[4], (2, 3, 4, 5)])
def test_pha_bin_field_must_match_after(vals):
    """The bin_lo/hi must match the data"""

    data_class, args = PHA_ARGS
    data = data_class(*args)
    data.bin_hi = [4, 5, 6]
    with pytest.raises(DataErr,
                       match=r"^size mismatch between independent axis and bin_lo: 3 vs [14]$"):
        data.bin_lo = vals


@pytest.mark.parametrize("related", ["staterror", "syserror", "grouping", "quality", "bin_lo", "bin_hi"])
@pytest.mark.parametrize("vals", [True, 0, np.asarray(1)])
def test_pha_related_field_can_not_be_a_scalar(related, vals):
    """The related fields (staterror/syserror/grouping/quality/...) can not be a scalar."""

    data_class, args = PHA_ARGS
    data = data_class(*args)
    with pytest.raises(DataErr,
                       match="Array must be a sequence or None"):
        setattr(data, related, vals)


@pytest.mark.parametrize("label", ["staterror", "syserror", "grouping", "quality", "bin_lo", "bin_hi"])
@pytest.mark.parametrize("vals", [[1, 1], np.ones(10)])
def test_pha_related_field_can_not_be_a_sequence_wrong_size(label, vals):
    """Check we error out if column=label has the wrong size: sequence"""

    data_class, args = PHA_ARGS
    data = data_class(*args)
    with pytest.raises(DataErr,
                       match=f"^size mismatch between independent axis and {label}: 3 vs (2|10)$"):
        setattr(data, label, vals)


def test_pha_independent_axis_can_not_be_a_set():
    """Check we error out if x is a set

    This is an attempt to see what happens if a user sends in a
    "surprising" value. We don't want to check too many of these
    cases, but having a few checks can be useful.

    """

    data_class, args = PHA_ARGS
    data = data_class(*args)

    with pytest.raises(DataErr,
                       match="Array must be a sequence or None"):
        data.set_indep(({"abc", False, 23.4}, ))


def test_pha_independent_axis_can_not_be_a_set_sequence():
    """Check we error out if x is a set

    This is an attempt to see what happens if a user sends in a
    "surprising" value. We don't want to check too many of these
    cases, but having a few checks can be useful.

    """

    data_class, args = PHA_ARGS
    data = data_class(*args)

    # this does not raise an error
    data.set_indep(([{"abc", False, 23.4}] * 3, ))


@pytest.mark.parametrize("field", ["y", "counts", "staterror", "syserror", pytest.param("grouping", marks=pytest.mark.xfail), pytest.param("quality", marks=pytest.mark.xfail)])
def test_pha_related_field_can_not_be_a_set(field):
    """Check we error out if y/counts/... is a set

    This is an attempt to see what happens if a user sends in a
    "surprising" value. We don't want to check too many of these
    cases, but having a few checks can be useful.

    """

    data_class, args = PHA_ARGS
    data = data_class(*args)

    with pytest.raises(DataErr,
                       match="Array must be a sequence or None"):
        # XFAIL does not raise an error for grouping/quality
        setattr(data, field, {"abc", False, 23.4})


@pytest.mark.parametrize("field", [pytest.param("y", marks=pytest.mark.xfail), pytest.param("counts", marks=pytest.mark.xfail), pytest.param("staterror", marks=pytest.mark.xfail), pytest.param("syserror", marks=pytest.mark.xfail), "grouping", "quality"])
def test_pha_related_field_can_not_be_a_set_sequence(field):
    """Check we error out if y/counts/... is a sequence of sets

    This is an attempt to see what happens if a user sends in a
    "surprising" value. We don't want to check too many of these
    cases, but having a few checks can be useful.

    """

    data_class, args = PHA_ARGS
    data = data_class(*args)

    # once non-quality/grouping fields fail the error check will need to be changed
    with pytest.raises(DataErr,
                       match="Array must be a sequence of integers or None"):
        # XFAIL does not raise an error except for grouping/quality
        setattr(data, field, [{"abc", False, 23.4}] * 3)


def test_pha_mask_size_must_match_ungrouped():
    """Check if the mask can be set to the wrong length: data not grouped"""

    chans = np.arange(1, 9)
    counts = np.ones(8)
    data = DataPHA("mask", chans, counts)
    data.grouping = [1, -1, -1, 1, -1, -1, 1, 1]
    assert not data.grouped

    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and mask: 8 vs 3"):
        data.mask = [1, 0, 1]


def test_pha_mask_size_must_match_grouped():
    """Check if the mask can be set to the wrong length: data grouped"""

    chans = np.arange(1, 9)
    counts = np.ones(8)
    data = DataPHA("mask", chans, counts)
    data.grouping = [1, -1, -1, 1, -1, -1, 1, 1]
    data.grouped = True

    # The length is chosen to match the un-grouped data.
    with pytest.raises(DataErr,
                       match="size mismatch between grouped data and mask: 4 vs 8"):
        data.mask = [1, 0, 1, 1, 0, 1, 0, 0]


@pytest.mark.parametrize("funcname", ["eval_model", "eval_model_to_fit"])
def test_pha_eval_model_checks_dimensionality_pha(funcname):
    """Does eval_model check the model dimensionality?"""

    data = DataPHA('tmp', np.arange(3), np.arange(3))
    model = Gauss2D()
    func = getattr(data, funcname)
    with pytest.raises(DataErr,
                       match="Data and model dimensionality do not match: 1D and 2D"):
        func(model)


def test_datapha_create_not_ndarray():
    """If sent non nd-array fields, does __init__ convert them?

    This is a regression test.
    """

    d = DataPHA('x', [1, 2, 3], (4, 5, 6),
                staterror=(8, 7, 6), syserror=[2, 3, 4],
                grouping=(1, 1, -1), quality=(0, 0, 0),
                backscal=[2, 3, 4], areascal=(0.1, 0.2, 0.9))

    assert isinstance(d.indep, tuple)
    assert len(d.indep) == 1
    assert isinstance(d.indep[0], np.ndarray)

    assert isinstance(d.channel, np.ndarray)
    assert isinstance(d.y, np.ndarray)
    assert isinstance(d.counts, np.ndarray)
    assert isinstance(d.staterror, np.ndarray)
    assert isinstance(d.syserror, np.ndarray)

    assert isinstance(d.grouping, np.ndarray)
    assert isinstance(d.quality, np.ndarray)

    assert isinstance(d.areascal, np.ndarray)
    assert isinstance(d.backscal, np.ndarray)


@pytest.mark.parametrize("field", ["staterror", "syserror",
                                   "grouping", "quality",
                                   "areascal", "backscal"])
def test_datapha_set_not_ndarray(field):
    """What happens if the field is set to a non-ndarray after creation?

    This is a regression test.
    """

    data_class, args = PHA_ARGS
    data = data_class(*args)

    setattr(data, field, tuple([1] * len(data.y)))
    got = getattr(data, field)

    assert isinstance(got, np.ndarray)


def test_datapha_mask_set_not_ndarray():
    """What happens if the mask field is set to a non-ndarray after creation?

    This is a regression test.
    """

    data_class, args = PHA_ARGS
    data = data_class(*args)

    data.mask = tuple([1] * len(data.y))

    assert isinstance(data.mask, np.ndarray)


def test_dataimg_create_not_ndarray():
    """If sent non nd-array fields, does __init__ convert them?

    This is a regression test.
    """

    x1, x0 = np.mgrid[2:5, 3:5]
    shape = x0.shape
    x0 = list(x0.flatten())
    x1 = tuple(x1.flatten())
    listval = [1] * len(x0)
    tupleval = tuple(listval)
    d = DataIMG('x', x0, x1, listval, shape=shape,
                staterror=tupleval, syserror=listval)

    assert isinstance(d.indep, tuple)
    assert len(d.indep) == 2
    assert isinstance(d.indep[0], np.ndarray)
    assert isinstance(d.indep[1], np.ndarray)

    assert isinstance(d.y, np.ndarray)
    assert isinstance(d.staterror, np.ndarray)
    assert isinstance(d.syserror, np.ndarray)


@pytest.mark.parametrize("field", ["staterror", "syserror"])
def test_dataimg_set_not_ndarray(field):
    """What happens if the field is set to a non-ndarray after creation?

    This is a regression test.
    """

    data_class, args = IMG_ARGS
    data = data_class(*args)

    setattr(data, field, tuple([1] * len(data.y)))
    got = getattr(data, field)

    assert isinstance(got, np.ndarray)


def test_dataimg_mask_set_not_ndarray():
    """What happens if the mask field is set to a non-ndarray after creation?

    This is a regression test.
    """

    data_class, args = IMG_ARGS
    data = data_class(*args)

    data.mask = tuple([1] * len(data.y))

    assert isinstance(data.mask, np.ndarray)


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS)
def test_data_is_empty(data_class, args):
    """There is no size attribute"""

    data = data_class("empty", *args)
    assert data.size is None


def test_datapha_size():
    """Check the size field."""

    data = PHA_ARGS[0](*PHA_ARGS[1])
    assert data.size == 3


@pytest.mark.parametrize("data_args", [IMG_ARGS, IMGINT_ARGS])
def test_data2d_size(data_args):
    """Check the size field."""

    data = data_args[0](*data_args[1])
    assert data.size == 12


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS)
def test_data_can_not_set_dep_to_scalar_when_empty(data_class, args):
    """Check out how we error out.

    This is a regression test.
    """

    data = data_class("empty", *args)
    with pytest.raises(DataErr,
                       match="The size of 'empty' has not been set"):
        data.set_dep(2)


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS_2D)
@pytest.mark.parametrize("index", ["x0", "x1"])
def test_data_empty_get_x_2d(data_class, args, index):
    """What happens when there's no data?

    This is a regression test.
    """

    data = data_class("empty", *args)
    getfunc = getattr(data, f"get_{index}")
    assert getfunc() is None


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS)
def test_data_empty_apply_filter(data_class, args):
    """What does apply_filter do when the data set is empty?

    We could error out or just return the supplied argument, so
    this is a regression test
    """

    data = data_class("empty", *args)
    orig = [2, 5]
    with pytest.raises(DataErr,
                       match="The size of 'empty' has not been set"):
        data.apply_filter(orig)


def test_pha_apply_grouping_empty():
    """What does apply_grouping do when the data set is empty?

    We could error out or just return the supplied argument, so
    this is a regression test
    """

    pha = DataPHA("example", None, None)
    orig = [2, 5]
    with pytest.raises(DataErr,
                       match="The size of 'example' has not been set"):
        pha.apply_grouping(orig)


def test_pha_change_independent_element():
    """Special case PHA as the handling is more-complex than the
    general data cases.
    """

    chans = np.arange(1, 10)
    counts = np.ones(len(chans))
    pha = DataPHA("change", chans, counts)

    assert len(pha.indep) == 1
    assert pha.indep[0][1] == 2

    # change the second element of the first component
    with pytest.raises(ValueError,
                       match="assignment destination is read-only"):
        pha.indep[0][1] = -100


@pytest.mark.parametrize("field", ["y", "dep", "counts"])
def test_pha_change_dependent_element(field):
    """Special case PHA as the handling is more-complex than the
    general data cases.
    """

    chans = np.arange(1, 10)
    counts = np.arange(11, 20)
    pha = DataPHA("change", chans, counts)

    attr = getattr(pha, field)
    assert attr[1] == 12

    # change the second element of the first component
    attr[1] = 100

    expected = counts.copy()
    expected[1] = 100
    assert attr == pytest.approx(expected)

    # explicitly check pha.y just to make sure we really are changing the
    # object
    assert pha.y == pytest.approx(expected)


def test_pha_do_we_copy_the_independent_axis():
    """Special case PHA as the handling is more-complex than the
    general data cases.

    This is a regression test.
    """

    chans = np.arange(1, 10)
    counts = np.ones(len(chans))
    pha = DataPHA("change", chans, counts)

    assert pha.indep[0] == pytest.approx(chans)

    # what happens if we change the chans array?
    orig = chans.copy()
    chans[1] = -100
    assert pha.indep[0] == pytest.approx(orig)


def test_pha_do_we_copy_the_dependent_axis():
    """Special case PHA as the handling is more-complex than the
    general data cases.

    This is a regression test.
    """

    chans = np.arange(1, 10)
    counts = np.ones(len(chans))
    pha = DataPHA("change", chans, counts)

    assert pha.y == pytest.approx(counts)

    # what happens if we change the counts array?
    counts[1] = 20
    assert pha.y == pytest.approx(counts)


def test_pha_compare_mask_and_filter():
    """We can use ignore/notice or change the mask to get the same result

    A response is added so we can check with energy filtering.
    """

    x = np.arange(1, 10)
    y = x * 10
    data = DataPHA("ex", x, y)

    elo = 0.4 + x * 0.1
    ehi = elo + 0.1
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    data.set_rmf(rmf)
    data.set_analysis("energy")

    assert data.units == "energy"
    assert data.mask

    # Use notice/ignore
    #
    data.notice(0.62, 1.12)
    data.ignore(0.74, 0.82)
    assert data.mask == pytest.approx([0, 1, 0, 0, 1, 1, 1, 0, 0])
    assert data.get_dep(filter=True) == pytest.approx([20, 50, 60, 70])

    data.notice()
    assert data.mask

    # Change the mask array directly
    data.mask = [0, 1, 0, 0, 1, 1, 1, 0, 0]

    assert data.mask == pytest.approx([0, 1, 0, 0, 1, 1, 1, 0, 0])
    assert data.get_dep(filter=True) == pytest.approx([20, 50, 60, 70])


def test_pha_notice_bkg_id_none():
    """Check bkg_id=None."""

    pha = DataPHA("src", [1, 2], [10, 10])
    b1 = DataPHA("1", [1, 2], [2, 2])
    bup = DataPHA("up", [1, 2], [3, 4])

    pha.set_background(b1, id=1)
    pha.set_background(bup, id="up")

    pha.notice(lo=2, bkg_id=None)  # the default

    assert pha.mask == pytest.approx([False, True])
    assert b1.mask == pytest.approx([False, True])
    assert bup.mask == pytest.approx([False, True])


@pytest.mark.parametrize("bkg_id", [1, "up"])
def test_pha_notice_bkg_id_scalar(bkg_id):
    """Check bkg_id=scalar."""

    pha = DataPHA("src", [1, 2], [10, 10])
    b1 = DataPHA("1", [1, 2], [2, 2])
    bup = DataPHA("up", [1, 2], [3, 4])

    pha.set_background(b1, id=1)
    pha.set_background(bup, id="up")

    pha.notice(lo=2, bkg_id=bkg_id)

    assert pha.mask is True
    if bkg_id == 1:
        assert b1.mask == pytest.approx([False, True])
        assert bup.mask is True
    else:
        assert b1.mask is True
        assert bup.mask == pytest.approx([False, True])


def test_pha_notice_bkg_id_array_all():
    """Check bkg_id=array of all ids."""

    pha = DataPHA("src", [1, 2], [10, 10])
    b1 = DataPHA("1", [1, 2], [2, 2])
    bup = DataPHA("up", [1, 2], [3, 4])

    pha.set_background(b1, id=1)
    pha.set_background(bup, id="up")

    pha.notice(lo=2, bkg_id=["up", 1])

    assert pha.mask is True
    assert b1.mask == pytest.approx([False, True])
    assert bup.mask == pytest.approx([False, True])


@pytest.mark.parametrize("bkg_id", [1, "up"])
def test_pha_notice_bkg_id_array_subset(bkg_id):
    """Check bkg_id=array of one."""

    pha = DataPHA("src", [1, 2], [10, 10])
    b1 = DataPHA("1", [1, 2], [2, 2])
    bup = DataPHA("up", [1, 2], [3, 4])

    pha.set_background(b1, id=1)
    pha.set_background(bup, id="up")

    pha.notice(lo=2, bkg_id=[bkg_id])

    assert pha.mask is True
    if bkg_id == 1:
        assert b1.mask == pytest.approx([False, True])
        assert bup.mask is True
    else:
        assert b1.mask is True
        assert bup.mask == pytest.approx([False, True])


def get_img_spatial_mask():
    """This is a regression test, but it does look sensible."""

    return np.asarray([[0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                       [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                       [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                       [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                       [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
                       [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                       [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                       [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                       [0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                      dtype=bool)


@requires_region
def test_img_spatial_filter():
    """This is really meant to check some of the assumptions of the
    next test: test_img_combine_spatial_filter_and_mask():

    """

    # x0 is 10 to 24
    # x1 is -5 to 4
    #
    x1, x0 = np.mgrid[-5:5, 10:25]
    shape = (10, 15)

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.arange(1, 151)
    data = DataIMG("img", x0, x1, y, shape)

    yimg = y.reshape(shape)

    assert data.mask
    assert data.get_filter() == ""
    assert data._region is None
    assert data.get_img() == pytest.approx(yimg)

    # Add a spatial filter.
    #
    data.notice2d("circle(15, -1, 5)")

    mask = get_img_spatial_mask()

    assert data.mask == pytest.approx(mask.flatten())
    assert data.get_filter() == "Circle(15,-1,5)"
    assert data._region is not None

    got = data.get_img()
    assert np.isfinite(got) == pytest.approx(mask)
    assert got[mask] == pytest.approx(yimg[mask])


@requires_region
def test_img_combine_spatial_filter_and_mask():
    """What happens when we have both a spatial filter and change the mask attribute?"""

    # x0 is 10 to 24
    # x1 is -5 to 4
    #
    x1, x0 = np.mgrid[-5:5, 10:25]
    shape = (10, 15)

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.arange(1, 151)
    data = DataIMG("img", x0, x1, y, shape)

    yimg = y.reshape(shape)

    assert data.mask
    assert data.get_filter() == ""
    assert data._region is None
    assert data.get_img() == pytest.approx(yimg)

    # Apply a simple mask filter. We exclude those pixels
    # with x1 equal to -3 or 2.
    #
    mask = np.invert((x1 == -3) | (x1 == 2))
    data.mask = mask

    yimg = yimg.astype(SherpaFloat)
    yimg[2, :] = np.nan
    yimg[7, :] = np.nan

    assert data.mask == pytest.approx(mask)
    assert data.get_filter() == ""
    assert data._region is None

    good = np.isfinite(yimg)
    got = data.get_img()

    assert np.isfinite(got) == pytest.approx(good)
    assert got[good] == pytest.approx(yimg[good])

    # Add a spatial filter. This overlaps the excluded lines so the
    # mask becomes complicated.
    #
    data.notice2d("circle(15, -1,5)")

    # Note that although the filter has now changed the mask
    # has not.
    #
    assert data.mask == pytest.approx(mask)
    assert data.get_filter() == "Circle(15,-1,5)"
    assert data._region is not None

    assert good.sum() == 120  # just to check
    good |= get_img_spatial_mask()
    assert good.sum() == 138  # just to check

    got = data.get_img()
    yimg = y.reshape(shape)

    assert np.isfinite(got) == pytest.approx(good)
    assert got[good] == pytest.approx(yimg[good])


@pytest.fixture
def make_image_sparse():
    """Create a grid with a missing point.

    This is intended for regression tests as it's not obvious what the
    meaning of DataIMG for this use case.
    """

    # Assume this is 3 x by 2 y and we are missing the x=2, y=2 data
    # ppoint which has index (1, 1), with the indexes starting at 0.
    #
    # Assume this is 3 x by 2 y and we are missing the x=2, y=2 data
    # ppoint which has index (1, 1), with the indexes starting at 0.
    #
    x0 = np.asarray([1, 2, 3, 1, 3])
    x1 = np.asarray([1, 1, 1, 2, 2])
    y = np.asarray([1, 2, 3, 4, 6])

    # The existing documentation does not say whether this is
    # (nx, ny) or (ny, nx). The sherpa.astro.io.read_image routine
    # suggests it's the NumPy shape order, but it depends on what
    # the backend.get_image_data call does for the "y" array.
    #
    # For now, assume it's ny, nx
    shape = (2, 3)

    return DataIMG("sparse", x0, x1, y, shape)


def test_image_sparse_show(make_image_sparse):
    """What happens if given a grid but missing points?

    This is a regression test as it's not obvious what the
    meaning of DataIMG for this use case.
    """

    out = str(make_image_sparse).split("\n")

    assert out[0] == "name      = sparse"
    assert out[1] == "x0        = Int64[5]"
    assert out[2] == "x1        = Int64[5]"
    assert out[3] == "y         = Int64[5]"
    assert out[4] == "shape     = (2, 3)"
    assert out[5] == "staterror = None"
    assert out[6] == "syserror  = None"
    assert out[7] == "sky       = None"
    assert out[8] == "eqpos     = None"
    assert out[9] == "coord     = logical"
    assert len(out) == 10


def test_image_sparse_get_indep(make_image_sparse):

    out = make_image_sparse.get_indep()
    assert len(out) == 2
    assert out[0] == pytest.approx([1, 2, 3, 1, 3])
    assert out[1] == pytest.approx([1, 1, 1, 2, 2])


def test_image_sparse_get_dep(make_image_sparse):

    out = make_image_sparse.get_dep()
    assert out == pytest.approx([1, 2, 3, 4, 6])


def test_image_sparse_get_img(make_image_sparse):

    with pytest.raises(ValueError,
                       match=r"cannot reshape array of size 5 into shape \(2,3\)"):
        make_image_sparse.get_img()


@requires_region
def test_image_sparse_region_filter(make_image_sparse):
    """Pick a region that overlaps the missing pixel"""

    data = make_image_sparse
    assert data.mask

    data.notice2d("rect(2,2,4,5)")
    assert data.mask == pytest.approx([0, 0, 0, 0, 1])

    out = data.get_indep(True)
    assert len(out) == 2
    assert out[0] == pytest.approx([3])
    assert out[1] == pytest.approx([2])

    out = data.get_dep(True)
    assert out == pytest.approx([6])


@requires_region
def test_image_sparse_region_filter_restore(make_image_sparse):
    """Check we can restore the region"""

    data = make_image_sparse
    assert data.mask

    data.notice2d("rect(2,2,4,5)")
    data.notice2d()
    assert data.mask

    out = data.get_indep(True)
    assert len(out) == 2
    assert out[0] == pytest.approx([1, 2, 3, 1, 3])
    assert out[1] == pytest.approx([1, 1, 1, 2, 2])

    out = data.get_dep(True)
    assert out == pytest.approx([1, 2, 3, 4, 6])


@requires_region
def test_image_sparse_region_filter_out_all(make_image_sparse):
    """Pick a region that removes all points"""

    data = make_image_sparse
    assert data.mask

    data.notice2d("rect(2,2,4,5)")
    data.notice2d("circle(0, 0, 100)", ignore=True)

    # The mask does not get converted to False
    assert data.mask == pytest.approx([0] * 5)

    out = data.get_indep(True)
    assert len(out) == 2
    assert out[0] == pytest.approx([])
    assert out[1] == pytest.approx([])

    out = data.get_dep(True)
    assert out == pytest.approx([])


@pytest.mark.parametrize("data_args",
                         [ARF_ARGS, RMF_ARGS, PHA_ARGS, IMG_ARGS, IMGINT_ARGS])
def test_len_data(data_args):
    data_class, args = data_args
    data = data_class(*args)

    # The last argument should be the "dependent" axis (it
    # is tricky for the RMF case).
    #
    assert len(data) == args[-1].size


def test_len_datapha_empty():
    data = DataPHA('empty', None, None)
    assert len(data) == 0


def test_len_datapha_grouped():
    """Special case a grouped PHA"""
    chans = np.arange(1, 51, 1, dtype=int)
    grps = [1] + [-1] * 49
    data = DataPHA('testdata', chans, np.zeros_like(chans),
                   grouping=grps)

    # Just make sure they behave as expected, as a reminder if we
    # change any routine.
    assert data.grouped
    assert len(data) == 50
    assert len(data.y) == 50
    assert len(data.get_y()) == 1
    assert len(data.get_dep()) == 50
    assert len(data.get_dep(filter=True)) == 1


def test_len_image_sparse(make_image_sparse):
    """Special case an image"""
    data = make_image_sparse

    assert len(data) == len(data.y)
    assert len(data) == 5


def test_pha_checks_background_is_pha():
    """What happens if we send in a non-PHA background?"""

    pha = DataPHA("x", [1, 2], [1, 2])
    bkg = Data1D("y", [1, 2], [1, 2])
    with pytest.raises(ArgumentTypeErr,
                       match="^'bkg' must be a PHA data set$"):
        pha.set_background(bkg)


def test_pha_checks_background_size():
    """What happens if we send in a background with a different number of channels?"""

    pha = DataPHA("x", [1, 2, 3], [1, 2, 3])
    bkg = DataPHA("y", [1, 2], [1, 2])
    with pytest.raises(DataErr,
                       match="^The source and background channels differ$"):
        pha.set_background(bkg)


def test_pha_checks_background_size_is_set_source():
    """Ensure the source PHA has channels when adding a background"""

    pha = DataPHA("x", None, None)
    bkg = DataPHA("y", [1, 2], [1, 2])
    with pytest.raises(DataErr,
                       match="^The channel field must be set before adding a background$"):
        pha.set_background(bkg)


def test_pha_checks_background_size_is_set_bkg():
    """Ensure the bkg PHA has channels when adding a background"""

    pha = DataPHA("x", [1, 2, 3], [1, 0, 1])
    bkg = DataPHA("y", None, None)
    with pytest.raises(DataErr,
                       match="^The channel field of the background must be set$"):
        pha.set_background(bkg)


@requires_group
def test_pha_group_xxx_with_background(caplog):
    """Do we get any warning from the background?

    Assume we can test just one of the group_xxx routines.

    This is issue #1881.
    """

    chans = [1, 2, 3, 4, 5]
    src = DataPHA("src", chans, [4, 2, 3, 2, 1])
    bkg = DataPHA("bkg", chans, [1, 2, 0, 2, 2])
    src.set_background(bkg)

    assert src.grouping is None
    assert src.quality is None
    assert bkg.grouping is None
    assert bkg.quality is None

    assert not src.grouped
    assert not bkg.grouped

    # Pick something which creates different grouping for source and
    # background.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):
        src.group_counts(4)

    assert src.grouping == pytest.approx([1, 1, -1, 1, -1])
    assert src.quality == pytest.approx([0, 0, 0, 2, 2])
    assert bkg.grouping == pytest.approx([1, -1, -1, -1, 1])
    assert bkg.quality == pytest.approx([0, 0, 0, 0, 2])

    assert src.grouped
    assert bkg.grouped

    assert len(caplog.records) == 0


def test_eval_model_when_empty_datapha():
    """This is a regression test."""

    def mdl(*args):
        assert len(args) == 1
        assert args[0] is None
        return [-9]  # easy to check for

    mdl.ndim = 1
    data = DataPHA("empty", None, None)
    resp = data.eval_model(mdl)
    assert resp == pytest.approx([-9])


def test_eval_model_to_fit_when_empty_datapha():
    """This is a regression test."""

    data = DataPHA("empty", None, None)
    with pytest.raises(DataErr,
                       match="The size of 'empty' has not been set"):
        _ = data.eval_model_to_fit(Gauss1D())


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS_2D)
@pytest.mark.parametrize("funcname", ["eval_model", "eval_model_to_fit"])
def test_eval_model_when_empty_2d(data_class, args, funcname):
    """This is a regression test."""

    def mdl(*args):
        assert len(args) > 1
        assert args[0] is None
        return [-9]  # easy to check for

    mdl.ndim = 2
    data = data_class("empty", *args)
    func = getattr(data, funcname)
    resp = func(mdl)
    assert resp == pytest.approx([-9])


def test_eval_model_when_all_ignored_datapha():
    """This is a regression test."""

    def mdl(*args):
        assert len(args) == 1
        assert args[0] is not None
        return [-9]  # easy to check for

    data = DataPHA("x", [1, 2, 3], [3, 2, 7])
    data.ignore()
    # TODO: should this error out?
    resp = data.eval_model(mdl)
    assert resp == pytest.approx([-9])


def test_eval_model_to_fit_when_all_ignored_datapha():
    """This is a regression test."""

    data = DataPHA("x", [1, 2, 3], [3, 2, 7])
    data.ignore()
    with pytest.raises(DataErr,
                       match="mask excludes all data"):
        _ = data.eval_model_to_fit(Gauss1D())


def test_eval_model_when_all_ignored_dataimg():
    """This is a regression test."""

    def mdl(*args):
        assert len(args) > 1
        assert args[0] is not None
        return [-9]  # easy to check for

    data = DataIMG("x", [1, 2, 1, 2], [1, 1, 2, 2], [3, 4, 5, 8])
    data.ignore()
    resp = data.eval_model(mdl)
    assert resp == pytest.approx([-9])


def test_eval_model_to_fit_when_all_ignored_dataimg():
    """This is a regression test."""

    data = DataIMG("x", [1, 2, 1, 2], [1, 1, 2, 2], [3, 4, 5, 8])
    data.ignore()
    with pytest.raises(DataErr,
                       match="mask excludes all data"):
        _ = data.eval_model_to_fit(Gauss2D())


def test_to_guess_when_empty_datapha():
    """This is a regression test."""

    data = DataPHA("empty", None, None)
    with pytest.raises(DataErr,
                       match="The size of 'empty' has not been set"):
        _ = data.to_guess()


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS_2D)
def test_to_guess_when_empty_2d(data_class, args):
    """This is a regression test."""

    data = data_class("empty", *args)
    resp = data.to_guess()

    # Ensure there are n None values, where n is the number of
    # independent + dependent axes - ie len(args)
    #
    assert len(resp) == len(args)
    for r in resp:
        assert r is None


def test_to_guess_when_all_ignored_datapha():
    """This is a regression test."""

    data = DataPHA("x", [1, 2, 3], [2, 1, 3])
    data.ignore()
    with pytest.raises(DataErr,
                       match="mask excludes all data"):
        _ = data.to_guess()


def test_to_guess_when_all_ignored_dataimg():
    """This is a regression test."""

    data = DataIMG("x", [1, 1, 2, 2], [1, 2, 1, 2], [3, 40, 7, 12])
    data.ignore()
    with pytest.raises(DataErr,
                       match="mask excludes all data"):
        _ = data.to_guess()


def test_get_x_when_empty_datapha():
    """This is a regression test."""

    data = DataPHA("empty", None, None)
    assert data.get_x() is None


def test_get_dims_when_empty_datapha():
    """This is a regression test."""

    data = DataPHA("empty", None, None)
    with pytest.raises(DataErr,
                       match="^The size of 'empty' has not been set$"):
        _ = data.get_dims()


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS_2D)
def test_get_dims_when_empty_2d(data_class, args):
    """This is a regression test."""

    data = data_class("empty", *args)
    with pytest.raises(DataErr,
                       match="^The size of 'empty' has not been set$"):
        _ = data.get_dims()


def test_get_filter_when_empty_datapha():
    """This is a regression test."""

    data = DataPHA("empty", None, None)
    with pytest.raises(DataErr,
                       match="^The size of 'empty' has not been set$"):
        _ = data.get_filter()


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS_2D)
def test_get_filter_when_empty_2d(data_class, args):
    """This is a regression test.

    Although Data2D hard-codes get_filter, the DataIMG class does
    change the result (technically only if region support is available
    but for this case it doesn't matter).

    """

    data = data_class("empty", *args)
    assert data.get_filter() == ''


@requires_data
@requires_fits
@pytest.mark.parametrize("analysis,label",
                         [("energy", "Energy (keV)"),
                          ("wavelength", "Wavelength (Angstrom)"),
                          ("channel", "Channel")])
def test_get_xlabel_pha(make_data_path, analysis, label):
    """get_xlabel for PHA file"""

    import sherpa.astro.io
    infile = make_data_path("3c273.pi")
    data = sherpa.astro.io.read_pha(infile)
    data.set_analysis(analysis)
    assert data.get_xlabel() == label


@requires_data
@requires_fits
@pytest.mark.parametrize("analysis", ["energy", "wavelength", "channel"])
@pytest.mark.parametrize("label", ["not a label", "", "Energy (keV)"])
def test_set_xlabel_pha(make_data_path, analysis, label):
    """set_xlabel for PHA file"""

    import sherpa.astro.io
    infile = make_data_path("3c273.pi")
    data = sherpa.astro.io.read_pha(infile)
    data.set_xlabel(label)
    data.set_analysis(analysis)
    assert data.get_xlabel() == label


@requires_data
@requires_fits
@pytest.mark.parametrize("analysis,units",
                         [("energy", "keV"),
                          ("wavelength", "Angstrom"),
                          ("channel", "channel")])
def test_get_ylabel_pha(make_data_path, analysis, units):
    """get_ylabel for PHA file

    This does not check all the options for the label.
    """

    import sherpa.astro.io
    infile = make_data_path("3c273.pi")
    data = sherpa.astro.io.read_pha(infile)
    data.set_analysis(analysis)
    assert data.get_ylabel() == f"Counts/sec/{units}"


@requires_data
@requires_fits
@pytest.mark.parametrize("analysis", ["energy", "wavelength", "channel"])
@pytest.mark.parametrize("label", ["not a label", "", "Energy (keV)"])
def test_set_ylabel_pha(make_data_path, analysis, label):
    """set_ylabel for PHA file

    This does not check all the options for the label.
    """

    import sherpa.astro.io
    infile = make_data_path("3c273.pi")
    data = sherpa.astro.io.read_pha(infile)
    data.set_ylabel(label)
    data.set_analysis(analysis)
    assert data.get_ylabel() == label


@requires_data
@requires_fits
def test_get_xlabel_arf(make_data_path):
    """get_xlabel for ARF"""

    import sherpa.astro.io
    infile = make_data_path("3c273.arf")

    data = sherpa.astro.io.read_arf(infile)
    assert data.get_xlabel() == 'Energy (keV)'


@requires_data
@requires_fits
def test_get_ylabel_arf(make_data_path):
    """get_ylabel for ARF"""

    import sherpa.astro.io
    infile = make_data_path("3c273.arf")

    data = sherpa.astro.io.read_arf(infile)
    # The label depends on the backend
    ylabel = data.get_ylabel()
    assert 'cm' in ylabel
    assert '2' in ylabel


@requires_data
@requires_fits
def test_get_xlabel_rmf(make_data_path):
    """get_xlabel for RMF"""

    import sherpa.astro.io
    infile = make_data_path("3c273.rmf")

    data = sherpa.astro.io.read_rmf(infile)
    assert data.get_xlabel() == 'Energy (keV)'


@requires_data
@requires_fits
def test_get_ylabel_rmf(make_data_path):
    """get_ylabel for RMF"""

    import sherpa.astro.io
    infile = make_data_path("3c273.rmf")

    data = sherpa.astro.io.read_rmf(infile)
    assert data.get_ylabel() == 'Counts'


@requires_data
@requires_fits
@pytest.mark.parametrize("label", ["not a label", "", "Energy (keV)"])
@pytest.mark.parametrize("resp", ["arf", "rmf"])
def test_set_xlabel_pha_response(make_data_path, label, resp):
    """set_xlabel for PHA related"""

    import sherpa.astro.io
    infile = make_data_path(f"3c273.{resp}")

    getfunc = getattr(sherpa.astro.io, f"read_{resp}")
    data = getfunc(infile)
    data.set_xlabel(label)
    assert data.get_xlabel() == label


@requires_data
@requires_fits
@pytest.mark.parametrize("label", ["not a label", "", "Energy (keV)"])
@pytest.mark.parametrize("resp", ["arf", "rmf"])
def test_set_ylabel_pha_response(make_data_path, label, resp):
    """set_xlabel for PHA related"""

    import sherpa.astro.io
    infile = make_data_path(f"3c273.{resp}")

    getfunc = getattr(sherpa.astro.io, f"read_{resp}")
    data = getfunc(infile)
    data.set_ylabel(label)
    assert data.get_ylabel() == label


@requires_data
@requires_fits
@pytest.mark.parametrize("coord", ["logical", "physical", "world"])
@pytest.mark.parametrize("axis", ["x0", "x1"])
def test_get_xlabel_img(make_data_path, coord, axis):
    """get_xlabel for DataIMG"""

    import sherpa.astro.io
    infile = make_data_path("acisf08478_000N001_r0043_regevt3_srcimg.fits")

    data = sherpa.astro.io.read_image(infile)
    getfunc = getattr(data, f"get_{axis}label")
    assert getfunc() == axis


@requires_data
@requires_fits
@pytest.mark.parametrize("coord", ["logical", "physical", "world"])
@pytest.mark.parametrize("axis", ["x0", "x1"])
@pytest.mark.parametrize("label", ["not a label", "", "Energy (keV)"])
def test_set_xlabel_img(make_data_path, coord, axis, label):
    """get_xlabel for DataIMG"""

    import sherpa.astro.io
    infile = make_data_path("acisf08478_000N001_r0043_regevt3_srcimg.fits")

    data = sherpa.astro.io.read_image(infile)
    getfunc = getattr(data, f"get_{axis}label")
    setfunc = getattr(data, f"set_{axis}label")

    setfunc(label)
    assert getfunc() == label


@requires_data
@requires_fits
@pytest.mark.parametrize("coord", ["logical", "physical", "world"])
@pytest.mark.parametrize("axis", ["x0", "x1"])
def test_get_ylabel_img(make_data_path, coord, axis):
    """get_ylabel for DataIMG"""

    import sherpa.astro.io
    infile = make_data_path("acisf08478_000N001_r0043_regevt3_srcimg.fits")

    data = sherpa.astro.io.read_image(infile)
    assert data.get_ylabel() == 'y'


@requires_data
@requires_fits
@pytest.mark.parametrize("coord", ["logical", "physical", "world"])
@pytest.mark.parametrize("label", ["not a label", "", "Energy (keV)"])
def test_set_xlabel_img(make_data_path, coord, label):
    """get_xlabel for DataIMG"""

    import sherpa.astro.io
    infile = make_data_path("acisf08478_000N001_r0043_regevt3_srcimg.fits")

    data = sherpa.astro.io.read_image(infile)
    data.set_ylabel(label)
    assert data.get_ylabel() == label
