#
#  Copyright (C) 2007, 2015, 2017, 2018, 2020, 2021, 2022
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

import warnings

import numpy as np

import pytest

from sherpa.astro.ui.utils import Session
from sherpa.astro.data import DataARF, DataPHA, DataRMF, DataIMG, DataIMGInt
from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa.models.basic import Gauss2D
from sherpa.utils import parse_expr, SherpaFloat
from sherpa.utils.err import DataErr
from sherpa.utils.testing import requires_data, requires_fits, requires_group


EMPTY_DATA_OBJECTS = [(DataPHA, [None] * 2),
                      (DataIMG, [None] * 3),
                      (DataIMGInt, [None] * 5)]


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

    with pytest.raises(ValueError) as exc:
        create_arf(energ_lo, energ_hi, specresp, ethresh=ethresh)

    emsg = "ethresh is None or > 0"
    assert str(exc.value) == emsg


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

    emsg = "The ARF 'user-arf' has an ENERG_LO value <= 0"

    with pytest.raises(DataErr) as exc:
        create_arf(energ_lo, energ_hi, specresp)

    assert str(exc.value) == emsg


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

    emsg = "The ARF 'user-arf' has an ENERG_HI value <= " + \
           "the replacement value of 1e-05"

    with pytest.raises(DataErr) as exc:
        create_arf(energ_lo, energ_hi, ethresh=1e-5)

    assert str(exc.value) == emsg


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

    emsg = "ethresh is None or > 0"
    with pytest.raises(ValueError) as exc:
        create_delta_rmf(energ_lo, energ_hi, ethresh=ethresh)

    assert str(exc.value) == emsg


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

    emsg = "The RMF 'delta-rmf' has an ENERG_LO value <= 0"

    with pytest.raises(DataErr) as exc:
        create_delta_rmf(energ_lo, energ_hi)

    assert str(exc.value) == emsg


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

    emsg = "The ARF 'user-arf' has an ENERG_LO value <= 0"
    with pytest.raises(DataErr) as exc:
        create_arf(energ_lo, energ_hi, specresp)

    assert str(exc.value) == emsg


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

    with pytest.raises(DataErr) as exc:
        create_delta_rmf(energ_lo, energ_hi)

    emsg = "The RMF 'delta-rmf' has an ENERG_LO value <= 0"
    assert str(exc.value) == emsg


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

    with pytest.raises(DataErr) as exc:
        create_delta_rmf(energ_lo, energ_hi, ethresh=1e-5)

    emsg = "The RMF 'delta-rmf' has an ENERG_HI value <= " + \
           "the replacement value of 1e-05"
    assert str(exc.value) == emsg


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
    testdata = DataPHA('testdata', np.arange(50, dtype=float) + 1.,
                       np.zeros(50))
    testbkg = DataPHA('testbkg', np.arange(50, dtype=float) + .5,
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
    with pytest.raises(DataErr) as exc:
        d.sum_background_data()

    assert str(exc.value) == "data set 'tmp' does not have any associated backgrounds"


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
    pha.ignore(None, 1)
    pha.ignore(7, None)

    pha.set_analysis(analysis)

    dep = np.array([15, 17, 16, 16, 17, 15, 19, 15, 16, 15,
                    16, 17, 15, 18, 16, 15, 15, 16, 15, 15,
                    15, 16, 16, 15, 15, 16, 16, 15, 16, 15])
    assert pha.get_dep(filter=True) == pytest.approx(dep)

    # The group mapping for group_width of 50 is listed in
    # test_grouping_filtering_binning.
    #
    pha.group_width(50)
    dep = np.array([213, 136,  79,  47,  47,  29,  27, 18])
    assert pha.get_dep(filter=True) == pytest.approx(dep)


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

    pha.ignore(None, 1)
    pha.ignore(7, None)

    pha.set_analysis(analysis)
    pha.group_width(50)

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

    with pytest.raises(DataErr) as de:
        func(lo, hi)

    assert str(de.value).startswith('unknown ')


@requires_data
@requires_fits
@pytest.mark.parametrize("notice", [True, False])
def test_pha_validates_range(notice, make_data_path):
    """Ensure lo <= hi check is made"""

    from sherpa.astro.io import read_pha

    pha = read_pha(make_data_path('3c273.pi'))
    func = pha.notice if notice else pha.ignore
    with pytest.raises(DataErr) as de:
        func(3, 2)

    assert str(de.value) == "unknown hi argument: 'must be >= lo'"


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
    with pytest.raises(DataErr) as de:
        d.notice(lo, hi, ignore=ignore)

    err = f'strings not allowed in {emsg} bound list'
    assert str(de.value) == err


def test_pha_creation_warns_about_non_numpy_channel():
    """What happens if the channel array is not a NumPy array?

    At the moment there is no warning for the independent axis.
    """

    chans = [1, 2, 3]
    counts = np.ones(3)
    with pytest.warns(UserWarning, match=r'Converting array \[1, 2, 3\] to numpy array\.'):
        d = DataPHA('tmp', chans, counts)

    assert isinstance(d.x, np.ndarray)
    assert isinstance(d.channel, np.ndarray)

    assert d.channel == pytest.approx(chans)


def test_pha_creation_warns_about_non_numpy_counts():
    """What happens if the counts array is not a NumPy array?"""

    chans = np.asarray([1, 2, 3])
    counts = [1, 0, 1]
    with pytest.warns(UserWarning, match=r'Converting array \[1, 0, 1\] to numpy array\.'):
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
    # for the channel/grup/energy range mapping for
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
    rmf = create_delta_rmf(chans + 10, chans + 11,
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
    assert pha.get_mask() is None

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
    assert pha.get_mask() is None

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

    for p in [p0, p1]:
        p.notice(2, 6)

    expected = np.zeros(9, dtype=bool)
    expected[1:6] = True
    assert p1.mask == pytest.approx(expected)

    expected = np.zeros(9, dtype=bool)
    expected[2:7] = True
    assert p0.mask == pytest.approx(expected)

    assert p1.get_dep(filter=True) == pytest.approx(counts[1:6])
    assert p0.get_dep(filter=True) == pytest.approx(counts[2:7])


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

    expected = np.array([1, 1, 1, -1, -1, 1, 1, -1, 1], dtype=np.int16)
    assert p1.grouping == pytest.approx(expected)
    assert p0.grouping == pytest.approx(expected)

    expected = np.zeros(9)
    expected[-1] = 2
    assert p1.quality == pytest.approx(expected)
    assert p0.quality == pytest.approx(expected)

    expected = np.array([3, 4, 3, 5, 4, 2], dtype=np.int16)
    assert p1.get_dep(filter=True) == pytest.approx(expected)
    assert p0.get_dep(filter=True) == pytest.approx(expected)

    for p in [p0, p1]:
        p.notice(2, 6)

    expected = np.zeros(9, dtype=bool)
    expected[1:6] = True
    assert p1.get_mask() == pytest.approx(expected)

    expected = np.zeros(9, dtype=bool)
    expected[2:8] = True  # Why is this not 2:7 since the p1 case used 1:6?
    assert p0.get_mask() == pytest.approx(expected)

    expected = np.array([3, 4, 3, 5, 4, 2], dtype=np.int16)
    assert p1.get_dep(filter=True) == pytest.approx(expected[1:4])
    assert p0.get_dep(filter=True) == pytest.approx(expected[2:5])


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
        p.notice(2, 6)

    assert p1.channel[0] == 1
    assert b1.channel[0] == 1
    assert p0.channel[0] == 0
    assert b0.channel[0] == 0

    assert p1.subtracted
    assert p0.subtracted

    expected = counts - bcounts
    assert p1.get_dep() == pytest.approx(expected)
    assert p0.get_dep() == pytest.approx(expected)

    # See also the confusion of test_pha_channel0_grouping
    #
    expected = np.array([3, 4, 3, 5, 4, 2], dtype=np.int16) - \
        np.array([0, 1, 2, 0, 3, 1], dtype=np.int16)
    assert p1.get_dep(filter=True) == pytest.approx(expected[1:4])
    assert p0.get_dep(filter=True) == pytest.approx(expected[2:5])


def test_set_channel_sets_independent_axis():
    """What happens if the channel attribute is set?

    This is meant to check if the independent axis is also
    changed to match the new channel setting.

    """
    chans = np.arange(1, 5)
    counts = chans + 1
    d = DataPHA("x", chans, counts)

    chans2 = np.arange(10, 20)
    d.channel = chans2
    assert len(d.indep) == 1
    assert d.indep[0] == pytest.approx(chans2)
    assert d.y == pytest.approx(counts)


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


def test_is_mask_reset_pha():
    """What happens to the mask attribute after the independent axis is changed? PHA"""

    data = PHA_ARGS[0](*PHA_ARGS[1])

    # Pick a value somewhere within the independent axis
    assert data.mask is True
    data.ignore(None, 2)
    assert data.mask == pytest.approx([False, False, True])

    # Change the independent axis, but to something of the same
    # length.
    data.indep = (data.channel + 20, )

    # This is a regression test as there is an argument that the mask
    # should be cleared.
    assert data.mask == pytest.approx([False, False, True])


def test_is_mask_reset_pha_channel():
    """What happens to the mask attribute after the channel is changed? PHA

    Extends test_is_mask_reset_pha
    """

    data = PHA_ARGS[0](*PHA_ARGS[1])
    data.ignore(None, 2)

    # Change the independent axis via the channel field
    data.channel += 40

    # This is a regression test as there is an argument that the mask
    # should be cleared.
    assert data.mask == pytest.approx([False, False, True])


@pytest.mark.parametrize("data_args",
                         [IMG_ARGS, IMGINT_ARGS])
def test_is_mask_reset_img(data_args):
    """What happens to the mask attribute after the independent axis is changed? img/imgint"""

    data_class, args = data_args
    data = data_class(*args)

    # Pick a value somewhere within the independent axis
    assert data.mask is True
    data.notice2d("rect(0,0,3,3)", ignore=True)

    # The filter depends on the intepretation of the bin edges so just
    # check something has happened (the values are for the img and
    # imgint variants).
    #
    omask = data.mask.copy()
    assert omask.sum() in [3, 8]

    # Change the independent axis, but to something of the same
    # length.
    indep = [x + 40 for x in data.indep]
    data.indep = tuple(indep)

    # This is a regression test as there is an argument that the mask
    # should be cleared.
    assert data.mask == pytest.approx(omask)


@pytest.mark.parametrize("data_args",
                         [ARF_ARGS, RMF_ARGS, PHA_ARGS, IMG_ARGS, IMGINT_ARGS])
def test_invalid_independent_axis(data_args):
    """What happens if we use the wrong number of independent axes?

    We just duplicate the current axes.
    """

    data_class, args = data_args
    data = data_class(*args)
    indep = data.indep
    with pytest.raises(TypeError):
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
    # At the moment this does not error out
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
    # This currently does not fail
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

    # This does not raise an error
    data_class("wrong", *args)


@pytest.mark.parametrize("data_args",
                         [ARF_ARGS, PHA_ARGS, IMG_ARGS, IMGINT_ARGS])
def test_set_independent_axis_to_none(data_args):
    """What happens if we clear the independent axis?"""

    data_class, args = data_args
    data = data_class(*args)

    assert all(d is not None for d in data.indep)

    indep = [None for d in data.indep]
    data.set_indep(tuple(indep))
    assert all(d is None for d in data.indep)


def test_set_independent_axis_to_none_pha_channel():
    """What happens if we clear the independent axis via channel?"""

    data = PHA_ARGS[0](*PHA_ARGS[1])
    assert data.channel is not None
    assert data.indep[0] is not None

    data.channel = None
    assert len(data.indep) == 1
    assert data.indep[0] is None


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

    # does not raise an error
    setattr(data, column, [1, 2])


@pytest.mark.parametrize("related", ["y", "counts"])
def test_pha_dependent_field_can_not_be_a_scalar(related):
    """This is to contrast with test_pha_related_fields_can_not_be_a_scalar.

    This is tested elsewhere but leave here to point out that the related
    fields are not all handled the same.
    """

    data_class, args = PHA_ARGS
    data = data_class(*args)
    # does not raise an error
    setattr(data, related, 2)


@pytest.mark.parametrize("vals", [[4], (2, 3, 4, 5)])
def test_pha_bin_field_must_match_initialization(vals):
    """The bin_lo/hi must match when creating the object"""

    data_class, args = PHA_ARGS
    # this does not raise an error
    data_class(*args, bin_lo=[1, 2, 3], bin_hi=vals)


@pytest.mark.parametrize("vals", [[4], (2, 3, 4, 5)])
def test_pha_bin_field_must_match_after(vals):
    """The bin_lo/hi must match when creating the object"""

    data_class, args = PHA_ARGS
    data = data_class(*args)
    data.bin_hi = [4, 5, 6]
    # this does not raise an error
    data.bin_lo = vals


@pytest.mark.parametrize("field", ["bin_lo", "bin_hi"])
@pytest.mark.parametrize("vals", [True, 0, np.asarray(1)])
def test_pha_bin_field_can_not_be_a_scalar(field, vals):
    """The bin_lo/hi fields can not be a scalar.

    Note that these fields are different to sys/staterror.
    """

    data_class, args = PHA_ARGS
    data = data_class(*args)
    # this does not raise an error
    setattr(data, field, vals)


@pytest.mark.parametrize("field", ["bin_lo", "bin_hi"])
@pytest.mark.parametrize("vals", [[1, 1], np.ones(10)])
def test_pha_bin_field_can_not_be_a_sequence_wrong_size(field, vals):
    """The bin_lo/hi fields can not be a scalar.

    Note that these fields are different to sys/staterror.
    """

    data_class, args = PHA_ARGS
    data = data_class(*args)
    # this does not raise an error
    setattr(data, field, vals)


@pytest.mark.parametrize("related", ["staterror", "syserror", "grouping", "quality"])
@pytest.mark.parametrize("vals", [True, 0, np.asarray(1)])
def test_pha_related_field_can_not_be_a_scalar(related, vals):
    """The related fields (staterror/syserror/grouping/quality) can not be a scalar."""

    data_class, args = PHA_ARGS
    data = data_class(*args)
    # this does not raise an error
    setattr(data, related, vals)


@pytest.mark.parametrize("label", ["staterror", "syserror", "grouping", "quality"])
@pytest.mark.parametrize("vals", [[1, 1], np.ones(10)])
def test_pha_related_field_can_not_be_a_sequence_wrong_size(label, vals):
    """Check we error out if column=label has the wrong size: sequence"""

    data_class, args = PHA_ARGS
    data = data_class(*args)
    # this does not raise an error
    setattr(data, label, vals)


@pytest.mark.parametrize("funcname", ["eval_model", "eval_model_to_fit"])
def test_pha_eval_model_checks_dimensionality_pha(funcname):
    """Does eval_model check the model dimensionality?"""

    data = DataPHA('tmp', np.arange(3), np.arange(3))
    model = Gauss2D()
    func = getattr(data, funcname)
    with pytest.raises(TypeError):
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
    assert isinstance(d.staterror, tuple)
    assert isinstance(d.syserror, list)

    assert isinstance(d.grouping, tuple)
    assert isinstance(d.quality, tuple)

    assert isinstance(d.areascal, tuple)
    assert isinstance(d.backscal, list)


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

    assert isinstance(got, tuple)


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
    assert isinstance(d.staterror, tuple)
    assert isinstance(d.syserror, list)


@pytest.mark.parametrize("field", ["staterror", "syserror"])
def test_dataimg_set_not_ndarray(field):
    """What happens if the field is set to a non-ndarray after creation?

    This is a regression test.
    """

    data_class, args = IMG_ARGS
    data = data_class(*args)

    setattr(data, field, tuple([1] * len(data.y)))
    got = getattr(data, field)

    assert isinstance(got, tuple)


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
    with pytest.raises(AttributeError):
        data.size


def test_datapha_size():
    """Check the size field."""

    data = PHA_ARGS[0](*PHA_ARGS[1])
    with pytest.raises(AttributeError):
        data.size


@pytest.mark.parametrize("data_args", [IMG_ARGS, IMGINT_ARGS])
def test_data2d_size(data_args):
    """Check the size field."""

    data = data_args[0](*data_args[1])
    with pytest.raises(AttributeError):
        data.size


@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS)
def test_data_can_not_set_dep_to_scalar_when_empty(data_class, args):
    """Check out how we error out.

    This is a regression test.
    """

    data = data_class("empty", *args)
    with pytest.raises(TypeError) as err:
        data.set_dep(2)

    assert str(err.value) == "object of type 'NoneType' has no len()"


#@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS[1:])
@pytest.mark.parametrize("data_class,args",
                         [(DataIMG, [None] * 3),
                          pytest.param(DataIMGInt, [None] * 5, marks=pytest.mark.xfail)])
@pytest.mark.parametrize("index", ["x0", "x1"])
def test_data_empty_get_x_2d(data_class, args, index):
    """What happens when there's no data?

    This is a regression test.
    """

    data = data_class("empty", *args)
    getfunc = getattr(data, f"get_{index}")
    # XFAIL: for Data2DInt there's a TypeError about adding None to None
    assert getfunc() is None


#@pytest.mark.parametrize("data_class,args", EMPTY_DATA_OBJECTS)
@pytest.mark.parametrize("data_class,args",
                         [pytest.param(DataPHA, [None] * 2, marks=pytest.mark.xfail),
                          (DataIMG, [None] * 3),
                          (DataIMGInt, [None] * 5)])
def test_data_empty_apply_filter(data_class, args):
    """What does apply_filter do when the data set is empty?

    We could error out or just return the supplied argument, so
    this is a regression test
    """

    data = data_class("empty", *args)
    orig = [2, 5]
    ans = data.apply_filter(orig)
    # XFAIL for DataPHA get a TypeError from calling len(None)
    assert ans == pytest.approx(orig)


def test_pha_apply_grouping_empty():
    """What does apply_grouping do when the data set is empty?

    We could error out or just return the supplied argument, so
    this is a regression test
    """

    pha = DataPHA("example", None, None)
    orig = [2, 5]
    ans = pha.apply_grouping(orig)
    assert ans == pytest.approx(orig)


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
    pha.indep[0][1] = -100

    expected = chans.copy()
    expected[1] = -100
    assert len(pha.indep) == 1
    assert pha.indep[0] == pytest.approx(expected)


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
    chans[1] = -100
    assert pha.indep[0] == pytest.approx(chans)


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
