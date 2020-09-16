#
#  Copyright (C) 2007, 2015, 2017, 2018, 2020
#        Smithsonian Astrophysical Observatory
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
from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.utils.err import DataErr
from sherpa.utils.testing import SherpaTestCase, requires_data, requires_fits

import logging
logger = logging.getLogger('sherpa')


def _monotonic_warning(response_type, filename):
    return UserWarning("The {} '{}' has a non-monotonic ENERG_LO array".format(response_type, filename))


def _bin_warning(response_type, filename):
    return UserWarning("The {} '{}' has at least one bin with ENERG_HI < ENERG_LO".format(response_type, filename))


def _assert_userwarning(expected_warnings, observed_warnings):

    expected_warnings_set = set([warning.args for warning in expected_warnings])
    observed_warnings_set = set([warning.message.args for warning in observed_warnings])

    assert  observed_warnings_set == expected_warnings_set


class test_filter_energy_grid(SherpaTestCase):

    _notice = np.ones(46, dtype=bool)
    _notice[44:46] = False

    _ignore = np.zeros(46, dtype=bool)
    _ignore[14:33] = True

    _emin = np.array([
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
        6.57000017e+00, 9.86960030e+00], np.float)

    _emax = np.array([
        0.2482    , 0.3066    , 0.46720001, 0.56940001, 0.64240003,
        0.7008    , 0.7446    , 0.78839999, 0.81760001, 0.86140001,
        0.89060003, 0.949     , 0.9928    , 1.03659999, 1.09500003,
        1.13880002, 1.19719994, 1.28480005, 1.4016    , 1.47459996,
        1.60599995, 1.69360006, 1.81040001, 1.898     , 1.9418    ,
        2.02940011, 2.08780003, 2.19000006, 2.27760005, 2.39439988,
        2.58419991, 2.71560001, 2.86159992, 3.08060002, 3.38720012,
        3.5624001 , 3.796     , 4.02960014, 4.24860001, 4.71579981,
        5.0223999 , 5.37279987, 5.89839983, 6.57000017, 9.8696003 ,
        14.95040035], np.float)

    def setUp(self):
        self.old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        self.pha = DataPHA('', np.arange(46, dtype=float) + 1.,
                           np.zeros(46),
                           bin_lo=self._emin,
                           bin_hi=self._emax)
        self.pha.units = "energy"

    def tearDown(self):
        logger.setLevel(self.old_level)

    def test_notice(self):
        self.pha.notice()
        self.pha.notice(0.0, 6.0)
        assert (self._notice == np.asarray(self.pha.mask)).all()

    def test_ignore(self):
        self.pha.notice()
        self.pha.ignore(0.0, 1.0)
        self.pha.ignore(3.0, 15.0)
        assert (self._ignore == np.asarray(self.pha.mask)).all()


class test_filter_energy_grid_reversed(SherpaTestCase):

    _notice = np.zeros(204, dtype=bool)
    _notice[0:42] = True

    _ignore = np.ones(204, dtype=bool)
    _ignore[66:70] = False
    _ignore[0:17] = False

    _emin = np.array([
        2.39196181, 2.35973215, 2.34076023, 2.30973101, 2.2884388,
        2.25861454, 2.22371697, 2.20662117, 2.18140674, 2.14317489,
        2.12185216, 2.09055495, 2.06256914, 2.04509854, 2.02788448,
        2.00133967, 1.97772908, 1.96379483, 1.93868744, 1.91855776,
        1.89444292, 1.87936974, 1.85819471, 1.84568763, 1.82923627,
        1.78920078, 1.77360916, 1.76206875, 1.74499893, 1.73006463,
        1.70084822, 1.6883322 , 1.67772949, 1.65171933, 1.63476169,
        1.59687376, 1.5745424 , 1.55736887, 1.54051399, 1.52546024,
        1.50043869, 1.48890531, 1.47329199, 1.46072423, 1.44289041,
        1.43344045, 1.41616774, 1.40441585, 1.3979584 , 1.38773119,
        1.37138033, 1.35170007, 1.33725214, 1.33249414, 1.31839108,
        1.30797839, 1.29657102, 1.28310275, 1.26550889, 1.25471842,
        1.24513853, 1.23672664, 1.22944438, 1.21509433, 1.21003771,
        1.20401597, 1.19705439, 1.18722582, 0.90194935, 0.89519638,
        0.88912934, 0.88492262, 0.87837797, 0.87366825, 0.8689999,
        0.86437255, 0.85693878, 0.84793305, 0.84404182, 0.83580172,
        0.82876647, 0.82395256, 0.81865752, 0.81185687, 0.80004948,
        0.79450154, 0.78852075, 0.77920061, 0.77340651, 0.76626247,
        0.76202762, 0.75783074, 0.75413191, 0.74727529, 0.74321008,
        0.73474538, 0.73166627, 0.72687   , 0.71785438, 0.71488959,
        0.71068853, 0.70199603, 0.69832331, 0.69387686, 0.68788701,
        0.68354762, 0.67847627, 0.67117327, 0.66512167, 0.66175646,
        0.65620857, 0.6518243 , 0.64605182, 0.64142239, 0.63754696,
        0.63128632, 0.62478495, 0.62006336, 0.61440694, 0.60915887,
        0.60591549, 0.60078359, 0.5938406 , 0.59103745, 0.58488411,
        0.58124125, 0.57883304, 0.57406437, 0.57023615, 0.56442606,
        0.56041539, 0.55701393, 0.55392498, 0.55030966, 0.54346251,
        0.53728294, 0.53515989, 0.5291304 , 0.52448714, 0.51990861,
        0.51589233, 0.50996011, 0.50509953, 0.49889025, 0.49512967,
        0.49003205, 0.48888513, 0.48524383, 0.48164544, 0.47720695,
        0.47283325, 0.46916556, 0.46660379, 0.46280268, 0.45925769,
        0.45514211, 0.45290345, 0.44987884, 0.44589564, 0.44333643,
        0.44099477, 0.43790293, 0.43446559, 0.43088335, 0.42605683,
        0.42131537, 0.41826019, 0.41506338, 0.41155648, 0.40895697,
        0.40502119, 0.40400422, 0.40164718, 0.39864835, 0.39584854,
        0.39389083, 0.39130434, 0.38890362, 0.38526753, 0.38292497,
        0.38075879, 0.37891743, 0.37648395, 0.37557775, 0.37347662,
        0.37154216, 0.36742872, 0.3641032 , 0.36167556, 0.35983625,
        0.35634032, 0.35248783, 0.35085678, 0.34843227, 0.34669766,
        0.34418666, 0.33912122, 0.33720407, 0.33505177, 0.33279634,
        0.33081138, 0.32847831, 0.32592943, 0.3111549 ], np.float)

    _emax = np.array([
        3.06803656, 2.39196181, 2.35973215, 2.34076023, 2.30973101,
        2.2884388 , 2.25861454, 2.22371697, 2.20662117, 2.18140674,
        2.14317489, 2.12185216, 2.09055495, 2.06256914, 2.04509854,
        2.02788448, 2.00133967, 1.97772908, 1.96379483, 1.93868744,
        1.91855776, 1.89444292, 1.87936974, 1.85819471, 1.84568763,
        1.82923627, 1.78920078, 1.77360916, 1.76206875, 1.74499893,
        1.73006463, 1.70084822, 1.6883322 , 1.67772949, 1.65171933,
        1.63476169, 1.59687376, 1.5745424 , 1.55736887, 1.54051399,
        1.52546024, 1.50043869, 1.48890531, 1.47329199, 1.46072423,
        1.44289041, 1.43344045, 1.41616774, 1.40441585, 1.3979584,
        1.38773119, 1.37138033, 1.35170007, 1.33725214, 1.33249414,
        1.31839108, 1.30797839, 1.29657102, 1.28310275, 1.26550889,
        1.25471842, 1.24513853, 1.23672664, 1.22944438, 1.21509433,
        1.21003771, 1.20401597, 1.19705439, 1.18722582, 0.90194935,
        0.89519638, 0.88912934, 0.88492262, 0.87837797, 0.87366825,
        0.8689999 , 0.86437255, 0.85693878, 0.84793305, 0.84404182,
        0.83580172, 0.82876647, 0.82395256, 0.81865752, 0.81185687,
        0.80004948, 0.79450154, 0.78852075, 0.77920061, 0.77340651,
        0.76626247, 0.76202762, 0.75783074, 0.75413191, 0.74727529,
        0.74321008, 0.73474538, 0.73166627, 0.72687   , 0.71785438,
        0.71488959, 0.71068853, 0.70199603, 0.69832331, 0.69387686,
        0.68788701, 0.68354762, 0.67847627, 0.67117327, 0.66512167,
        0.66175646, 0.65620857, 0.6518243 , 0.64605182, 0.64142239,
        0.63754696, 0.63128632, 0.62478495, 0.62006336, 0.61440694,
        0.60915887, 0.60591549, 0.60078359, 0.5938406 , 0.59103745,
        0.58488411, 0.58124125, 0.57883304, 0.57406437, 0.57023615,
        0.56442606, 0.56041539, 0.55701393, 0.55392498, 0.55030966,
        0.54346251, 0.53728294, 0.53515989, 0.5291304 , 0.52448714,
        0.51990861, 0.51589233, 0.50996011, 0.50509953, 0.49889025,
        0.49512967, 0.49003205, 0.48888513, 0.48524383, 0.48164544,
        0.47720695, 0.47283325, 0.46916556, 0.46660379, 0.46280268,
        0.45925769, 0.45514211, 0.45290345, 0.44987884, 0.44589564,
        0.44333643, 0.44099477, 0.43790293, 0.43446559, 0.43088335,
        0.42605683, 0.42131537, 0.41826019, 0.41506338, 0.41155648,
        0.40895697, 0.40502119, 0.40400422, 0.40164718, 0.39864835,
        0.39584854, 0.39389083, 0.39130434, 0.38890362, 0.38526753,
        0.38292497, 0.38075879, 0.37891743, 0.37648395, 0.37557775,
        0.37347662, 0.37154216, 0.36742872, 0.3641032 , 0.36167556,
        0.35983625, 0.35634032, 0.35248783, 0.35085678, 0.34843227,
        0.34669766, 0.34418666, 0.33912122, 0.33720407, 0.33505177,
        0.33279634, 0.33081138, 0.32847831, 0.32592943], np.float)

    def setUp(self):
        self.pha = DataPHA('', np.arange(204, dtype=float) + 1.,
                           np.zeros(204),
                           bin_lo=self._emin,
                           bin_hi=self._emax)
        self.pha.units = "energy"

    def tearDown(self):
        pass

    def test_notice(self):
        self.pha.notice()
        self.pha.notice(4., 8.3)
        assert (self._notice == np.asarray(self.pha.mask)).all()

    def test_ignore(self):
        self.pha.notice()
        self.pha.ignore(10.3, 13.8)
        self.pha.ignore(4.6, 6.2)
        assert (self._ignore == np.asarray(self.pha.mask)).all()


class test_filter_wave_grid(SherpaTestCase):

    _notice = np.ones(16384, dtype=bool)
    _notice[8465:16384] = False

    _ignore = np.zeros(16384, dtype=bool)
    _ignore[14064:15984] = True

    _emin = np.arange(205.7875, 0.9875, -0.0125)

    _emax = _emin + 0.0125

    def setUp(self):
        self.old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        self.pha = DataPHA('', np.arange(16384, dtype=float) + 1,
                           np.zeros(16384),
                           bin_lo=self._emin,
                           bin_hi=self._emax)

    def tearDown(self):
        logger.setLevel(self.old_level)

    def test_notice(self):
        self.pha.units = 'wavelength'
        self.pha.notice()
        self.pha.notice(100.0, 225.0)
        assert (self._notice == np.asarray(self.pha.mask)).all()

    def test_ignore(self):
        self.pha.units = 'wavelength'
        self.pha.notice()
        self.pha.ignore(30.01, 225.0)
        self.pha.ignore(0.1, 6.0)
        assert (self._ignore == np.asarray(self.pha.mask)).all()


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
# The create_arf/create_delta_rmf routines are similar to those in
# test_instrument.py
#
def create_arf(elo, ehi, specresp=None, exposure=None, ethresh=None):
    """Create an ARF.

    Parameters
    ----------
    elo, ehi : array
        The energy bins (low and high, in keV) for the ARF. It is
        assumed that ehi_i > elo_i, elo_j > 0, the energy bins are
        either ascending - so elo_i+1 > elo_i - or descending
        (elo_i+1 < elo_i), and that there are no overlaps.
    specresp : None or array, optional
        The spectral response (in cm^2) for the ARF. It is assumed
        to be >= 0. If not given a flat response of 1.0 is used.
    exposure : number or None, optional
        If not None, the exposure of the ARF in seconds.
    ethresh : number or None, optional
        Passed through to the DataARF call. It controls whether
        zero-energy bins are replaced.

    Returns
    -------
    arf : DataARF instance

    """

    assert elo.size == ehi.size
    assert (exposure is None) or (exposure > 0.0)

    if specresp is None:
        specresp = np.ones(elo.size, dtype=np.float32)

    return DataARF('test-arf', energ_lo=elo, energ_hi=ehi,
                   specresp=specresp, exposure=exposure, ethresh=ethresh)


def create_delta_rmf(rmflo, rmfhi, startchan=1,
                     e_min=None, e_max=None, ethresh=None):
    """Create a RMF for a delta-function response.

    This is a "perfect" (delta-function) response.

    Parameters
    ----------
    rmflo, rmfhi : array
        The energy bins (low and high, in keV) for the RMF.
        It is assumed that emfhi_i > rmflo_i, rmflo_j > 0, that the energy
        bins are either ascending, so rmflo_i+1 > rmflo_i or descending
        (rmflo_i+1 < rmflo_i), and that there are no overlaps.
        These correspond to the Elow and Ehigh columns (represented
        by the ENERG_LO and ENERG_HI columns of the MATRIX block) of
        the OGIP standard.
    startchan : int, optional
        The starting channel number: expected to be 0 or 1 but this is
        not enforced.
    e_min, e_max : None or array, optional
        The E_MIN and E_MAX columns of the EBOUNDS block of the
        RMF.
    ethresh : number or None, optional
        Passed through to the DataARF call. It controls whether
        zero-energy bins are replaced.

    Returns
    -------
    rmf : DataRMF instance

    Notes
    -----
    I do not think I have the startchan=0 case correct (does the
    f_chan array have to change?).
    """

    assert rmflo.size == rmfhi.size
    assert startchan >= 0

    # Set up the delta-function response.
    # TODO: should f_chan start at startchan?
    #
    nchans = rmflo.size
    matrix = np.ones(nchans, dtype=np.float32)
    dummy = np.ones(nchans, dtype=np.int16)
    f_chan = np.arange(1, nchans + 1, dtype=np.int16)

    return DataRMF('delta-rmf', detchans=nchans,
                   energ_lo=rmflo, energ_hi=rmfhi,
                   n_grp=dummy, n_chan=dummy,
                   f_chan=f_chan, matrix=matrix,
                   offset=startchan,
                   e_min=e_min, e_max=e_max,
                   ethresh=ethresh)


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
        expected_warnings = [_bin_warning('ARF', 'test-arf'), _monotonic_warning('ARF', 'test-arf')]
    else:
        expected_warnings = [_bin_warning('ARF', 'test-arf')]

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

    expected_warnings = [_monotonic_warning('ARF', 'test-arf')]

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

    emsg = "The ARF 'test-arf' has an ENERG_LO value <= 0"

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

    expected_warnings = [UserWarning("The minimum ENERG_LO in the ARF 'test-arf' was 0 " + \
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

    emsg = "The ARF 'test-arf' has an ENERG_HI value <= " + \
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

    expected_warnings = [UserWarning("The minimum ENERG_LO in the RMF 'delta-rmf' was 0 " + \
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

    emsg = "The ARF 'test-arf' has an ENERG_LO value <= 0"
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

    expected_warnings = [UserWarning("The ARF 'test-arf' has an ENERG_LO value < 0")]

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
def test_ungroup():
    '''Make sure that ungrouped data can be ungrouped.

    This test just groups and ungroups a few times.
    '''
    session = Session()
    testdata = DataPHA('testdata', np.arange(50, dtype=float) + 1.,
                                   np.zeros(50),
                                   bin_lo=1,
                                   bin_hi=10)
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
                       np.zeros(50),
                       bin_lo=1, bin_hi=10)
    testbkg = DataPHA('testbkg', np.arange(50, dtype=float) + .5,
                      np.zeros(50),
                      bin_lo=1, bin_hi=10)
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
@pytest.mark.parametrize("infile,expected",
                         [("9774.pi", '0.0080:14.9431'),
                          ("3c273.pi", '0.1248:12.4100')])
def test_pha_get_filter_none(infile, expected, make_data_path):
    """Check get_filter with no filter

    It would be nice to do this with faked data, but easier this
    way.
    """

    import sherpa.astro.io

    pha = sherpa.astro.io.read_pha(make_data_path(infile))

    assert pha.get_filter(format='%.4f') == expected


@requires_data
@requires_fits
@pytest.mark.parametrize("infile,expected",
                         [("9774.pi", '0.5037:1.9783,3.0149:7.0007'),
                          ("3c273.pi", '0.5183:1.9199,3.2339:8.2198')])
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
                         [("9774.pi", '0.5037:0.6059'),
                          ("3c273.pi", '0.5183:0.6059')])
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
                         [("9774.pi", '0.5037:7.0007'),
                          ("3c273.pi", '0.4745:9.8623')])
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

    assert pha.mask.dtype == np.bool
    assert pha.mask.size == size
    assert pha.mask.sum() == nset
