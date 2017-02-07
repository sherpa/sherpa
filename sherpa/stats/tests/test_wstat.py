#
#  Copyright (C) 2016, 2017  Smithsonian Astrophysical Observatory
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

# Explicit tests of wstat, related to issue 248
# https://github.com/sherpa/sherpa/issues/248
# The intention is the functionality of these tests will be moved
# elsewhere, but for now it is useful to have some testing of
# the various options.
#
# Note: sherpa/tests/test_fit_unit.py has been added which replicates
#       some of these tests. A review has not been done to see if
#       these tests can be removed.
#
# At present the statistic tests are not broken up into unit and
# integration tests, as most of the functionality is tested either
# in integration tests (directly or implicitly).
#
# These tests are regression tests, designed to check
# that the wstat statistic can be used with PHA data sets:
#
# a) one data set (backscal is a scalar)
#   a.1) ungrouped, no filter
#   a.2) ungrouped, energy filter
#   a.3) ungrouped, wavelength filter
#   a.4) grouped, no filter
#   a.5) grouped, energy filter
#   a.6) grouped, wavelength filter
#
#   The filters should be complex (e.g. multiple ignore regions)
#   but this part is not exercised.
#
# b) multiple data sets (backscal is a scalar)
#
# Repeat a but for multiple spectra.
#
# c) as a but with backscal being a vector
#
# d) as b but with backscal being a vector
#
# Note that for test a a previously-grouped PHA file is used, whereas
# for test b one of the files is grouped on the fly.
#
# TODO:
#   - replace the data used for test c with actual grating data
#
#   - replace unittest with py.test
#
#   - replace with a synthetic data set so the test can be run
#     without data or the IO module; however, this version does test
#     the calculation when the responses (i.e. ARF and RMF) are
#     included.
#
#   - validate some/all of the statistic values to check the
#     answer is meaningful
#
#   - should there be any tests that include quality flags?
#
#   - to test out if issue 227 is fixed
#     https://github.com/sherpa/sherpa/issues/227
#     Actually, DougBurke believes this is due to the behavior
#     of group_counts, which we already work around below,
#     search for "Note: this is related to issue 227".
#

import numpy as np

from sherpa.utils import SherpaTestCase, requires_data, requires_fits, \
    requires_group
from sherpa.astro import ui

from unittest import expectedFailure

import logging

logger = logging.getLogger('sherpa')


@requires_data
@requires_fits
class test_wstat_single_scalar(SherpaTestCase):
    """Single PHA file with a scalar backscal value"""

    def setUp(self):

        # defensive programming (one of the tests has been seen to fail
        # when the whole test suite is run without this)
        ui.clean()

        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        ui.set_stat('wstat')

        infile = self.make_path('3c273.pi')
        ui.load_pha(1, infile)

        ui.set_source(1, ui.powlaw1d.pl)

        # The powerlaw slope and normalization are
        # intended to be "a reasonable approximation"
        # to the data, just to make sure that any statistic
        # calculation doesn't blow-up too much.
        #
        ui.set_par("pl.gamma", 1.782)
        ui.set_par("pl.ampl", 1.622e-4)

    def tearDown(self):
        ui.clean()

        try:
            logger.setLevel(self._old_logger_level)
        except AttributeError:
            pass

    def _filter_data(self):
        """Filter the data.

        Use a slightly-complex filter - e.g. not just the
        end points - to exercise the system.
        """

        ui.ignore(None, 0.5)
        ui.ignore(3, 4)
        ui.ignore(7, None)

    def _check_stat(self, nbins, expected):

        # check the filter sizes (mainly so that these tests
        # get flagged as in need of a look if anything changes
        # in other parts of the code, such as filtering and binning
        #
        self.assertEqual(nbins, ui.get_data().get_dep(True).size)

        stat = ui.calc_stat()
        self.assertAlmostEqual(expected, stat, places=7)

    def test_wstat_grouped_all(self):

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        self._check_stat(46, 62.12275005769554)

    def test_wstat_grouped_filtered(self):
        self._filter_data()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        self._check_stat(35, 39.84715083397662)

    def test_wstat_ungrouped_all(self):
        ui.ungroup()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        self._check_stat(1024, 652.9968167212116)

    def test_wstat_ungrouped_filtered(self):
        ui.ungroup()
        self._filter_data()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        self._check_stat(375, 416.0601496345599)


@requires_group
@requires_data
@requires_fits
class test_wstat_two_scalar(SherpaTestCase):
    """Two PHA files with a scalar backscal value"""

    def setUp(self):

        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        ui.set_stat('wstat')

        infile1 = self.make_path('3c273.pi')
        infile2 = self.make_path('9774.pi')
        ui.load_pha(1, infile1)
        ui.load_pha(2, infile2)

        # Since 9774.pi isn't grouped, group it. Note that this
        # call groups the background to 20 counts per bin. In this
        # case we do not want that; instead we want to use the same
        # grouping scheme as the source file.
        #
        # Note: this is related to issue 227
        #
        ui.group_counts(2, 20)
        ui.set_grouping(2, bkg_id=1, val=ui.get_grouping(2))

        # There's no need to have the same model in both datasets,
        # but assume the same source model can be used, with a
        # normalization difference.
        #
        ui.set_source(1, ui.powlaw1d.pl1)
        ui.set_source(2, ui.const1d.c2 * ui.get_source(1))

        # The powerlaw slope and normalization are
        # intended to be "a reasonable approximation"
        # to the data, just to make sure that any statistic
        # calculation doesn't blow-up too much.
        #
        # Note: the model values for 3c273 are slighly different
        #       to the single-PHA-file case, so stat results are
        #       slightly different
        #
        ui.set_par("pl1.gamma", 1.7)
        ui.set_par("pl1.ampl", 1.6e-4)
        ui.set_par("c2.c0", 45)

    def tearDown(self):
        ui.clean()

        try:
            logger.setLevel(self._old_logger_level)
        except AttributeError:
            pass

    def _filter_data(self):
        """Filter the data.

        Use a slightly-complex filter - e.g. not just the
        end points - to exercise the system.
        """

        ui.ignore(None, 0.5)
        ui.ignore(3, 4)
        ui.ignore(7, None)

    def _check_stat(self, idval, nbins, expected):

        # check the filter sizes (mainly so that these tests
        # get flagged as in need of a look if anything changes
        # in other parts of the code, such as filtering and binning
        #
        self.assertEqual(nbins, ui.get_data(idval).get_dep(True).size)

        stat = ui.calc_stat(idval)
        self.assertAlmostEqual(expected, stat, places=7)

    def _check_stat2(self, expected):
        """Is calc_stat() == expected?"""

        stat = ui.calc_stat()
        self.assertAlmostEqual(expected, stat, places=7)

    def test_wstat_grouped_all(self):

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        exp1 = 66.20114213871206
        exp2 = 401.75572944361613
        self._check_stat(1, 46, exp1)
        self._check_stat(2, 148, exp2)
        self._check_stat2(exp1 + exp2)

    def test_wstat_grouped_filtered(self):
        self._filter_data()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        exp1 = 45.93447168433694
        exp2 = 127.35556915677182
        self._check_stat(1, 35, exp1)
        self._check_stat(2, 120, exp2)
        self._check_stat2(exp1 + exp2)

    def test_wstat_ungrouped_all(self):
        ui.ungroup(1)
        ui.ungroup(2)

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        exp1 = 657.2275371837611
        exp2 = 880.8442022201893
        self._check_stat(1, 1024, exp1)
        self._check_stat(2, 1024, exp2)
        self._check_stat2(exp1 + exp2)

    def test_wstat_ungrouped_filtered(self):
        ui.ungroup(1)
        ui.ungroup(2)
        self._filter_data()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        exp1 = 420.50951039235736
        exp2 = 397.4089421041855
        self._check_stat(1, 375, exp1)
        self._check_stat(2, 375, exp2)
        self._check_stat2(exp1 + exp2)


@requires_group
@requires_data
@requires_fits
class test_wstat_group_counts(SherpaTestCase):
    """PHA file that has been grouped.

    This is the error reported in issue #227, whereby the
    grouping in the source and background files is different,
    which causes an error. It is not clear yet whether this
    should be supported or not, but a regression test is added
    to check this behavior.

    The settings are the same as the test_wstat_two_scalar
    case, for the 9774.pi file. It would be easier to move this
    into test_wstat_two_scalar (since it avoids having multiple
    copies of the settings, and lets the wstat be tested on
    > 2 data sets), but for now keep as a separate set of tests
    to make sure it is obvious what is and isn't failing.
    """

    def setUp(self):

        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        ui.set_stat('wstat')

        infile = self.make_path('9774.pi')
        ui.load_pha(1, infile)

        ui.group_counts(1, 20)

        # Unlike the test_wstat_two_scalar case, the grouping
        # is not copied over.
        # ui.set_grouping(1, bkg_id=1, val=ui.get_grouping(1))

        ui.set_source(1, ui.const1d.c1 * ui.powlaw1d.pl1)

        # These should be the same as test_wstat_two_scalar
        #
        ui.set_par("pl1.gamma", 1.7)
        ui.set_par("pl1.ampl", 1.6e-4)
        ui.set_par("c1.c0", 45)

    def tearDown(self):
        ui.clean()

        try:
            logger.setLevel(self._old_logger_level)
        except AttributeError:
            pass

    def _filter_data(self):
        """Filter the data.

        Use a slightly-complex filter - e.g. not just the
        end points - to exercise the system.
        """

        ui.ignore(None, 0.5)
        ui.ignore(3, 4)
        ui.ignore(7, None)

    def _check_stat(self, idval, nbins, expected):

        # check the filter sizes (mainly so that these tests
        # get flagged as in need of a look if anything changes
        # in other parts of the code, such as filtering and binning
        #
        self.assertEqual(nbins, ui.get_data(idval).get_dep(True).size)

        stat = ui.calc_stat(idval)
        self.assertAlmostEqual(expected, stat, places=7)

    def test_wstat_grouped_all(self):

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        expval = 401.75572944361613
        self._check_stat(1, 148, expval)

    def test_wstat_grouped_filtered(self):
        self._filter_data()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        expval = 127.35556915677182
        self._check_stat(1, 120, expval)

    def test_wstat_ungrouped_all(self):
        ui.ungroup(1)

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        expval = 880.8442022201893
        self._check_stat(1, 1024, expval)

    def test_wstat_ungrouped_filtered(self):
        ui.ungroup(1)
        self._filter_data()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        expval = 397.4089421041855
        self._check_stat(1, 375, expval)


@requires_data
@requires_fits
class test_wstat_single_array(SherpaTestCase):
    """Single PHA file with an array of backscal values.

    This really should use a grating PHA dataset, but it's not
    obvious we have one (along with the necessary responses) in
    the sherpa-test-data/ repository. So for now "hack" in
    one. The statistic values were calculated by changing the
    backscal by 0.9 but leaving it as a scalar. As the tests
    currently fail, they have not been validated.
    """

    def setUp(self):

        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        ui.set_stat('wstat')

        infile = self.make_path('3c273.pi')
        ui.load_pha(1, infile)

        # Change the backscale value slightly so that the
        # results are different to other runs with this file.
        #
        nbins = ui.get_data(1).get_dep(False).size
        bscal = 0.9 * np.ones(nbins) * ui.get_backscal(1)
        ui.set_backscal(1, backscale=bscal)

        ui.set_source(1, ui.powlaw1d.pl)

        # The powerlaw slope and normalization are
        # intended to be "a reasonable approximation"
        # to the data, just to make sure that any statistic
        # calculation doesn't blow-up too much.
        #
        ui.set_par("pl.gamma", 1.7)
        ui.set_par("pl.ampl", 1.7e-4)

    def tearDown(self):
        ui.clean()

        try:
            logger.setLevel(self._old_logger_level)
        except AttributeError:
            pass

    def _filter_data(self):
        """Filter the data.

        Use a slightly-complex filter - e.g. not just the
        end points - to exercise the system.
        """

        ui.ignore(None, 0.5)
        ui.ignore(3, 4)
        ui.ignore(7, None)

    def _check_stat(self, nbins, expected):

        # check the filter sizes (mainly so that these tests
        # get flagged as in need of a look if anything changes
        # in other parts of the code, such as filtering and binning
        #
        self.assertEqual(nbins, ui.get_data().get_dep(True).size)

        stat = ui.calc_stat()
        self.assertAlmostEqual(expected, stat, places=7)

    def test_wstat_grouped_all(self):

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        self._check_stat(46, 71.21845954979574)

    def test_wstat_grouped_filtered(self):
        self._filter_data()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        self._check_stat(35, 45.6311990089982)

    def test_wstat_ungrouped_all(self):
        ui.ungroup()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        self._check_stat(1024, 663.0160968458746)

    def test_wstat_ungrouped_filtered(self):
        ui.ungroup()
        self._filter_data()

        # Used git commit 770359b5004374b969ebb63c173f293419397b4c
        # to create the oracle value, on a linux 64-bit machine.
        self._check_stat(375, 420.8390856766203)
