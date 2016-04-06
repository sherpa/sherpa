#
#  Copyright (C) 2016  Smithsonian Astrophysical Observatory
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

import unittest
from sherpa.astro import ui
from numpy.testing import assert_almost_equal
import logging
import sys

logger = logging.getLogger("sherpa")


def run(verbosity, require_failure):
    test_suite = SmokeTestSuite(require_failure=require_failure)
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(test_suite)
    if result is None or result.failures or result.errors or result.unexpectedSuccesses:
        sys.exit("Test failures were detected")


class SmokeTest(unittest.TestCase):
    def setUp(self):
        self._old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

    def tearDown(self):
        if hasattr(self, "old_level"):
            logger.setLevel(self._old_level)

    def test_fit(self):
        ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
        ui.set_source("polynom1d.p")
        ui.thaw("p.c1")
        ui.set_method("levmar")
        ui.fit()
        model = ui.get_model_component("p")
        expected = [0, 1]
        observed = [model.c0.val, model.c1.val]
        assert_almost_equal(observed, expected)

    def test_failure(self):
        self.fail("Requested to fail")


class SmokeTestSuite(unittest.TestSuite):
    test_classes = (SmokeTest,)
    failure_methods = ('test_failure',)

    def __init__(self, *args, **kwargs):
        self.require_failure = kwargs.pop("require_failure", False)
        unittest.TestSuite.__init__(self, *args, **kwargs)

        for test_class in self.test_classes:
            self.addTest(unittest.makeSuite(test_class))

        self._set_skips()

    def _set_skips(self):
        if not self.require_failure:
            for case in self:
                for test in case:
                    if test.id().split(".")[-1] in self.failure_methods:
                        setattr(test, 'setUp', lambda: test.skipTest('Smoke Test not required to fail, skipping'))