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
from tempfile import NamedTemporaryFile

from sherpa.utils import requires_fits, requires_xspec
from sherpa.astro import ui
from numpy.testing import assert_almost_equal
import logging
import sys
import os

logger = logging.getLogger("sherpa")


def run(verbosity=0, require_failure=False, fits=None, xspec=False):
    test_suite = SmokeTestSuite(require_failure=require_failure)
    runner = unittest.TextTestRunner(verbosity=verbosity)
    result = runner.run(test_suite)

    missing_requirements = []

    if fits:
        try:
            __import__(fits, globals(), locals())
        except ImportError:
            missing_requirements.append(_import_error(fits, "fits"))

    if xspec:
        try:
            import sherpa.astro.xspec
        except ImportError:
            missing_requirements.append(_import_error("xspec", "xspec"))

    for missing in missing_requirements:
        result.addFailure(None, missing)

    if result is None or result.failures or result.errors or result.unexpectedSuccesses:
        sys.exit("Test failures were detected")


class SmokeTest(unittest.TestCase):
    def setUp(self):
        self._old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        from sherpa.astro import datastack
        folder = os.path.dirname(datastack.__file__)
        self.fits = os.path.join(folder, "tests", "data", "acisf07867_000N001_r0002_pha3.fits")

        self.x = [1, 2, 3]
        self.y = [1, 2, 3]

    def tearDown(self):
        if hasattr(self, "old_level"):
            logger.setLevel(self._old_level)

    def test_fit(self):
        ui.load_arrays(1, self.x, self.y)
        ui.set_source("polynom1d.p")
        ui.thaw("p.c1")
        ui.set_method("levmar")
        ui.fit()
        model = ui.get_model_component("p")
        expected = [0, 1]
        observed = [model.c0.val, model.c1.val]
        assert_almost_equal(observed, expected)

    @requires_fits
    def test_fits_io(self):
        ui.load_pha(self.fits)
        with NamedTemporaryFile() as f:
            ui.save_pha(f.name, ascii=False, clobber=True)

    @requires_xspec
    def test_xspec(self):
        ui.load_arrays(1, self.x, self.y)
        ui.set_source("xspowerlaw.p")
        ui.set_method("moncar")
        ui.set_stat("chi2xspecvar")
        ui.fit()
        model = ui.get_model_component("p")
        expected = [-1.3686404, 0.5687635]
        observed = [model.PhoIndex.val, model.norm.val]
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


def _import_error(module, name):
    sys.exit("ERROR: Requested {} as {} but module not found".format(module, name))