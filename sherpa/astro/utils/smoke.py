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

from sherpa.astro import ui
from sherpa.utils import has_package_from_list
from numpy.testing import assert_almost_equal
import logging
import sys
import os

logger = logging.getLogger("sherpa")


def run(verbosity=0, require_failure=False, fits=None, xspec=False, ds9=False):
    """
    Run the smoke tests.

    A smoke test is a small set of tests designed to verify that the installation does not have
    any obvious issues.

    Parameters
    ----------
    verbosity : non-negative int
        Set verbosity level: the higher the level, the more verbose the test. This
        parameter is passed to the unittest TextTestRunner class.
        require_failure boolean. For debugging purposes, the test may be required to always fail.

    fits : basestring
        The string representing the FITS backend module to require. If the module cannot be imported,
        the smoke test will fail. However, the other tests will keep running, although they may be skipped if they
        require the FITS backend.

    xspec : boolean
        Whether the xspec module is required. If the xspec module cannot be imported, the smoke test will
        fail. However, the other tests will keep running, although they may be skipped if they require xspec.

    Returns
    -------
    The function will exit with a non-zero exit status if any errors are detected.
    """

    missing_requirements = []

    class DependencyTest(unittest.TestCase):

        @unittest.skipIf(not fits, reason="FITS dependency not requested")
        def test_fits_dependency(self):
            try:
                __import__(fits, globals(), locals())
            except ImportError:
                missing_requirements.append(_import_error(fits, "fits"))
                raise

        @unittest.skipIf(not xspec, reason="XSPEC dependency not requested")
        def test_xspec_dependency(self):
            try:
                import sherpa.astro.xspec
            except ImportError:
                missing_requirements.append(_import_error("xspec", "xspec"))
                raise

        @unittest.skipIf(not ds9, reason="DS9 dependency not requested")
        def test_ds9_dependency(self):
            try:
                import sherpa.image.ds9_backend
            except:
                missing_requirements.append("ERROR: DS9 not usable or not in path")
                raise

    test_suite = SmokeTestSuite(require_failure=require_failure)
    test_suite.addTest(unittest.makeSuite(DependencyTest))
    runner = unittest.TextTestRunner(verbosity=int(verbosity))
    result = runner.run(test_suite)

    # This makes sure the exit status of the process is not OK, and prints out the
    # missing dependencies, if any.
    if result is None or result.failures or result.errors or result.unexpectedSuccesses:
        sys.exit("Test failures were detected:\n{}".format("\n".join(missing_requirements)))


class SmokeTest(unittest.TestCase):
    """
    Smoke test test case.

    Historically, CIAO software has been accompanied by "smoke tests", tests that are shipped with the
    main code distribution and that can be run to assess whether the software was properly installed.

    The Smoke Test is not a full tests suite, but simply tries to detect obvious failures that signal an
    improper setup.

    This class implements the Smoke Test as a standalone unittest test case. This means that there are no
    external dependencies required for the Smoke Test to run, and the only data files it may use are shipped
    with the main distribution as well.

    It is also run as part fo the full Sherpa tests suite.

    Some methods require external dependencies. If such dependencies are not found the tests are skipped.

    There is more logic (see the run method in this module) that checks the requirements are actually installed
    and fails the Smoke Test if they are not.
    """
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
        """
        Perform a very simple fit with built-in models, and check that the results make sense.
        """
        ui.load_arrays(1, self.x, self.y)
        ui.set_source("polynom1d.p")
        ui.thaw("p.c1")
        ui.set_method("levmar")
        ui.fit()
        model = ui.get_model_component("p")
        expected = [0, 1]
        observed = [model.c0.val, model.c1.val]
        assert_almost_equal(observed, expected)

    # Using skipIf directly as we need these tests to run when there is no pytest installed, so we can't use
    # the decorators in sherpa.utils.testing
    @unittest.skipIf(not has_package_from_list('astropy.io.fits', 'pycrates'), reason="Requires fits backend")
    def test_fits_io(self):
        """
        Test that basic FITS I/O functions work.

        This test ensures that the FITS backend can be used to perform basic
        I/O functions.
        """
        ui.load_pha(self.fits)
        with NamedTemporaryFile() as f:
            ui.save_pha(f.name, ascii=False, clobber=True)

    @unittest.skipIf(not has_package_from_list('sherpa.astro.xspec'), reason="Requires xspec")
    def test_xspec(self):
        """
        Perform a very simple fit with an xspec model.

        Also check that the results make sense.

        This test proves that the xspec extension properly works, and that there are no obvious building, linking, or
        environment issues that would prevent the xspec model from running.
        """
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
        """
        This is a debug test that always fails.

        It may be used to check that error conditions upstream are correctly
        handled by scripts calling the Smoke Test. Whether or not this test is run depends on the option passed to the
        run function in this module.
        """
        self.fail("Requested to fail")


class SmokeTestSuite(unittest.TestSuite):
    """
    Smoke test suite.

    At this point the Smoke Test is comprised of a single unittest test case. However, it may be extended in the future,
    so we set up a test suite.

    Also, the suite disables the always-failing debugging test depending on the options with which the Smoke Test
    is running.
    """

    # classes in the suite
    test_classes = (SmokeTest,)

    # tests that always fail
    failure_methods = ('test_failure',)

    def __init__(self, *args, **kwargs):
        """
        The init method passes all the arguments along to the unittest.TestSuite init method.

        It also accepts a "require_failure" keyword argument used to decide whether the debugging test method must
        fail.
        """
        self.require_failure = kwargs.pop("require_failure", False)
        unittest.TestSuite.__init__(self, *args, **kwargs)

        for test_class in self.test_classes:
            self.addTest(unittest.makeSuite(test_class))

        self._set_skips()

    def _set_skips(self):
        """
        Make sure the test method that always fails is skipped, unless "require_failure" is set to True.
        """
        if not self.require_failure:
            for case in self:
                for test in case:
                    if test.id().split(".")[-1] in self.failure_methods:
                        setattr(test, 'setUp', lambda: test.skipTest('Smoke Test not required to fail, skipping'))


def _import_error(module, name):
    return "ERROR: Requested {} as {} but module not found".format(module, name)
