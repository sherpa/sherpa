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

import numpy
import unittest
import os
import importlib

from sherpa.utils._utils import sao_fcmp

try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False


def _get_datadir():
    import os
    try:
        import sherpatest
        datadir = os.path.dirname(sherpatest.__file__)
    except ImportError:
        try:
            import sherpa
            datadir = os.path.join(os.path.dirname(sherpa.__file__), os.pardir,
                                   'sherpa-test-data', 'sherpatest')
            if not os.path.exists(datadir) or not os.listdir(datadir):
                # The dir is empty, maybe the submodule was not initialized
                datadir = None
        except ImportError:
            # neither sherpatest nor sherpa can be found, falling back to None
            datadir = None
    return datadir


class SherpaTestCase(unittest.TestCase):
    """
    Base class for Sherpa unit tests. The use of this class is deprecated in favor of pytest functions.
    """

    # The location of the Sherpa test data (it is optional)
    datadir = _get_datadir()

    def make_path(self, *segments):
        """Add the segments onto the test data location.

        Parameters
        ----------
        *segments
           Path segments to combine together with the location of the
           test data.

        Returns
        -------
        fullpath : None or string
           The full path to the repository, or None if the
           data directory is not set.

        """
        if self.datadir is None:
            return None
        return os.path.join(self.datadir, *segments)

    # What is the benefit of this over numpy.testing.assert_allclose(),
    # which was added in version 1.5 of NumPy?
    def assertEqualWithinTol(self, first, second, tol=1e-7, msg=None):
        """Check that the values are equal within an absolute tolerance.

        Parameters
        ----------
        first : number or array_like
           The expected value, or values.
        second : number or array_like
           The value, or values, to check. If first is an array, then
           second must be an array of the same size. If first is
           a scalar then second can be a scalar or an array.
        tol : number
           The absolute tolerance used for comparison.
        msg : string
           The message to display if the check fails.

        """

        self.assertFalse(numpy.any(sao_fcmp(first, second, tol)), msg)

    def assertNotEqualWithinTol(self, first, second, tol=1e-7, msg=None):
        """Check that the values are not equal within an absolute tolerance.

        Parameters
        ----------
        first : number or array_like
           The expected value, or values.
        second : number or array_like
           The value, or values, to check. If first is an array, then
           second must be an array of the same size. If first is
           a scalar then second can be a scalar or an array.
        tol : number
           The absolute tolerance used for comparison.
        msg : string
           The message to display if the check fails.

        """

        self.assertTrue(numpy.all(sao_fcmp(first, second, tol)), msg)

    # for running regression tests from sherpa-test-data
    def run_thread(self, name, scriptname='fit.py'):
        """Run a regression test from the sherpa-test-data submodule.

        Parameters
        ----------
        name : string
           The name of the science thread to run (e.g., pha_read,
           radpro). The name should match the corresponding thread
           name in the sherpa-test-data submodule. See examples below.
        scriptname : string
           The suffix of the test script file name, usually "fit.py."

        Examples
        --------
        Regression test script file names have the structure
        "name-scriptname.py." By default, scriptname is set to "fit.py."
        For example, if one wants to run the regression test
        "pha_read-fit.py," they would write

        >>> run_thread("pha_read")

        If the regression test name is "lev3fft-bar.py," they would do

        >>> run_thread("lev3fft", scriptname="bar.py")

        """

        scriptname = name + "-" + scriptname
        self.locals = {}
        cwd = os.getcwd()
        os.chdir(self.datadir)
        try:
            with open(scriptname, "rb") as fh:
                cts = fh.read()
            exec(compile(cts, scriptname, 'exec'), {}, self.locals)
        finally:
            os.chdir(cwd)


def has_package_from_list(*packages):
    """
    Returns True if at least one of the ``packages`` args is importable.
    """
    for package in packages:
        try:
            importlib.import_module(package)
            return True
        except:
            pass
    return False


def requires_package(msg=None, *packages):
    """
    Decorator for test functions requiring specific packages.
    """
    condition = has_package_from_list(*packages)
    msg = msg or "required module missing among {}.".format(
        ", ".join(packages))

    def decorator(test_function):
        return pytest.mark.skipif(not condition, reason=msg)(test_function)

    return decorator


if HAS_PYTEST:
    #  Pytest cannot be assumed to be installed by the regular user, unlike unittest, which is part of Python's
    #  standard library. The decorator will be defined if pytest is missing, but if the tests are run they throw
    #  and exception prompting users to install pytest, in those cases where pytest is not installed automatically.

    requires_data = pytest.mark.skipif(SherpaTestCase.datadir is None,
                                       reason="required test data missing")


    def requires_plotting(test_function):
        """
        Decorator for test functions requiring a plotting library.
        """
        packages = ('pylab', 'pychips')
        msg = "plotting backend required"
        return requires_package(msg, *packages)(test_function)


    def requires_pylab(test_function):
        """
        Returns True if the pylab module is available (pylab).
        Used to skip tests requiring matplotlib
        """
        packages = ('pylab',
                    )
        msg = "matplotlib backend required"
        return requires_package(msg, *packages)(test_function)


    def requires_fits(test_function):
        """
        Returns True if there is an importable backend for FITS I/O.
        Used to skip tests requiring fits_io
        """
        packages = ('pyfits',
                    'astropy.io.fits',
                    'pycrates',
                    )
        msg = "FITS backend required"
        return requires_package(msg, *packages)(test_function)


    def requires_group(test_function):
        """Decorator for test functions requiring group library"""
        return requires_package("group library required", 'group')(test_function)


    def requires_stk(test_function):
        """Decorator for test functions requiring stk library"""
        return requires_package("stk library required", 'stk')(test_function)


    def requires_ds9(test_function):
        """Decorator for test functions requiring ds9"""
        return requires_package('ds9 required', 'sherpa.image.ds9_backend')(test_function)


    def requires_xspec(test_function):
        return requires_package("xspec required", "sherpa.astro.xspec")(test_function)

else:

    def wrapped():
        raise ImportError(PYTEST_MISSING_MESSAGE)


    def make_fake():
        def wrapper(*arg, **kwargs):
            return wrapped
        return wrapper

    requires_data = make_fake()

    requires_plotting = make_fake()

    requires_pylab = make_fake()

    requires_fits = make_fake()

    requires_group = make_fake()

    requires_stk = make_fake()

    requires_ds9 = make_fake()

    requires_xspec = make_fake()


PYTEST_MISSING_MESSAGE = "Package `pytest` is missing. Please install `pytest` before running tests or using the test" \
                         "decorators"
