#
#  Copyright (C) 2017, 2020  Smithsonian Astrophysical Observatory
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
import os
import importlib
import logging
import pkg_resources

import numpy

from sherpa.utils._utils import sao_fcmp

try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False


warning = logging.getLogger(__name__).warning


# I am not convinced the conversion from __file__ access to
# pkg_resources is correct here. Using __file__ is frowned upon,
# hence the change, and as we can not assume Python 3.7 or later
# at the moment we can not use importlib_resources.
#
def _get_datadir():
    """Setup the location of the Sherpa test data, if installed."""

    try:
        import sherpatest
        datadir = pkg_resources.resource_filename('sherpatest', '')

    except ImportError:
        try:
            import sherpa
            datadir = pkg_resources.resource_filename('sherpa', '')
            datadir = os.path.join(datadir, os.pardir,
                                   'sherpa-test-data', 'sherpatest')
        except ImportError:
            # neither sherpatest nor sherpa can be found, falling back to None
            datadir = None

    # Check the directory exists
    if datadir is None or not os.path.exists(datadir) or not os.listdir(datadir):
        return None

    return datadir


DATADIR = _get_datadir()

def get_datadir():
    """Return the location of the Sherpa test data.

    The data directory is determined either by the existance of the
    sherpatest module, sherpa-test-directory submodule, or
    by being set manually with set_datadir.
    """

    return DATADIR


def set_datadir(dname):
    """Set the location of the Sherpa test data.

    """

    global DATADIR
    if not os.path.exists(dname) or not os.listdir(dname):
        warning("Unable to set the datadir to {}".format(dname))
        dname = None

    DATADIR = dname
    SherpaTestCase.datadir = dname


class SherpaTestCase(unittest.TestCase):
    """
    Base class for Sherpa unit tests. The use of this class is deprecated in favor of pytest functions.
    """

    # The location of the Sherpa test data (it is optional)
    datadir = get_datadir()

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


if HAS_PYTEST:
    #  Pytest cannot be assumed to be installed by the regular user, unlike unittest, which is part of Python's
    #  standard library. The decorator will be defined if pytest is missing, but if the tests are run they throw
    #  and exception prompting users to install pytest, in those cases where pytest is not installed automatically.

    def requires_data(test_function):
        """
         Decorator for functions requiring external data (i.e. data not distributed with Sherpa
         itself) is missing.  This is used to skip tests that require such data.

         See PR #391 for why this is a function: https://github.com/sherpa/sherpa/pull/391
         """
        condition = SherpaTestCase.datadir is None
        msg = "required test data missing"
        return pytest.mark.skipif(condition, reason=msg)(test_function)


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


    def requires_plotting(test_function):
        """
        Decorator for test functions requiring a plotting library.
        """
        packages = ('pylab', )
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
        packages = ('astropy.io.fits',
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
        def wrapper(*args, **kwargs):
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

    def requires_package(*args):
        return make_fake()


PYTEST_MISSING_MESSAGE = "Package `pytest` is missing. Please install `pytest` before running tests or using the test" \
                         "decorators"
