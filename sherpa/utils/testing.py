#
#  Copyright (C) 2017, 2020 - 2025
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

import importlib
import os

from sherpa.plot import TemporaryPlottingBackend

try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False

from sherpa.utils.err import RuntimeErr


def _get_datadir():
    """Return the location of the test data files, if installed.

    Returns
    -------
    path : str or None
        The path to the Sherpa test data directory or None.
    """

    try:
        import sherpatest
        datadir = os.path.dirname(sherpatest.__file__)
    except ImportError:
        try:
            import sherpa
            datadir = os.path.join(os.path.dirname(sherpa.__file__), os.pardir,
                                   'sherpa-test-data', 'sherpatest')
            if not os.path.exists(datadir) or not os.path.isdir(datadir) \
               or not os.listdir(datadir):
                # The dir is empty, maybe the submodule was not initialized
                datadir = None
        except ImportError:
            # neither sherpatest nor sherpa can be found, falling back to None
            datadir = None
    return datadir


DATADIR = _get_datadir()


def set_datadir(datadir):
    """Set the data directory.

    Parameters
    ----------
    datadir : str
        The name to the data directory. It must exist.

    Raises
    ------
    OSError
        If datadir is not a directory or is empty.

    """

    if not os.path.exists(datadir) or not os.path.isdir(datadir) \
       or not os.listdir(datadir):
        raise OSError(f"datadir={datadir} is empty or not a directory")

    global DATADIR
    DATADIR = datadir


def get_datadir():
    """Return the data directory.

    Returns
    -------
    datadir : str or None
        The name to the data directory, if it exists.

    """

    return DATADIR


def has_package_from_list(*packages):
    """
    Returns True if at least one of the ``packages`` args is importable.
    """
    for package in packages:
        try:
            importlib.import_module(package)
            return True
        except (ImportError, RuntimeErr):
            # Apparently we need to also catch RuntimeErr (the Sherpa
            # version) as sherpa.image.DS9 can raise it (e.g. when
            # DS9 is not installed).
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
        condition = DATADIR is None
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

    def requires_fits(test_function):
        """Returns True if a working backend for FITS I/O is importable.

        Used to skip tests requiring fits_io. The dummy backend itself
        is not a "working backend" in the sense that it cannot be used to
        read or write files.

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

    def requires_region(test_function):
        """Decorator for test functions requiring 2D region support"""
        return requires_package("Region required", 'sherpa.astro.utils._region')(test_function)

    def requires_wcs(test_function):
        """Decorator for test functions requiring WCS support"""
        return requires_package("WCS required", 'sherpa.astro.utils._wcs')(test_function)

    def requires_ds9(test_function):
        """Decorator for test functions requiring ds9

        This also ensures that the test is run in the same xdist
        'group' as the other DS9 tests, which means they will all be
        run sequentially, avoiding issue #1704 when using pytest-xdist
        to run multiple tests in parallel.

        """

        pkg = requires_package('ds9 required', 'sherpa.image.ds9_backend')
        grp = pytest.mark.xdist_group("ds9-tests")

        return grp(pkg(test_function))

    def requires_xspec(test_function):
        return requires_package("xspec required", "sherpa.astro.xspec")(test_function)

    def requires_psf(test_function):
        """Decorator for test functions requiring PSF support"""
        # Note that 'from sherpa.utils import _psf' returns None
        # but 'import sherpa.utils._psf' fails.
        #
        return requires_package("FFT/PSF rquired", 'sherpa.utils._psf')(test_function)

else:

    def wrapped():
        raise ImportError(PYTEST_MISSING_MESSAGE)

    def make_fake():
        def wrapper(*args, **kwargs):
            return wrapped
        return wrapper

    requires_data = make_fake()

    requires_fits = make_fake()

    requires_group = make_fake()

    requires_stk = make_fake()

    requires_region = make_fake()

    requires_wcs = make_fake()

    requires_ds9 = make_fake()

    requires_xspec = make_fake()

    requires_psf = make_fake()

    def requires_package(*args):
        return make_fake()


PYTEST_MISSING_MESSAGE = "Package `pytest` is missing. Please install `pytest` before running tests or using the test" \
                         "decorators"
