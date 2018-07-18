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

import pytest
import os
import sys
import re
from sherpa.utils import SherpaTestCase

from six.moves import reload_module

try:  # Python 3
    from unittest import mock
except ImportError:  # Python 2
    import mock

try:
    from astropy.io.fits.verify import VerifyWarning
    have_astropy = True
except ImportError:
    have_astropy = False


TEST_DATA_OPTION = "--test-data"


def pytest_addoption(parser):
    parser.addoption("-D", TEST_DATA_OPTION, action="store",
                     help="Alternative location of test data files")


# Whilelist of known warnings. One can associate different warning messages
# to the same warning class
known_warnings = {
    DeprecationWarning:
        [
            r"unorderable dtypes.*",
            r"Non-string object detected for the array ordering.*",
            r"using a non-integer number instead of an integer will result in an error in the future",
            r"Use load_xstable_model to load XSPEC table models",
            #  This does not have to do with Sherpa and is coming from some versions of
            #  jupyter_client
            r"metadata .* was set from the constructor.*",
        ],
    UserWarning:
        [
            r"File '/data/regression_test/master/in/sherpa/aref_sample.fits' does not have write permission.  Changing to read-only mode.",
            r"File '/data/regression_test/master/in/sherpa/aref_Cedge.fits' does not have write permission.  Changing to read-only mode."
        ],
    RuntimeWarning:
        [r"invalid value encountered in sqrt",
         # See https://github.com/ContinuumIO/anaconda-issues/issues/6678
         r"numpy.dtype size changed, may indicate binary " +
         r"incompatibility. Expected 96, got 88"
         ],
}

if sys.version_info >= (3, 2):
    python3_warnings = {
        ResourceWarning:
            [
                r"unclosed file .*king_kernel.txt.* closefd=True>",
                r"unclosed file .*phas.dat.* closefd=True>",
                r"unclosed file .*data.txt.* closefd=True>",
                r"unclosed file .*cstat.dat.* closefd=True>",
                r"unclosed file .*data1.dat.* closefd=True>",
                r"unclosed file .*aref_Cedge.fits.* closefd=True>",
                r"unclosed file .*aref_sample.fits.* closefd=True>",
                r"unclosed file .*/tmp.* closefd=True>",
                # added for sherpa/astro/ui/tests/test_astro_ui_utils_unit.py
                r"unclosed file .*/dev/null.* closefd=True>",
                r"unclosed file .*table.txt.* closefd=True>",
            ],
        RuntimeWarning:
            [r"invalid value encountered in sqrt",
             # See https://github.com/ContinuumIO/anaconda-issues/issues/6678
             r"numpy.dtype size changed, may indicate binary " +
             r"incompatibility. Expected 96, got 88"
             ],
    }
    known_warnings.update(python3_warnings)


if have_astropy:
    astropy_warnings = {
        # See bug #372 on GitHub
        ImportWarning:
        [
            r"can't resolve.*__spec__.*__package__.*",
        ],
        VerifyWarning:
        [
            r"Invalid keyword for column.*",
        ],
    }
    known_warnings.update(astropy_warnings)


# Can this be replaced by the warning support added in pytest 3.1?
# See https://docs.pytest.org/en/latest/warnings.html#warnings
#
@pytest.fixture(scope="function", autouse=True)
def capture_all_warnings(request, recwarn, pytestconfig):
    """
    This fixture will run automatically before and after every test function is executed.
    It uses pytest's infrastructure to get all recorded warnings and match them against the while list. If an
    unknown warning is found, then the fixture finalizer will fail the specific test function.

    In the verbose pytest report the test function will show twice if an unknown warning is captured: one with the
    actual result of the test and one with an ERROR. The warning will be shown as part of the stderr stream.

    Parameters
    ----------
    request standard injected service for pytest fixtures
    recwarn injected pytest service for accessing recorded warnings
    pytestconfig injected service for accessing the configuration data

    """
    def known(warning):
        message = warning.message
        for known_warning in known_warnings[type(message)]:
            pattern = re.compile(known_warning)
            if pattern.match(str(message)):
                return True
        return False

    def fin():
        warnings = [w for w in recwarn.list
                    if type(w.message) not in known_warnings or not known(w)]

        nwarnings = len(warnings)
        if nwarnings > 0:
            # Only print out the first time a warning is seen; we could
            # report this to the user explicitly, but rely on the
            # counter information (i.e. i+1/nwarnings) to show that
            # there were repeats.
            #
            swarnings = set([])
            print("*** Warnings created: {}".format(nwarnings))
            for i, w in enumerate(warnings):
                sw = str(w)
                if sw in swarnings:
                    continue
                print("{}/{} {}".format(i + 1, nwarnings, w))
                swarnings.add(sw)

        assert 0 == nwarnings

    request.addfinalizer(fin)


def pytest_configure(config):
    """
    This configuration hook overrides the default mechanism for test data self-discovery, if the --test-data command line
    option is provided

    Parameters
    ----------
    config standard service injected by pytest
    """
    try:
        path = config.getoption(TEST_DATA_OPTION)
        if path:
            SherpaTestCase.datadir = path
    except ValueError:  # option not defined from command line, no-op
        pass


@pytest.fixture(scope="session")
def make_data_path():
    """
    Fixture for tests requiring the test data dir. It returns a function that can be used to make paths by using
    path elements relative to the test data folder (which is flat, so in principle only the first element is required)

    Returns
    -------
    make_data_path : func
        A function that accepts a list of path elements to be joined with
        the base data dir path. This function exits with a RuntimeError
        if the data directory is None, pointing out the requires_data
        decorator is needed.
    """
    path = SherpaTestCase.datadir

    def wrapped(arg):
        if path is None:
            raise RuntimeError("Test needs the requires_data decorator")
        return os.path.join(path, arg)

    return wrapped


@pytest.fixture
def mock_chips(monkeypatch, tmpdir, request):
    """
    Fixture for tests mocking chips

    Returns
    -------
    The tuple (backend, mock_chips)
    """

    # First, inject a mock chips module in the backend.
    chips = mock.MagicMock()
    monkeypatch.setitem(sys.modules, name="pychips", value=chips)

    # figure out what IO module we can use
    try:
        import pycrates
        io = "crates"
    except ImportError:
        io = "pyfits"  # Even if this is not available, config code will fall back to dummy

    # Now, write a fake configuration file to a temporary location
    config = tmpdir.mkdir("config").join("sherpa.rc")
    config.write("""
[options]
plot_pkg : chips
io_pkg : {}
    """.format(io))

    # Then, inject a function that returns the fake file
    def get_config():
        return str(config)
    import sherpa
    monkeypatch.setattr(sherpa, name="get_config", value=get_config)

    # Force reload of sherpa modules that might have already read the configuration
    from sherpa import plot
    from sherpa.astro import plot as astro_plot

    reload_module(plot)
    reload_module(astro_plot)

    # Force a reload, to make sure we always return a fresh instance, so we track the correct mock object
    from sherpa.plot import chips_backend
    reload_module(chips_backend)

    def fin():
        monkeypatch.undo()
        reload_module(sherpa)
        reload_module(plot)
        reload_module(astro_plot)
        reload_module(sherpa.all)
        reload_module(sherpa.astro.all)  # These are required because otherwise Python will not match imported classes.

    request.addfinalizer(fin)

    return chips_backend, chips
