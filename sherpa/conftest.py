#
#  Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021
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

import os
import re
import logging

import numpy as np
from numpy import VisibleDeprecationWarning

import pytest

from sherpa.utils.testing import get_datadir, set_datadir

try:
    from astropy.io.fits.verify import VerifyWarning
    have_astropy = True
except ImportError:
    have_astropy = False

try:
    from sherpa.astro import xspec
    has_xspec = True
except ImportError:
    has_xspec = False


import sherpa.plot
import sherpa.astro.plot
from sherpa.plot import dummy_backend
try:
    from sherpa.plot import pylab_backend
except ImportError:
    pylab_backend = None


# In some instances, if some xvfb processes did not stop cleanly
# pytest-xfvb starts complaining that a virtual screen is already on
# the following code works around that. The issue is hard to reproduce
# but I (OL) tested it locally as it reappeared on my workstation.
try:
    import random
    from pyvirtualdisplay import abstractdisplay
    abstractdisplay.RANDOMIZE_DISPLAY_NR = True
    abstractdisplay.random = random
    random.seed()
except ImportError:
    pass


TEST_DATA_OPTION = "--test-data"


# Follow https://docs.pytest.org/en/latest/example/simple.html#control-skipping-of-tests-according-to-command-line-option
# for adding a command-line option to let slow-running tests be run
# (not by default), which combines with existing code and options.
#
def pytest_addoption(parser):
    parser.addoption("-D", TEST_DATA_OPTION, action="store",
                     help="Alternative location of test data files")

    parser.addoption("--runslow", action="store_true", default=False,
                     help="run slow tests")

    parser.addoption("--runzenodo", action="store_true", default=False,
                     help="run tests that query Zenodo (requires internet)")


def pytest_collection_modifyitems(config, items):

    # Skip tests unless --runxxx given in cli
    #
    for label in ["slow", "zenodo"]:
        opt = "--run{}".format(label)
        if config.getoption(opt):
            continue

        skip = pytest.mark.skip(reason="need {} option to run".format(opt))
        for item in items:
            if label in item.keywords:
                item.add_marker(skip)


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

            # Matplotlib version 2 warnings (from HTML notebook represention)
            #
            r'np.asscalar\(a\) is deprecated since NumPy v1.16, use a.item\(\) instead',
            r"Using or importing the ABCs from 'collections' instead of from 'collections.abc' is deprecated since Python 3.3,and in 3.9 it will stop working",

        ],
    UserWarning:
        [
            r"File '/data/regression_test/master/in/sherpa/aref_sample.fits' does not have write permission.  Changing to read-only mode.",
            r"File '/data/regression_test/master/in/sherpa/aref_Cedge.fits' does not have write permission.  Changing to read-only mode.",
            r"Converting array .* to numpy array",

            # Matplotlib version 2 warnings (from HTML notebook represention)
            #
            r'Attempting to set identical bottom==top results\nin singular transformations; automatically expanding.\nbottom=1.0, top=1.0',

            # See issue #571 and tests
            # test_ui_plot.py
            #    test_plot_fit_xxx_pylab
            #    test_plot_fit_xxx_overplot_pylab
            r"Attempted to set non-positive xlimits for log-scale axis; invalid limits will be ignored.",
            r"Attempted to set non-positive left xlim on a log-scaled axis.*",
        ],
    RuntimeWarning:
        [r"invalid value encountered in sqrt",
         # See https://github.com/ContinuumIO/anaconda-issues/issues/6678
         r"numpy.dtype size changed, may indicate binary " +
         r"incompatibility. Expected 96, got 88",
         # See https://github.com/numpy/numpy/pull/432
         r"numpy.ufunc size changed",
         # numpy 1.20 shows this in some tests
         r"numpy.ndarray size changed, may indicate binary "
         ],
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
            # tests in sherpa/astro/ui/tests/test_astro_ui_unit.py
            # seen in some macOS test runs (e.g. pip, not conda) and python 3.8
            r"unclosed file .*/data.dat'.* closefd=True>",
            r"unclosed file .*/model.dat'.* closefd=True>",
            r"unclosed file .*/resid.out'.* closefd=True>",
        ],
    VisibleDeprecationWarning:
        [],
}

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


"""
# Currently no matplotlib warnings to ignore, but leave code so it is
# easy to add one back in

try:
    from matplotlib import MatplotlibDeprecationWarning

    matplotlib_warnings = {
        MatplotlibDeprecationWarning:
            []
    }
    known_warnings.update(matplotlib_warnings)
except ImportError:
    pass
"""


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

    known = check_known_warning

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


def check_known_warning(warning):
    """Return True if this is an "allowed" warning."""

    message = warning.message
    for known_warning in known_warnings[type(message)]:
        pattern = re.compile(known_warning)
        if pattern.match(str(message)):
            return True

    return False


@pytest.fixture
def is_known_warning():
    """Returns a function that returns True if this is an "allowed" warning.

    It is not expected that this will see much use.
    """

    return check_known_warning


def pytest_configure(config):
    """
    This configuration hook overrides the default mechanism for test data self-discovery, if the --test-data command line
    option is provided

    It also adds support for the "slow" and "zenodo" test markers.

    Parameters
    ----------
    config standard service injected by pytest
    """
    try:
        path = config.getoption(TEST_DATA_OPTION)
        if path:
            set_datadir(path)
    except ValueError:  # option not defined from command line, no-op
        pass

    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line("markers", "zenodo: indicates requires internet access to Zenodo")


@pytest.fixture(scope="session")
def test_data_path():
    """
    Fixture for getting the path to the test data

    Returns
    -------

    path : basestring
        The string with the path to the main test data folder.
    """
    path = get_datadir()

    if path is None:
        raise RuntimeError("Test needs the requires_data decorator")

    return path


@pytest.fixture(scope="session")
def make_data_path(test_data_path):
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

    def wrapped(arg):
        return os.path.join(test_data_path, arg)

    return wrapped


def run_thread_function(name, scriptname, test_data_path):
    """Run a regression test from the sherpa-test-data submodule.

    Parameters
    ----------
    name : string
       The name of the science thread to run (e.g., pha_read,
       radpro). The name should match the corresponding thread
       name in the sherpa-test-data submodule. See examples below.
    scriptname : string
       The suffix of the test script file name, usually "fit.py."
    test_data_path : string
       The path to the test data folder

    Returns
    -------
    localsyms : dict
        Any variables created by the script (includes model
        components).

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
    from sherpa.astro import ui

    scriptname = name + "-" + scriptname
    cwd = os.getcwd()
    os.chdir(test_data_path)

    localsyms = {}

    def assign_model(name, val):
        localsyms[name] = val

    old_assign_model = ui.get_model_autoassign_func()

    try:
        with open(scriptname, "rb") as fh:
            cts = fh.read()
        ui.set_model_autoassign_func(assign_model)
        exec(compile(cts, scriptname, 'exec'), {}, localsyms)
    finally:
        ui.set_model_autoassign_func(old_assign_model)
        os.chdir(cwd)

    return localsyms


@pytest.fixture
def run_thread(test_data_path):
    """
    Fixture that returns a function that can be used to run a thread.

    Returns
    -------
    run_thread : function
        The returned function is basically the run_thread_function function, using the test data path found at runtime
    """
    def run_thread_wrapper(name, scriptname='fit.py'):
        return run_thread_function(name, scriptname, test_data_path)

    return run_thread_wrapper


@pytest.fixture
def clean_astro_ui():
    """Ensure sherpa.astro.ui.clean is called before AND after the test.

    This also resets the XSPEC settings (if XSPEC support is provided).

    See Also
    --------
    clean_ui

    Notes
    -----
    It does NOT change the logging level; perhaps it should, but the
    screen output is useful for debugging at this time.
    """
    from sherpa.astro import ui

    if has_xspec:
        old_xspec = xspec.get_xsstate()
    else:
        old_xspec = None

    ui.clean()
    yield

    ui.clean()
    if old_xspec is not None:
        xspec.set_xsstate(old_xspec)


@pytest.fixture
def clean_ui():
    """Ensure sherpa.ui.clean is called before AND after the test.

    See Also
    --------
    clean_astro_ui

    Notes
    -----
    It does NOT change the logging level; perhaps it should, but the
    screen output is useful for debugging at this time.
    """
    from sherpa import ui

    ui.clean()
    yield

    ui.clean()


@pytest.fixture
def restore_xspec_settings():
    """Fixture to ensure that XSPEC settings are restored after the test.

    The aim is to allow the test to change XSPEC settings but to
    ensure they are restored to their default values after the test.

    The logic for what should and should not be reset is left to
    xspec.get_state/set_state. There are known issues with trying to
    reset "XSET" values (e.g. you can only set them to "", not
    "delete" them, and it is up to the model implementation to note
    the value has changed).
    """

    try:
        from sherpa.astro import xspec
    except ImportError:
        # can I return from here safely?
        return

    # grab initial values
    #
    state = xspec.get_xsstate()

    # return control to test
    yield

    # clean up after test
    xspec.set_xsstate(state)


@pytest.fixture
def reset_seed(request):
    """Force a random seed after the test.

    The random seed is set to np.random.seed() after the
    test is done. It is expected that the test sets the
    seed to an explicit value. Ideally we would use the
    new NumPy RNG but we still need to support older NumPy
    versions.

    """

    yield
    np.random.seed()


@pytest.fixture
def hide_logging():
    """Set Sherpa's logging to ERROR for the test.

    This code is somehwat redundant with the decorator
    in utils/logging:SherpaVerbosity and, should this
    fixture ever need to be redone, it might be worth
    investigating if the same decorator can be used.
    """

    logger = logging.getLogger('sherpa')
    olvl = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield
    logger.setLevel(olvl)


@pytest.fixture
def old_numpy_printing():
    """Force NumPy to use old-style printing for the test.

    This is only needed whilst we still support NumPy 1.13
    (I think). Calling this fixture will ensure we have
    consistent printing of NumPy arrays (to some degree
    anyway).

    The original printoptions is reset after the test is
    run.
    """

    oldopts = np.get_printoptions()
    if 'legacy' in oldopts:
        np.set_printoptions(legacy='1.13')

    yield

    if 'legacy' in oldopts:
        np.set_printoptions(legacy=oldopts['legacy'])


PLOT_BACKENDS = [dummy_backend]
if pylab_backend is not None:
    PLOT_BACKENDS.append(pylab_backend)


@pytest.fixture(params=PLOT_BACKENDS)
def override_plot_backend(request):
    """Override the plot backend for this test

    Runs the test with the given plot backend and then restores the
    original value.

    Note that this does not reload the ui layers, which will have
    retain reference to the original backend throughout.

    """

    old = sherpa.plot.backend

    # safety check
    #
    assert old.name == sherpa.astro.plot.backend.name

    changed = old.name != request.param.name
    if changed:
        sherpa.plot.backend = request.param
        sherpa.astro.plot.backend = request.param

    yield
    if changed:
        sherpa.plot.backend = old
        sherpa.astro.plot.backend = old
