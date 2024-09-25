#
#  Copyright (C) 2016 - 2024
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
try:
    from numpy.exceptions import VisibleDeprecationWarning
except ImportError:
    # Must be Earlier than NumPy 1.25
    from numpy import VisibleDeprecationWarning

import pytest

from sherpa.utils.err import RuntimeErr
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


from sherpa.plot.backends import BaseBackend, IndepOnlyBackend, PLOT_BACKENDS
from sherpa.plot import TemporaryPlottingBackend

FUNCTIONAL_PLOT_BACKENDS = []
'''Unfortunately, there is no metadata that tell us which plotting backends are
dummies and which ones actually produce useful plots. That's because
"one person's dummy might be another person's use case", e.g. a backend might be a
dummy for image, but not for a histogram.

The goal of this list is to be able to run a large number of plotting tests
with any installed, functional backend. These tests mostly execute
`plot`, `histogram`, `vline`, `hline`, etc, but not image, 3D plots etc.
So, a better name might be BACKENDS_FUNCTIONAL_FOR_1D_PLOT, but this can be
refined when we actually have more plotting backends and more multi-D plot tests.

Thus, we can't just autogenerate this list, but need to make it by hand here.
'''

try:
    from sherpa.plot.pylab_backend import PylabBackend
    FUNCTIONAL_PLOT_BACKENDS.append('pylab')

    # Override the matplotlib backend from TkAgg to Agg to try
    # and avoid errors from multiprocessing thanks to the error
    # "Tcl_AsyncDelete: async handler deleted by the wrong thread"
    # See issue #1590.
    #
    # This has been moved over from the existing code, which was
    # added after this commit was made. Is this correct?
    #
    import matplotlib
    matplotlib.use("agg")

except ImportError:
    pass


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
            r"Use load_xstable_model to load XSPEC table models",

            # NumPy 1.25 warnings that are raised by (mid-2023) crates code.
            # Hopefully this can be removed by December 2023.
            #
            r"Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future. Ensure you extract a single element from your array before performing this operation. \(Deprecated NumPy 1.25.\)",

            # See https://github.com/sherpa/sherpa/issues/1953
            # Technically this should only be needed if bokeh is installed
            # but it doesn't seem worth setting up that machinery.
            #
            "\nPyarrow will become a required dependency of pandas ",

            # See issues
            # https://github.com/sherpa/sherpa/issues/2158
            # https://github.com/sherpa/sherpa/issues/2007
            #
            # For now we hide them.
            #
            r".* is multi-threaded, use of fork\(\) may lead to deadlocks in the child\."

        ],
    UserWarning:
        [
            r"File '/data/regression_test/master/in/sherpa/aref_sample.fits' does not have write permission.  Changing to read-only mode.",
            r"File '/data/regression_test/master/in/sherpa/aref_Cedge.fits' does not have write permission.  Changing to read-only mode.",
            r"Converting array .* to numpy array",

            # Matplotlib version 2 warnings (from HTML notebook representation)
            #
            r'Attempting to set identical bottom==top results\nin singular transformations; automatically expanding.\nbottom=1.0, top=1.0',

            # See issues #571, #1584, and matplotlib issue
            # https://github.com/matplotlib/matplotlib/issues/9970
            #
            # The exact warning message appears to depend on the
            # matplotlib version and is related to doing
            #
            #     plot_fit_xxx(xlog=True)
            #     plot_fit_yyy(...)
            #
            # with the second call causing the problem (there may be
            # other cases, but this is the one that fits the mpl
            # issue).
            #
            r"Attempted to set non-positive xlimits for log-scale axis; invalid limits will be ignored.",
            r"Attempted to set non-positive left xlim on a log-scaled axis.*",
            r"Attempt to set non-positive xlim on a log-scaled axis will be ignored."

        ],
    RuntimeWarning:
        [r"invalid value encountered in sqrt",
         ],
    ResourceWarning:
        [
        ],
    VisibleDeprecationWarning:
        [],
    ImportWarning:
    [
        # It is not clear why this happens - seen when testing
        # https://github.com/sherpa/sherpa/pull/1752 - and the problem appears
        # to be a setuptools/pip interaction, based on
        #   https://github.com/pypa/setuptools/issues/2104
        #   https://github.com/pypa/setuptools/issues/2052
        #   https://github.com/pypa/setuptools/issues/1383
        # At present Sherpa forces an old version of setuptoools (<60) so it is
        # unlikely to be fixed. When the setuptools restriction is removed we can
        # hopefully remove this.
        r"VendorImporter.find_spec\(\) not found; falling back to find_module\(\)",
    ]
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
        warnings = [w for w in recwarn.list if not known(w)]

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
    for known_warning in known_warnings.get(type(message), []):
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
def hide_logging():
    """Set Sherpa's logging to ERROR for the test.

    This code is somewhat redundant with the decorator
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


@pytest.fixture(params=PLOT_BACKENDS.keys())
def all_plot_backends(request):
    """Override the plot backend for this test

    Runs the test with the given plot backend and then restores the
    original value.

    Note that this does not reload the ui layers, which will
    retain a reference to the original backend throughout.
    """
    with TemporaryPlottingBackend(request.param):
        yield
        if request.param in ('pylab', 'PylabErrorArea'):
            from matplotlib import pyplot as plt
            plt.close(fig="all")


@pytest.fixture(params=PLOT_BACKENDS.keys())
def all_plot_backends_astro_ui(request):
    """Override the plot backend for the astro ui session.

    Runs the test with the given plot backend and then restores the
    original value.
    """

    # Unfortunately can not easily use TemporaryPlottingBackend here
    #
    from sherpa.plot import backend
    old = backend

    from sherpa.astro.ui import set_plot_backend
    set_plot_backend(request.param)

    yield

    # Backend-specific clean up.
    if request.param == 'pylab':
        from matplotlib import pyplot as plt
        plt.close(fig="all")

    # Restore the original backend.
    #
    set_plot_backend(old)


@pytest.fixture(params=FUNCTIONAL_PLOT_BACKENDS)
def plot_backends(request):
    """Override the plot backend for this test (only for functional backends)

    This fixtures runs tests for all backends that actually produce plots, i.e.
    not for any type of dummy.

    Runs the test with the given plot backend and then restores the
    original value.

    Note that this does not reload the ui layers, which will
    retain a reference to the original backend throughout.
    """

    with TemporaryPlottingBackend(request.param):
        yield
        if request.param == 'pylab':
            from matplotlib import pyplot as plt
            plt.close(fig="all")


@pytest.fixture(autouse=True, scope="session")
def cleanup_pylab_backend():
     """Ensure that the pylab backend has closed down all windows.

     This is related to https://github.com/The-Compiler/pytest-xvfb/issues/11
     and the idea is to ensure that all matplotlib windows are closed at the
     end of the tests.
     """

     yield

     # Technically the system could have matplotlib installed but not
     # selected. However, this is not easily checked, so just check
     # if we can install matplotlib and, if so, close down any windows.
     #
     try:
         from matplotlib import pyplot as plt
     except ImportError:
         return

     plt.close(fig="all")


@pytest.fixture(autouse=True, scope="session")
def cleanup_ds9_backend():
    """Ensure that the DS9 window has closed down.

    The DS9 window can be left after the tests are run. This fixture
    attempts to close down any DS9 instance that has been started.

    """

    yield

    # If DS9 is not available this will lead to a warning message,
    # but there's no way to stop this. The code could import ds9_backend
    # instead of backend but then it would have to deal with failure,
    # and it does not seem an improvement.
    #
    from sherpa.image import backend
    backend.close()


@pytest.fixture(autouse=True)
def add_sherpa_test_data_dir(doctest_namespace):
    '''Define `data_dir` for doctests

    We make this an autouse=True fixture, because that means that
    the variable `data_dir` is available in every file that we may
    want to test without any special markup to that file. There is a
    risk that this is too magical and could confuse developers
    in the future, but it seems more important that this keeps any
    extra markup out of those files and keep them clean.
    '''
    doctest_namespace["data_3c273"] = 'sherpa/astro/datastack/tests/data/'

    # We could error out here, but we only want the tests that
    # use it to fail.
    #
    path = get_datadir()
    if path is None:
        return None

    if not path.endswith('/'):
        path += '/'

    doctest_namespace["data_dir"] = path


# Try to come up with a regexp for a numeric literal.  This is
# informed by pytest-doctestplus but this is a lot easier as we only
# care about a single number.
#
eterm = r"(?:e[+-]?\d+)"
term = "|".join([fr"[+-]?\d+\.\d*{eterm}?",
                 fr"[+-]?\.\d+{eterm}?"
                 fr"[+-]?\d+{eterm}"])
PATTERN = re.compile(fr"(.*)(?<![\d+-])({term})")


class NumberChecker:
    """Allow a string to be compared, allowing for tolerances.

    It is a limited version of
    pytest_doctestplus.output_checker.OutputChecker: it only checks a
    single number value in the text. This could be updated to be
    closer to doctestplus if found to be useful.

    """

    def __init__(self, expected, rtol=1e-4, atol=None):
        lhs, value, rhs = self.split(expected)
        self.lhs = lhs
        self.number = value
        self.rhs = rhs

        self.rtol = rtol
        self.atol = atol

    def split(self, text):
        """Split a string on a number

        Parameters
        ----------
        text : str
            The string containing a single number.
        """

        match = PATTERN.match(text)
        if match is None:
            # Make sure that we note a case where the regexp has
            # failed to find a number.
            #
            raise ValueError(f"Error in test - no number found in '{text}'")

        lhs = match[1]
        number = float(match[2])
        rhs = text[match.end():]

        return lhs, number, rhs

    def check(self, got):
        """check that got matches the expected text"""
        lhs, value, rhs = self.split(got)
        assert lhs == self.lhs
        assert rhs == self.rhs
        assert value == pytest.approx(self.number,
                                      rel=self.rtol, abs=self.atol)


def check_str_fixture(out, expecteds):
    """Check that out, when split on a newline, matches expecteds.

    It would be nice to use use doctestplus - that is
    pytest_doctestplus.output_checker.OutputChecker - for the checking
    (to avoid the need to mark up expected text as a regexp) but this
    is problematic at the moment, so it is not done. For now we
    require users to explicitly mark up those lines which we want to
    check with numeric constraints by adding the text (taken from
    doctestplus)

        "# doctest: +FLOAT_CMP"

    to the end of each line that needs it.

    Parameters
    ----------
    out : str
        The output to check (each line is compared to the expecteds
        array).
    expecteds : list of str or Pattern
        The expected output, split on new-line. Each line is one of: a
        Pattern, a string ending in "# doctest: +FLOAT_CMP", or a
        string (which is just checked for equality).

    """

    toks = str(out).split("\n")
    for tok, expected in zip(toks, expecteds):
        # expected is one of:
        #  - a pattern
        #  - a string ending in #doctest: +FLOAT_CMP
        #  - a normal string
        #
        try:
            assert expected.match(tok), (expected.pattern, tok)
        except AttributeError:
            idx = expected.find("# doctest: +FLOAT_CMP")
            if idx > -1:
                pat = NumberChecker(expected[:idx].rstrip())
                pat.check(tok)
                continue

            assert tok == expected

    assert len(toks) == len(expecteds)


@pytest.fixture
def check_str():
    """Returns the check_str_fixture routine.

    See check_str_fixture for documentation. This is a bit of a hack
    to treat a fixture as a way of importing check_str_fixture.

    """

    return check_str_fixture


@pytest.fixture
def xsmodel():
    """fixture that returns an XS<name> model instance.

    The test needs to be marked @requires_xspec when using this
    fixture.

    """

    try:
        from sherpa.astro import xspec
    except ImportError:
        raise RuntimeError("Test needs the requires_xspec decorator")

    def func(name, mname=None):
        cls = getattr(xspec, f"XS{name}")
        if mname is None:
            return cls()

        return cls(mname)

    return func


# Fixtures that control access to the tests, based on the availability
# of external "features" (normally this is the presence of optional
# modules). These used to be decorators in sherpa.utils.testing.
#
@pytest.fixture
def requires_pylab():
    """Runs the test with the pylab plotting backend if available."""

    pylab_backend = pytest.importorskip("sherpa.plot.pylab_backend")
    plt = pytest.importorskip("matplotlib.pyplot")

    with TemporaryPlottingBackend(pylab_backend.PylabBackend()):
        yield

    plt.close(fig="all")
