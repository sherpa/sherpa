#
# Copyright (C) 2015  Smithsonian Astrophysical Observatory
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
import sys
import imp
from sherpa import plot
from sherpa.plot import dummy_backend
import logging

try:
    from mock import MagicMock, NonCallableMagicMock, patch
except ImportError:
    from unittest.mock import MagicMock, NonCallableMagicMock, patch

# Some tests require that we precicely set up a clean
# environment so that some module are not importable,
# These constants makes sure that we can revert the system
# modules back to a vanilla state
BLACKLIST = ['matplotlib.pyplot', 'pychips']
INIT_MODULES = [
    module for module in sys.modules.keys() if module not in BLACKLIST]


def clean_modules():
    """
    reset system modules to a clean slate
    :return:
    """
    for m in [x for x in sys.modules.keys() if x not in INIT_MODULES]:
        del(sys.modules[m])


@pytest.fixture(params=["chips", "pylab", "dummy"])
def patch_module(request):
    """
    Fixture for basic tests making sure that the correct backend is selected
    depending on the user's selection. This is self contained, and no modules are required
    to be preset for this fixture to set them up.

    However, this fixture relies on sherpa.plot.dummy_backend to provide a faithful
    stub of the plotting backend implementations.

    This fixture is parametrized so to select each known backend
    :param request
    :return: the tuple (backend_name, backend_under_test)
    """
    clean_modules()
    sys.modules['chips_backend'] = MagicMock(spec=plot.dummy_backend)
    sys.modules['sherpa.plot.pylab_backend'] = MagicMock(
        spec=plot.dummy_backend)
    sys.modules['pychips'] = MagicMock()
    sys.modules['matplotlib'] = MagicMock()
    sys.modules['matplotlib.pyplot'] = MagicMock()
    plot_backend = _load_plotter_backend(request.param)

    def fin():
        del sys.modules['chips_backend']
        del sys.modules['sherpa.plot.pylab_backend']
        del sys.modules['pychips']
        del sys.modules['matplotlib']
        del sys.modules['matplotlib.pyplot']

    request.addfinalizer(fin)
    return request.param, plot_backend


def _load_plotter_backend(name):
    """
    We need to make sure that we can get a freshly imported backend module,
    rather than relying on the ones in memory. This utility function
    makes sure modules are appropriately reloaded between tests.
    :return:
    """
    plot.plotter.name = name
    plot_backend = _import_from_dotted_path(
        "sherpa.astro.datastack.plot_backend")
    reload(plot_backend)
    return plot_backend


def _import_from_dotted_path(dotted_names, path=None):
    """
    utility function used by _load_plotter_backend to recursively use imp to load
    modules.

    Ex: import_from_dotted_path('foo.bar') -> from foo import bar; return bar
    :param dotted_names:
    :param path:
    :return:
    """
    next_module, remaining_names = dotted_names.split('.', 1)
    fp, pathname, description = imp.find_module(next_module, path)
    module = imp.load_module(next_module, fp, pathname, description)
    if hasattr(module, remaining_names):
        return getattr(module, remaining_names)
    if '.' not in remaining_names:
        return module
    return _import_from_dotted_path(remaining_names, path=module.__path__)


@pytest.fixture(params=["chips", "pylab"])
def non_existent(request):
    """
    Fixture for error handling tests making sure that a backup
    implementation is used when the known backends are not available.

    This fixture ensures that the environment correctly emulates the lack
    of backends, even when they are actually available.
    :param request:
    :return:
    """
    clean_modules()
    plot_backend = _load_plotter_backend(request.param)
    return request.param, plot_backend


class PyChipsSpec(object):

    """
    Spec for PyChips Mock replacement
    We need to do this because pychips is an optional
    dependency, so we cannot rely on the autospec mechanism
    provided by mock.
    """

    def stub(self, *args, **kwargs):
        pass

    print_window = stub
    set_plot_xlabel = stub
    set_plot_ylabel = stub
    set_plot_title = stub
    add_window = stub
    current_window = stub
    set_current_window = stub
    limits = stub
    get_plot_range = stub


@pytest.fixture
def pychips_fixture():
    """
    Fixture for tests exercising the PyChips API

    The fixture sets up and returns:
      - the chips backend implementation in
    datastack.plot_backend, making sure it can be imported even if chips
    is not present on the system (and mock it if it is available)
      - the pychips mock module
      - a dataset stub used by several tests served by this fixture.
    :return: the tuple (plot_chips_backend, pychips_mock, dataset_stub)
    """
    # Force selection of backend, overriding configuration
    backend_name = "chips"

    # Generate mock for pychips. This module may not be present,
    # so we need to mock it before we try to import it,
    # otherwise the import would fail.
    # Without this constraints we would still mock it,
    # but after the code requiring it has been loaded into memory,
    # and with a different pattern.
    # This object also serves as a stub in some tests
    pychips = MagicMock(spec=PyChipsSpec)
    pychips.AUTO = NonCallableMagicMock()
    pychips.X_AXIS = NonCallableMagicMock()
    pychips.Y_AXIS = NonCallableMagicMock()
    pychips.chips_plus = NonCallableMagicMock()

    # clear Python's memory from any cached modules that it might want
    # to reuse.
    clean_modules()

    # force loading our mock pychips implementation instead of the (possibly
    # missing) system one.
    sys.modules['datastack.plot_backend.pychips'] = pychips
    sys.modules['pychips'] = pychips

    # Now the environment is finally ready to load the actual backend.
    # The backend will now import the mock pychips module, whether the
    # module is present or not.
    plot_backend = _load_plotter_backend(backend_name)

    # Many tests served by this fixture require a stub of a dataset
    dataset = MagicMock()
    dataset.__getitem__.return_value = "TESTID"

    return plot_backend, pychips, dataset


def test_backend(patch_module):
    """
    test that the correct backend is imported by datastack.plot_backend,
    and ensure it's not empty. Only the presence of one function is exercised.

    :param patch_module: fixture required by this test
    :return:
    """
    backend_name, plot_backend = patch_module
    assert backend_name == plot_backend.name
    assert backend_name + "_backend" == plot_backend.backend_module.name
    assert plot_backend.initialize_backend is not None


def test_backend_import_error(non_existent):
    """
    test that dummy_backend is selected and imported when an error
    occurs importing an existing one.

    :param non_existent: the fixture required by this test
    :return:
    """
    backend_name, plot_backend = non_existent
    assert backend_name == plot_backend.name
    assert "dummy_backend" == plot_backend.backend_module.name
    assert plot_backend.initialize_backend is not None


def test_dummy_warnings(non_existent):
    """
    test that warnings are issued when the dummy backend is selected
    and its functions called.

    :param non_existent: the fixture required by this test
    :return:
    """
    plot_backend = _load_plotter_backend("dummy")
    logger = logging.getLogger("datastack.plot_backend.plot_dummy")
    with patch.object(logger, "warning") as warning:
        for function in [
            plot_backend.initialize_backend,
            plot_backend.initialize_plot,
            plot_backend.select_plot,
            plot_backend.save_plot,
            plot_backend.plot_savefig,
            plot_backend.plot_xlabel,
            plot_backend.plot_ylabel,
            plot_backend.plot_title,
            plot_backend.plot_xlim,
            plot_backend.plot_ylim,
            plot_backend.plot_set_xscale,
            plot_backend.plot_set_yscale
        ]:
            function()
            warning.assert_called_with("using dummy plotting backend")


def test_pychips_initialize(pychips_fixture):
    """
    test pychips API calls are sound: initialize
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mock Context
    # Since we can always assume that the backend implementation is importable,
    # as long as we mock its dependency on pychips (even if pychips is missing),
    # we can use the patch context manager to simply replace the chips_backend
    # *after it is loaded* by the backend_module. We rely on autospec to avoid
    # having to decribe the specification ourselves.
    # For these tests we are not interested that much in the backend implementation.
    # However, it is important that we ensure that the backend initialization
    # is performed (or performable, rather, as we will need to test this again
    # in higher level tests). Without initialization of the backend,
    # chips would not be launched.
    with patch("datastack.plot_backend.backend_module.chips_backend", autospec=True) as chips_backend:
        # Test Calls
        plot_backend.initialize_backend()

        # Behavior Test
        chips_backend.begin.assert_called_once_with()
        chips_backend.end.assert_called_once_with()


def test_pychips_new_window(pychips_fixture):
    """
    test pychips API calls are sound: new_window
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    dataset = pychips_fixture[2]

    # Test Calls
    plot_backend.initialize_plot(dataset, MagicMock())

    # Behavior Test
    pychips.add_window.assert_called_once_with(['id', 'TESTID'])
    pychips.current_window.assert_called_once_with('TESTID')


def test_pychips_select_plot(pychips_fixture):
    """
    test pychips API calls are sound: select_plot
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    dataset = pychips_fixture[2]

    # Test Calls
    plot_backend.initialize_plot(dataset, MagicMock())

    # Behavior Test
    pychips.add_window.assert_called_once_with(['id', 'TESTID'])
    pychips.current_window.assert_called_once_with('TESTID')


def test_pychips_save_plot(pychips_fixture):
    """
    test pychips API calls are sound: save_plot
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    dataset = pychips_fixture[2]

    # Test Calls
    options = {
        'clobber': True,
        'fittopage': True,
    }
    did = dataset["id"]
    filename = "test.png"
    plot_backend.save_plot(did, filename, options)

    # Behavior Test
    pychips.print_window.assert_called_once_with(did, filename, options)


def test_pychips_savefig(pychips_fixture):
    """
    test pychips API calls are sound: save_fig
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    dataset = pychips_fixture[2]

    # Test Calls
    options = {
        'clobber': True,
        'fittopage': True,
    }
    did = dataset["id"]
    filename = "test.png"
    plot_backend.save_plot(did, filename, options)

    # Behavior Test
    pychips.print_window.assert_called_once_with(did, filename, options)


def test_pychips_plot_xlabel(pychips_fixture):
    """
    test pychips API calls are sound: plot_xlabel
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    dataset = pychips_fixture[2]

    # Test Calls
    did = dataset["id"]
    label = "XTEXT"
    plot_backend.plot_xlabel(did, label)

    # Behavior Test
    pychips.set_plot_xlabel.assert_called_once_with(did, label)


def test_pychips_plot_ylabel(pychips_fixture):
    """
    test pychips API calls are sound: plot_ylabel
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    dataset = pychips_fixture[2]

    # Test Calls
    did = dataset["id"]
    label = "YTEXT"
    plot_backend.plot_ylabel(did, label)

    # Behavior Test
    pychips.set_plot_ylabel.assert_called_once_with(did, label)


def test_pychips_plot_title(pychips_fixture):
    """
    test pychips API calls are sound: plot_title
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    dataset = pychips_fixture[2]

    # Test Calls
    did = dataset["id"]
    label = "TITLE"
    plot_backend.plot_title(did, label)

    # Behavior Test
    pychips.set_plot_title.assert_called_once_with(did, label)


def test_pychips_plot_xlim_none(pychips_fixture):
    """
    test pychips API calls are sound: plot_xlim
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    pychips.get_plot_range.return_value = [0, 1, 2, 3]

    # Test Calls
    xlimits = plot_backend.plot_xlim()

    # Behavior Test
    pychips.get_plot_range.assert_called_once_with()

    # State Test
    assert [0, 1] == xlimits


def test_pychips_plot_xlim_min(pychips_fixture):
    """
    test pychips API calls are sound: plot_xlim
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    pychips.get_plot_range.return_value = [-10, 1, 2, 3]

    # Test Calls
    xlimits = plot_backend.plot_xlim(xmin=-10)

    # Behavior Test
    pychips.limits.assert_called_once_with(pychips.X_AXIS, xmin=-10, xmax=1)

    # State Test
    assert [-10, 1] == xlimits


def test_pychips_plot_xlim_max(pychips_fixture):
    """
    test pychips API calls are sound: plot_xlim
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    pychips.get_plot_range.return_value = [0, 10, 2, 3]

    # Test Calls
    xlimits = plot_backend.plot_xlim(xmax=10)

    # Behavior Test
    pychips.limits.assert_called_once_with(pychips.X_AXIS, xmin=0, xmax=10)

    # State Test
    assert [0, 10] == xlimits


def test_pychips_plot_xlim_both(pychips_fixture):
    """
    test pychips API calls are sound: plot_xlim
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    pychips.get_plot_range.return_value = [-10, 10, 2, 3]

    # Test Calls
    xlimits = plot_backend.plot_xlim(xmin=-10, xmax=10)

    # Behavior Test
    pychips.limits.assert_called_once_with(pychips.X_AXIS, xmin=-10, xmax=10)

    # State Test
    assert [-10, 10] == xlimits


def test_pychips_plot_ylim_none(pychips_fixture):
    """
    test pychips API calls are sound: plot_ylim
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    pychips.get_plot_range.return_value = [0, 1, 2, 3]

    # Test Calls
    ylimits = plot_backend.plot_ylim()

    # Behavior Test
    pychips.get_plot_range.assert_called_once_with()

    # State Test
    assert [2, 3] == ylimits


def test_pychips_plot_ylim_min(pychips_fixture):
    """
    test pychips API calls are sound: plot_ylim
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    pychips.get_plot_range.return_value = [0, 1, -20, 3]

    # Test Calls
    ylimits = plot_backend.plot_ylim(ymin=-20)

    # Behavior Test
    pychips.limits.assert_called_once_with(pychips.Y_AXIS, ymin=-20, ymax=3)

    # State Test
    assert [-20, 3] == ylimits


def test_pychips_plot_ylim_max(pychips_fixture):
    """
    test pychips API calls are sound: plot_ylim
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    pychips.get_plot_range.return_value = [0, 1, 2, 30]

    # Test Calls
    ylimits = plot_backend.plot_ylim(ymax=30)

    # Behavior Test
    pychips.limits.assert_called_once_with(pychips.Y_AXIS, ymin=2, ymax=30)

    # State Test
    assert [2, 30] == ylimits


def test_pychips_plot_ylim_both(pychips_fixture):
    """
    test pychips API calls are sound: plot_ylim
    :param pychips_fixture: the fixture required by this test
    :return:
    """
    # Unit under test
    plot_backend = pychips_fixture[0]

    # Mocks
    pychips = pychips_fixture[1]

    # Stubs
    pychips.get_plot_range.return_value = [0, 1, -20, 30]

    # Test Calls
    ylimits = plot_backend.plot_ylim(ymin=-20, ymax=30)

    # Behavior Test
    pychips.limits.assert_called_once_with(pychips.Y_AXIS, ymin=-20, ymax=30)

    # State Test
    assert [-20, 30] == ylimits
