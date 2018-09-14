#
#  Copyright (C) 2016, 2018  Smithsonian Astrophysical Observatory
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
from sherpa.utils.testing import requires_data, requires_fits
from six.moves import reload_module


def assert_chips_called(chips):
    """
    Just make sure the mock is being called once. At this point we are not too
    concerned about the correctness of the calls. We are assuming the code works.

    Parameters
    ----------
    chips : the mock object backing the tests
    """
    chips.lock.assert_any_call()
    chips.unlock.assert_any_call()


@pytest.fixture
def chips(mock_chips):
    return mock_chips[1]


@pytest.fixture(autouse=True)
def ui(mock_chips):
    """
    Load and return the sherpa.astro.ui module
    """
    from sherpa.astro import ui

    reload_module(ui)
    return ui


@pytest.fixture(autouse=True)
def clean_up(request, ui):
    def fin():
        ui.clean()
    request.addfinalizer(fin)


@pytest.fixture(autouse=True)
def load_data(ui, make_data_path):
    """
    Load dataset before every test.
    """
    ui.load_data(make_data_path("3c273.pi"))
    ui.set_source("powlaw1d.p")
    ui.set_bkg_model("const1d.c")
    ui.fit()


@requires_data
@requires_fits
@pytest.mark.parametrize("call, args", [
    ("model", ()),
    ("arf", ()),
    ("source", ()),
    ("bkg", ()),
    ("bkg_model", ()),
    ("bkg_resid", ()),
    ("bkg_ratio", ()),
    ("bkg_delchi", ()),
    ("bkg_chisqr", ()),
    ("bkg_fit", ()),
    ("bkg_source", ()),
    ("energy_flux", ()),
    ("photon_flux", ()),
    ("bkg_fit_resid", ()),
    ("bkg_fit_delchi", ()),
    ("source_component", (1, "p")),
    ("model_component", (1, "p")),
    ("order", (1,)),
])
def test_plot(call, args, chips, ui):
    function = getattr(ui, "plot_" + call)
    function(*args)
    assert_chips_called(chips)
