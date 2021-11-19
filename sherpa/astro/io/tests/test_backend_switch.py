#
#  Copyright (C) 2021
#  MIT
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
from unittest.mock import patch

import pytest

from sherpa.utils.testing import requires_data, requires_fits
from sherpa.astro import io
from sherpa.astro import ui


@requires_fits
@requires_data
def test_backend_switch(make_data_path):
    '''Test that we can switch backends.

    Of course, switching only makes sense if there is more than one backend
    to try, so this is listed as "requires_fits".
    '''
    pha = make_data_path("3c273.pi")
    ui.clean()
    ui.load_pha(pha)
    assert ui.get_data().name.endswith('3c273.pi')

    orig_backend = io.backend

    import sherpa.astro.io.dummy_backend
    io.backend = sherpa.astro.io.dummy_backend

    with pytest.raises(NotImplementedError,
                       match="No usable I/O backend was imported."):
         ui.load_pha(pha)

    io.backend = orig_backend
    infile = make_data_path("9774.pi")
    ui.load_pha(infile)
    assert ui.get_data().name.endswith('9774.pi')


def test_io_load_config_invalid_io_pkg(tmp_path):
    """What happens when io_pkg is invalid?

    Technically not a backend switch, but fits in with
    the "reload" feature.
    """

    conf = tmp_path / 'test-config.rc'
    conf.write_text("""[options]
io_pkg : not_a_module
""", encoding='ascii')

    with patch("sherpa.get_config", lambda: str(conf)):
        mod = importlib.reload(io)

    assert mod.backend is None

    # Restore the import
    mod = importlib.reload(io)


def test_io_load_config_ogip_min_energy_negative(tmp_path):
    """What happens when ogip.minimum_energy is < 0?

    Technically not a backend switch, but fits in with
    the "reload" feature.
    """

    conf = tmp_path / 'test-config.rc'
    conf.write_text("""[options]
io_pkg : dummy

[ogip]
minimum_energy = -2
""", encoding='ascii')

    with patch("sherpa.get_config", lambda: str(conf)):
        with pytest.raises(ValueError) as ve:
            importlib.reload(io)

    assert str(ve.value) == 'Invalid value for [ogip] minimum_energy config value; it must be None or a float > 0'

    # Restore the import
    importlib.reload(io)


def test_io_load_config_ogip_min_energy_invalid(tmp_path):
    """What happens when ogip.minimum_energy is not numeric?

    Technically not a backend switch, but fits in with
    the "reload" feature.
    """

    conf = tmp_path / 'test-config.rc'
    conf.write_text("""[options]
io_pkg : dummy

[ogip]
minimum_energy = foo
""", encoding='ascii')

    with patch("sherpa.get_config", lambda: str(conf)):
        with pytest.raises(ValueError) as ve:
            importlib.reload(io)

    assert str(ve.value) == 'Invalid value for [ogip] minimum_energy config value; it must be None or a float > 0'

    # Restore the import
    importlib.reload(io)


def test_io_load_config_ogip_min_energy_none(tmp_path):
    """What happens when ogip.minimum_energy is none?

    Technically not a backend switch, but fits in with
    the "reload" feature.
    """

    conf = tmp_path / 'test-config.rc'
    conf.write_text("""[options]
io_pkg : dummy

[ogip]
minimum_energy = NONE
""", encoding='ascii')

    orig = io.ogip_emin
    assert orig is not None

    with patch("sherpa.get_config", lambda: str(conf)):
        mod = importlib.reload(io)

    assert mod.ogip_emin is None

    # Restore the import
    mod = importlib.reload(io)
    assert mod.ogip_emin == orig


def test_io_load_config_ogip_min_energy_set(tmp_path):
    """What happens when ogip.minimum_energy is valid?

    Technically not a backend switch, but fits in with
    the "reload" feature.
    """

    conf = tmp_path / 'test-config.rc'
    conf.write_text("""[options]
io_pkg : dummy

[ogip]
minimum_energy = 0.2
""", encoding='ascii')

    orig = io.ogip_emin
    assert orig is not None

    with patch("sherpa.get_config", lambda: str(conf)):
        mod = importlib.reload(io)

    assert mod.ogip_emin == pytest.approx(0.2)

    # Restore the import
    mod = importlib.reload(io)
    assert mod.ogip_emin == orig


def test_io_load_config_ogip_min_energy_unset(tmp_path):
    """What happens when ogip.minimum_energy is not set?

    Technically not a backend switch, but fits in with
    the "reload" feature.
    """

    conf = tmp_path / 'test-config.rc'
    conf.write_text("""[options]
io_pkg : dummy
""", encoding='ascii')

    orig = io.ogip_emin
    assert orig is not None

    with patch("sherpa.get_config", lambda: str(conf)):
        mod = importlib.reload(io)

    # the default is 1e-10, which should match orig
    assert mod.ogip_emin == pytest.approx(1e-10)
    assert orig == pytest.approx(1e-10)

    # Restore the import
    mod = importlib.reload(io)
