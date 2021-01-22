#
#  Copyright (C) 2020, 2021  Smithsonian Astrophysical Observatory
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

"""Support simulation-related routines in sherpa.astro.ui.utils.

At present it is *very* limited.
"""

import numpy as np

import pytest

from sherpa.astro import ui
from sherpa.utils.err import DataErr,IOErr


def xfail(arg):
    return pytest.param(arg, marks=pytest.mark.xfail)


@pytest.mark.parametrize("id", [None, 1, "faked"])
def test_fake_pha_no_rmf(id, clean_astro_ui):
    """Check we error out if RMF is None."""

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays(id, channels, counts, ui.DataPHA)

    # RMF is checked first
    #
    with pytest.raises(DataErr) as exc:
        ui.fake_pha(id, arf=None, rmf=None, exposure=1000.0)

    emsg = f'An RMF has not been found or supplied for data set {id}'
    assert str(exc.value) == emsg


@pytest.mark.parametrize("id", [None, 1, "faked"])
def test_fake_pha_missing_rmf(id, clean_astro_ui, tmp_path):
    """Check we error out if RMF is not valid."""

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays(id, channels, counts, ui.DataPHA)

    rmf = tmp_path / 'rmf'
    with pytest.raises(IOErr) as exc:
        ui.fake_pha(id, None, str(rmf), 1000.0)

    assert str(exc.value) == f"file '{rmf}' not found"


@pytest.mark.parametrize("id", [None, 1, "faked"])
def test_fake_pha_missing_arf(id, clean_astro_ui, tmp_path):
    """Check we error out if ARF is not valid."""

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays(id, channels, counts, ui.DataPHA)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    arf = tmp_path / 'arf'

    with pytest.raises(IOErr) as exc:
        ui.fake_pha(id, str(arf), rmf, 1000.0)

    assert str(exc.value) == f"file '{arf}' not found"


@pytest.mark.parametrize("id", [xfail(None), 1, "faked"])
def test_fake_pha_incompatible_rmf(id, clean_astro_ui):
    """Check we error out if RMF is wrong size."""

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays(id, channels, counts, ui.DataPHA)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6, 1.8, 2.0])
    elo = ebins[:-1]
    ehi = ebins[1:]
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    with pytest.raises(DataErr) as exc:
        ui.fake_pha(id, None, rmf, 1000.0)

    id = 1 if id is None else id
    emsg = f"RMF 'delta-rmf' is incompatible with PHA dataset '{id}'"
    assert str(exc.value) == emsg


@pytest.mark.parametrize("id", [None, 1, "faked"])
@pytest.mark.parametrize("has_bkg", [True, False])
def test_fake_pha_basic(id, has_bkg, clean_astro_ui):
    """No background.

    See also test_fake_pha_add_background

    For simplicity we use perfect responses.

    A background dataset can be added, but it should
    not be used in the simulation.
    """

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)
    bcounts = 100 * counts

    ui.load_arrays(id, channels, counts, ui.DataPHA)
    ui.set_exposure(id, 100)

    if has_bkg:
        bkg = ui.DataPHA('bkg', channels, bcounts,
                         exposure=200, backscal=0.4)
        ui.set_bkg(id, bkg, bkg_id='faked-bkg')

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    mdl = ui.create_model_component('const1d', 'mdl')
    mdl.c0 = 2
    ui.set_source(id, mdl)

    ui.fake_pha(id, arf, rmf, 1000.0)

    faked = ui.get_data(id)
    assert faked.exposure == pytest.approx(1000.0)
    assert (faked.channel == channels).all()

    assert faked.name == 'faked'
    assert faked.get_arf().name == 'test-arf'
    assert faked.get_rmf().name == 'delta-rmf'

    if has_bkg and id is not None:
        assert faked.background_ids == ['faked-bkg']
        bkg = ui.get_bkg(id, 'faked-bkg')
        assert bkg.name == 'bkg'
        assert bkg.counts == pytest.approx(bcounts)
        assert bkg.exposure == pytest.approx(200)

    else:
        assert faked.background_ids == []

    # check we've faked counts (the scaling is such that it is
    # very improbable that this condition will fail)
    assert (faked.counts > counts).all()

    # For reference the predicted source signal is
    #    [200, 400, 400]
    #
    # What we'd like to say is that the predicted counts are
    # similar, but this is not easy to do. What we can try
    # is summing the counts (to average over the randomness)
    # and then a simple check
    #
    assert faked.counts.sum() > 200


@pytest.mark.parametrize("id", [None, 1, "faked"])
def test_fake_pha_add_background(id, clean_astro_ui):
    """Check we can add a background component.

    See also test_fake_pha_basic.

    For simplicity we use perfect responses.
    """

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)
    bcounts = 100 * counts

    ui.load_arrays(id, channels, counts, ui.DataPHA)
    ui.set_exposure(id, 100)
    ui.set_backscal(id, 0.1)

    bkg = ui.DataPHA('bkg', channels, bcounts,
                     exposure=200, backscal=0.4)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    mdl = ui.create_model_component('const1d', 'mdl')
    mdl.c0 = 2
    ui.set_source(id, mdl)

    ui.fake_pha(id, arf, rmf, 1000.0, bkg=bkg)

    faked = ui.get_data(id)
    assert faked.exposure == pytest.approx(1000.0)
    assert (faked.channel == channels).all()

    # check we've faked counts (the scaling is such that it is
    # very improbable that this condition will fail)
    assert (faked.counts > counts).all()

    # For reference the predicted source signal is
    #    [200, 400, 400]
    # and the background signal is
    #    [125, 125, 125]
    # so, even with randomly drawn values, the following
    # checks should be robust.
    #
    predicted_by_source = 1000 * mdl(elo, ehi)
    predicted_by_bkg = (1000/200) * (0.1/0.4) * bcounts
    assert (faked.counts > predicted_by_source).all()
    assert (faked.counts > predicted_by_bkg).all()


@pytest.mark.parametrize("id", [None, 1, "faked"])
def test_fake_pha_no_data(id, clean_astro_ui):
    """What happens if there is no data loaded at the id?
    """

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    mdl = ui.create_model_component('const1d', 'mdl')
    mdl.c0 = 2
    ui.set_source(id, mdl)

    ui.fake_pha(id, arf, rmf, 1000.0)

    # We don't really check anything sensible with the counts.
    # It is unlikely the simulated counts will be <= 1.
    #
    # For reference the predicted source signal is
    #    [200, 400, 400]
    #
    channels = np.arange(1, 4)
    counts = [1, 1, 1]

    faked = ui.get_data(id)
    assert faked.exposure == pytest.approx(1000.0)
    assert (faked.channel == channels).all()

    assert faked.name == 'faked'
    assert faked.get_arf().name == 'test-arf'
    assert faked.get_rmf().name == 'delta-rmf'

    assert faked.background_ids == []

    # check we've faked counts (the scaling is such that it is
    # very improbable that this condition will fail)
    assert (faked.counts > counts).all()

    # What we'd like to say is that the predicted counts are
    # similar, but this is not easy to do. What we can try
    # is summing the counts (to average over the randomness)
    # and then a simple check
    #
    assert (faked.counts.sum() > 200) and (faked.counts.sum() < 10000)
