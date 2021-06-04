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

"""Test the OO layer fake_pha command.

At present it is *very* limited.
"""

import numpy as np

import pytest

from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa.astro.data import DataPHA
from sherpa.astro.fake import fake_pha
from sherpa.models import Const1D

channels = np.arange(1, 4, dtype=np.int16)
counts = np.ones(3, dtype=np.int16)
bcounts = 1000 * counts

ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
elo = ebins[:-1]
ehi = ebins[1:]
arf = create_arf(elo, ehi)
rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)


@pytest.mark.parametrize("has_bkg", [True, False])
@pytest.mark.parametrize("is_source", [True, False])
def test_fake_pha_basic(has_bkg, is_source):
    """No background.

    See also test_fake_pha_add_background

    For simplicity we use perfect responses.

    A background dataset can be added, but it should
    not be used in the simulation woth default settings
    """
    data = DataPHA('any', channels, counts, exposure=1000.)

    if has_bkg:
        bkg = DataPHA('bkg', channels, bcounts,
                      exposure=2000, backscal=0.4)
        data.set_background(bkg, id='unused-bkg')

    data.set_arf(arf)
    data.set_rmf(rmf)

    mdl = Const1D('mdl')
    mdl.c0 = 2

    fake_pha(data, mdl, is_source=is_source, add_bkgs=False)

    assert data.exposure == pytest.approx(1000.0)
    assert (data.channel == channels).all()

    assert data.name == 'any'
    assert data.get_arf().name == 'user-arf'
    assert data.get_rmf().name == 'delta-rmf'

    if has_bkg:
        assert data.background_ids == ['unused-bkg']
        bkg = data.get_background('unused-bkg')
        assert bkg.name == 'bkg'
        assert bkg.counts == pytest.approx(bcounts)
        assert bkg.exposure == pytest.approx(2000)

    else:
        assert data.background_ids == []

    if is_source:
        # check we've faked counts (the scaling is such that it is
        # very improbable that this condition will fail)
        assert (data.counts > counts).all()

        # For reference the predicted source signal is
        #    [200, 400, 400]
        #
        # What we'd like to say is that the predicted counts are
        # similar, but this is not easy to do. What we can try
        # is summing the counts (to average over the randomness)
        # and then a simple check
        #
        assert data.counts.sum() > 500
        assert data.counts.sum() < 1500
        # This is more likely to fail by chance, but still very unlikely
        assert data.counts[1] > data.counts[0]
    else:
        # No multiplication with exposure time, arf binning, etc.
        # so we just expect very few counts
        assert data.counts.sum() < 20

    # Essentially double the exposure by having two identical arfs
    data.set_arf(arf, 2)
    data.set_rmf(rmf, 2)
    fake_pha(data, mdl, is_source=is_source, add_bkgs=False)
    if is_source:
        assert data.counts.sum() > 1200
        assert data.counts.sum() < 3000
        assert data.counts[1] > data.counts[0]
    else:
        assert data.counts.sum() < 50
        assert data.counts.sum() > 1

def test_fake_pha_background_pha():
    '''Sample from background pha'''
    data = DataPHA('any', channels, counts, exposure=1000.)
    bkg = DataPHA('bkg', channels, bcounts, exposure=2000, backscal=2.5)
    data.set_background(bkg, id='used-bkg')

    data.set_arf(arf)
    data.set_rmf(rmf)

    mdl = Const1D('mdl')
    mdl.c0 = 0
    # Just make sure that the model does not contribute
    fake_pha(data, mdl, is_source=True, add_bkgs=False)
    assert data.counts.sum() == 0

    fake_pha(data, mdl, is_source=True, add_bkgs=True)
    # expected is [200, 200, 200]
    assert data.counts.sum() > 400
    assert data.counts.sum() < 1000

    # Add several more backgrounds. Actual background should be average.
    # We add 5 background with half the exposure time as this first oned
    # and essentially 0 counts. So, we should find 1/11 of the counts
    # we found in the last run.
    for i in range(5):
        bkg = DataPHA('bkg', channels, np.ones(3, dtype=np.int16),
                      exposure=1000, backscal=2.5)
        data.set_background(bkg, id=i)

    fake_pha(data, mdl, is_source=True, add_bkgs=True)
    # expected is about [18, 18, 18]
    assert data.counts.sum() > 10
    assert data.counts.sum() < 200


def test_fake_pha_bkg_model():
    """Test background model
    """
    data = DataPHA('any', channels, counts, exposure=1000.)

    bkg = DataPHA('bkg', channels, bcounts,
                  exposure=2000, backscal=1.)
    data.set_background(bkg, id='used-bkg')

    data.set_arf(arf)
    data.set_rmf(rmf)

    bkg.set_arf(arf)
    bkg.set_rmf(rmf)

    mdl = Const1D('mdl')
    mdl.c0 = 0

    bmdl = Const1D('bmdl')
    bmdl.c0 = 2

    fake_pha(data, mdl, is_source=True, add_bkgs=True,
             bkg_models={'used-bkg': bmdl})

    assert data.exposure == pytest.approx(1000.0)
    assert (data.channel == channels).all()

    assert data.name == 'any'
    assert data.get_arf().name == 'user-arf'
    assert data.get_rmf().name == 'delta-rmf'

    # The background itself is unchanged
    assert data.background_ids == ['used-bkg']
    bkg = data.get_background('used-bkg')
    assert bkg.name == 'bkg'
    assert bkg.counts == pytest.approx(bcounts)
    assert bkg.exposure == pytest.approx(2000)

    # check we've faked counts (the scaling is such that it is
    # very improbable that this condition will fail)
    assert (data.counts > counts).all()

    # For reference the predicted signal is
    #    [200, 400, 400]
    # but, unlike in the test above, this time it's all coming
    # from the background.
    #
    # What we'd like to say is that the predicted counts are
    # similar, but this is not easy to do. What we can try
    # is summing the counts (to average over the randomness)
    # and then a simple check
    #
    assert data.counts.sum() > 500
    assert data.counts.sum() < 1500
    # This is more likely to fail by chance, but still very unlikely
    assert data.counts[1] > 1.5 * data.counts[0]

    # Now add a second set of arf/rmf for the data.
    # However, all the signal is background, so this does not change
    # any of the results.
    data.set_arf(arf, 2)
    data.set_rmf(rmf, 2)
    fake_pha(data, mdl, is_source=True, add_bkgs=True,
             bkg_models={'used-bkg': bmdl})
    assert data.counts.sum() > 500
    assert data.counts.sum() < 1500
    assert data.counts[1] > 1.5 * data.counts[0]
