#
#  Copyright (C) 2020, 2021, 2023
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

"""Test the OO layer fake_pha command.

At present it is *very* limited.
"""

import numpy as np

import pytest

from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa.astro.data import DataPHA
from sherpa.astro.fake import fake_pha
from sherpa.astro import io
from sherpa.models import Box1D, Const1D
from sherpa.utils.err import ArgumentTypeErr, DataErr
from sherpa.utils.testing import requires_data, requires_fits


channels = np.arange(1, 4, dtype=np.int16)
counts = np.ones(3, dtype=np.int16)
bcounts = 1000 * counts

ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
elo = ebins[:-1]
ehi = ebins[1:]
arf = create_arf(elo, ehi)
rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)


def test_fake_pha_requires_pha():
    """Check what happens when the first argument is not a PHA.

    This is a regression test, in that the error may not be
    directly caught (i.e. from an explicit error check) but
    from using an object incorrectly.
    """

    pha = DataPHA("x", [1, 2, 3], [0, 0, 0])
    mdl = Const1D()
    with pytest.raises(AttributeError,
                       match="^'Const1D' object has no attribute 'response_ids'$"):
        fake_pha(mdl, pha)


def test_fake_pha_requires_model():
    """Check what happens when the second argument is not a model.

    This is a regression test, in that the error may not be
    directly caught (i.e. from an explicit error check) but
    from using an object incorrectly.
    """

    # We require that the PHA has a response
    pha = DataPHA("x", channels, channels * 0)
    pha.set_rmf(rmf)
    pha.set_arf(arf)
    with pytest.raises(AttributeError,
                       match="^'list' object has no attribute 'name'$"):
        fake_pha(pha, [1, 2, 3])


def test_fake_pha_requires_callable():
    """Check what happens when method is not a callable."""

    # We require that the PHA has a response
    pha = DataPHA("x", channels, channels * 0)
    pha.set_rmf(rmf)
    pha.set_arf(arf)
    with pytest.raises(ArgumentTypeErr,
                       match="^'method' must be a callable$"):
        fake_pha(pha, Const1D(), method=True)


def test_fake_pha_requires_response():
    """Check we error out when there's no response"""

    pha = DataPHA("x", [1, 2, 3], [0, 0, 0])
    mdl = Const1D()
    with pytest.raises(DataErr,
                       match="^An RMF has not been found or supplied for data set x$"):
        fake_pha(pha, mdl)


def identity(x):
    """Return the input data"""
    return x


@pytest.mark.parametrize("method,is_source,expected1,expected2",
                         [(None, True, [184, 423, 427], [428, 773, 800]),
                          (identity, True, [200, 400, 400], [400, 800, 800]),
                          (None, False, [1, 1, 0], [2, 2, 0]),
                          (identity, False, [2, 2, 2], [2, 2, 2])
                          ])
@pytest.mark.parametrize("has_bkg", [True, False])
def test_fake_pha_basic(method, is_source, expected1, expected2, has_bkg, reset_seed):
    """No background.

    For simplicity we use perfect responses.

    A background dataset can be added, but it should
    not be used in the simulation with default settings
    """
    np.random.seed(4276)
    data = DataPHA("any", channels, counts, exposure=1000.)

    if has_bkg:
        bkg = DataPHA("bkg", channels, bcounts,
                      exposure=2000, backscal=0.4)
        data.set_background(bkg, id="unused-bkg")

    data.set_arf(arf)
    data.set_rmf(rmf)

    mdl = Const1D("mdl")
    mdl.c0 = 2

    fake_pha(data, mdl, is_source=is_source, add_bkgs=False, method=method)

    assert data.exposure == pytest.approx(1000.0)
    assert (data.channel == channels).all()

    assert data.name == "any"
    assert data.get_arf().name == "user-arf"
    assert data.get_rmf().name == "delta-rmf"

    if has_bkg:
        assert data.background_ids == ["unused-bkg"]
        bkg = data.get_background("unused-bkg")
        assert bkg.name == "bkg"
        assert bkg.counts == pytest.approx(bcounts)
        assert bkg.exposure == pytest.approx(2000)

    else:
        assert data.background_ids == []

    # For reference the predicted source signal is
    #    [200, 400, 400]
    #
    # When is_source is not set then the values are significantly
    # smaller as there's no response, in particular no exposure term
    # or changes due to bin-widths varying with position.
    #
    assert data.counts == pytest.approx(expected1)

    # Essentially double the exposure by having two identical arfs
    #
    data.set_arf(arf, 2)
    data.set_rmf(rmf, 2)
    fake_pha(data, mdl, is_source=is_source, add_bkgs=False, method=method)

    assert data.counts == pytest.approx(expected2)


@pytest.mark.parametrize("method,expected1,expected2",
                         [(None, [186, 197, 212], [9, 4, 6]),  # expected2 is wrong
                          pytest.param(identity, [200, 200, 200], [202 / 6, 202 / 6, 202 / 6], marks=pytest.mark.xfail)
                          ])
def test_fake_pha_background_pha(method, expected1, expected2, reset_seed):
    """Sample from background pha (no background model)"""
    np.random.seed(1234)

    data = DataPHA("any", channels, counts, exposure=1000.)
    bkg = DataPHA("bkg", channels, bcounts, exposure=2000, backscal=2.5)
    data.set_background(bkg, id="used-bkg")

    data.set_arf(arf)
    data.set_rmf(rmf)

    mdl = Const1D("mdl")
    mdl.c0 = 0

    # Just make sure that the model does not contribute
    fake_pha(data, mdl, is_source=True, add_bkgs=False)
    assert data.counts == pytest.approx([0, 0, 0])

    # The background signal (as each bin is the same) is
    #
    #   background = 1000 counts
    #   exposure_background / exposure_source = 2
    #   backscal_background / backscal_source = 2.5
    #
    # and there's no other scaling (no areascal), so the
    # predicted background count in the source is
    #
    #   1000 / 2 / 2.5 = 1000 / 5 = 200
    #
    fake_pha(data, mdl, is_source=True, add_bkgs=True, method=method)
    assert data.counts == pytest.approx(expected1)

    # Add several more backgrounds with, compared to the first
    # background,
    #
    #   - same backscal
    #   - half the exposure time (so same as the source)
    #   - much smaller counts (1 rather than 1000)
    #
    # For each of the new background components the background
    # contribution to the source is
    #
    #    1 / 1 / 2.5 = 1 / 2.5 = 0.4
    #
    # So the full contribution is
    #
    #    200 + 0.4 * 5 = 202
    #
    # and the average of this is 202 / 6 = 33 2/3.
    #
    for i in range(5):
        bkg = DataPHA("bkg", channels, counts,
                      exposure=1000, backscal=2.5)
        data.set_background(bkg, id=i)

    fake_pha(data, mdl, is_source=True, add_bkgs=True, method=method)
    assert data.counts == pytest.approx(expected2)


@pytest.mark.parametrize("method,expected1,expected2",
                         [(None, [186, 411, 405], [197, 396, 389]),
                          (identity, [200, 400, 400], [200, 400, 400])
                          ])
def test_fake_pha_bkg_model(method, expected1, expected2, reset_seed):
    """Test background model
    """

    np.random.seed(5329853)

    data = DataPHA("any", channels, counts, exposure=1000.)

    bkg = DataPHA("bkg", channels, bcounts,
                  exposure=2000, backscal=1.)
    data.set_background(bkg, id="used-bkg")

    data.set_arf(arf)
    data.set_rmf(rmf)

    bkg.set_arf(arf)
    bkg.set_rmf(rmf)

    mdl = Const1D("mdl")
    mdl.c0 = 0

    bmdl = Const1D("bmdl")
    bmdl.c0 = 2

    # With no background model the simulated source counts
    # are 0.
    #
    bmodels = {"used-bkg": bmdl}
    fake_pha(data, mdl, is_source=True, add_bkgs=False,
             bkg_models=bmodels, method=method)

    assert data.counts == pytest.approx([0, 0, 0])

    # Check we have created source counts this time.
    #
    fake_pha(data, mdl, is_source=True, add_bkgs=True,
             bkg_models=bmodels, method=method)

    assert data.exposure == pytest.approx(1000.0)
    assert (data.channel == channels).all()

    assert data.name == "any"
    assert data.get_arf().name == "user-arf"
    assert data.get_rmf().name == "delta-rmf"

    # The background itself is unchanged by the simulation.
    assert data.background_ids == ["used-bkg"]
    bkg = data.get_background("used-bkg")
    assert bkg.name == "bkg"
    assert bkg.counts == pytest.approx(bcounts)
    assert bkg.exposure == pytest.approx(2000)

    # The expected number of counts here comes from the background
    # model, and not the background dataset. This means that, unlike
    # test_fake_pha_background_pha where the background is flat, we
    # have the potential for a different expected count per bin. The
    # ARF and RMF are "flat", as is the model (when plotted per keV),
    # *but* the bin widths are not constant in energy, so when the
    # model is integrated across each bin you get
    #
    #    bin 0: 1.1 - 1.2 keV, width = 0.1 keV
    #        1  1.2 - 1.4              0.2
    #        2  1.4 - 1.6              0.2
    #
    # So the signal in the first bin will be half the other two.
    #
    # With a model of 2 photon/cm^2/s/keV and a perfect ARF/RMF
    # the predicted signal is, for a source exposure time of 1000s,
    #
    #    2 * 1000 * [0.1, 0.2, 0.2] = [200, 400, 400]
    #
    assert data.counts == pytest.approx(expected1)

    # Now add a second set of arf/rmf for the data. However, all the
    # signal is background, so this does not change any of the
    # results.
    #
    data.set_arf(arf, 2)
    data.set_rmf(rmf, 2)
    fake_pha(data, mdl, is_source=True, add_bkgs=True,
             bkg_models=bmodels, method=method)

    assert data.counts == pytest.approx(expected2)


@pytest.mark.parametrize("method,expected1,expected2",
                         [(None, [388, 777, 755], [479, 957, 785]),  # expected2 values look wrong
                          pytest.param(identity, [400, 800, 800], [400 + 1000/3, 800 + 1000/3, 800 - 400/3], marks=pytest.mark.xfail)
                          ])
def test_fake_pha_bkg_model_multiple(method, expected1, expected2, reset_seed):
    """Test background model with multiple backgounds

    This is to check the different scaling factors, in particular
    the exposure to background exposure differences.

    This test takes advantage of the method argument to allow
    the actual model values to be checked as well as a random
    realisation using poisson_noise).
    """

    np.random.seed(873)

    data = DataPHA("any", channels, counts, exposure=1000.)

    bkg1 = DataPHA("bkg1", channels, bcounts,
                   exposure=2000, backscal=1.)
    bkg2 = DataPHA("bkg2", channels, bcounts,
                   exposure=6000, backscal=2.)
    bkg3 = DataPHA("bkg3", channels, bcounts,
                   exposure=500, backscal=4.)

    data.set_arf(arf)
    data.set_rmf(rmf)

    for i, bkg in enumerate([bkg1, bkg2, bkg3], 1):
        data.set_background(bkg, id=i)
        bkg.set_arf(arf)
        bkg.set_rmf(rmf)

    # Add a flat source component that is not 0 unlike the
    # other background checks.
    #
    mdl = Const1D("mdl")
    mdl.c0 = 4

    # Have the background models depend on energy to make it
    # obvious how the scaling works.
    #
    # The energy bins are 1.1-1.2, 1.2-1.4, and 1.4-1.6 keV.
    #
    bmdl1 = Box1D("bmdl1")
    bmdl1.xlow = 1.1
    bmdl1.xhi = 1.2
    bmdl1.ampl.set(10, max=10)

    bmdl2 = Box1D("bmdl2")
    bmdl2.xlow = 1.2
    bmdl2.xhi = 1.4
    bmdl2.ampl.set(10, max=10)

    bmdl3 = Box1D("bmdl3")
    bmdl3.xlow = 1.4
    bmdl3.xhi = 1.6
    bmdl3.ampl.set(-8, min=-8)

    bmodels = {1: bmdl1, 2: bmdl2, 3: bmdl3}

    # Source-model only.
    #
    # The predicted counts are
    #
    #   1000 * 4 * [0.1, 0.2, 0.2] = [400, 800, 800]
    #
    fake_pha(data, mdl, is_source=True, add_bkgs=False,
             bkg_models=bmodels, method=method)

    assert data.exposure == pytest.approx(1000.0)
    assert data.counts == pytest.approx(expected1)

    # Check we create the expected signal, which is all from the
    # backgrounds and the source.
    #
    fake_pha(data, mdl, is_source=True, add_bkgs=True,
             bkg_models=bmodels, method=method)

    # The background identifiers haven't changed.
    #
    assert data.background_ids == [1, 2, 3]

    # Check the exposure times have not changed.
    #
    assert data.exposure == pytest.approx(1000.0)
    assert data.get_background(1).exposure == pytest.approx(2000)
    assert data.get_background(2).exposure == pytest.approx(6000)
    assert data.get_background(3).exposure == pytest.approx(500)

    # The signal is a sum of the source and weighted background
    # components. The source exposure time should be used, not the
    # background, so all we care about is the ratio of the backscal
    # values (as the areascal is None in all cases here).
    #
    # In the following we show
    #
    #     exposure * model * bin-width / (backscal relative to source)
    #
    #   Source:
    #     1000 * 4 * [0.1, 0.2, 0.2] / 1           = [400, 800, 800]
    #
    #   Bkg 1:
    #     1000 * [10, 0, 0] * [0.1, 0.2, 0.2] / 1  = [1000, 0, 0]
    #
    #   Bkg 2:
    #     1000 * [0, 10, 0] * [0.1, 0.2, 0.2] / 2   = [0, 1000, 0]
    #
    #   Bkg 1:
    #     1000 * [0, 0, -8] * [0.1, 0.2, 0.2] / 4  = [0, 0, -400]
    #
    # So the expected total is, because we average the background
    # contributions:
    #
    #     [400, 800, 800] + [1000, 1000, -400] / 3
    #   = [733.333, 1133.333, 666.66]
    #
    assert data.counts == pytest.approx(expected2)


@pytest.mark.parametrize("method,expected",
                         [(None, [433, 793, 806]),
                          (identity, [400, 800, 800])
                          ])
@pytest.mark.parametrize("add_bkgs", [False, True])
def test_fake_pha_with_no_background(method, expected, add_bkgs, reset_seed):
    """Check add_bkgs keyword does nothing if no background."""

    np.random.seed(3567)

    # For fun we only set a RMF, not a ARF.
    data = DataPHA("x", channels, counts, exposure=1000)
    data.set_rmf(rmf)
    data.set_analysis("energy")

    mdl = Const1D()
    mdl.c0 = 4

    assert data.exposure == pytest.approx(1000.0)
    assert data.channel == pytest.approx(channels)
    assert data.counts == pytest.approx(counts)

    fake_pha(data, mdl, add_bkgs=add_bkgs, method=method)

    assert data.exposure == pytest.approx(1000.0)
    assert data.channel == pytest.approx(channels)
    assert data.counts == pytest.approx(expected)


@requires_fits
def test_fake_pha_has_valid_ogip_keywords_all_fake(tmp_path, reset_seed):
    """See #1209

    When everything is faked, what happens?
    """

    np.random.seed(5)

    data = DataPHA("any", channels, counts, exposure=1000.)

    bkg = DataPHA("bkg", channels, bcounts,
                  exposure=2000, backscal=1.)
    data.set_background(bkg, id="used-bkg")

    data.set_arf(arf)
    data.set_rmf(rmf)

    bkg.set_arf(arf)
    bkg.set_rmf(rmf)

    mdl = Const1D("mdl")
    mdl.c0 = 0

    bmdl = Const1D("bmdl")
    bmdl.c0 = 2

    fake_pha(data, mdl, is_source=True, add_bkgs=True,
             bkg_models={"used-bkg": bmdl})

    outfile = tmp_path / "sim.pha"
    io.write_pha(str(outfile), data, ascii=False)

    inpha = io.read_pha(str(outfile))
    assert inpha.channel == pytest.approx(channels)

    # it is not required that we check counts (that is, we can drop this
    # if it turns out not to be repeatable across platforms), but for
    # now keep the check.
    #
    assert inpha.counts == pytest.approx([188, 399, 416])

    for field in ["staterror", "syserror", "bin_lo", "bin_hi",
                  "grouping", "quality"]:
        assert getattr(inpha, field) is None

    assert inpha.exposure == pytest.approx(1000.0)
    assert inpha.backscal == pytest.approx(1.0)
    assert inpha.areascal == pytest.approx(1.0)
    assert not inpha.grouped
    assert not inpha.subtracted
    assert inpha.response_ids == []
    assert inpha.background_ids == []

    hdr = inpha.header
    assert hdr["TELESCOP"] == "none"
    assert hdr["INSTRUME"] == "none"
    assert hdr["FILTER"] == "none"

    for key in ["EXPOSURE", "AREASCAL", "BACKSCAL",
                "ANCRFILE", "BACKFILE", "RESPFILE"]:
        assert key not in hdr


@requires_fits
@requires_data
def test_fake_pha_has_valid_ogip_keywords_from_real(make_data_path, tmp_path, reset_seed):
    """See #1209

    In this version we use a "real" PHA file as the base.

    See sherpa/astro/ui/tests/test_astro_ui_utils_simulation.py

        test_fake_pha_issue_1209

    which is closer to the reported case in #1209
    """

    np.random.seed(5)

    infile = make_data_path("acisf01575_001N001_r0085_pha3.fits.gz")
    data = io.read_pha(infile)

    mdl = Const1D("mdl")
    mdl.c0 = 0

    bmdl = Const1D("bmdl")
    bmdl.c0 = 2

    fake_pha(data, mdl, is_source=True, add_bkgs=True,
             bkg_models={"used-bkg": bmdl})

    outfile = tmp_path / "sim.pha"
    io.write_pha(str(outfile), data, ascii=False)

    inpha = io.read_pha(str(outfile))
    assert inpha.channel == pytest.approx(np.arange(1, 1025))

    # it is not required that we check counts (that is, we can drop this
    # if it turns out not to be repeatable across platforms), but for
    # now keep the check.
    #
    expected = np.zeros(1024)
    for idx in [70, 92, 104, 138, 171, 457, 623, 670, 685, 706, 776,
                792, 860, 894, 950, 1019]:
        expected[idx] = 1

    expected[1023] = 3
    assert inpha.counts == pytest.approx(expected)

    for field in ["staterror", "syserror", "bin_lo", "bin_hi",
                  "grouping", "quality"]:
        assert getattr(inpha, field) is None

    assert inpha.exposure == pytest.approx(37664.157219191)
    assert inpha.backscal == pytest.approx(2.2426552620567e-06)
    assert inpha.areascal == pytest.approx(1.0)
    assert not inpha.grouped
    assert not inpha.subtracted
    assert inpha.response_ids == []
    assert inpha.background_ids == []

    hdr = inpha.header
    assert hdr["TELESCOP"] == "CHANDRA"
    assert hdr["INSTRUME"] == "ACIS"
    assert hdr["FILTER"] == "none"

    for key in ["EXPOSURE", "AREASCAL", "BACKSCAL",
                "ANCRFILE", "BACKFILE", "RESPFILE"]:
        assert key not in hdr
