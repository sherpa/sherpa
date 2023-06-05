#
#  Copyright (C) 2020, 2021, 2022, 2023
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

"""Support simulation-related routines in sherpa.astro.ui.utils.

At present it is *very* limited.
"""

import numpy as np

import pytest

from sherpa.utils.testing import requires_data, requires_fits
from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa.astro import ui
from sherpa.utils.err import DataErr, IOErr


@pytest.mark.parametrize("idval", [None, 1, "faked"])
def test_fake_pha_no_rmf(idval, clean_astro_ui):
    """Check we error out if RMF is None."""

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays(idval, channels, counts, ui.DataPHA)

    # RMF is checked first
    #
    idval = 1 if idval is None else idval
    emsg = f'An RMF has not been found or supplied for data set {idval}'
    with pytest.raises(DataErr, match=emsg):
        ui.fake_pha(idval, arf=None, rmf=None, exposure=1000.0)


@pytest.mark.parametrize("idval", [None, 1, "faked"])
def test_fake_pha_missing_rmf(idval, clean_astro_ui, tmp_path):
    """Check we error out if RMF is not valid."""

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays(idval, channels, counts, ui.DataPHA)

    rmf = tmp_path / 'rmf'
    with pytest.raises(IOErr,
                       match=f"file '{rmf}' not found"):
        ui.fake_pha(idval, None, str(rmf), 1000.0)


@pytest.mark.parametrize("idval", [None, 1, "faked"])
def test_fake_pha_missing_arf(idval, clean_astro_ui, tmp_path):
    """Check we error out if ARF is not valid."""

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays(idval, channels, counts, ui.DataPHA)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    arf = tmp_path / 'arf'

    with pytest.raises(IOErr,
                       match=f"file '{arf}' not found"):
        ui.fake_pha(idval, str(arf), rmf, 1000.0)


@pytest.mark.parametrize("idval", [None, 1, "faked"])
def test_fake_pha_incompatible_rmf(idval, clean_astro_ui):
    """Check we error out if RMF is wrong size."""

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays(idval, channels, counts, ui.DataPHA)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6, 1.8, 2.0])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    idval = 1 if idval is None else idval
    emsg = f"RMF 'delta-rmf' is incompatible with PHA dataset '{idval}'"
    with pytest.raises(DataErr, match=emsg):
        ui.fake_pha(idval, arf, rmf, 1000.0)


def identity(x):
    """Return the input data"""
    return x


@pytest.mark.parametrize("method,expected",
                         [(None, [13, 16, 16]),
                          (identity, [12, 12, 12])
                          ])
@pytest.mark.parametrize("idval", [None, 1, "faked"])
def test_fake_pha_exposure_is_none_data_is_none(method, expected, idval, clean_astro_ui, reset_seed):
    """What happens if the exposure is None in the fake_pha call and the DataPHA?"""

    np.random.seed(7239652)

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    # 0.5 keV bin widths mean that source evaluation is half the
    # c0 value for a const1d model.
    #
    egrid = np.arange(2, 6) * 0.5
    elo = egrid[:-1]
    ehi = egrid[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    ui.load_arrays(idval, channels, counts, ui.DataPHA)
    ui.set_rmf(idval, rmf)
    ui.set_arf(idval, arf)

    cpt = ui.create_model_component('const1d', 'cpt')
    cpt.c0 = 24
    ui.set_source(idval, cpt)

    ui.fake_pha(idval, arf=None, rmf=None, exposure=None, method=method)

    # The exposure time is ignored.
    d = ui.get_data(idval)
    assert d.counts == pytest.approx(expected)
    assert d.exposure is None


@pytest.mark.parametrize("method,expected",
                         [(None, [9, 12, 18]),
                          (identity, [12, 12, 12])
                          ])
@pytest.mark.parametrize("idval", [None, 1, "faked"])
def test_fake_pha_exposure_is_none_data_is_set(method, expected, idval, clean_astro_ui, reset_seed):
    """What happens if the exposure is None in the fake_pha call but is set in DataPHA?"""

    np.random.seed(4596)

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    egrid = np.arange(2, 6) * 0.5
    elo = egrid[:-1]
    ehi = egrid[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    ui.load_arrays(idval, channels, counts, ui.DataPHA)
    ui.set_exposure(idval, 5)
    ui.set_rmf(idval, rmf)
    ui.set_arf(idval, arf)

    cpt = ui.create_model_component('const1d', 'cpt')
    cpt.c0 = 24
    ui.set_source(idval, cpt)

    ui.fake_pha(idval, arf=None, rmf=None, exposure=None, method=method)

    # The model evaluation should now be texp * c0 * bin-width
    # but the exposure time is being over-written so it is
    # not included.
    #
    d = ui.get_data(idval)
    assert d.counts == pytest.approx(expected)
    assert d.exposure is None


@pytest.mark.parametrize("method,expected",
                         [(None, [215, 390, 425]),
                          (identity, [200, 400, 400])
                          ])
@pytest.mark.parametrize("idval", [None, 1, "faked"])
@pytest.mark.parametrize("has_bkg", [True, False])
def test_fake_pha_basic(method, expected, idval, has_bkg, clean_astro_ui, reset_seed):
    """No background.

    See also test_fake_pha_add_background

    For simplicity we use perfect responses.

    A background dataset can be added, but it should
    not be used in the simulation.
    """

    np.random.seed(20347)

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)
    bcounts = 100 * counts

    ui.load_arrays(idval, channels, counts, ui.DataPHA)
    ui.set_exposure(idval, 100)

    if has_bkg:
        bkg = ui.DataPHA('bkg', channels, bcounts,
                         exposure=200, backscal=0.4)
        ui.set_bkg(idval, bkg, bkg_id='faked-bkg')

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    mdl = ui.create_model_component('const1d', 'mdl')
    mdl.c0 = 2
    ui.set_source(idval, mdl)

    ui.fake_pha(idval, arf, rmf, 1000.0, method=method)

    faked = ui.get_data(idval)
    assert faked.exposure == pytest.approx(1000.0)
    assert (faked.channel == channels).all()

    assert faked.name == 'faked'
    assert faked.get_arf().name == 'test-arf'
    assert faked.get_rmf().name == 'delta-rmf'

    if has_bkg:
        assert faked.background_ids == ['faked-bkg']
        bkg = ui.get_bkg(idval, 'faked-bkg')
        assert bkg.name == 'bkg'
        assert bkg.counts == pytest.approx(bcounts)
        assert bkg.exposure == pytest.approx(200)

    else:
        assert faked.background_ids == []

    # For reference the predicted source signal is
    #    [200, 400, 400]
    #
    assert faked.counts == pytest.approx(expected)


@pytest.mark.parametrize("method,expected",
                         [(None, [198, 384, 409]),
                          (identity, [200, 400, 400])
                          ])
def test_fake_pha_basic_arfrmf_set_in_advance(method, expected, clean_astro_ui, reset_seed):
    """Similar to test_fake_pha_basic but instead of passing in
    the RMF, we set it before. The result should be the same, so we
    don't have to go through all the parameterization of that test.
    """

    np.random.seed(20348)

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)

    ui.load_arrays('id123', channels, counts, ui.DataPHA)
    ui.set_exposure('id123', 100)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf('id123', rmf)
    ui.set_arf('id123', arf)

    mdl = ui.create_model_component('const1d', 'mdl')
    mdl.c0 = 2
    ui.set_source('id123', mdl)

    ui.fake_pha('id123', None, None, 1000.0, method=method)

    faked = ui.get_data('id123')
    assert faked.exposure == pytest.approx(1000.0)
    assert (faked.channel == channels).all()

    assert faked.name == 'faked'
    assert faked.background_ids == []

    # For reference the predicted source signal is
    #    [200, 400, 400]
    #
    assert faked.counts == pytest.approx(expected)


@pytest.mark.parametrize("method,expected",
                         [(None, [321, 506, 504]),
                          (identity, [325, 525, 525])
                          ])
@pytest.mark.parametrize("idval", [None, 1, "faked"])
def test_fake_pha_add_background(method, expected, idval, clean_astro_ui, reset_seed):
    """Check we can add a background component.

    See also test_fake_pha_basic.

    For simplicity we use perfect responses.
    """

    np.random.seed(20349)

    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)
    bcounts = 100 * counts

    ui.load_arrays(idval, channels, counts, ui.DataPHA)
    ui.set_exposure(idval, 100)
    ui.set_backscal(idval, 0.1)

    bkg = ui.DataPHA('bkg', channels, bcounts,
                     exposure=200, backscal=0.4)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    mdl = ui.create_model_component('const1d', 'mdl')
    mdl.c0 = 2
    ui.set_source(idval, mdl)

    ui.fake_pha(idval, arf, rmf, 1000.0, bkg=bkg, method=method)

    faked = ui.get_data(idval)
    assert faked.exposure == pytest.approx(1000.0)
    assert (faked.channel == channels).all()

    # For reference the predicted source signal is
    #    [200, 400, 400]
    # and the background signal is
    #    [125, 125, 125]
    #
    assert faked.counts == pytest.approx(expected)


@pytest.mark.parametrize("method,expected",
                         [(None, [189, 382, 400]),
                          (identity, [200, 400, 400])
                          ])
@pytest.mark.parametrize("idval", [None, 1, "faked"])
def test_fake_pha_no_data(method, expected, idval, clean_astro_ui, reset_seed):
    """What happens if there is no data loaded at the id?
    """

    np.random.seed(21347)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    mdl = ui.create_model_component('const1d', 'mdl')
    mdl.c0 = 2
    ui.set_source(idval, mdl)

    ui.fake_pha(idval, arf, rmf, 1000.0, method=method)

    channels = np.arange(1, 4)
    counts = [1, 1, 1]

    faked = ui.get_data(idval)
    assert faked.exposure == pytest.approx(1000.0)
    assert (faked.channel == channels).all()

    assert faked.name == 'faked'
    assert faked.get_arf().name == 'test-arf'
    assert faked.get_rmf().name == 'delta-rmf'

    assert faked.background_ids == []

    # For reference the predicted source signal is
    #    [200, 400, 400]
    #
    assert faked.counts == pytest.approx(expected)


@requires_fits
@requires_data
def test_fake_pha_file(make_data_path, clean_astro_ui, reset_seed):
    '''Test fake_pha using real input file.

    Note that HEG orders -1 and +1 should really be treated spearately,
    but for this test we just need two files to load.
    '''

    np.random.seed(22347)

    ui.set_source("gauss1d.g1")
    g1 = ui.get_source()
    g1.pos = 3
    g1.FWHM = .5

    ui.fake_pha(None,
                make_data_path('3c120_heg_-1.arf.gz'),
                make_data_path('3c120_heg_-1.rmf.gz'),
                1000.)
    data = ui.get_data()

    # Even with noise, maximum should be close to 3 keV. With the
    # fixed seed this check could be made exact, but leave as is.
    #
    assert np.isclose(data.get_x()[np.argmax(data.counts)], 3., atol=.2)

    # The model is localised around 3 keV so we could check that most
    # of the data is zero, but that seems excessive. Unlike other tests
    # we can not just check data.counts.sum() equals a certain value
    # since the value depends on OS: linux gets 6845 but macOS gets 6842.
    #
    total = data.counts.sum()
    assert 6840 <= total <= 6850


@requires_fits
@requires_data
def test_fake_pha_file_as_list(make_data_path, clean_astro_ui, reset_seed):
    '''Test fake_pha using real input file. See #1722

    This is test_fake_pha_file with the responses given as a list
    as that was found to be problematic. The data should match
    test_fake_pha_file with the same random seed.

    '''

    np.random.seed(22347)

    ui.set_source("gauss1d.g1")
    g1 = ui.get_source()
    g1.pos = 3
    g1.FWHM = .5

    # TODO: this fails with a No instrument response found
    ui.fake_pha(None,
                [make_data_path('3c120_heg_-1.arf.gz')],
                [make_data_path('3c120_heg_-1.rmf.gz')],
                1000.)
    data = ui.get_data()

    # Even with noise, maximum should be close to 3 keV. With the
    # fixed seed this check could be made exact, but leave as is.
    #
    assert np.isclose(data.get_x()[np.argmax(data.counts)], 3., atol=.2)

    # The model is localised around 3 keV so we could check that most
    # of the data is zero, but that seems excessive. Unlike other tests
    # we can not just check data.counts.sum() equals a certain value
    # since the value depends on OS: linux gets 6845 but macOS gets 6842.
    #
    total = data.counts.sum()
    assert 6840 <= total <= 6850


@requires_fits
@requires_data
def test_fake_pha_multi_file(make_data_path, clean_astro_ui, reset_seed):
    '''Test fake_pha using multiple real input files.

    Note that HEG orders -1 and +1 should really be treated spearately,
    but for this test we just need two files to load.
    '''

    np.random.seed(22349)

    ui.set_source("gauss1d.g1")
    g1 = ui.get_source()
    g1.pos = 3
    g1.FWHM = .5

    ui.fake_pha(None,
                [make_data_path('3c120_heg_-1.arf.gz'),
                 make_data_path('3c120_heg_1.arf.gz')],
                [make_data_path('3c120_heg_-1.rmf.gz'),
                 make_data_path('3c120_heg_1.rmf.gz')],
                500.)
    data = ui.get_data()

    # Even with noise, maximum should be close to 3 keV. With the
    # fixed seed this check could be made exact, but leave as is.
    #
    assert np.isclose(data.get_x()[np.argmax(data.counts)], 3., atol=.2)

    # The model is localised around 3 keV so we could check that most
    # of the data is zero, but that seems excessive. Unlike other tests
    # we can not just check data.counts.sum() equals a certain value
    # since the value depends on OS: linux gets 7053 but macOS gets 7133.
    #
    total = data.counts.sum()
    assert 7050 <= total <= 7150


@pytest.mark.parametrize("method,expected",
                         [(None, [42, 137, 53]),
                          (identity, [50, 151.25, 50])  # TODO: should be [50, 150, 50]
                          ])
def test_fake_pha_background_model(method, expected, clean_astro_ui, reset_seed):
    """Check we can add a background component.

    See also test_fake_pha_basic.

    For simplicity we use perfect responses.
    """

    np.random.seed(27347)

    idval = 'qwerty'
    channels = np.arange(1, 4, dtype=np.int16)
    counts = np.ones(3, dtype=np.int16)
    bcounts = 100 * counts

    ui.load_arrays(idval, channels, counts, ui.DataPHA)
    ui.set_exposure(idval, 100)
    ui.set_backscal(idval, 0.1)

    bkg = ui.DataPHA('bkg', channels, bcounts,
                     exposure=200, backscal=0.4)

    ebins = np.asarray([1.1, 1.2, 1.4, 1.6])
    elo = ebins[:-1]
    ehi = ebins[1:]
    arf = ui.create_arf(elo, ehi)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)

    mdl = ui.create_model_component('box1d', 'mdl')
    mdl.xlow = 1
    mdl.xhi = 1.4
    mdl.ampl = 0.5
    bkgmdl = ui.create_model_component('box1d', 'bkgmdl')
    bkgmdl.xlow = 1.2
    bkgmdl.xhi = 2
    bkgmdl.ampl = 1
    ui.set_source(idval, mdl)
    ui.set_bkg(idval, bkg)
    ui.set_bkg_source(idval, bkgmdl)
    ui.set_arf(idval, arf, bkg_id=1)
    ui.set_rmf(idval, rmf, bkg_id=1)

    ui.fake_pha(idval, arf, rmf, 1000.0, bkg='model', method=method)

    faked = ui.get_data(idval)
    assert faked.exposure == pytest.approx(1000.0)
    assert (faked.channel == channels).all()

    # For reference the predicted source signal is
    #    1000 * [0.5, 0.5, 0] * [0.1, 0.2, 0.2] = [50, 100, 0]
    # and the background signal for the source is, before backscal
    #    1000 * [0, 1, 1]     * [0.1, 0.2, 0.2] = [0, 200, 200]
    # and after source_backscal / bkg_backscal
    #    [0, 200, 200] * 0.1 / 0.4              = [0, 50, 50]
    #
    # So the predicted signal is [50, 150, 50]
    #
    # The idea here is that the observed background data should
    # not contribute.
    #
    # The current predicted signal is actually
    #
    #    [50, 150, 50] + [0, 1.25, 0]
    #
    # where the extra term is because the background model has
    # been evaluated again, this time with no response, to calculate
    # [0, 1, 0] and then scaled by 100 * (100 / 200) * (0.1 / 0.4).
    #
    assert faked.counts == pytest.approx(expected)


@requires_fits
@requires_data
def test_fake_pha_issue_1209(make_data_path, clean_astro_ui, tmp_path):
    """Check issue #1209.

    See also sherpa/astro/tests/test_fake_pha.py for

        test_fake_pha_has_valid_ogip_keywords_all_fake
        test_fake_pha_has_valid_ogip_keywords_from_real

    The session fake_pha includes quite a lot of logic which
    makes that the test case for #1209 should be done at this
    level, to complement the tests mentioned above.

    """

    infile = make_data_path("acisf01575_001N001_r0085_pha3.fits.gz")
    ui.load_pha(infile)
    ui.set_source(ui.powlaw1d.pl)
    pl.gamma = 1.8
    pl.ampl = 1e-4

    arf = ui.get_arf()
    rmf = ui.get_rmf()

    # check the TOTCTS setting in the input file
    d1 = ui.get_data()
    assert d1.header["TOTCTS"] == 855
    assert d1.counts.sum() == 855

    assert not d1.grouped
    assert not d1.subtracted

    ui.set_source("newid", pl)
    ui.fake_pha("newid",
                exposure=ui.get_exposure(),
                bkg=ui.get_bkg(),
                rmf=rmf,
                arf=arf,
                backscal=ui.get_backscal())
    stat = ui.calc_stat("newid")

    outfile = tmp_path / "sim.pha"
    ui.save_pha("newid", str(outfile))

    ui.load_pha(3, str(outfile))
    d3 = ui.get_data(3)
    assert isinstance(d3, ui.DataPHA)

    assert d3.exposure == pytest.approx(37664.157219191)
    assert d3.areascal == pytest.approx(1.0)
    assert d3.backscal == pytest.approx(2.2426552620567e-06)

    assert d3.background_ids == []
    assert d3.response_ids == []

    assert not d3.grouped
    assert not d3.subtracted

    # check the header
    hdr = d3.header
    assert hdr["TELESCOP"] == "CHANDRA"
    assert hdr["INSTRUME"] == "ACIS"
    assert hdr["FILTER"] == "none"

    # check some other values related to #1209 and #488 (OGIP)
    #
    assert "TOTCTS" not in hdr
    assert hdr["GROUPING"] == 0
    assert hdr["QUALITY"] == 0
    assert hdr["SYS_ERR"] == 0

    # We should get the same value - the responses are not written
    # to the temporary directory and so we need to load them
    # directly.
    #
    ui.set_rmf(3, rmf)
    ui.set_arf(3, arf)
    ui.set_source(3, pl)
    assert ui.calc_stat(3) == pytest.approx(stat)


@requires_fits
@requires_data
def test_fake_pha_issue_1568(make_data_path, clean_astro_ui, reset_seed):
    """Check issue #1568.

    In some cases, in particular XMM/RGS we only have an RMF, but no ARF.
    Make sure faking still works.
    """
    infile = make_data_path("xmmrgs/P0112880201R1S004RSPMAT1003.FTZ")
    np.random.seed(22347)

    ui.set_source("gauss1d.g1")
    g1 = ui.get_source()
    g1.pos = 2.
    g1.FWHM = .5

    ui.fake_pha(None, arf=None, rmf=infile, exposure=1000.)
    data = ui.get_data()

    # The maximum count value should be near 2 keV, thanks to the response,
    # and there's a reasonable number of counts. With this random seed
    # we can check the values.
    #
    assert data.counts.sum() == 7968
    assert data.get_x()[np.argmax(data.counts)] == pytest.approx(1.957132)


def check_analysis_settings(expected, bexpected=None,
                            mexpected=None, mbexpected=None,
                            swidth=None, bwidth=None):
    """Check plot results for type=counts then type=rates

    expected - source data
    bexpected - background data
    mexpected - source model data
    mbexpected - source background data

    This is for the tests that use setup_fake_pha_test.

    It checks that the type=rate plots are related to the
    type=counts plots with the expected ratios (given by
    the swidth/bwidth arguments).

    """

    if bexpected is None:
        bexpected = [300, 300, 300]

    if mexpected is None:
        mexpected = [80, 80, 300]

    if mbexpected is None:
        mbexpected = [0, 0, 1200]

    if swidth is None:
        swidth = [20, 40, 100]

    if bwidth is None:
        bwidth = [100, 200, 500]

    expected = np.asarray(expected)
    bexpected = np.asarray(bexpected)
    mexpected = np.asarray(mexpected)
    mbexpected = np.asarray(mbexpected)

    swidth = np.asarray(swidth)
    bwidth = np.asarray(bwidth)

    # These values don't change.
    splot = [2, 1, 0]
    bsplot = [0, 0, 3]

    ui.set_analysis(1, "energy", type="counts")
    assert ui.get_data_plot().y == pytest.approx(expected)
    assert ui.get_bkg_plot().y == pytest.approx(bexpected)
    assert ui.get_model_plot().y == pytest.approx(mexpected)
    assert ui.get_bkg_model_plot().y == pytest.approx(mbexpected)

    assert ui.get_source_plot().y == pytest.approx(splot)
    assert ui.get_bkg_source_plot().y == pytest.approx(bsplot)

    # Correct by "exposure"
    expected = expected / swidth
    bexpected = bexpected / bwidth
    mexpected = mexpected / swidth
    mbexpected = mbexpected / bwidth

    ui.set_analysis(1, "energy", type="rate")
    assert ui.get_data_plot().y == pytest.approx(expected)
    assert ui.get_bkg_plot().y == pytest.approx(bexpected)
    assert ui.get_model_plot().y == pytest.approx(mexpected)
    assert ui.get_bkg_model_plot().y == pytest.approx(mbexpected)

    assert ui.get_source_plot().y == pytest.approx(splot)
    assert ui.get_bkg_source_plot().y == pytest.approx(bsplot)


def setup_fake_pha_test():
    """Set up for test_fake_pha_xxx

    The test calling this should use the reset_seed fixture as this
    sets the random seed to a known value.

    """

    # Fake up a 3-channel source and background dataset with different
    # exposure times, backscal, and ARFs (the RMF will be "perfect" to
    # make this simple).
    #
    #               Source  Bgnd
    #    exposure =   100    500
    #    backscal =   0.2    0.4
    #         ARF =     2    0.8
    #
    # The RMF does not have equal-width energy bins just to show that
    # this has an affect.
    #
    #     1.0 - 1.2, 1.2 - 1.6, 1.6 - 2.6 keV
    #
    egrid = np.asarray([1, 1.2, 1.6, 2.6])
    elo = egrid[:-1]
    ehi = egrid[1:]

    arf_src = create_arf(elo, ehi, [2, 2, 2])
    arf_bkg = create_arf(elo, ehi, [0.8, 0.8, 0.8])

    rmf_src = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    rmf_bkg = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

    ui.load_arrays(1, [1, 2, 3], [20, 20, 20], ui.DataPHA)
    ui.set_exposure(100)
    ui.set_backscal(0.2)

    ui.set_rmf(rmf_src)
    ui.set_arf(arf_src)

    ui.set_analysis("energy")

    # We do not have a "ui" way to create the background
    #
    bkg = ui.DataPHA("bkg", [1, 2, 3], [300, 300, 300], exposure=500, backscal=0.4)
    bkg.set_arf(arf_bkg)
    bkg.set_rmf(rmf_bkg)

    ui.set_bkg(bkg)

    ui.set_source(ui.box1d.sbox)
    ui.set_bkg_source(ui.box1d.bbox)

    sbox.xlow = 0
    sbox.xhi = 1.4
    sbox.ampl.set(2, max=2)

    bbox.xlow = 1.6
    bbox.xhi = 3
    bbox.ampl.set(3, max=3)

    assert ui.get_data().counts == pytest.approx([20, 20, 20])

    # Explaining get_model_plot return (for counts):
    #
    # - expected source signal, from sbox
    #
    #    texp * arf * model * bin-widths
    #    100 * [2, 2, 2] * [2, 1, 0] * [0.2, 0.4, 1.0]
    #
    #    = [80, 80, 0]
    #
    # - expected background signal, from bbox
    #
    #    texp * arf * model * bin-widths * src_backscal / bkg_backscal
    #    100 * [2, 2, 2] * [0, 0, 3] * [0.2, 0.4, 1.0] * 0.2 / 0.4
    #
    #    = [0, 0, 300]
    #
    # For get_bkg_model_plot we use the background exposure and ARF values:
    #
    #    500 * [0.8, 0.8, 0.8] * [0, 0, 3] * [0.2, 0.4, 1.0]
    #
    #    = [0, 0, 1200]
    #
    check_analysis_settings(expected=[20, 20, 20])

    np.random.seed(392347)
    return (arf_src, rmf_src)


@pytest.mark.parametrize("method,expected",
                         [(None, [85, 75, 297]),
                          (identity, [80, 80, 300])
                          ])
def test_fake_pha_without_bgnd(method, expected, clean_astro_ui, reset_seed):
    """Extend test_fake_pha_with_bgnd_model to check with no background"""

    arf, rmf = setup_fake_pha_test()
    ui.fake_pha(1, arf, rmf, 100, bkg=None, method=method)

    check_analysis_settings(expected=expected)


@pytest.mark.parametrize("method,expected",
                         [(None, [85, 76, 297]),
                          (identity, [80, 80.3, 300.3])  # TODO: should be [80, 80, 300]
                          ])
def test_fake_pha_with_bgnd_model(method, expected, clean_astro_ui, reset_seed):
    """Add a test found when investigating #1685

    At some point working on #1684 this included the background twice,
    so add the test in to ensure this never happens.

    """

    arf, rmf = setup_fake_pha_test()
    ui.fake_pha(1, arf, rmf, 100, bkg="model", method=method)

    # The simulated data should match the model plot - that is
    # [80, 80, 300], but the background model gets included twice,
    # with the second time without a response, so it ends up adding
    # signal to the second two bins from 0.1 * [0, 3, 3], which
    # accounts for the [0, 0.3, 0.3] addition used above for the
    # "expected" value when the identity function is used.
    #
    check_analysis_settings(expected=expected)


@pytest.mark.parametrize("method,expected",
                         [(None, [90, 89, 1196]),
                          (identity, [85, 87.5, 1202.5])
                          ])
def test_fake_pha_with_bgnd_data(method, expected, clean_astro_ui, reset_seed):
    """Extend test_fake_pha_with_bgnd_model to check with a dataset."""

    bexpected = [10, 15, 5]

    # Choose a different background to the actual background,
    # particularly with different exposure/backscal.
    #
    bkg = ui.DataPHA("fakey", [1, 2, 3], bexpected)
    bkg.exposure = 400
    bkg.backscal = 0.1

    arf, rmf = setup_fake_pha_test()
    ui.fake_pha(1, arf, rmf, 100, bkg=bkg, method=method)

    # Check what's changed and what hasn't
    #
    # The expected signal should be the sum of the expected model
    # and the scaled background counts. That is
    #
    #     [80, 80, 1200] +
    #     (100 / 400) * (0.2 / 0.1) * [10, 15, 5]
    #
    #   = [80, 80, 1200] + 0.5 * [10, 15, 5]
    #   = [85, 87.5, 1202.5]
    #
    # Note that this includes the background model as well as the
    # background data.
    #
    check_analysis_settings(expected=expected,
                            bexpected=bexpected,
                            mexpected=[80, 80, 1200],
                            mbexpected=[0, 0, 2400],
                            swidth=[20, 40, 100],
                            bwidth=[80, 160, 400])


@pytest.mark.parametrize("grouped", [False, True])
def test_fake_pha_grouping(grouped, clean_astro_ui):
    """What happens with grouped data?"""

    ui.dataspace1d(1, 5, dstype=ui.DataPHA)
    ui.set_exposure(10)

    grping = [1, -1, 1, -1, -1]
    ui.set_grouping(1, grping)
    ui.group(1)

    egrid = np.arange(1, 7) * 0.5
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf(1, rmf)

    assert ui.get_data().grouped

    mdl = ui.create_model_component("polynom1d", "mdl")
    ui.set_source(1, mdl)

    ui.fake_pha(1, arf=None, rmf=None, exposure=10, grouped=grouped)

    assert ui.get_grouping(1) == pytest.approx(grping)
    assert ui.get_data().grouped is grouped


def test_fake_pha_subtracted(clean_astro_ui):
    """What happens with subtracted data?"""

    ui.dataspace1d(1, 5, dstype=ui.DataPHA)
    ui.set_exposure(10)

    bkg = ui.DataPHA("b", [1, 2, 3, 4, 5], [0, 0, 0, 0, 0])
    bkg.exposure = 5
    ui.set_bkg(1, bkg)

    ui.subtract()

    egrid = np.arange(1, 7) * 0.5
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    ui.set_rmf(1, rmf)

    assert ui.get_data().subtracted

    mdl = ui.create_model_component("polynom1d", "mdl")
    ui.set_source(1, mdl)

    ui.fake_pha(1, arf=None, rmf=None, exposure=10)

    assert ui.get_data().subtracted
