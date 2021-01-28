#
#  Copyright (C) 2020, 2021
#     Smithsonian Astrophysical Observatory
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

"""
Flux-related tests of the sherpa.astro.ui module.
"""

import numpy as np

import pytest

from sherpa.astro import ui
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_plotting, requires_xspec
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, FitErr, IdentifierErr, IOErr, \
    ModelErr
import sherpa.astro.utils


def fail(*arg):
    return pytest.param(*arg, marks=pytest.mark.xfail)


# TODO: use wavelength and channels analysis

@pytest.mark.parametrize("id", [None, 1, "foo"])
@pytest.mark.parametrize("func", [ui.calc_photon_flux, ui.calc_energy_flux])
def test_calc_flux_pha_invalid_range(id, func, clean_astro_ui):
    """Ensure an error is raised if lo > hi"""

    x = np.arange(3, 6)
    y = np.ones(x.size - 1)

    ui.load_arrays(id, x[:-1], x[1:], y, ui.Data1DInt)
    mdl = ui.create_model_component('const1d', 'm')

    if id is None:
        ui.set_source(mdl)
    else:
        ui.set_source(id, mdl)

    # Note: really the error message should not include energy since in
    # this case (Data1DInt) there's no energy, and if this were
    # a DataPHA case the message says energy even if analysis=wave
    # or channel.
    #
    emsg = 'the energy range is not consistent, 12 !< 5'
    with pytest.raises(IOErr, match=emsg):
        if id is None:
            func(12, 5)
        else:
            func(12, 5, id=id)


@requires_data
@requires_fits
@pytest.mark.parametrize("func", [ui.calc_photon_flux, ui.calc_energy_flux])
def test_calc_flux_pha_invalid_model(func, make_data_path, clean_astro_ui):
    """Don't allow strings for model parameter"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)
    ui.set_source('powerlaw.pl')

    emsg = "'model' must be a model object"
    with pytest.raises(ArgumentTypeErr, match=emsg):
        func(0.5, 7, model='pl')


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.calc_energy_flux,
                                    ui.calc_photon_flux])
def test_calc_foo_flux_no_bkg(method, make_data_path, clean_astro_ui):
    """No background model.
    """

    setup_sample(1, make_data_path, fit=False)
    with pytest.raises(ModelErr) as exc:
        method(lo=0.5, hi=7, bkg_id=1)

    assert str(exc.value) == 'background model 1 for data set 1 has not been set'


def test_calc_flux_pha_bin_edges(clean_astro_ui):
    """What happens when filter edges partially overlap bins?

    Later tests may also cover this condition, but here we use
    faked data that is made to make the behavior "obvious".
    """

    chans = np.arange(1, 11, 1, dtype=np.int)
    counts = np.zeros(chans.size, dtype=np.int)

    # "perfect" response
    energies = np.arange(1, 12, 1)
    elo, ehi = energies[:-1], energies[1:]
    flat = np.ones(chans.size, dtype=np.int)

    d = ui.DataPHA('example', chans, counts)
    arf = ui.create_arf(elo, ehi, flat)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=elo, startchan=1,
                        fname=None)

    d.set_arf(arf)
    d.set_rmf(rmf)
    ui.set_data(1, d)

    ui.set_source(ui.powlaw1d.pl)
    pl.ampl = 1e-4
    pl.gamma = 1.7

    # Evaluate the model on the energy grid
    ymdl = pl(elo, ehi)

    pflux = ui.calc_photon_flux(2.6, 7.8)
    eflux = ui.calc_energy_flux(2.6, 7.8)

    enscale = sherpa.astro.utils._charge_e

    # Left in as notes:
    #
    # This are the true values. Given that the edge bins (should be)
    # using linear interpolation, how close do we get to this?
    #
    # gterm1 = 1.0 - 1.7
    # gterm2 = 2.0 - 1.7
    # true_pflux = 1e-4 * (7.8**gterm1 - 2.6**gterm1) / gterm1
    # true_eflux = enscale * 1e-4 * (7.8**gterm2 - 2.6**gterm2) / gterm2
    #
    # Comparing the linear interpolation scheme to that used by
    # Sherpa prior to fixing #619, I find
    #
    # Photon fluxes:
    #     linear_interp / true_pflux = 1.042251
    #     sherpa        / true_pflux = 0.837522
    #
    # Energy fluxes:
    #     linear_interp / true_eflux = 1.017872
    #     sherpa        / true_eflux = 0.920759
    #
    scale = np.asarray([0.0, 0.4, 1.0, 1.0, 1.0, 1.0, 0.8, 0, 0.0, 0.0])
    expected_pflux = (ymdl * scale).sum()

    emid = enscale * (elo + ehi) / 2
    expected_eflux = (emid * ymdl * scale).sum()

    assert pflux == pytest.approx(expected_pflux)

    # check against log as values ~ 3e-13
    eflux = np.log10(eflux)
    expected_eflux = np.log10(expected_eflux)
    assert eflux == pytest.approx(expected_eflux)


def test_calc_flux_pha_density_bin_edges(clean_astro_ui):
    """What happens when filter edges partially overlap bins? flux density

    Later tests may also cover this condition, but here we use
    faked data that is made to make the behavior "obvious".
    """

    chans = np.arange(1, 11, 1, dtype=np.int)
    counts = np.zeros(chans.size, dtype=np.int)

    # "perfect" response
    energies = np.arange(1, 12, 1)
    elo, ehi = energies[:-1], energies[1:]
    flat = np.ones(chans.size, dtype=np.int)

    d = ui.DataPHA('example', chans, counts)
    arf = ui.create_arf(elo, ehi, flat)
    rmf = ui.create_rmf(elo, ehi, e_min=elo, e_max=elo, startchan=1,
                        fname=None)

    d.set_arf(arf)
    d.set_rmf(rmf)
    ui.set_data(1, d)

    ui.set_source(ui.powlaw1d.pl)
    pl.ampl = 1e-4
    pl.gamma = 1.7

    # choose an energy that is not equal to the center of the bin
    # just to check how this is handled
    #
    pdens = ui.calc_photon_flux(2.6)
    edens = ui.calc_energy_flux(2.6)

    enscale = sherpa.astro.utils._charge_e

    # Evaluate the model over the bin 2-3 keV; since the grid
    # has a width of 1 keV we do not need to divide by the bin
    # width when calculating the density.
    #
    ymdl = pl([2], [3])
    expected_pdens = ymdl.sum()
    expected_edens = enscale * 2.5 * expected_pdens

    # Prior to fixing #619, Sherpa returns 0 for both densities
    #
    assert pdens == pytest.approx(expected_pdens)

    # check against log as values ~ 5e-13
    edens = np.log10(edens)
    expected_edens = np.log10(expected_edens)
    assert edens == pytest.approx(expected_edens)


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "foo"])
@pytest.mark.parametrize("lo, hi", [(None, None),
                                    (0.1, 11),
                                    (0.1, 7),
                                    (0.5, 11),
                                    (0.5, 7),
                                    (0.05, 7),
                                    (0.5, 15),
                                    (0.05, 15),
                                    (0.01, 0.05),
                                    (12, 14)])
def test_calc_flux_pha(id, lo, hi, make_data_path, clean_astro_ui):
    """Do calc_photon/energy_flux return the expected results: fluxes?

    This skips those combinations where only one of lo or hi
    is None, since this is handle by test_calc_flux_density_pha.

    Flues are in units of <value>/cm^2/s  {value=photon, erg}

    The checks are made against ranges that are chosen to cover
    matching the grid, a subset of the grid (with and without
    matching the start/end), partial overlaps, and no overlap.
    """

    infile = make_data_path('3c273.pi')

    # The 3c273 RMF was generated over the range 0.1 to 11 keV
    # (inclusive). By picking gamma = 1 the photon-flux
    # integral is just ampl * (log(ehi) - log(elo))
    # and the energy-flux ampl is ampl * (ehi - elo) * scale
    # where scale converts 1 keV to erg.
    #
    ampl = 1e-4

    if lo is None:
        loval = 0.1
    elif lo < 0.1:
        loval = 0.1
    else:
        loval = lo

    if hi is None:
        hival = 11.0
    elif hi > 11.0:
        hival = 11.0
    else:
        hival = hi

    # expected fluxes; special case the handling of there being no
    # overlap between the user grid and the data grid.
    #
    if lo is not None and (lo > 11.0 or hi < 0.1):
        pflux_exp = 0.0
        eflux_exp = 0.0
    else:
        pflux_exp = ampl * (np.log(hival) - np.log(loval))
        eflux_exp = 1.602e-9 * ampl * (hival - loval)

    pl = ui.create_model_component('powlaw1d', 'pl')
    pl.ampl = ampl
    pl.gamma = 1

    if id is None:
        ui.load_pha(infile)
        ui.set_source(pl)
    else:
        ui.load_pha(id, infile)
        ui.set_source(id, pl)

    # Use a subset of the data range (to check that the calc routines
    # ignores them, ie uses the full 0.1 to 11 keV range.
    #
    ui.ignore(None, 0.5)
    ui.ignore(7, None)

    # Do not use named arguments, but assume positional arguments
    if id is None:
        pflux = ui.calc_photon_flux(lo, hi)
        eflux = ui.calc_energy_flux(lo, hi)
    else:
        pflux = ui.calc_photon_flux(lo, hi, id)
        eflux = ui.calc_energy_flux(lo, hi, id)

    # Since the energy fluxes are ~1e-12 we want to rescale the
    # value before comparison. Here we use a log transform.
    #
    eflux_exp = np.log10(eflux_exp)
    eflux = np.log10(eflux)

    assert pflux == pytest.approx(pflux_exp, rel=1e-3)
    assert eflux == pytest.approx(eflux_exp, rel=1e-4)


# The lo/hi range which match in the different settings; using _hc
# to convert from keV to Angstroms is a bit low-level.
#
# The energy-to-channel conversion was done by filtering in energy
# space and asking Sherpa what channels this selected.
#
@requires_data
@requires_fits
@pytest.mark.parametrize("elo, ehi, setting, lo, hi",
                         [(None, None, 'wave', None, None),
                          (None, None, 'channel', None, None),
                          (0.5, 7.0, 'wave',
                           ui.DataPHA._hc / 7.0, ui.DataPHA._hc / 0.5),
                          fail(0.5, 7.0, 'channel', 35, 480)  # issue 619 see also 308
                         ])
def test_calc_flux_pha_analysis(elo, ehi, setting, lo, hi, make_data_path, clean_astro_ui):
    """Do calc_photon/energy_flux return the expected results: fluxes + analysis setting

    Basic test for different analysis settings: the
    same range (modulo precision of conversion) gives the
    same results.
    """

    infile = make_data_path('3c273.pi')
    pl = ui.create_model_component('powlaw1d', 'pl')

    ui.load_pha(infile)
    ui.set_source(pl)

    pflux = ui.calc_photon_flux(elo, ehi)
    eflux = ui.calc_energy_flux(elo, ehi)

    ui.set_analysis(setting)
    pflux2 = ui.calc_photon_flux(lo, hi)
    eflux2 = ui.calc_energy_flux(lo, hi)

    # use approx here since the bin edges are not guaranteed
    # to line up, and use a large tolerance.
    #
    assert pflux2 == pytest.approx(pflux, rel=1e-2)

    eflux = np.log10(eflux)
    eflux2 = np.log10(eflux2)
    assert eflux2 == pytest.approx(eflux, rel=1e-3)


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "foo"])
@pytest.mark.parametrize("energy", [0.05,
                                    0.1,
                                    0.109,
                                    0.11,
                                    1,
                                    1.015,
                                    5,
                                    10.99,
                                    10.991,
                                    11,
                                    15])
def test_calc_flux_density_pha(id, energy, make_data_path, clean_astro_ui):
    """Do calc_photon/energy_flux return the expected results: densities

    The answer should be the same when lo is set and hi
    is None or vice versa. The flux densities are going to
    be in units of <value>/cm^2/s/keV  {value=photon, erg}

    Note: this tests the "edge" condition when lo=hi; this
    is not documented, but left in as a check (and perhaps
    it should be documented).

    """

    infile = make_data_path('3c273.pi')

    # The 3c273 RMF was generated over the range 0.1 to 11 keV
    # (inclusive). By picking gamma = 1 the photon flux
    # density is ampl / e and the energy flux is scale * ampl.
    # However, this is the exact calculation, but the one done by
    # Sherpa involves calculating the model over a bin and then
    # dividing by that bin width, which is different enough to
    # the analytic formula that we use this approach here when
    # calculating the expected values. The bin width is 0.01 keV,
    # with bins starting at 0.1.
    #
    ampl = 1e-4

    # flux densities: exact
    # pflux_exp = ampl / energy
    # eflux_exp = 1.602e-9 * ampl

    # Note that you can calculate an answer at the left edge of the grid, but
    # not at the right (this is either a < vs <= comparison, or numeric
    # issues with the maximum grid value).
    #
    # Note that the RMF emin is just above 0.1 keV, ie
    # 0.10000000149011612, which is why an energy of 0.1 gives
    # an answer of 0. Now, sao_fcmp(0.1, 0.10000000149011612, sherpa.utils.eps)
    # returns 0, so we could consider these two values equal, but
    # that would complicate the flux calculation so is (currently) not
    # done.
    #
    de = 0.01
    if energy <= 0.1 or energy >= 11:
        pflux_exp = 0.0
        eflux_exp = 0.0
    else:
        # assuming a bin centered on the energy; the actual grid
        # is not this, but this should be close
        hwidth = de / 2
        pflux_exp = ampl * (np.log(energy + hwidth) - np.log(energy - hwidth)) / de
        eflux_exp = 1.602e-9 * energy * pflux_exp

    pl = ui.create_model_component('powlaw1d', 'pl')
    pl.ampl = ampl
    pl.gamma = 1

    if id is None:
        ui.load_pha(infile)
        ui.set_source(pl)
    else:
        ui.load_pha(id, infile)
        ui.set_source(id, pl)

    # Use a subset of the data range (to check that the calc routines
    # ignores them, ie uses the full 0.1 to 11 keV range.
    #
    ui.ignore(None, 0.5)
    ui.ignore(7, None)

    # Do not use named arguments, but assume positional arguments
    if id is None:
        pflux1 = ui.calc_photon_flux(energy)
        pflux2 = ui.calc_photon_flux(None, energy)
        pflux3 = ui.calc_photon_flux(energy, energy)

        eflux1 = ui.calc_energy_flux(energy)
        eflux2 = ui.calc_energy_flux(None, energy)
        eflux3 = ui.calc_energy_flux(energy, energy)
    else:
        pflux1 = ui.calc_photon_flux(energy, None, id)
        pflux2 = ui.calc_photon_flux(None, energy, id)
        pflux3 = ui.calc_photon_flux(energy, energy, id)

        eflux1 = ui.calc_energy_flux(energy, None, id)
        eflux2 = ui.calc_energy_flux(None, energy, id)
        eflux3 = ui.calc_energy_flux(energy, energy, id)

    eflux1 = np.log10(eflux1)
    eflux2 = np.log10(eflux2)
    eflux3 = np.log10(eflux3)

    # Use equality here since the numbers should be the same
    assert pflux1 == pflux2
    assert pflux1 == pflux3
    assert eflux1 == eflux2
    assert eflux1 == eflux3

    # Note the "large" tolerance here
    eflux_exp = np.log10(eflux_exp)
    assert pflux1 == pytest.approx(pflux_exp, rel=5e-2)
    assert eflux1 == pytest.approx(eflux_exp, rel=1e-3)


@requires_data
@requires_fits
def test_calc_flux_pha_unabsorbed(make_data_path, clean_astro_ui):
    """Can we calculate an unabsorbed flux?"""

    # The idea is that with a model expression of
    #    const1d.scale * powlaw1d.pl
    # when scale is not 1 (and not integrated) then we can
    # just look to see if the "absorbed" flux is scale * the
    # "unabsorbed" flux.
    #
    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    scale = ui.create_model_component('const1d', 'scale')
    pl = ui.create_model_component('powlaw1d', 'pl')

    scale.c0 = 0.8
    scale.integrate = False
    pl.gamma = 1.5
    pl.ampl = 1e-4

    ui.set_source(scale * pl)

    pflux_abs = ui.calc_photon_flux(0.5, 7)
    pflux_unabs = ui.calc_photon_flux(0.5, 7, model=pl)

    eflux_abs = ui.calc_energy_flux(0.5, 7)
    eflux_unabs = ui.calc_energy_flux(0.5, 7, model=pl)

    pflux_scale = pflux_abs / pflux_unabs
    eflux_scale = eflux_abs / eflux_unabs

    assert pflux_scale == pytest.approx(0.8)
    assert eflux_scale == pytest.approx(0.8)


@requires_data
@requires_fits
@pytest.mark.parametrize("method,sflux,bflux",
                         [(ui.calc_energy_flux, 8.296745792997814e-13, 1.342566520894761e-13),
                          (ui.calc_photon_flux, 0.00034801650283229483, 2.9675854718842685e-05)])
def test_calc_flux_bkg(method, sflux, bflux, make_data_path, clean_astro_ui):
    """Basic test of a background dataset.

    This does not test all combinations, as it is expected that
    they can be checked with existing tests.
    """

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    # ignore low-energy so do not need XSPEC
    ui.ungroup()
    ui.notice(0.8, 7)

    ui.set_stat('cstat')

    spl = ui.create_model_component('powlaw1d', 'spl')
    bpl = ui.create_model_component('powlaw1d', 'bpl')

    spl.gamma = 1.91
    spl.ampl = 1.85e-4
    bpl.gamma = 0.73
    bpl.ampl = 9.27e-6

    ui.set_source(spl)
    ui.set_bkg_source(bpl)

    ui.fit()

    # just to check the fit hasn't significantly moved
    check = ui.calc_stat()
    assert check == pytest.approx(775.6453986960231)

    # These should be the same, by construction
    sflux_all = method(0.5, 7)
    sflux_spl = method(0.5, 7, model=spl)
    assert sflux_all == sflux_spl  # do not use pytest.approx

    bflux_all = method(0.5, 7, bkg_id=1)
    bflux_bpl = method(0.5, 7, bkg_id=1, model=bpl)
    assert bflux_all == bflux_bpl  # do not use pytest.approx

    # check expected flux (regression test).
    #
    got = np.log10(sflux_spl)
    exp = np.log10(sflux)
    assert got == pytest.approx(exp)

    got = np.log10(bflux_bpl)
    exp = np.log10(bflux)
    assert got == pytest.approx(exp)


def setup_sample(id, make_data_path, fit=True):
    """Set up the given dataset for a sample*flux call

    The calling function needs @requires_data, @requires_fits, @requires_xspec
    decorators.
    """

    infile = make_data_path('3c273.pi')
    if id is None:
        ui.load_pha(infile)
        ui.subtract()
    else:
        ui.load_pha(id, infile)
        ui.subtract(id)

    ui.notice(0.5, 7)
    ui.set_stat('chi2datavar')
    ui.set_method('levmar')

    # get the starting point close to the best fit
    #
    # Use xswabs rather than xsphabs as it's faster and
    # we don't care here about scientific validity of the
    # results.
    #
    # gal = ui.create_model_component('xsphabs', 'gal')
    gal = ui.create_model_component('xswabs', 'gal')
    pl = ui.create_model_component('powlaw1d', 'pl')
    mdl = gal * pl
    gal.nh = 0.0393
    pl.gamma = 2.027
    pl.ampl = 1.96e-4

    if id is None:
        ui.set_source(mdl)
        if fit:
            ui.fit()
    else:
        ui.set_source(id, mdl)
        if fit:
            ui.fit(id)

    return gal, pl


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.sample_energy_flux,
                                    ui.sample_photon_flux])
@pytest.mark.parametrize("niter", [0, -1])
@pytest.mark.parametrize("id", [None, 1, "foo"])
def test_sample_foo_flux_invalid_niter(method, niter, id,
                                       make_data_path, clean_astro_ui):
    """What happens for sample_energy/photon_flux when num is <= 0

    See also test_sample_flux_invalid_niter
    """

    setup_sample(id, make_data_path, fit=False)
    with pytest.raises(ArgumentErr):
        method(lo=0.5, hi=7, id=id, num=niter)


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.sample_energy_flux,
                                    ui.sample_photon_flux])
@pytest.mark.parametrize("etype,correlated,scales", [(ModelErr, False, []),
                                                     (ArgumentErr, False, [[]]),
                                                     (ModelErr, False, [1, 2]),
                                                     (ModelErr, False, [1, 2, 3, 4]),
                                                     (ModelErr, False, [[1, 2], [3, 4]]),
                                                     (ArgumentErr, False, [[1, 2, 3], [3, 4, 5]]),
                                                     (ModelErr, False, np.asarray([])),
                                                     (ArgumentErr, False, np.asarray([[]])),
                                                     (ModelErr, False, np.asarray([1, 2])),
                                                     (ModelErr, False, np.asarray([1, 2, 3, 4])),
                                                     (ModelErr, False, np.asarray([[1, 2], [3, 4]])),
                                                     (ArgumentErr, False, np.asarray([[1, 2, 3], [3, 4, 5]])),
                                                     (ArgumentErr, True, []),
                                                     (ArgumentErr, True, [[]]),
                                                     (ArgumentErr, True, np.asarray([])),
                                                     (ArgumentErr, True, np.asarray([[]])),
                                                     (ArgumentErr, True, [1, 2, 3]),
                                                     (ArgumentErr, True, [[1, 2, 3], [1, 2, 3]]),
                                                     (ArgumentErr, True, np.asarray([1, 2, 3])),
                                                     (ModelErr, True, np.ones((2, 2))),
                                                     (ArgumentErr, True, np.ones((3, 3, 3))),
                                                     (ArgumentErr, False, [1, np.inf, 2]),
                                                     (ArgumentErr, False, [1, 2, None]),
                                                     (ArgumentErr, True, [[0.1, 0.01, 0.02], [0.01, np.nan, 0.05], [0.02, 0.01, 0.08]]),
                                                     (ArgumentErr, False, np.ones(3).reshape(1, 3, 1)),
                                                     (ArgumentErr, True, np.ones(9).reshape(1, 3, 3))])
def test_sample_foo_flux_invalid_scales(method, etype, correlated, scales,
                                        make_data_path, clean_astro_ui):
    """What happens for sample_energy/photon_flux when scales is
    the wrong shape, or contains invalid values

    The scales parameter should be (for this fit with 3 free parameters):
       correlated=True   3 by 3 covariance matrix
                  False  3 by 3 covariance matrix or 3-element errors vector

    Since these checks are done at a low level, we do not have to
    loop over every element (e.g. method, id setting) as done in some
    other checks.
    """

    setup_sample('x', make_data_path, fit=False)
    with pytest.raises(etype):
        method(lo=0.5, hi=7, id='x', num=10,
               correlated=correlated, scales=scales)


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.sample_energy_flux,
                                    ui.sample_photon_flux])
@pytest.mark.parametrize("correlated,scales", [(False, [0.1]),
                                               (False, np.ones((4, 4)).diagonal()),
                                               (False, np.ones((4, 4))),
                                               (True, np.ones((1, 1))),
                                               (True, np.ones((4, 4)))])
def test_sample_foo_flux_invalid_scales2(method, correlated, scales,
                                         make_data_path, clean_astro_ui):
    """A repeat of test_sample_foo_flux_invalid_scales for explicit model components.

    Unlike test_sample_foo_flux_invalid_scales, do not repeat as
    many times. In this case we are testing that the grid is
    either 3 by 3 or 2 by 2, and not checking for "invalid"
    values in the inputs.
    """

    cpts = setup_sample(1, make_data_path, fit=False)

    # Test correlated=False with 1D and 2D
    #                = True      2D
    #
    # The full model has 3 free parameters and cpts[1]
    # has 2 free parameters, so try 1 and 4 parameters.
    #
    with pytest.raises(ModelErr):
        method(lo=0.5, hi=7, id=1, num=10, model=cpts[1],
               correlated=correlated, scales=scales)


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.sample_energy_flux,
                                    ui.sample_photon_flux])
def test_sample_foo_flux_no_free_params(method, make_data_path,
                                        hide_logging, clean_astro_ui):
    """sample_energy/photon_flux when no free parameters.
    """

    cpts = setup_sample(1, make_data_path, fit=False)
    for cpt in cpts:
        ui.freeze(cpt)

    with pytest.raises(FitErr) as exc:
        method(lo=0.5, hi=7, num=1)

    assert str(exc.value) == 'model has no thawed parameters'

    # try with explict model setting
    #
    with pytest.raises(FitErr) as exc:
        method(lo=0.5, hi=7, num=1, model=cpts[1])

    assert str(exc.value) == 'model has no thawed parameters'


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.sample_energy_flux,
                                    ui.sample_photon_flux])
def test_sample_foo_flux_not_a_model(method, make_data_path, clean_astro_ui):
    """sample_energy/photon_flux when src model is not a model.

    This is currently not well handled (you get an error but
    it's not a "nice" error).
    """

    setup_sample(1, make_data_path, fit=False)
    with pytest.raises(AttributeError) as exc:
        method(lo=0.5, hi=7, num=1, model='x')

    assert str(exc.value) == "'str' object has no attribute 'thawedpars'"


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.sample_energy_flux,
                                    ui.sample_photon_flux])
def test_sample_foo_flux_invalid_model(method, make_data_path,
                                       hide_logging, clean_astro_ui):
    """sample_energy/photon_flux when src model is not part of the fit.
    """

    setup_sample(1, make_data_path, fit=False)
    mdl = ui.create_model_component('powlaw1d', 'p1')
    with pytest.raises(ArgumentErr) as exc:
        method(lo=0.5, hi=7, num=1, model=mdl)

    assert str(exc.value) == "Invalid src: 'model contains term not in fit'"


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("multi,single", [(ui.sample_energy_flux,
                                           ui.calc_energy_flux),
                                          (ui.sample_photon_flux,
                                           ui.calc_photon_flux)])
@pytest.mark.parametrize("id", [None, 1, "foo"])
@pytest.mark.parametrize("niter", [1, 2, 10])
@pytest.mark.parametrize("correlated", [False, True])
def test_sample_foo_flux_niter(multi, single, id, niter, correlated,
                               make_data_path, clean_astro_ui):
    """Do the sample_energy/photon_flux do what we expect?

    Iterate n times, check that each iteration is different
    (at least, doesn't match the previous version), and that
    the returned flux is as expected.
    """

    gal, pl = setup_sample(id, make_data_path)

    nh0 = gal.nh.val
    gamma0 = pl.gamma.val
    ampl0 = pl.ampl.val

    ans = multi(lo=0.5, hi=7, id=id, num=niter, correlated=correlated)
    assert ans.shape == (niter, 5)

    # the routine hasn't changed the parameter values
    assert gal.nh.val == nh0
    assert pl.gamma.val == gamma0
    assert pl.ampl.val == ampl0

    # we expect that at least one of the parameter variables
    # is different (technically the random sampler could pick
    # the exact positions we started with, but this is unlikely)
    # pytest.approx is used in case there is any loss in
    # precision in storing the parameter values in the returned
    # NumPy array.
    #
    # Since we loop over iterations the first line checks against
    # the best fit, and then we check that the next iteration is
    # different from the previous iteration.
    #
    for i in range(niter):
        diffs = [ans[i, j] != pytest.approx(p.val)
                 for j, p in enumerate([gal.nh, pl.gamma, pl.ampl], 1)]
        assert any(diffs)

        # Can we re-create this flux?
        gal.nh = ans[i, 1]
        pl.gamma = ans[i, 2]
        pl.ampl = ans[i, 3]

        flux = single(lo=0.5, hi=7, id=id)
        assert ans[i, 0] == flux

    # check clipped is 0/1; ideally this is a boolean and not a
    # float representation
    #
    clipped = ans[:, 4]
    v0 = clipped == 0
    v1 = clipped == 1
    v = v0 | v1
    assert v.all()


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("multi", [ui.sample_energy_flux,
                                   ui.sample_photon_flux])
@pytest.mark.parametrize("correlated,lnh0,gamma0,lampl0,slnh0,sgamma0,slampl0",
                         [(False,
                           -1.4327586455259567, 2.023901305737941, -3.708831467739439,
                           -1.5187931648736508, 0.10547355541516636, -4.643460099308331),
                          (True, -1.412181958091415, 2.0328696702831586, -3.7063310517413326,
                           -1.4914185211470155, 0.10389408045222855, -4.637303782755868)])
def test_sample_foo_flux_params(multi, correlated, lnh0, gamma0, lampl0,
                                slnh0, sgamma0, slampl0,
                                make_data_path, clean_astro_ui,
                                hide_logging, reset_seed):
    """Is the parameter sampling in sample_energy/photon_flux sensible?

    Do the parameter values used in the sampling make sense?

    It is not obvious why we need a relatively high tolerance
    (~1e-3) to test the numbers, given the seed is fixed.
    """

    np.random.seed(4276)

    # Rather than loop over ids, like earlier tests, just pick a non-default
    # one.
    #
    id = 2
    gal, pl = setup_sample(id, make_data_path)
    ui.covar(id)
    errs = ui.get_covar_results()

    # The parameter sampling uses the hard-limit boundaries,
    # not the soft-limit boundaries. This can be seen by
    # artificially constricting the soft limits of the
    # gamma parameter. With the default limits, gamma is 2.03 +/- 0.11,
    # so 3-sigma limits are ~ 1.70 to 2.36. In 1000 iterations
    # we would expect ~ 3 bins to fall outside this range.
    # A 2-sigma limit range is 1.81 - 2.25, and we'd expect
    # ~ 45 bins in 1000 to fall outside this range.
    #
    pl.gamma.min = 1.81
    pl.gamma.max = 2.25

    ans = multi(lo=0.5, hi=7, id=id, num=1000, correlated=correlated)

    # do not expect any IEEE special values here
    assert np.isfinite(ans).all()

    nh = ans[:, 1]
    gamma = ans[:, 2]
    ampl = ans[:, 3]

    # Checks we can run?
    #
    #  - mean/median should be close to best-fit
    #  - sigma should be close to error
    #    (at least for correlated=False)
    #  - min/max should be restricted to parameter limits
    #    (this check is currently not enforced as the sampler
    #    does not enforce this).
    #
    # The checks use relatively-large tolerances, as there is
    # no reason the values will match closely
    #
    # The nH value, as it bumps against the lower bound of 0, has
    # been seen to require a larger tolerance than the other parameters.
    #
    assert np.log10(np.median(nh)) == pytest.approx(lnh0, rel=1e-3)
    assert np.median(gamma) == pytest.approx(gamma0, rel=1e-3)
    assert np.log10(np.median(ampl)) == pytest.approx(lampl0, rel=1e-3)

    assert np.log10(np.std(nh)) == pytest.approx(slnh0, rel=1e-3)
    assert np.std(gamma) == pytest.approx(sgamma0, rel=1e-3)
    assert np.log10(np.std(ampl)) == pytest.approx(slampl0, rel=1e-3)

    # Check against the hard limits, although this is not particularly
    # informative since the hard limits for these parameters are
    # all +/- 3.4e38
    #
    assert nh.min() >= gal.nh.hard_min
    assert nh.max() <= gal.nh.hard_max

    assert gamma.min() >= pl.gamma.hard_min
    assert gamma.max() <= pl.gamma.hard_max

    assert ampl.min() >= pl.ampl.hard_min
    assert ampl.max() <= pl.ampl.hard_max

    # A probabilistic check that the gamma range lies outside
    # the soft limits. It is possible for this check to fail
    # because of the RNG, but with ~45 expected values outside
    # the limits this is unlikely. Now we have a set seed we can
    # check both are triggered.
    #
    gmin = gamma.min() < pl.gamma.min
    gmax = gamma.max() > pl.gamma.max
    assert gmin
    assert gmax

    # a simple check on the flux, it should be > 0
    #
    assert ans[:, 0].min() > 0


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("multi", [ui.sample_energy_flux, ui.sample_photon_flux])
@pytest.mark.parametrize("id", [None, "foo"])
def test_sample_foo_flux_clip_soft(multi, id,
                                   make_data_path, clean_astro_ui):
    """Can we clip with soft limits?
    """

    gal, pl = setup_sample(id, make_data_path)

    nh0 = gal.nh.val
    gamma0 = pl.gamma.val
    ampl0 = pl.ampl.val

    niter = 100

    # make this even more obvious than test_sample_foo_flux_params
    pl.gamma.min = 1.9
    pl.gamma.max = 2.1

    ans = multi(lo=0.5, hi=7, id=id, num=niter, clip='soft')
    assert ans.shape == (niter, 5)

    # check clipped is 0/1; expect just under 50% clipped but just check
    # we have some clipped, not the number
    #
    clipped = ans[:, 4]
    v0 = clipped == 0
    v1 = clipped == 1
    v = v0 | v1
    assert v.all()
    assert v0.any()
    assert v1.any()


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("multi", [ui.sample_energy_flux, ui.sample_photon_flux])
@pytest.mark.parametrize("id", [None, "foo"])
def test_sample_foo_flux_clip_none(multi, id,
                                   make_data_path, clean_astro_ui):
    """We get no clipping.

    Same setup as test_sample_foo_flux_clip_soft
    """

    gal, pl = setup_sample(id, make_data_path)

    nh0 = gal.nh.val
    gamma0 = pl.gamma.val
    ampl0 = pl.ampl.val

    niter = 100

    # make this even more obvious than test_sample_foo_flux_params
    pl.gamma.min = 1.9
    pl.gamma.max = 2.1

    ans = multi(lo=0.5, hi=7, id=id, num=niter, clip='none')
    assert ans.shape == (niter, 5)

    # check clipped is 0/1
    #
    clipped = ans[:, 4]
    v0 = clipped == 0
    v1 = clipped == 1
    v = v0 | v1
    assert v.all()
    assert v0.all()
    assert not(v1.any())


# The covariance matrix should be close to the following
# (found from running the code):
#
#   1.28e-3  3.09e-3  7.41e-7
#   3.09e-3  1.10e-2  2.19e-6
#   7.41e-7  2.19e-6  5.38e-10
#
COVMAT = np.asarray([[1.28e-3, 3.09e-3, 7.41e-7],
                     [3.09e-3, 1.10e-2, 2.19e-6],
                     [7.41e-7, 2.19e-6, 5.38e-10]])


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("multi", [ui.sample_energy_flux,
                                   ui.sample_photon_flux])
@pytest.mark.parametrize("correlated,scales,lnh0,gamma0,lampl0,slnh0,sgamma0,slampl0",
                         [(False, COVMAT,
                           -1.3910107005746297, 2.02960216422235, -3.7055883565005048,
                           -1.5045463684643257, 0.10508987949127284, -4.63785707737552),
                          # since np.sqrt(COVMAT.diagonal()) should == COVMAT
                          # this is a bit pointess
                          (False, np.sqrt(COVMAT.diagonal()),
                           -1.3910107005746297, 2.02960216422235, -3.7055883565005048,
                           -1.5045463684643257, 0.10508987949127284, -4.63785707737552),
                          (True, COVMAT,
                           -1.4143474006245673, 2.0232833551455367, -3.707730785910153,
                           -1.5014997538936377, 0.10428938426161602, -4.6312159319100346)])
def test_sample_foo_flux_scales(multi, correlated, scales,
                                lnh0, gamma0, lampl0,
                                slnh0, sgamma0, slampl0,
                                make_data_path, clean_astro_ui,
                                hide_logging, reset_seed):
    """What happens when we specify the scale parameters?

    Test out sending in the errors to the sample_*_flux routines.
    The form depends on the correlated parameter (1D vs 2D).

    The tests are based on those in test_sample_foo_flux_params
    """

    np.random.seed(888)

    id = 2
    gal, pl = setup_sample(id, make_data_path)

    pl.gamma.min = 1.81
    pl.gamma.max = 2.25

    ans = multi(lo=0.5, hi=7, id=id, num=1000,
                correlated=correlated, scales=scales)

    assert np.isfinite(ans).all()

    nh = ans[:, 1]
    gamma = ans[:, 2]
    ampl = ans[:, 3]

    assert np.log10(np.median(nh)) == pytest.approx(lnh0, rel=1e-3)
    assert np.median(gamma) == pytest.approx(gamma0, rel=1e-3)
    assert np.log10(np.median(ampl)) == pytest.approx(lampl0, rel=1e-3)

    if scales.ndim == 2:
        errs = np.sqrt(scales.diagonal())
    else:
        errs = scales

    snh = np.std(nh)
    sgamma = np.std(gamma)
    sampl = np.std(ampl)

    assert np.log10(snh) == pytest.approx(slnh0, rel=1e-3)
    assert sgamma == pytest.approx(sgamma0, rel=1e-3)
    assert np.log10(sampl) == pytest.approx(slampl0, rel=1e-3)

    assert nh.min() >= gal.nh.hard_min
    assert nh.max() <= gal.nh.hard_max

    assert gamma.min() >= pl.gamma.hard_min
    assert gamma.max() <= pl.gamma.hard_max

    assert ampl.min() >= pl.ampl.hard_min
    assert ampl.max() <= pl.ampl.hard_max

    gmin = gamma.min() < pl.gamma.min
    gmax = gamma.max() > pl.gamma.max
    # assert gmin or gmax
    assert gmin
    assert gmax

    assert ans[:, 0].min() > 0


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("multi", [ui.sample_energy_flux,
                                   ui.sample_photon_flux])
def test_sample_foo_flux_scales_example(multi, make_data_path, clean_astro_ui,
                                        hide_logging, reset_seed):
    """Ensure that one of the examples works as expected.

    It is a simplified version of test_sample_foo_flux_scales
    but parameter errors sent in from the parmaxes field
    of the covariance output.
    """

    np.random.seed(3975529)

    id = None
    gal, pl = setup_sample(id, make_data_path)

    ui.covar()
    scales = ui.get_covar_results().parmaxes

    ans = multi(lo=0.5, hi=7, id=id, num=1000,
                correlated=False, scales=scales)

    assert np.isfinite(ans).all()

    nh = ans[:, 1]
    gamma = ans[:, 2]
    ampl = ans[:, 3]

    assert np.log10(np.median(nh)) == pytest.approx(-1.434730948849745, rel=1e-3)
    assert np.median(gamma) == pytest.approx(2.032662859390957, rel=1e-3)
    assert np.log10(np.median(ampl)) == pytest.approx(-3.707907209935886, rel=1e-3)

    assert np.log10(np.std(nh)) == pytest.approx(-1.5118443869902682, rel=1e-3)
    assert np.std(gamma) == pytest.approx(0.10622782336666241, rel=1e-3)
    assert np.log10(np.std(ampl)) == pytest.approx(-4.620006818656268, rel=1e-3)

    assert ans[:, 0].min() > 0


@requires_data
@requires_fits
def test_bug_276(make_data_path):
    ui.load_pha(make_data_path('3c273.pi'))
    ui.set_model('polynom1d.p1')
    ui.fit()
    ui.covar()
    scal = ui.get_covar_results().parmaxes
    ui.sample_flux(ui.get_model_component('p1'), 0.5, 1, num=5,
                   correlated=False, scales=scal)


@pytest.mark.xfail
@pytest.mark.parametrize("niter", [0, -1])
def test_sample_flux_invalid_niter(niter):
    """What happens when num is 0 or negative?

    To simplify things, run with a non PHA dataset.

    This unfortunately fails at this time due to (probably two)
    different bugs, so the test isn't actually effective.
    """

    xbins = np.arange(5, 10, 1)
    y = np.asarray([10, 12, 11, 12])

    ui.load_arrays(1, xbins[:-1], xbins[1:], y, ui.Data1DInt)
    ui.set_source(1, ui.const1d.mdl)

    ui.set_stat('cash')
    ui.fit()
    ui.covar()

    with pytest.raises(ArgumentErr) as exc:
        ui.sample_flux(lo=6, hi=8, Xrays=False, num=niter)

    assert str(exc.value) == "Invalid num: 'must be a positive integer'"


@requires_data
@requires_fits
@pytest.mark.parametrize("niter", [pytest.param(0, marks=pytest.mark.xfail), -1])
def test_sample_flux_invalid_niter_pha(niter, make_data_path, clean_astro_ui):
    """What happens when num is 0 or negative? PHA dataset

    Since test_sample_flux_invalid_niter (currently) fails,
    also try with a PHA dataset.

    Note that niter=0 fails because internally we add 1 to it
    inside sample_flux before we check its value.
    """

    ui.load_pha(make_data_path('3c273.pi'))
    ui.set_model('polynom1d.p1')
    ui.notice(1, 7)

    ui.set_stat('cash')
    ui.fit()
    ui.covar()

    with pytest.raises(ArgumentErr) as exc:
        ui.sample_flux(lo=2, hi=8, num=niter)

    assert str(exc.value) == "Invalid num: 'must be a positive integer'"


@requires_data
@requires_fits
def test_sample_flux_not_a_model(make_data_path, clean_astro_ui):
    """What happens when the model component is not a model?

    Unfortunately we need a PHA dataset as xrays=False errors
    out at the moment.
    """

    ui.load_pha(make_data_path('3c273.pi'))
    ui.set_model('polynom1d.p1')
    ui.notice(1, 7)

    with pytest.raises(ArgumentTypeErr) as exc:
        ui.sample_flux(lo=6, hi=8, modelcomponent="p1")

    assert str(exc.value) == "'modelcomponent' must be a model"


@requires_data
@requires_fits
def test_sample_flux_invalid_model(hide_logging, make_data_path, clean_astro_ui):
    """What happens when the model component is not part of the source?

    Unfortunately we need a PHA dataset as xrays=False errors
    out at the moment.

    At the moment there is no requirement that modelcomponent is
    part of the model expression, so all this does is check the
    call succeeds. It really should fail!
    """

    ui.load_pha(make_data_path('3c273.pi'))
    ui.set_model('polynom1d.p1')
    ui.notice(1, 7)

    ui.set_stat('cash')
    ui.fit()
    ui.covar()

    m2 = ui.create_model_component('const1d', 'm2')

    # This should error out! For now we just check the return values
    # so we can check when the code changes.
    res = ui.sample_flux(lo=6, hi=8, modelcomponent=m2)
    assert len(res) == 3
    assert res[2].shape == (2, 4)


@pytest.mark.parametrize("conf", [0, 100.1])
def test_sample_flux_invalid_confidence(conf):
    """What happens when confidence is outside the range (0, 100]?
    """

    xbins = np.arange(5, 10, 1)
    y = np.asarray([10, 12, 11, 12])

    ui.load_arrays(1, xbins[:-1], xbins[1:], y, ui.Data1DInt)
    ui.set_source(1, ui.const1d.mdl)

    # confidence is checked before anything else, so this works even
    # if the data isn't sensible for X-ray data
    with pytest.raises(ArgumentErr):
        ui.sample_flux(confidence=conf)


@pytest.mark.xfail
def test_sample_flux_1d_int():
    """Basic test of sample_flux with Data1DInt not DataPHA

    Tests out behavior of
      - Xrays=False
      - a non-chi-square statistic
      - model expression with a single parameter

    At the moment this fails for the same reason test_sample_flux_invalid_niter
    fails.
    """

    xbins = np.arange(5, 10, 1)
    y = np.asarray([10, 12, 11, 12])

    ui.load_arrays(1, xbins[:-1], xbins[1:], y, ui.Data1DInt)
    ui.set_source(1, ui.const1d.mdl)

    ui.set_stat('cash')
    ui.fit()
    ui.covar()

    f1, f2, vals = ui.sample_flux(lo=6, hi=8, Xrays=False, num=10)

    assert f1.shape == (3, )
    assert np.all(f1 == f2)
    assert vals.shape == (11, 3)

    assert f1[0] < f1[1]
    assert f1[0] > f1[2]

    # TODO: write more tests, once bug has been fixed to allow sample_flux
    # to get this far


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("multi,fac", [(ui.sample_energy_flux, 1.1435100396354445),
                                       (ui.sample_photon_flux, 1.1901038918815168)])
@pytest.mark.parametrize("correlated", [False, True])
def test_sample_foo_flux_component(multi, fac, correlated,
                                   make_data_path, clean_astro_ui,
                                   hide_logging, reset_seed):
    """Can we sample just a component?

    The idea is to check that the flux for the unabsorbed
    component is larger (on average) than the absorbed
    component. We use the model argument for both calls,
    which also tests we can send in a multiple-component
    model expression.

    The expected ratio between absorbed and unabsorbed flux
    (the fac argument) was calculated from Sherpa's
    calc_energy/photon_flux. This could have been included in
    this test, but an external value acts as a regression test.

    The assumption is that *for this dataset* the use of correlated
    errors does not significantly affect the distributions,
    so the same scale factor can be used (this is based on the
    best-fit location, so shouldn't depend on the errors).
    """

    np.random.seed(39401)

    id = 'xx'
    gal, pl = setup_sample(id, make_data_path)
    mdl = gal * pl

    nh0 = gal.nh.val
    gamma0 = pl.gamma.val
    ampl0 = pl.ampl.val

    ui.covar()
    errs = ui.get_covar_results().parmaxes

    # restrict to 0.5-2 where the absorption makes a bigger
    # difference, relatively speaking, than the 0.5-7 keV
    # band used in a number of other tests.
    #
    absorbed_def = multi(lo=0.5, hi=2, id=id, num=1000,
                         correlated=correlated)
    absorbed = multi(lo=0.5, hi=2, id=id, num=1000,
                     model=mdl, correlated=correlated)

    unabsorbed = multi(lo=0.5, hi=2, id=id, num=1000,
                       model=pl, correlated=correlated)

    assert absorbed.shape == (1000, 5)
    assert unabsorbed.shape == (1000, 5)

    assert np.isfinite(absorbed).all()
    assert np.isfinite(unabsorbed).all()

    flux_absorbed = absorbed[:, 0]
    flux_unabsorbed = unabsorbed[:, 0]
    assert flux_absorbed.min() > 0
    assert flux_unabsorbed.min() > 0

    # quick check that absorbed_def and absorbed are "the same"
    # (modulo the RNG), by comparing the flux distribution. We
    # do not perform a more-quantitative analysis here as this
    # is assumed to be subsumed by the tests run below on the
    # absorbed values, coupled with the previous tests that have
    # been run.
    #
    flux_def = np.median(absorbed_def[:, 0])
    flux_abs = np.median(flux_absorbed)
    assert flux_abs == pytest.approx(flux_def, rel=0.1)

    # Is the ratio similar to the expected values (which was
    # calculated at the best-fit location)
    #
    ratio = np.median(flux_unabsorbed) / flux_abs
    assert ratio == pytest.approx(fac, rel=0.004)

    # The distributions of the two sets of parameters should be
    # similar, since they are drawn from the same distributions.
    #
    # Compare the medians to the best-fit values and the
    # standard deviations to the covariance estimates.
    #
    nh = absorbed[:, 1]
    gamma = absorbed[:, 2]
    ampl = absorbed[:, 3]

    # We could convert this to a specific check now we have
    # a fixed seed.
    #
    assert np.median(nh) == pytest.approx(nh0, rel=0.2)
    assert np.median(gamma) == pytest.approx(gamma0, rel=0.1)
    assert np.median(ampl) == pytest.approx(ampl0, rel=0.1)

    assert np.std(nh) == pytest.approx(errs[0], rel=0.2)
    assert np.std(gamma) == pytest.approx(errs[1], rel=0.1)
    assert np.std(ampl) == pytest.approx(errs[2], rel=0.1)

    gamma = unabsorbed[:, 2]
    ampl = unabsorbed[:, 3]

    assert np.median(gamma) == pytest.approx(gamma0, rel=0.1)
    assert np.median(ampl) == pytest.approx(ampl0, rel=0.1)

    assert np.std(gamma) == pytest.approx(errs[1], rel=0.1)
    assert np.std(ampl) == pytest.approx(errs[2], rel=0.1)


@requires_data
@requires_fits
@pytest.mark.parametrize("idval", [1, 2])
def test_sample_flux_pha_num1(idval, make_data_path, clean_astro_ui,
                              hide_logging, reset_seed):
    """What happens with 1 iteration?

    This is the same data as test_sample_flux_751_752 but
    with some extra tests (just to check that the fit is
    where we expect it to be; if it changes we should
    review the tests, so make them asserts).
    """

    # Use a different value to test_sample_flux_751_752
    np.random.seed(3704)

    gamma0 = 1.95014
    ampl0 = 1.77506e-4
    stat0 = 16.270233678440196

    # This data set is not heavily absorbed, so can get away with
    # just a powerlaw fit for this example (which means we do not
    # need XSPEC to run the test).
    #
    ui.load_pha(idval, make_data_path('3c273.pi'))
    ui.subtract(idval)
    ui.ignore(None, 1)
    ui.ignore(7, None)
    ui.set_source(idval, ui.powlaw1d.p1)
    ui.fit(idval)
    ui.covar(idval)

    # just to check we are in the right place; that is, if these fail
    # it is not a problem with sample_flux but it (may be) worth
    # erroring out quickly as it means we want to review this test
    #
    sinfo = ui.get_stat_info()
    assert len(sinfo) == 1
    assert sinfo[0].ids == [idval]
    assert sinfo[0].statval == pytest.approx(stat0)
    assert sinfo[0].dof == 28

    pmdl = ui.get_model_component('p1')
    assert pmdl.gamma.val == pytest.approx(gamma0)
    assert pmdl.ampl.val == pytest.approx(ampl0)

    scal = ui.get_covar_results().parmaxes
    flux1, flux2, vals = ui.sample_flux(id=idval, lo=1, hi=5, num=1,
                                        correlated=False, scales=scal)

    assert np.all(flux1 == flux2)
    assert len(flux1) == 3
    assert vals.shape == (2, 5)

    # With a fixed seed we can "guarantee" the values are different.
    #
    for row in vals:
        assert row[1] != pytest.approx(gamma0)
        assert row[2] != pytest.approx(ampl0)
        assert row[3] == 0
        assert row[4] > stat0

    # check that the parameter values have not changed, nor the
    # statitic value, by calling sample_flux
    #
    assert pmdl.gamma.val == pytest.approx(gamma0)
    assert pmdl.ampl.val == pytest.approx(ampl0)
    sinfo = ui.get_stat_info()
    assert len(sinfo) == 1
    assert sinfo[0].ids == [idval]
    assert sinfo[0].statval == pytest.approx(stat0)
    assert sinfo[0].dof == 28


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.sample_energy_flux, ui.sample_photon_flux])
@pytest.mark.parametrize("correlated,scales3",
                         [(False, [0.04, 0.12, 2.5e-5]),
                          (False, COVMAT),
                          (True, COVMAT)])
def test_sample_foo_flux_component_scales(method, correlated, scales3,
                                          make_data_path, clean_astro_ui,
                                          hide_logging, reset_seed):
    """Can we sample just a component and send in errors?

    Since the full model (gal * pl) has 3 free parameters
    and the power-law 2, do we get the "same" results
    when using the full errors (nh, gamma, ampl) as
    just (gamma, ampl)?

    Checks for
       - correlated=False  scales=1D array
       - correlated=False  scales=2D covmat
       - correlated=True   scales=2D covmat

    """

    np.random.seed(283491)

    id = 2
    cpts = setup_sample(id, make_data_path)
    pl = cpts[1]

    gamma0 = pl.gamma.val
    ampl0 = pl.ampl.val

    unabsorbed3 = method(lo=0.5, hi=2, id=id, num=1000,
                         model=pl, correlated=correlated,
                         scales=scales3)

    assert unabsorbed3.shape == (1000, 5)

    assert np.isfinite(unabsorbed3).all()

    flux_unabsorbed3 = unabsorbed3[:, 0]
    assert flux_unabsorbed3.min() > 0

    # The distributions of the two sets of parameters should be
    # similar, since they are drawn from the same distributions.
    #
    # Compare the medians to the best-fit values and the
    # standard deviations to the covariance estimates.
    #
    # Now we have a fixed seed we could check the actual values
    #
    ans = np.asarray(scales3)
    if ans.ndim == 2:
        errs3 = np.sqrt(ans.diagonal())
    else:
        errs3 = ans

    gamma = unabsorbed3[:, 2]
    ampl = unabsorbed3[:, 3]

    assert np.median(gamma) == pytest.approx(gamma0, rel=0.1)
    assert np.median(ampl) == pytest.approx(ampl0, rel=0.1)

    assert np.std(gamma) == pytest.approx(errs3[1], rel=0.1)
    assert np.std(ampl) == pytest.approx(errs3[2], rel=0.1)


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method", [ui.sample_energy_flux, ui.sample_photon_flux])
@pytest.mark.parametrize("id", [None, 1, "foo"])
def test_sample_foo_flux_component_scales_fitpars(method, id,
                                                  make_data_path, clean_astro_ui):
    """Is the fit unchanged when a component + errors are used?

    sample_energy/photon_flux can change the model parameters,
    including frozen status, in particular when a component
    of the fit is used and errors for just that component are
    given. This may no-longer be true, but worth a check.

    This test checks that running the sample does not change
    the fit statistic, number of degrees of freedom, and
    parameter values. This repeats some of the checks in
    test_sample_foo_flux_niter but for a different set of
    inputs.
    """

    # save some time by not fitting
    gal, pl = setup_sample(id, make_data_path, fit=False)
    mdl = gal * pl

    pnames0 = [p.fullname for p in mdl.pars]
    pvals0 = [p.val for p in mdl.pars]
    pstate0 = [p.frozen for p in mdl.pars]

    stats0 = ui.get_stat_info()
    assert len(stats0) == 1
    stats0 = stats0[0]
    if id is None:
        assert stats0.ids == [1]
    else:
        assert stats0.ids == [id]

    # just to check that the "fit" hasn't changed from the expected value
    assert stats0.numpoints == 42
    assert stats0.dof == 39

    def validate():
        """Check the fit is unchanged"""

        stats = ui.get_stat_info()
        assert len(stats) == 1
        stats = stats[0]
        for field in ["name", "ids", "bkg_ids", "statname",
                      "statval", "numpoints", "dof", "qval", "rstat"]:
            v = getattr(stats, field)
            v0 = getattr(stats0, field)
            assert v == v0

        for par, name, val, state in zip(mdl.pars,
                                         pnames0,
                                         pvals0,
                                         pstate0):
            assert par.fullname == name
            assert par.val == val
            assert par.frozen == state

    # add in fake errors for the nh values
    errs = [0.1, 0.12, 3e-5]
    cmat = [[0.1, 0, 0], [0, 0.12, 3e-6], [0, 3e-6, 6e-10]]

    # uncorrelated, give errors
    ans = method(lo=0.2, hi=10, id=id, num=2, model=pl, correlated=False,
                 scales=errs)
    assert ans.shape == (2, 5)
    validate()

    # uncorrelated, give covariance matrix
    ans = method(lo=0.2, hi=10, id=id, num=2, model=pl, correlated=False,
                 scales=cmat)
    assert ans.shape == (2, 5)
    validate()

    # correlated, give covariance matrix
    ans = method(lo=0.2, hi=10, id=id, num=2, model=pl, correlated=True,
                 scales=cmat)
    assert ans.shape == (2, 5)
    validate()


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("method,fluxval1,fluxval2",
                         [(ui.sample_energy_flux,
                           8.416796789450763e-13, 1.5474830654516002e-13),
                          (ui.sample_photon_flux,
                           0.00034170234129722176, 3.3578827014096626e-05)])
@pytest.mark.parametrize("id", [None, "foo"])
def test_sample_foo_flux_bkg(method, fluxval1, fluxval2,
                             id, make_data_path, clean_astro_ui,
                             hide_logging, reset_seed):
    """Basic test when calculating flux with a background model.

    fluxval is a value used to check that the median values for the
    source and background runs are "very roughly" correct: for
    the energy flux, we expect the source and background rates to
    be ~8e13 and ~1e-13 (with variance larger on the background),
    so we use 6e-13 as the separation. Ditto for photon flux with
    ~4e-4 and ~3e-5 fluxes.

    """

    np.random.seed(22843)

    infile = make_data_path('3c273.pi')
    if id is None:
        ui.load_pha(infile)
        ui.ungroup()
    else:
        ui.load_pha(id, infile)
        ui.ungroup(id)

    ui.notice(0.5, 7)
    ui.set_stat('cstat')
    ui.set_method('levmar')  # use levmar even with cstat for this case

    gal = ui.create_model_component('xswabs', 'gal')
    spl = ui.create_model_component('powlaw1d', 'spl')
    bpl = ui.create_model_component('powlaw1d', 'bpl')
    mdl = gal * spl

    gal.nh = 0.029
    spl.gamma = 1.96
    spl.ampl = 1.9e-4
    bpl.gamma = 0.70
    bpl.ampl = 9.0e-6

    if id is None:
        ui.set_source(mdl)
        ui.set_bkg_source(bpl)
        ui.fit()
    else:
        ui.set_source(id, mdl)
        ui.set_bkg_source(id, bpl)
        ui.fit(id)

    # check bkg fit is sensible
    res = ui.get_fit_results()
    assert res.numpoints == 892
    assert res.dof == 887
    assert res.statval == pytest.approx(803.208946605444)
    assert res.parnames == ('gal.nH', 'spl.gamma', 'spl.ampl', 'bpl.gamma', 'bpl.ampl')
    assert res.succeeded

    # Check how many parameters are returned:
    #   aflux: nh, source gamma, source ampl
    #   bflux: bgnd gamma, bgnd ampl
    #
    niter = 10
    aflux = method(0.5, 7, id=id, num=niter)
    assert aflux.shape == (niter, 7)

    bflux = method(0.5, 7, id=id, bkg_id=1, num=niter)
    assert bflux.shape == (niter, 7)

    # Compare the median values to the input values.
    #
    assert np.log10(np.median(aflux[:, 0])) == pytest.approx(np.log10(fluxval1), rel=1e-4)
    assert np.log10(np.median(bflux[:, 0])) == pytest.approx(np.log10(fluxval2), rel=1e-4)

    # check the gamma values: source~2, bgnd~0.7
    #
    assert np.median(aflux[:, 2]) == pytest.approx(1.9575001493165511, rel=1e-4)
    assert np.median(bflux[:, 4]) == pytest.approx(0.6022031224500035, rel=1e-4)

    # This assumes there's at least one clipped value; this is not
    # guaranteed but it looks like it happens consistently
    print(aflux[:, 6])
    assert aflux[:, 6].min() == 0
    assert aflux[:, 6].max() == 1


@requires_data
@requires_fits
@requires_xspec
def test_sample_foo_flux_multi(make_data_path, clean_astro_ui,
                               hide_logging, reset_seed):
    """Basic test when calculating flux with multiple datasets
    and, for completeness, fitting the background.
    """

    np.random.seed(8290573)

    # use xswabs as it is simple and unlikely to change
    # (rather than one of the newwe absorption models)
    #
    gal = ui.create_model_component('xswabs', 'gal')
    spl = ui.create_model_component('powlaw1d', 'spl')
    bpl = ui.create_model_component('powlaw1d', 'bpl')
    mdl = gal * spl

    gal.nh = 0.04
    gal.nh.freeze()

    spl.gamma = 1.77
    spl.ampl = 7.3e-7
    bpl.gamma = 1.42
    bpl.ampl = 1.0e-6

    # read in obs1, obs3, obs4
    for did in range(1, 5):
        if did == 2:
            continue

        # mix and match integer and string ids
        if did == 3:
            did = str(did)

        ui.load_pha(did, make_data_path('obs{}.pi'.format(did)))
        ui.set_source(did, mdl)
        ui.set_bkg_source(did, bpl)

    ui.notice(0.5, 7)

    ui.set_stat('cstat')
    ui.set_method('levmar')  # use levmar even with cstat for this case

    ui.fit()

    # check fit is sensible
    res = ui.get_fit_results()
    assert res.numpoints == 2676
    assert res.dof == 2672
    assert res.statval == pytest.approx(1038.2827482111202)
    assert res.parnames == ('spl.gamma', 'spl.ampl', 'bpl.gamma', 'bpl.ampl')
    assert res.datasets == (1, '3', 4)
    assert res.succeeded

    niter = 100

    # Fitting to all datasets should return similar errors (not identical
    # because random, and the energy grids could be different, but aren't
    # in this dataset).
    #
    s1 = ui.sample_energy_flux(lo=0.5, hi=7, id=1, otherids=('3', 4),
                               model=spl, num=niter)
    b1 = ui.sample_energy_flux(lo=0.5, hi=7, id=1, otherids=('3', 4),
                               bkg_id=1, model=bpl, num=niter)

    s3 = ui.sample_energy_flux(lo=0.5, hi=7, id='3', otherids=(1, 4),
                               model=spl, num=niter)
    b3 = ui.sample_energy_flux(lo=0.5, hi=7, id='3', otherids=(1, 4),
                               bkg_id=1, model=bpl, num=niter)

    assert s1.shape == (niter, 6)
    assert b1.shape == (niter, 6)

    assert s3.shape == (niter, 6)
    assert b3.shape == (niter, 6)

    # Compare the median and std dev of the gamma parameter (as of order 1)
    # for both the source and background measurements, as they should be
    # similar.
    #
    y1 = np.median(s1[:, 1])
    y3 = np.median(s3[:, 1])
    assert y1 == pytest.approx(1.7798978430394854), 'source gamma: median'
    assert y3 == pytest.approx(1.83540454021743)

    y1 = np.median(b1[:, 3])
    y3 = np.median(b3[:, 3])
    assert y1 == pytest.approx(1.388683605319035), 'background gamma: median'
    assert y3 == pytest.approx(1.3891467682303402)

    y1 = np.std(s1[:, 1])
    y3 = np.std(s3[:, 1])
    assert y1 == pytest.approx(0.35157485056777504), 'source gamma: std'
    assert y3 == pytest.approx(0.370558381331603)

    y1 = np.std(b1[:, 3])
    y3 = np.std(b3[:, 3])
    assert y1 == pytest.approx(0.1402803370971452), 'background gamma: std'
    assert y3 == pytest.approx(0.14460146969815185)

    # If we compare to a single run we should see larger errors
    # for the single-run case.
    #
    s4 = ui.sample_energy_flux(lo=0.5, hi=7, id=4, otherids=(),
                               model=spl, num=niter)
    b4 = ui.sample_energy_flux(lo=0.5, hi=7, id=4, otherids=(),
                               bkg_id=1, model=bpl, num=niter)

    assert s4.shape == (niter, 6)
    assert b4.shape == (niter, 6)

    # The point is that y4 > y3 in both these (compare to previously)
    y4 = np.std(s4[:, 1])
    assert y4 == pytest.approx(0.5597874155740422)

    y4 = np.std(b4[:, 3])
    assert y4 == pytest.approx(0.23594231624955062)

    # would like to assume there's at least one clipped value for s4
    # but not clear
    assert s4[:, 5].min() == 0
    assert s4[:, 5].max() <= 1


@requires_data
@requires_fits
@pytest.mark.parametrize("idval", [1, 2])
def test_sample_flux_751_752(idval, make_data_path, clean_astro_ui,
                             hide_logging, reset_seed):
    """Very basic test of sample_flux.

    Based around issue #751 (not all iterations have a statistic
    value) and #752 (fails when the id is not the default value).

    Based on test_sample_flux_pha_num1.

    """

    # Try to ensure we get some clipping.
    #
    np.random.seed(4073)

    # we want niter to be even so that the returned number of
    # values (niter+1) is odd, as that makes the median check
    # easier.
    #
    niter = 100

    ui.load_pha(idval, make_data_path('3c273.pi'))
    ui.subtract(idval)
    ui.ignore(None, 1)
    ui.ignore(7, None)
    ui.set_source(idval, ui.powlaw1d.p1)
    ui.fit(idval)
    ui.covar(idval)
    scal = ui.get_covar_results().parmaxes

    # Restrict gamma so that (hopefully, as it depends on the RNG, but seems
    # to do so with niter >~ 10) it triggers clipping in
    # sample_flux, which is the cause of #751. Using niter=100 I
    # see clipping fractions of ~20-30% (before fixing #751).
    #
    p1 = ui.get_model_component('p1')
    ui.set_par(p1.gamma, min=1.8, max=2.1)

    flux1, flux2, vals = ui.sample_flux(id=idval, lo=1, hi=5, num=niter,
                                        correlated=False, scales=scal)

    # as modelcomponent is None, first two arguments should be the same
    assert np.all(flux1 == flux2)
    assert flux1.shape == (3, )
    assert vals.shape == (niter + 1, 5)

    # Basic checks on the flux errors:
    #   - are they in the right order
    #   - manual check that they are close to the expected value
    #     (the default is confidence=68).
    #     as we have an odd number of iterations
    #     THIS IS NOT DONE YET - see #286, in particular
    #     https://github.com/sherpa/sherpa/issues/286#issuecomment-596207243
    #
    fmed, fusig, flsig = flux1
    assert fusig > fmed
    assert flsig < fmed

    # It isn't clear what should be tested here (see #286) so just
    # test something very basic
    #
    fluxes = vals[:, 0]
    assert flsig > fluxes.min()
    assert fusig < fluxes.max()

    # The clip values are not set (since these are set to the hard limits)
    #
    clips = vals[:, -2]
    assert (clips == 0).all()

    # All the statistic values should be set (> 0). This is #751.
    #
    stats = vals[:, -1]
    # assert (stats > 0).all()
    assert (stats > 0).sum() == 81


@requires_data
@requires_fits
@pytest.mark.parametrize("getfunc,medflux",
                         [(ui.get_photon_flux_hist, 1.276591979716474e-4),
                          (ui.get_energy_flux_hist, 4.550271338687814e-13)])
def test_get_xxx_flux_hist_unabsorbed(getfunc, medflux, make_data_path, clean_astro_ui,
                                      reset_seed, hide_logging):
    """Can we get the histogram data for fluxes (for an unabsorbed flux?)"""

    np.random.seed(427247)

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    scale = ui.create_model_component('const1d', 'scale')
    pl = ui.create_model_component('powlaw1d', 'pl')

    scale.c0 = 0.8
    scale.integrate = False
    pl.gamma = 1.5
    pl.ampl = 1e-4

    ui.set_source(scale * pl)

    # Would like to use model=pl but needs #803 (I think)
    flux = getfunc(0.7, 7, num=500, bins=10)
    assert flux.flux.shape == (500, )
    assert flux.modelvals.shape == (500, 3)
    assert flux.xlo.size == 11
    assert flux.xhi.size == 11
    assert flux.y.size == 11

    assert (flux.xlo[1:] == flux.xhi[:-1]).all()

    assert np.median(flux.flux) == pytest.approx(medflux)


@requires_plotting
@requires_data
@requires_fits
@pytest.mark.parametrize("plotfunc",
                         [ui.plot_photon_flux, ui.plot_energy_flux])
def test_plot_xxx_flux_unabsorbed(plotfunc, make_data_path, clean_astro_ui,
                                  reset_seed, hide_logging):
    """Can we plot the histogram data for fluxes (for an unabsorbed flux?)

    There is essentially no check that we do the right thing
    """

    np.random.seed(427248)

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    scale = ui.create_model_component('const1d', 'scale')
    pl = ui.create_model_component('powlaw1d', 'pl')

    scale.c0 = 0.8
    scale.integrate = False
    pl.gamma = 1.5
    pl.ampl = 1e-4

    ui.set_source(scale * pl)

    # Would like to use model=pl but needs #803 (I think)
    plotfunc(0.7, 7, num=500, bins=10)

    # What do we do to test this?


@requires_data
@requires_fits
def test_sample_flux_errors(make_data_path, clean_astro_ui,
                            hide_logging, reset_seed):
    """Does sample_flux return sensible sigma limits?

    This can be thought of as test_sample_flux_751_752 but just
    looking at the flux median/upper/lower limits (to have a
    simple test for addressing #286).

    There is also a test of the statistics "column" information.

    """

    np.random.seed(28737)

    # now we have a seed can we get tighter constraints

    # want to get good stats
    niter = 1000

    ui.load_pha(make_data_path('3c273.pi'))
    ui.subtract()
    ui.ignore(None, 1)
    ui.ignore(7, None)
    ui.set_source(ui.powlaw1d.p1)
    ui.fit()

    stat0 = ui.calc_stat()

    res = ui.sample_flux(lo=0.5, hi=7, num=niter, correlated=False)

    # check the source model hasn't been changed (note: do not use
    # pytest.approx as expect the same value, but can switch if
    # needed).
    #
    assert ui.calc_stat() == stat0

    # Get values close to 1 to make comparison easier
    res[0] *= 1e14
    fmed = res[0][0]
    fusig = res[0][1]
    flsig = res[0][2]

    # Regression test
    #
    elo = fmed - flsig
    ehi = fusig - fmed
    assert fmed == pytest.approx(77.27966775437665)
    assert elo == pytest.approx(9.773389936230018)
    assert ehi == pytest.approx(11.417501513386611)

    # The last column of res[2] is the statistic value for a row -
    # although due to filtering we don't know which row without some
    # significant work. Check that the values make sense (ie if we
    # have found the best-fit then all values should be >= stat0).
    #
    assert res[2].shape == (niter + 1, 5)
    stats = res[2][:, 4]

    # For this dataset no rows are excluded.
    #
    assert (stats > 0).all()

    # Ideally the best-fit location is the best fit!
    assert stats.min() >= stat0

    # Check the clip column too
    #
    clips = res[2][:, 3]
    assert (clips == 0).all()


@requires_data
@requires_fits
@pytest.mark.parametrize("idval", [1, 2])
def test_sample_flux_pha_component(idval, make_data_path, clean_astro_ui,
                                   hide_logging, reset_seed):
    """Does the component analysis work?

    Apply a model which has a normalization term in it (which is fixed)
    so that we can easily compare the full and component fluxes.

    The parameter checks are "generic" (ie don't need to be
    done with a model component given), but are included here because
    this uses multiple model components in the source expression.

    """

    np.random.seed(97284)

    # Now we have a fixed seed we don't need to run as many iterations.
    #
    niter = 100

    ui.load_pha(idval, make_data_path('3c273.pi'))
    ui.subtract(idval)
    ui.ignore(None, 1)
    ui.ignore(7, None)
    ui.set_source(idval, ui.const1d.scal * ui.powlaw1d.p1)

    p1 = ui.get_model_component('p1')
    scal = ui.get_model_component('scal')
    scal.c0 = 0.8
    scal.integrate = False
    scal.c0.frozen = True

    ui.fit(idval)
    ui.covar(idval)

    stat0 = ui.calc_stat(idval)

    scal = ui.get_covar_results().parmaxes
    flux1, flux2, vals = ui.sample_flux(modelcomponent=p1,
                                        id=idval, lo=1, hi=5, num=niter,
                                        correlated=False, scales=scal)

    # check the source model hasn't been changed (note: do not use
    # pytest.approx as expect the same value, but can switch if
    # needed).
    #
    assert ui.calc_stat(idval) == stat0

    assert flux1.shape == (3, )
    assert flux2.shape == (3, )
    assert np.any(flux1 != flux2)
    assert vals.shape == (niter + 1, 5)

    assert flux1[0] < flux1[1]
    assert flux1[0] > flux1[2]

    # Now, flux2 should be 1/0.8 = 1.25 times larger,
    # and this should hold for all three terms (since drawn from the
    # same distributions modulo the different model components there
    # is no "randomness" in this ratio).
    #
    for i in range(3):
        rlim = flux2[i] / flux1[i]
        assert rlim == pytest.approx(1 / 0.8)

    # Regression test
    assert np.log10(flux1[0]) == pytest.approx(-12.320132853894787)

    # Regression test of the output (if these change, compare against
    # p1.gamma.val, p1.ampl.val).
    #
    mgamma = np.median(vals[:, 1])
    mampl = np.log10(np.median(vals[:, 2]))
    assert mgamma == pytest.approx(1.9476566223146197)
    assert mampl == pytest.approx(-3.6498020707240153)

    # No clipping
    #
    assert (vals[:, 3] == 0).all()

    # No missing values (statistic is set for all bins)
    #
    assert (vals[:, 4] > 0).all()


@requires_data
@requires_fits
@pytest.mark.parametrize("idval", [1, 2])
def test_sample_flux_pha_bkg_no_source(idval, make_data_path, clean_astro_ui,
                                       hide_logging, reset_seed):
    """Can we get the flux for the background fit when there's no source model?"""

    np.random.seed(287241)

    niter = 100

    ui.load_pha(idval, make_data_path('3c273.pi'))
    ui.ignore(None, 0.8)
    ui.ignore(7, None)

    ui.set_bkg_source(idval, ui.powlaw1d.bpl)

    bpl = ui.get_model_component('bpl')
    bpl.gamma = 0.54
    bpl.ampl = 6.4e-6

    ui.fit_bkg(idval)

    # At the moment we need a source model (I guess to get the
    # covariance), so this fails.
    #
    with pytest.raises(IdentifierErr) as exc:
        ui.sample_flux(id=idval, bkg_id=1,
                       lo=1, hi=5, num=niter,
                       correlated=False)

    assert str(exc.value) == 'model stack is empty'


@requires_data
@requires_fits
@requires_xspec
@pytest.mark.parametrize("idval", [1, 2])
def test_sample_flux_pha_bkg(idval, make_data_path, clean_astro_ui,
                             hide_logging, reset_seed):
    """Can we get the flux for the background fit.

    In this case we have a source fit that has also been run,
    so the error analysis should be fine.
    """

    ui.set_stat('chi2datavar')
    ui.set_method('levmar')

    niter = 100

    ui.load_pha(idval, make_data_path('3c273.pi'))

    # Need to ungroup so we can use the mask attribute
    ui.ungroup(idval)
    ui.ignore(None, 0.8)
    ui.ignore(7, None)

    d = ui.get_data(idval)
    b = ui.get_bkg(idval)

    # We need to do it this way (background first)
    ui.group_counts(num=20, id=idval, bkg_id=1, tabStops=~b.mask)
    ui.group_counts(num=20, id=idval, tabStops=~d.mask)

    ui.set_bkg_source(idval, ui.powlaw1d.bpl)
    ui.set_source(idval, ui.xswabs.gabs * ui.powlaw1d.pl)

    bpl = ui.get_model_component('bpl')
    bpl.gamma = 0.54
    bpl.ampl = 6.4e-6

    pl = ui.get_model_component('pl')
    pl.gamma = 2.0
    pl.ampl = 2e-4

    gabs = ui.get_model_component('gabs')
    gabs.nh = 0.05
    gabs.nh.freeze()

    ui.fit_bkg(idval)
    ui.fit(idval)

    # We don't try to run covar or get the statistic values
    # (the API doesn't make that easy for background fits).
    #

    # Use the same seed for both runs.
    #
    np.random.seed(287247)
    bflux1, bflux2, bvals = ui.sample_flux(id=idval, bkg_id=1,
                                           lo=1, hi=5, num=niter,
                                           correlated=False)

    np.random.seed(287247)
    flux1, flux2, vals = ui.sample_flux(id=idval,
                                        lo=1, hi=5, num=niter,
                                        correlated=False)

    assert flux2 == pytest.approx(flux1)
    assert bflux2 == pytest.approx(bflux1)

    assert np.log10(flux1) == pytest.approx([-12.31348652, -12.28200881, -12.34610975])
    assert np.log10(bflux1) == pytest.approx([-13.15241991, -12.97986006, -13.31382087])

    # Just check that the flux we've got is close to the one calculated by
    # calc_energy_flux. This is just to test the values we have checked for
    # the median value of flux1 and bflux1.
    #
    # Note that I am slightly concerned that the background flux has been
    # calculated using the source response, not the background response.
    # This is too hard to identify with this dataset. We need a test which
    # has drastically different background response (and maybe exposure
    # time).
    #
    sflux = ui.calc_energy_flux(id=idval, lo=1, hi=5)
    bflux = ui.calc_energy_flux(id=idval, bkg_id=1, lo=1, hi=5)

    assert np.log10(sflux) == pytest.approx(-12.310733810154623)
    assert np.log10(bflux) == pytest.approx(-13.110679628882489)

    # The random parameters are assumed to be the same. In reality the
    # background doesn't need to include the source dataset but that's an
    # issue to look at later.
    #
    assert bvals.shape == (101, 7)
    assert vals.shape == (101, 7)

    # For the value checks we need to drop the first column, which has the
    # fluxes.
    #
    assert (bvals[:, 1:] == vals[:, 1:]).all()
