#
#  Copyright (C) 2020
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

import pytest

import numpy as np

from sherpa.astro import ui
from sherpa.utils.testing import requires_data, requires_fits
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, IOErr
import sherpa.astro.utils


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


def fails_619(*x):
    """See issue 619

    This isn't strictly needed, but useful for documentation
    """
    return pytest.param(*x, marks=pytest.mark.xfail)


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
                          fails_619(0.5, 7.0, 'channel', 35, 480)  # see also 308
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


def setup_sample(id, make_data_path):
    """Set up the given dataset for a sample*flux call

    The calling function needs @requires_data and @requires_fits
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
    gal = ui.create_model_component('xsphabs', 'gal')
    pl = ui.create_model_component('powlaw1d', 'pl')
    mdl = gal * pl
    gal.nh = 0.0393
    pl.gamma = 2.027
    pl.ampl = 1.96e-4

    if id is None:
        ui.set_source(mdl)
        ui.fit()
    else:
        ui.set_source(id, mdl)
        ui.fit(id)

    return gal, pl


@requires_data
@requires_fits
@pytest.mark.parametrize("method", [ui.sample_energy_flux, ui.sample_photon_flux])
@pytest.mark.parametrize("niter", [0, -1])
@pytest.mark.parametrize("id", [None, 1, 2, "foo"])
def test_sample_foo_flux_invalid_niter(method, niter, id, make_data_path, clean_astro_ui):
    """What happens for sample_energy/photon_flux when num is <= 0

    See also test_sample_flux_invalid_niter
    """

    setup_sample(id, make_data_path)
    with pytest.raises(ArgumentErr):
        method(lo=0.5, hi=7, id=id, num=niter)
