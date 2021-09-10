#
#  Copyright (C) 2020, 2021
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

import numpy as np
from numpy.testing import assert_allclose

import pytest

from sherpa.utils.testing import requires_xspec, requires_data, requires_fits
from sherpa.data import Data1DInt
from sherpa.models.basic import Box1D, Const1D, PowLaw1D
from sherpa.models.parameter import Parameter
from sherpa.models.model import RegridWrappedModel

from sherpa.astro import ui

# Unlike test_xspec.py which tries to evaluate all the additive and
# multiplicative models, here we just try several ones for which
# we can easily evaluate the results
#
#   - cflux: calculates the flux over an energy range as a parameter
#   - [zv][am]shift: shift a model
#

# Identify when new models are added (in case tests should be written),
# as well as the number of thawed and frozen parameters.
#
XSPEC_CON_MODELS = [('cflux', 1, 2),
                    ('clumin', 1, 3),
                    ('cpflux', 1, 2),
                    ('gsmooth', 1, 1),
                    ('ireflect', 1, 6),
                    ('kdblur', 1, 3),
                    ('kdblur2', 1, 5),
                    ('kerrconv', 2, 5),
                    ('kyconv', 2, 10),
                    ('lsmooth', 1, 1),
                    ('partcov', 1, 0),
                    ('rdblur', 1, 3),
                    ('reflect', 1, 4),
                    ('rfxconv', 2, 3),
                    ('rgsxsrc', 0, 1),
                    ('simpl', 2, 1),
                    ('thcomp', 3, 1),
                    ('vashift', 0, 1),
                    ('vmshift', 0, 1),
                    ('xilconv', 2, 4),
                    ('zashift', 0, 1),
                    ('zmshift', 0, 1)]


@requires_xspec
def test_count_xspec_convolution_models():
    """Do we have all the models we expect"""

    from sherpa.astro import xspec

    for n, _, _ in XSPEC_CON_MODELS:
        clsname = 'XS{}'.format(n)
        try:
            cls = getattr(xspec, clsname)
        except AttributeError:
            assert False, 'missing convolution model {}'.format(clsname)

        mdl = cls()
        assert issubclass(cls, xspec.XSConvolutionKernel), mdl


@requires_xspec
def test_check_no_extra_xspec_convolution_models():
    """Are we missing any?

    This is for the test suite (check we know about everything)
    rather than catching any errors in the code.
    """

    from sherpa.astro import xspec

    exclude = ['XSModel', 'XSAdditiveModel', 'XSMultiplicativeModel',
               'XSConvolutionModel', 'XSConvolutionKernel',
               'XSBaseParameter', 'XSParameter']

    names = [(n, getattr(xspec, n))
             for n in dir(xspec)
             if n.startswith('XS') and n not in exclude]

    classes = set([])
    for name, cls in names:
        try:
            mdl = cls()
        except TypeError:
            pass

        if isinstance(mdl, xspec.XSConvolutionKernel):
            npars = len(mdl.pars)
            nfrozen = len([p for p in mdl.pars if p.frozen])
            classes.add((name[2:], npars - nfrozen, nfrozen))

    assert classes == set(XSPEC_CON_MODELS)


def setup_data(elo=0.1, ehi=10.0, ebin=0.01):
    """Return a data set.

    Parameters
    ----------
    elo, ehi : number, optional
        The start and end energy of the grid, in keV.
    ebin : number, optional
        The bin width, in keV.

    Returns
    -------
    data : Data1DInt
        The data object. The Y values are not expected to be used,
        so are set to 1.

    """

    if elo >= ehi:
        raise ValueError("elo >= ehi")
    if elo <= 0:
        raise ValueError("elo <= 0")
    if ebin <= 0:
        raise ValueError("ebin <= 0")

    x = np.arange(elo, ehi + ebin, ebin)
    if x.size < 2:
        raise ValueError("elo, ehi, ebin not sensible")

    y = np.ones(x.size - 1)
    return Data1DInt('dummy', x[:-1], x[1:], y)


def _check_pars(label, mdl, parvals):
    """Check that we have the expected parameters.

    Parameters
    ----------
    label : str
        The label to use in the assertions
    mdl
        The model object. It must have a pars attribute.
    parvals : sequence of tuples
        Each entry gives the parameter name, value, frozen flag,
        and units field.
    """

    nexp = len(parvals)
    ngot = len(mdl.pars)
    assert nexp == ngot, '{}: number of parameters'.format(label)

    for i, vals in enumerate(parvals):

        par = mdl.pars[i]
        plbl = '{} param {}: '.format(label, i + 1)
        assert isinstance(par, Parameter), plbl + 'is a parameter'
        assert par.name == vals[0], plbl + 'name'
        assert par.val == vals[1], plbl + 'value'
        assert par.frozen == vals[2], plbl + 'frozen'
        assert par.units == vals[3], plbl + 'units'


@requires_xspec
def test_cflux_settings():
    """Do the expected things happen when a model is calculated?"""

    from sherpa.astro import xspec

    kern = xspec.XScflux('cflux')
    assert isinstance(kern, xspec.XSConvolutionKernel), \
        "cflux creates XSConvolutionKernel"

    cfluxpars = [('Emin', 0.5, True, 'keV'),
                 ('Emax', 10.0, True, 'keV'),
                 ('lg10Flux', -12, False, 'cgs')]
    _check_pars('cflux', kern, cfluxpars)

    mdl = kern(PowLaw1D('pl'))
    assert isinstance(mdl, xspec.XSConvolutionModel), \
        "cflux(mdl) creates XSConvolutionModel"

    plpars = [('gamma', 1.0, False, ''),
              ('ref', 1.0, True, ''),
              ('ampl', 1.0, False, '')]
    _check_pars('model', mdl, cfluxpars + plpars)


def _test_cflux_calc(mdl, slope, ampl):
    """Test the CFLUX convolution model calculation.

    This is a test of the convolution interface, as the results of
    the convolution can be easily checked.

    Parameters
    ----------
    mdl
        The unconvolved model. It is assumed to be a power law (so
        either a XSpowerlaw or PowLaw1D instance).
    slope, ampl : number
        The slope (gamma or PhoIndex) and amplitude
        (ampl or norm) of the power law.

    See Also
    --------
    test_cflux_calc_sherpa, test_cflux_calc_xspec

    """

    from sherpa.astro import xspec

    d = setup_data()

    kern = xspec.XScflux('cflux')

    mdl_unconvolved = mdl
    mdl_convolved = kern(mdl)

    # As it's just a power law, we can evaluate analytically
    #
    pterm = 1.0 - slope
    emin = d.xlo[0]
    emax = d.xhi[-1]
    counts_expected = ampl * (emax**pterm - emin**pterm) / \
        pterm

    # Could evaluate the model directly, but check that it's working
    # with the Sherpa setup. It is not clear that going through
    # the eval_model method really buys us much testing since the
    # main check we would really want to do is with the DataPHA
    # class, but that requires a lot of set up; the idea is that
    # the interface is abstract enough that going through
    # Data1DInt's eval_model API is sufficient.
    #
    y_unconvolved = d.eval_model(mdl_unconvolved)
    nbins_unconvolved = y_unconvolved.size
    assert nbins_unconvolved == d.xlo.size, \
        'unconvolved model has correct # bins'

    counts_unconvolved = y_unconvolved.sum()

    # This is mainly a sanity check of everything, as it is not
    # directly related to the convolution model. The counts_expected
    # term is > 0.01 with the chosen settings, so 1e-15
    # is a hand-wavy term to say "essentially zero". It may need
    # tweaking depending on platform/setup.
    #
    assert np.abs(counts_expected - counts_unconvolved) < 1e-15, \
        'can integrate a power law'

    # How to verify the convolution model? Well, one way is to
    # give it a flux, get the model values, and then check that
    # these values are very close to the expected values for
    # that flux.
    #
    kern.emin = 0.5
    kern.emax = 7.0

    desired_flux = 3.2e-14
    kern.lg10flux = np.log10(desired_flux)

    y_convolved = d.eval_model(mdl_convolved)
    nbins_convolved = y_convolved.size
    assert nbins_convolved == d.xlo.size, \
        'convolved model has correct # bins'

    # The model evaluated by mdl_convolved should have a
    # log_10(flux), over the kern.Emin to kern.Emax energy range,
    # of kern.lg10Flux, when the flux is in erg/cm^2/s (assuming
    # the model it is convolving has units of photon/cm^2/s).
    #
    # d.xlo and d.xhi are in keV, so need to convert them to
    # erg.
    #
    # Use 1 keV = 1.60218e-9 erg for the conversion, which limits
    # the accuracy (in part because I don't know how this compares
    # to the value that XSPEC uses, which could depend on the
    # version of XSPEC).
    #
    econv = 1.60218e-9
    emid = econv * (d.xlo + d.xhi) / 2.0

    y_signal = emid * y_convolved

    # Do not bother with trying to correct for any partially-filled
    # bins at the end of this range (there shouldn't be any).
    #
    idx = np.where((d.xlo >= kern.emin.val) &
                   (d.xhi <= kern.emax.val))
    flux = y_signal[idx].sum()

    # Initial tests had desired_flux = 3.2e-14 and the difference
    # ~ 1.6e-16. This has been used to determine the tolerance.
    # Note that 2e-16/3.2e-14 ~ 0.006 ie ~ 0.6%
    #
    assert np.abs(flux - desired_flux) < 2e-16, \
        'flux is not as expected'

    # The cflux model should handle any change in the model amplitude
    # (I believe; as there's no one definition of what the model
    # amplitude means in XSPEC). So, increasing the amplitude - assumed
    # to the last parameter of the input model - should result in
    # the same flux values. The unconvolved model is checked to make
    # sure it does scale (as a sanity check).
    #
    rescale = 1000.0
    mdl.pars[-1].val *= rescale
    y_unconvolved_rescaled = d.eval_model(mdl_unconvolved)
    y_convolved_rescaled = d.eval_model(mdl_convolved)

    assert_allclose(y_unconvolved, y_unconvolved_rescaled / rescale,
                    atol=0, rtol=1e-7)
    assert_allclose(y_convolved, y_convolved_rescaled,
                    atol=0, rtol=1e-7)


@requires_xspec
def test_cflux_calc_xspec():
    """Test the CFLUX convolution model calculations (XSPEC model)

    This is a test of the convolution interface, as the results of
    the convolution can be easily checked. The model being convolved
    is an XSPEC model.

    See Also
    --------
    test_cflux_calc_sherpa

    """

    from sherpa.astro import xspec

    mdl = xspec.XSpowerlaw('xspec')
    mdl.phoindex = 1.7
    mdl.norm = 0.025
    _test_cflux_calc(mdl, mdl.phoindex.val, mdl.norm.val)


@requires_xspec
def test_cflux_calc_sherpa():
    """Test the CFLUX convolution model calculations (sherpa model)

    This is a test of the convolution interface, as the results of
    the convolution can be easily checked. The model being convolved
    is a Sherpa model.

    See Also
    --------
    test_cflux_calc_xspec

    """

    mdl = PowLaw1D('sherpa')
    mdl.gamma = 1.7
    mdl.ampl = 0.025
    _test_cflux_calc(mdl, mdl.gamma.val, mdl.ampl.val)


@requires_xspec
def test_cflux_nbins():
    """Check that the number of bins created by cflux is correct.

    The test_cflux_calc_xxx routines do include a test of the number
    of bins, but that is just for a Data1DInt dataset, so the model
    only ever gets called with explicit lo and hi edges. Now that we
    no-longer support evaluating the model with a single grid this
    test may not add much power, but leave in for now.

    Notes
    -----
    There's no check of a non-contiguous grid.

    """

    from sherpa.astro import xspec

    spl = PowLaw1D('sherpa')
    xpl = xspec.XSpowerlaw('xspec')

    spl.gamma = 0.7
    xpl.phoindex = 0.7

    egrid = np.arange(0.1, 2, 0.01)
    elo = egrid[:-1]
    ehi = egrid[1:]

    nbins = elo.size

    def check_bins(lbl, mdl):
        y = mdl(elo, ehi)
        assert y.size == nbins, f'{lbl}: elo/ehi'

    # verify assumptions
    #
    check_bins('Sherpa model', spl)
    check_bins('XSPEC model', xpl)

    # Now try the convolved versions
    #
    cflux = xspec.XScflux("conv")
    check_bins('Convolved Sherpa', cflux(spl))
    check_bins('Convolved XSPEC', cflux(xpl))


@requires_xspec
def test_calc_xspec_regrid():
    """Test the CFLUX convolution model calculations (XSPEC model)

    Can we regrid the model and get a similar result to the
    direct version? This uses zashift but with 0 redshift.

    See Also
    --------
    test_calc_sherpa_regrid

    """

    from sherpa.astro import xspec

    mdl = xspec.XSpowerlaw('xspec')
    mdl.phoindex = 1.7
    mdl.norm = 0.025

    # rather than use the 'cflux' convolution model, try
    # the redshift model but with 0 redshift.
    #
    kern = xspec.XSzashift('zshift')
    kern.redshift = 0
    mdl_convolved = kern(mdl)

    d = setup_data(elo=0.2, ehi=5, ebin=0.01)

    dr = np.arange(0.1, 7, 0.005)
    mdl_regrid = mdl_convolved.regrid(dr[:-1], dr[1:])

    yconvolved = d.eval_model(mdl_convolved)
    yregrid = d.eval_model(mdl_regrid)

    # Do a per-pixel comparison (this will automatically catch any
    # difference in the output sizes).
    #
    ydiff = np.abs(yconvolved - yregrid)
    mdiff = ydiff.max()
    # rdiff = ydiff / yconvolved

    # in testing see the max difference being ~ 3e-17
    assert mdiff < 1e-15, \
        'can rebin an XSPEC powerlaw: max diff={}'.format(mdiff)


@requires_xspec
def test_calc_sherpa_regrid():
    """Test the CFLUX convolution model calculations (Sherpa model)

    Can we redshift a sherpa model and get the expected
    result, with a regrid applied? In this case a feature
    is added outside the default energy range, but in the
    regridded range, to check things are working.

    See Also
    --------
    test_calc_xspec_regrid

    """

    from sherpa.astro import xspec

    mdl = Box1D('box')
    mdl.xlow = 4
    mdl.xhi = 12
    # why is the box amplitude restricted like this?
    mdl.ampl.max = 2
    mdl.ampl = 2

    kern = xspec.XSzashift('zshift')
    kern.redshift = 1.0
    mdl_convolved = kern(mdl)

    d = setup_data(elo=0.1, ehi=10, ebin=0.01)
    dr = setup_data(elo=0.1, ehi=13, ebin=0.005)

    mdl_regrid = mdl_convolved.regrid(dr.xlo, dr.xhi)

    yconvolved = d.eval_model(mdl_convolved)
    yregrid = d.eval_model(mdl_regrid)

    # We expect values < 2 keV to be 0
    #   yconvolved:  > 2 keV to < 5 keV to be 0.02 (ampl / bin width)
    #   yregrid:     > 2 keV to < 6 keV to be 0.02
    # above this to be 0, with some fun at the edges.
    #
    ehi = d.xhi

    idx = np.where(ehi < 2)
    ymax = np.abs(yconvolved[idx]).max()
    assert ymax == 0.0, 'yconvolved < 2 keV: max={}'.format(ymax)

    idx = np.where((ehi >= 2) & (ehi < 5))
    ymax = np.abs(yconvolved[idx] - 0.02).max()
    assert ymax < 1e-14, 'yconvolved: 2-5 keV max={}'.format(ymax)

    idx = np.where(ehi >= 5)
    ymax = np.abs(yconvolved[idx]).max()
    assert ymax == 0.0, 'yconvolved: > 5 keV max={}'.format(ymax)

    # expect last bin to be ~ 0 but preceeding ones to be zero
    idx = np.where(ehi < 2)
    ymax = np.abs(yregrid[idx][:-1]).max()
    assert ymax == 0, 'yregrid < 2 keV: max={}'.format(ymax)

    ymax = np.abs(yregrid[idx])[-1]
    assert ymax < 1e-14, 'yregrid < 2 keV: max={}'.format(ymax)

    idx = np.where((ehi >= 2) & (ehi < 6))
    ymax = np.abs(yregrid[idx] - 0.02).max()
    assert ymax < 2e-14, 'yregrid: 2-6 keV {}'.format(ymax)

    # expect first bin to be ~ 0 but following ones to be zero
    idx = np.where(ehi >= 6)
    ymax = np.abs(yregrid[idx])[0]
    assert ymax < 2e-14, 'yregrid: > 6 keV max={}'.format(ymax)

    ymax = np.abs(yregrid[idx][1:]).max()
    assert ymax == 0.0, 'yregrid: > 6 keV max={}'.format(ymax)


@requires_xspec
def test_xspec_con_ui_registered():
    """Are the convolution models registered for use by the UI layer?"""

    mdls = set(ui.list_models('xspec'))
    for n, _, _ in XSPEC_CON_MODELS:
        xsname = "xs{}".format(n)
        assert xsname in mdls, xsname


@requires_xspec
@requires_data
@requires_fits
def test_xspec_con_ui_cflux(make_data_path, clean_astro_ui, restore_xspec_settings):
    """Check cflux from the UI layer with a response."""

    from sherpa.astro import xspec

    infile = make_data_path('3c273.pi')
    ui.load_pha('random', infile)
    ui.subtract('random')
    ui.ignore(None, 0.5)
    ui.ignore(7, None)

    ui.set_source('random', 'xsphabs.gal * xscflux.sflux(powlaw1d.pl)')
    mdl = ui.get_source('random')

    assert mdl.name == '(xsphabs.gal * xscflux.sflux(powlaw1d.pl))'
    assert len(mdl.pars) == 7
    assert mdl.pars[0].fullname == 'gal.nH'
    assert mdl.pars[1].fullname == 'sflux.Emin'
    assert mdl.pars[2].fullname == 'sflux.Emax'
    assert mdl.pars[3].fullname == 'sflux.lg10Flux'
    assert mdl.pars[4].fullname == 'pl.gamma'
    assert mdl.pars[5].fullname == 'pl.ref'
    assert mdl.pars[6].fullname == 'pl.ampl'

    assert isinstance(mdl.lhs, xspec.XSphabs)
    assert isinstance(mdl.rhs, xspec.XSConvolutionModel)

    gal = ui.get_model_component('gal')
    sflux = ui.get_model_component('sflux')
    pl = ui.get_model_component('pl')
    assert isinstance(gal, xspec.XSphabs)
    assert isinstance(sflux, xspec.XScflux)
    assert isinstance(pl, PowLaw1D)

    # the convolution model needs the normalization to be fixed
    # (not for this example, as we are not fitting, but do this
    # anyway for reference)
    pl.ampl.frozen = True

    sflux.emin = 1
    sflux.emax = 5
    sflux.lg10Flux = -12.3027

    pl.gamma = 2.03
    gal.nh = 0.039

    ui.set_xsabund('angr')
    ui.set_xsxsect('vern')

    # check we get the "expected" statistic (so this is a regression
    # test).
    #
    ui.set_stat('chi2gehrels')
    sinfo = ui.get_stat_info()

    assert len(sinfo) == 1
    sinfo = sinfo[0]
    assert sinfo.numpoints == 40
    assert sinfo.dof == 37
    assert sinfo.statval == pytest.approx(21.25762265234619)

    # Do we get the same flux from Sherpa's calc_energy_flux?
    #
    cflux = ui.calc_energy_flux(id='random', model=sflux(pl), lo=1, hi=5)
    lcflux = np.log10(cflux)
    assert lcflux == pytest.approx(sflux.lg10Flux.val)


@requires_xspec
@requires_data
@requires_fits
def test_xspec_con_ui_shift(make_data_path, clean_astro_ui, restore_xspec_settings):
    """Check shifted models from the UI layer with a response.

    There is no regrid here, so we see the issue with the upper edge of the
    RMF cutting off the source.
    """

    from sherpa.astro import xspec

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)
    ui.subtract()
    ui.ignore(None, 0.5)
    ui.ignore(7, None)

    msource = ui.box1d.box + ui.const1d.bgnd
    ui.set_source(ui.xsphabs.gal * ui.xszashift.zsh(msource))
    mdl = ui.get_source()

    assert mdl.name == '(xsphabs.gal * xszashift.zsh((box1d.box + const1d.bgnd)))'
    assert len(mdl.pars) == 6
    assert mdl.pars[0].fullname == 'gal.nH'
    assert mdl.pars[1].fullname == 'zsh.Redshift'
    assert mdl.pars[2].fullname == 'box.xlow'
    assert mdl.pars[3].fullname == 'box.xhi'
    assert mdl.pars[4].fullname == 'box.ampl'
    assert mdl.pars[5].fullname == 'bgnd.c0'

    assert isinstance(mdl.lhs, xspec.XSphabs)
    assert isinstance(mdl.rhs, xspec.XSConvolutionModel)

    gal = ui.get_model_component('gal')
    zsh = ui.get_model_component('zsh')
    box = ui.get_model_component('box')
    bgnd = ui.get_model_component('bgnd')
    assert isinstance(gal, xspec.XSphabs)
    assert isinstance(zsh, xspec.XSzashift)
    assert isinstance(box, Box1D)
    assert isinstance(bgnd, Const1D)

    zsh.redshift = 1

    # turn off the absorption to make the comparison easier
    gal.nh = 0

    # pick an energy range that exceeds the RMF maximum energy (11 keV)
    box.xlow = 10
    box.xhi = 13
    box.ampl = 0.5

    bgnd.c0 = 0.001
    bgnd.integrate = False

    mplot = ui.get_source_plot()

    # Expect, as z = 1
    #
    #    0.1 for E < 5 keV     (10 / (1+z))
    #    0       E > 5.5 keV   (11 / (1+z))  due to RMF cut off
    #    0.6     5 - 5.5 keV
    #
    idx1 = mplot.xhi <= 5
    idx2 = mplot.xlo >= 5.5

    # ensure we pick the "expected" range
    assert idx1.sum() == 490
    assert idx2.sum() == 550

    assert mplot.xhi[idx1].max() == pytest.approx(5)
    assert mplot.xlo[idx2].min() == pytest.approx(5.5)

    # use the inverse of the two index arrays to ensure we are not
    # missing any bins.
    #
    idx3 = ~(idx1 | idx2)

    # The tolerance has to be relatively large otherwise things fail
    assert mplot.y[idx1] == pytest.approx(0.1, rel=3e-5)
    assert mplot.y[idx2] == pytest.approx(0)
    assert mplot.y[idx3] == pytest.approx(0.6, rel=1e-5)


@requires_xspec
@requires_data
@requires_fits
def test_xspec_con_ui_shift_regrid(make_data_path, clean_astro_ui, restore_xspec_settings):
    """Check shifted models from the UI layer with a response and regrid.

    Unlike test_xspec_con_ui_shift, the convolution model is run on an extended
    grid compared to the RMF.
    """

    from sherpa.astro import xspec

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)
    ui.subtract()
    ui.ignore(None, 0.5)
    ui.ignore(7, None)

    # Ensure the grid contains the RMF grid (0.1-11 keV).
    # Really we should have emax to be > 11 * (1+z) but
    # I purposefully pick a smaller maximum to check we
    # get 0 values in the output
    #
    rgrid = np.arange(0.1, 20, 0.01)
    rlo = rgrid[:-1]
    rhi = rgrid[1:]

    msource = ui.box1d.box + ui.const1d.bgnd
    csource = ui.xszashift.zsh(msource)
    ui.set_source(ui.xsphabs.gal * csource.regrid(rlo, rhi))
    mdl = ui.get_source()

    # What should the string representation be?
    #
    assert mdl.name == '(xsphabs.gal * regrid1d(xszashift.zsh((box1d.box + const1d.bgnd))))'

    assert len(mdl.pars) == 6
    assert mdl.pars[0].fullname == 'gal.nH'
    assert mdl.pars[1].fullname == 'zsh.Redshift'
    assert mdl.pars[2].fullname == 'box.xlow'
    assert mdl.pars[3].fullname == 'box.xhi'
    assert mdl.pars[4].fullname == 'box.ampl'
    assert mdl.pars[5].fullname == 'bgnd.c0'

    assert isinstance(mdl.lhs, xspec.XSphabs)
    assert isinstance(mdl.rhs, RegridWrappedModel)

    gal = ui.get_model_component('gal')
    zsh = ui.get_model_component('zsh')
    box = ui.get_model_component('box')
    bgnd = ui.get_model_component('bgnd')
    assert isinstance(gal, xspec.XSphabs)
    assert isinstance(zsh, xspec.XSzashift)
    assert isinstance(box, Box1D)
    assert isinstance(bgnd, Const1D)

    zsh.redshift = 1

    # turn off the absorption to make the comparison easier
    gal.nh = 0

    # pick an energy range that exceeds the RMF maximum energy (11 keV)
    box.xlow = 10
    box.xhi = 13
    box.ampl = 0.5

    bgnd.c0 = 0.001
    bgnd.integrate = False

    mplot = ui.get_source_plot()

    # Expect, as z = 1
    #
    #    0.1 for E < 5 keV  or  6.5 - 10 keV
    #    0.6     5 - 6.5 keV
    #    0       > 10 keV
    #
    idx1 = (mplot.xhi <= 5) | ((mplot.xlo >= 6.5) & (mplot.xhi <= 10))
    idx2 = (mplot.xlo >= 5) & (mplot.xhi <= 6.5)
    idx3 = mplot.xlo >= 10

    # ensure we pick the "expected" range (there are 1090 bins in the
    # RMF energy grid)
    assert idx1.sum() == 840
    assert idx2.sum() == 150
    assert idx3.sum() == 100

    # The tolerance has to be relatively large otherwise things fail
    #
    # It appears that the very-last bin of idx1 is actually ~ 0.05,
    # so separate that out here. It is the last bin of the "valid"
    # array, and so at z=1 may only have been half-filled by the
    # convolution.
    #
    assert mplot.y[idx1][:-1] == pytest.approx(0.1, rel=1e-5)
    assert mplot.y[idx1][-1] == pytest.approx(0.05, rel=3e-5)

    assert mplot.y[idx2] == pytest.approx(0.6, rel=3e-5)
    assert mplot.y[idx3] == pytest.approx(0)
