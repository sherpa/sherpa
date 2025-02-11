#
#  Copyright (C) 2007, 2015 - 2024
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

import numpy

import pytest

from sherpa.astro import ui
from sherpa.models import Parameter
from sherpa.utils.testing import requires_data, \
    requires_fits, requires_xspec
from sherpa.utils.err import ParameterErr

# How many models should there be?
# This number includes all additive, multiplicative, and convolition models,
# even the ones that would be disabled by a decoration from .utils.
# The number can be calculated by counting the occurrences of the strings
#    '(XSAdditiveModel)'
#    '(XSMultiplicativeModel)'
#    '(XSConvolutionKernel)'
# in `xspec/__init__.py`
#
XSPEC_MODELS_COUNT = 281

# Conversion between wavelength (Angstrom) and energy (keV)
# The values used are from sherpa/include/constants.hh
#
_hc = 6.6260693e-27 * 2.99792458e+18 / 1.60217653e-9


def is_proper_subclass(obj, cls):
    if type(cls) is not tuple:
        cls = (cls,)
    if obj in cls:
        return False
    return issubclass(obj, cls)


# Make it easier to run these tests on subsets of the XSPEC
# model library, or in case model names get changed.
#
def remove_item(xs, x):
    """Remove x from xs or do nothing (if x is not a member of xs)."""
    try:
        xs.remove(x)
    except ValueError:
        pass


# There is an argument to be made that these tests only need
# to exercise a small number of models (e.g. one each of the
# different templates used in xspec_extension.hh - so single
# precision, double precision, and table models) but there
# is also benefit in testing them all (in particular, to check
# that interface changes do not significantly change the
# results for any model). For now stick with testing most of
# the models.
#
# This function will contain the convolution models, so
# downstream code needs to handle these differently. Note that
# explicit testing of these models is handled by
# test_xspec_con.py.
#
def get_xspec_models():
    """What are the XSpec model names to test.

    """

    try:
        import sherpa.astro.xspec as xs
    except ImportError:
        return []

    version = xs.get_xsversion()

    # The alternate approach is to use that taken by
    # test_create_model_instances
    #
    model_names = [model_name for model_name in dir(xs)
                   if model_name.startswith('XS')]

    # Could just exclude any names that end in 'Model', but this
    # could remove valid model names, so be explicit.
    for n in ['XSModel', 'XSMultiplicativeModel', 'XSAdditiveModel',
              'XSTableModel', 'XSConvolutionModel', 'XSConvolutionKernel',
              'XSBaseParameter', 'XSParameter']:
        remove_item(model_names, n)

    # the grbjet model with XSPEC 12.12.0 (and presumably 12.12.0.a) can
    # occasionally evaluate to all 0's. This has been reported, but for now
    # skip this model.
    #
    remove_item(model_names, 'XSgrbjet')

    # The bsedov model causes a crashe with XSPEC 12.14.0 to 12.14.0e
    # (it should be fixed in 12.4.0f and later). The model is not
    # present before 12.14.0.
    #
    if version in ["12.14.0", "12.14.0a", "12.14.0b",
                   "12.14.0c", "12.14.0d", "12.14.0e"]:
        remove_item(model_names, 'XSbsedov')

    # The bvvnei model causes a crash with XSPEC 12.14.0 to 12.14.0h
    # (it should be fixed in 12.14.0i and later). The model is not
    # present before 12.14.0.
    #
    if version in ["12.14.0", "12.14.0a", "12.14.0b",
                   "12.14.0c", "12.14.0d", "12.14.0e",
                   "12.14.0f", "12.14.0g", "12.14.0h"]:
        remove_item(model_names, 'XSbvvnei')

    models = [getattr(xs, model_name) for model_name in model_names]
    models = list(filter(lambda mod: mod.version_enabled, models))

    return models


def make_grid():
    """Return the 'standard' contiguous grid used in these tests.

    Returns elo, ehi, wlo, whi where the e-xxx values are in keV and
    the w-xxx values in Angstrom. The *lo/*hi values are this grid
    separated out into bin edges.

    The following condition holds: *hi[i] > *lo[i], and the
    values in the two arrays are either monotonically increasing
    or decreasing.

    """

    # This is close to the standard energy range used by Chandra,
    # but it can take time to run when iterating through over
    # 100 models. So chose a coarser grid and a smaller energy
    # range.
    #
    # egrid = numpy.arange(0.1, 11.01, 0.01)
    # egrid = numpy.arange(0.1, 5.01, 0.1)

    # XSgaussian has peak at 6.5 keV and sigma ~ 0.1 keV
    # so need to get to around 7 keV
    egrid = numpy.arange(0.1, 7.01, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    wgrid = _hc / egrid
    whi = wgrid[:-1]
    wlo = wgrid[1:]

    return elo, ehi, wlo, whi


def make_grid_noncontig2():
    """Return the 'standard' non-contiguous grid used in these tests,
    with 2 gaps.

    The grid from make_grid is modified to make one gap.
    """

    elo, ehi, wlo, whi = make_grid()

    # remove two sections
    idx1 = (elo <= 1.1) | (elo >= 1.4)
    idx2 = (elo <= 2.3) | (elo >= 2.9)
    idx = idx1 & idx2

    # at the moment do not return the original grid and filter
    return elo[idx], ehi[idx], wlo[idx], whi[idx]


# In Sherpa 4.8.0 development, the XSPEC model class included
# an explicit check that all the bins were finite. This has
# been removed as testing has shown there are complications when
# using this with some real-world responses. So add in an
# explicit - hopefully temporary - test here.
#
# Also ensure that at least one bin is > 0
# (technically we could have a model with all negative
#  values, or all 0 with the default values, but this
#  seems unlikely; it is more likely a model evaluates
#  to 0 - or at least not positive - because a data file
#  is missing).
#
def assert_is_finite(vals, modelcls, label):
    """modelcls is the model class (e.g. xspec.XSphabs)"""

    import sherpa.astro.xspec as xs

    emsg = f"model {modelcls} is finite [{label}]"
    assert numpy.isfinite(vals).all(), emsg

    # Some models seem to return 0's, so skip them for now:
    # these models have a default redshift parameter of 0 but
    # the code complains if z <= 0 and returns 0's.
    #
    if modelcls in [xs.XSbcph, xs.XSbvcph, xs.XScph, xs.XSvcph]:
        assert (vals == 0.0).all(), \
            f'Expected {modelcls} to evaluate to all zeros [{label}]'
        return

    emsg = f"Expected model {modelcls} to have a value > 0 [{label}]"
    assert (vals > 0.0).any(), emsg


@requires_xspec
def test_create_model_instances(clean_astro_ui):
    import sherpa.astro.xspec as xs
    count = 0

    for cls in dir(xs):
        if not cls.startswith('XS'):
            continue

        cls = getattr(xs, cls)

        if is_proper_subclass(cls, (xs.XSAdditiveModel,
                                    xs.XSMultiplicativeModel,
                                    xs.XSConvolutionKernel)):
            # Ensure that we can create an instance, but do
            # nothing with it.
            cls()
            count += 1

    assert count == XSPEC_MODELS_COUNT


@requires_xspec
def test_check_default_name():
    import sherpa.astro.xspec as xs

    for clname in dir(xs):
        if not clname.startswith('XS'):
            continue

        cls = getattr(xs, clname)
        if is_proper_subclass(cls, (xs.XSAdditiveModel,
                                    xs.XSMultiplicativeModel,
                                    xs.XSConvolutionKernel)):

            # At the moment we have some defaulting to xs... and some just ...
            # (the forner are convolution cases which should probably be
            # switched to drop the leading xs).
            #
            mdl = cls()
            expected = clname.lower()
            assert mdl.name in [expected, expected[2:]]


@requires_xspec
def test_norm_works():
    # Check that the norm parameter for additive models
    # works, as it is handled separately from the other
    # parameters.
    import sherpa.astro.xspec as xs

    # need an additive model
    mdl = xs.XSpowerlaw()
    mdl.PhoIndex = 2
    egrid1 = [0.1, 0.2, 0.3, 0.4]
    egrid2 = [0.2, 0.3, 0.4, 0.5]

    mdl.norm = 1.2
    y1 = mdl(egrid1, egrid2)

    mfactor = 2.1
    mdl.norm = mdl.norm.val * mfactor
    y2 = mdl(egrid1, egrid2)

    # check that the sum is not 0 and that it
    # scales as expected.
    s1 = y1.sum()
    s2 = y2.sum()

    assert s1 > 0.0, 'powerlaw is positive'
    assert s2 == pytest.approx(mfactor * s1)


@requires_xspec
def test_evaluate_model():
    import sherpa.astro.xspec as xs
    mdl = xs.XSbbody()
    out = mdl([1, 2, 3, 4], [2, 3, 4, 5])
    if mdl.calc.__name__.startswith('C_'):
        otype = numpy.float64
    else:
        otype = numpy.float32

    assert out.dtype.type == otype
    # check all values are > 0
    assert (out > 0).all()


# Select a few models which are likely to cover
# additive/multiplicative and language (e.g. FORTRAN vs C/C++).
#
BASIC_MODELS = ['powerlaw', 'gaussian',
                'vapec',  # pick this as scientifically "useful"
                'constant', 'wabs']


@requires_xspec
@pytest.mark.parametrize('model', BASIC_MODELS)
def test_lowlevel(model):
    """The XSPEC class interface requires lo,hi but the low-level allows just x

    Pick a few additive and multiplicative models.
    """

    import sherpa.astro.xspec as xs

    cls = getattr(xs, 'XS{}'.format(model))
    mdl = cls()

    pars = [p.val for p in mdl.pars]

    # grid chosen to match XSgaussian's default parameter setting
    # (to make sure evaluates to > 0).
    #
    egrid = numpy.arange(6, 7, 0.1)
    e1 = egrid[:-1]
    e2 = egrid[1:]

    y1 = mdl._calc(pars, egrid)
    y2 = mdl._calc(pars, e1, e2)

    # should be able to use equality rather than approx, but
    # this is easier to use
    assert y1[:-1] == pytest.approx(y2)
    assert y1[-1] == 0.0
    assert (y2 > 0).all()


@requires_xspec
@pytest.mark.parametrize('model', BASIC_MODELS)
def test_lowlevel_checks_too_many_arguments(model):
    """Check we get a sensible error when called with no arguments.

    Note that this tests the interface to the actual XSPEC model
    not the Python class.
    """

    import sherpa.astro.xspec as xs

    cls = getattr(xs, 'XS{}'.format(model))
    mdl = cls()

    with pytest.raises(TypeError) as exc:
        mdl._calc()

    assert str(exc.value) == "function missing required argument 'pars' (pos 1)"


@requires_xspec
def test_checks_input_length():
    import sherpa.astro.xspec as xs
    mdl = xs.XSpowerlaw()

    # Check when input array is too small (< 2 elements)
    with pytest.raises(TypeError) as exc1:
        mdl([0.1], [0.2])

    assert str(exc1.value) == "input array must have at least 2 elements, found 1"

    # Check when input arrays are not the same size.
    with pytest.raises(TypeError) as exc2:
        mdl([0.1, 0.2, 0.3], [0.2, 0.3])

    assert str(exc2.value) == "input arrays are not the same size: 3 and 2"

    with pytest.raises(TypeError) as exc3:
        mdl([0.1, 0.2], [0.2, 0.3, 0.4])

    assert str(exc3.value) == "input arrays are not the same size: 2 and 3"


@requires_xspec
@requires_data
@requires_fits
@pytest.mark.parametrize('loadfunc', [ui.load_xstable_model, ui.load_table_model])
def test_xstablemodel_checks_input_length(loadfunc, clean_astro_ui, make_data_path):

    loadfunc('mdl', make_data_path('xspec-tablemodel-RCS.mod'))
    mdl = ui.get_model_component('mdl')

    # Check when input array is too small (< 2 elements)
    with pytest.raises(TypeError) as exc1:
        mdl([0.1], [0.2])

    emsg = "input array must have at least 2 elements, found 1"
    assert str(exc1.value) == emsg

    # Check when input arrays are not the same size (when the
    # low and high bin edges are given)
    with pytest.raises(TypeError) as exc2:
        mdl([0.1, 0.2, 0.3], [0.2, 0.3])

    emsg = "input arrays are not the same size: 3 and 2"
    assert str(exc2.value) == emsg

    with pytest.raises(TypeError) as exc3:
        mdl([0.1, 0.2], [0.2, 0.3, 0.4])

    emsg = "input arrays are not the same size: 2 and 3"
    assert str(exc3.value) == emsg


@requires_xspec
@requires_data
@requires_fits
@pytest.mark.parametrize('loadfunc', [ui.load_xstable_model, ui.load_table_model])
def test_xspec_xstablemodel(loadfunc, clean_astro_ui, make_data_path):
    # Just test one table model; use the same scheme as
    # test_xspec_models_noncontiguous().
    #
    # The table model is from
    # https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/rcs.html
    # retrieved July 9 2015. The exact model is irrelevant for this
    # test, so this was chosen as it's relatively small.
    loadfunc('tmod', make_data_path('xspec-tablemodel-RCS.mod'))

    # when used in the test suite it appears that the tmod
    # global symbol is not created, so need to access the component
    tmod = ui.get_model_component('tmod')

    assert tmod.name == 'xstablemodel.tmod'

    elo, ehi, wlo, whi = make_grid()
    evals = tmod(elo, ehi)
    wvals = tmod(wlo, whi)

    assert_is_finite(evals, tmod, "energy")
    assert_is_finite(wvals, tmod, "wavelength")
    assert wvals == pytest.approx(evals)


@requires_xspec
@requires_data
@requires_fits
@pytest.mark.parametrize('loadfunc', [ui.load_xstable_model, ui.load_table_model])
def test_xspec_xstablemodel_noncontiguous2(loadfunc, clean_astro_ui, make_data_path):
    loadfunc('tmod', make_data_path('xspec-tablemodel-RCS.mod'))
    tmod = ui.get_model_component('tmod')

    elo, ehi, wlo, whi = make_grid_noncontig2()

    evals = tmod(elo, ehi)
    wvals = tmod(wlo, whi)

    assert_is_finite(evals, tmod, "energy")
    assert_is_finite(wvals, tmod, "wavelength")
    assert wvals == pytest.approx(evals)
    assert (wvals > 0).all()


@requires_xspec
@requires_data
@requires_fits
def test_xpec_tablemodel_outofbound(clean_astro_ui, make_data_path):
    ui.load_xstable_model('tmod', make_data_path('xspec-tablemodel-RCS.mod'))
    # when used in the test suite it appears that the tmod
    # global symbol is not created, so need to access the component
    tmod = ui.get_model_component('tmod')
    elo = numpy.arange(1, 5)
    ehi = elo + 1

    print(tmod)

    with pytest.raises(ParameterErr) as e:
        tmod.calc([0., .2, 1., 1.], elo, ehi)

    assert 'minimum' in str(e)


@requires_xspec
def test_convolution_model_cflux():
    """This tests the low-level interfce of the convolution model"""

    # Use the cflux convolution model, since this gives
    # an easily-checked result.
    #
    import sherpa.astro.xspec as xs

    if not hasattr(xs._xspec, 'C_cflux'):
        pytest.skip('cflux convolution model is missing')

    # The energy grid should extend beyond the energy grid
    # used to evaluate the model, to avoid any edge effects.
    # It also makes things easier if the elo/ehi values align
    # with the egrid bins.
    elo = 0.55
    ehi = 1.45
    egrid = numpy.linspace(0.5, 1.5, 101)
    eg1 = egrid[:-1]
    eg2 = egrid[1:]

    mdl1 = xs.XSpowerlaw()
    mdl1.PhoIndex = 2

    # flux of mdl1 over the energy range of interest; converting
    # from a flux in photon/cm^2/s to erg/cm^2/s, when the
    # energy grid is in keV.
    y1 = mdl1(eg1, eg2)
    idx, = numpy.where((egrid >= elo) & (egrid < ehi))

    # To match XSpec, need to multiply by (Ehi^2-Elo^2)/(Ehi-Elo)
    # which means that we need the next bin to get Ehi. Due to
    # the form of the where statement, we should be missing the
    # Ehi value of the last bin
    e1 = egrid[idx]
    e2 = egrid[idx + 1]

    f1 = 8.01096e-10 * ((e2 * e2 - e1 * e1) * y1[idx] / (e2 - e1)).sum()

    # The cflux parameters are elo, ehi, and the log of the
    # flux within this range (this is log base 10 of the
    # flux in erg/cm^2/s). The parameters chosen for the
    # powerlaw, and energy range, should have f1 ~ 1.5e-9
    # (log 10 of this is -8.8).
    lflux = -5.0
    pars = [elo, ehi, lflux]

    y1_a = numpy.zeros(y1.size + 1)
    y1_a[:-1] = y1
    y2_a = xs._xspec.C_cflux(pars, y1_a, egrid)
    y2_b = xs._xspec.C_cflux(pars, y1, eg1, eg2)

    assert y2_a[:-1] == pytest.approx(y2_b)
    assert y2_a[-1] == 0.0

    assert_is_finite(y2_b, "cflux", "energy")

    elo = egrid[:-1]
    ehi = egrid[1:]
    wgrid = _hc / egrid
    whi = wgrid[:-1]
    wlo = wgrid[1:]

    expected = y1 * 10**lflux / f1
    assert y2_b == pytest.approx(expected)

    y1 = mdl1(elo, ehi)
    y2 = xs._xspec.C_cflux(pars, y1, elo, ehi)
    assert y2 == pytest.approx(expected)

    y1 = mdl1(wlo, whi)
    y2 = xs._xspec.C_cflux(pars, y1, wlo, whi)
    assert y2 == pytest.approx(expected)


@requires_xspec
def test_convolution_model_cpflux_noncontiguous():
    """convolution models require a contiguous grid"""

    import sherpa.astro.xspec as xs

    if not hasattr(xs._xspec, 'C_cpflux'):
        pytest.skip('cpflux convolution model is missing')

    elo, ehi, wlo, whi = make_grid_noncontig2()

    lflux = -5.0
    pars = [0.2, 0.8, lflux]
    y1 = numpy.zeros(elo.size)

    emsg = "XSPEC convolution model requires a contiguous grid"

    with pytest.raises(ValueError) as exc1:
        xs._xspec.C_cpflux(pars, y1, elo, ehi)

    assert str(exc1.value) == emsg

    with pytest.raises(ValueError) as exc2:
        xs._xspec.C_cpflux(pars, y1, wlo, whi)

    assert str(exc2.value) == emsg


@requires_xspec
@requires_data
@requires_fits
def test_set_analysis_wave_fabrizio(clean_astro_ui, make_data_path):
    rmf = make_data_path('3c273.rmf')
    arf = make_data_path('3c273.arf')

    ui.set_model("fabrizio", "xspowerlaw.p1")
    ui.fake_pha("fabrizio", arf, rmf, 10000)

    parvals = [1, 1]

    model = ui.get_model("fabrizio")
    bare_model, _ = ui._session._get_model_status("fabrizio")
    y = bare_model.calc(parvals, model.xlo, model.xhi)
    y_m = numpy.mean(y)

    ui.set_analysis("fabrizio", "wave")

    model2 = ui.get_model("fabrizio")
    bare_model2, _ = ui._session._get_model_status("fabrizio")
    y2 = bare_model2.calc(parvals, model2.xlo, model2.xhi)
    y2_m = numpy.mean(y2)

    assert y2_m == pytest.approx(y_m)


@requires_xspec
def test_xsxset_get(clean_astro_ui):
    import sherpa.astro.xspec as xs
    # TEST CASE #1 Case insentitive keys
    xs.set_xsxset('fooBar', 'somevalue')
    assert xs.get_xsxset('Foobar') == 'somevalue'


@requires_xspec
def test_additive_single_norm_model():
    """Check that we can not sneak in a separate norm parameter"""

    from sherpa.astro.xspec import XSAdditiveModel

    class XSNotAClassName(XSAdditiveModel):
        __function__ = "foo"

        def __init__(self, name='foo'):
            self.kT = Parameter(name, 'kT', 1.0)
            self.norm = False
            XSAdditiveModel.__init__(self, name, (self.kT, ))

    # Setting self.norm and not including it in the pars array
    # is an error.
    with pytest.raises(ParameterErr,
                       match="norm is set but not included in pars"):
        XSNotAClassName()


@requires_xspec
def test_nonexistent_model():
    from sherpa.astro.xspec.utils import include_if
    from sherpa.astro.xspec import XSAdditiveModel

    @include_if(False)
    class XSbtapec(XSAdditiveModel):
        __function__ = "foo"

        def __init__(self, name='foo'):
            self.kT = Parameter(name, 'kT', 1.0)
            XSAdditiveModel.__init__(self, name, (self.kT, ))

    m = XSbtapec()

    with pytest.raises(AttributeError) as exc:
        m([], [])

    assert include_if.DISABLED_MODEL_MESSAGE.format("XSbtapec") == str(exc.value)


@requires_xspec
def test_not_compiled_model():
    """
    Test the error handling case where a model is included according to the conditional decorator, but it wraps a
    function that had not been compiled.
    """
    from sherpa.astro.xspec.utils import include_if, ModelMeta
    from sherpa.astro.xspec import XSAdditiveModel

    @include_if(True)
    class XSfoo(XSAdditiveModel):
        __function__ = "C_foo"

        def __init__(self, name='foo'):
            self.kT = Parameter(name, 'kT', 1.0)
            XSAdditiveModel.__init__(self, name, (self.kT, ))

    m = XSfoo()

    with pytest.raises(AttributeError) as exc:
        m([], [])

    assert ModelMeta.NOT_COMPILED_FUNCTION_MESSAGE == str(exc.value)


@requires_xspec
def test_old_style_xspec_class():
    """
    We changed the way xspec models are declared, but just in case let's make sure old-style declarations still work.
    """
    from sherpa.astro.xspec import XSzbabs, XSMultiplicativeModel, _xspec

    class XSfoo(XSMultiplicativeModel):
        _calc = _xspec.xszbabs

        def __init__(self, name='zbabs'):
            self.nH = Parameter(name, 'nH', 1.e-4, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
            self.nHeI = Parameter(name, 'nHeI', 1.e-5, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
            self.nHeII = Parameter(name, 'nHeII', 1.e-6, 0.0, 1.0e5, 0.0, 1.0e6, '10^22 atoms / cm^2')
            self.redshift = Parameter(name, 'redshift', 0.0, 0.0, 1.0e5, 0.0, 1.0e6)
            XSMultiplicativeModel.__init__(self, name, (self.nH, self.nHeI, self.nHeII, self.redshift))

    m = XSfoo()

    actual = m([1, 2, 3], [2, 3, 4])
    expected = XSzbabs()([1, 2, 3], [2, 3, 4])
    assert actual == pytest.approx(expected)


@requires_xspec
@pytest.mark.parametrize("modelcls", get_xspec_models())
def test_evaluate_xspec_model(modelcls):
    """Can we call a model with its default parameters?

    There is limited validation of the results. A regression-style test could
    be made (e.g. with use of pytest-arraydiff) but is it worth it?

    Rather than loop over each model within the test, use pytest to loop over
    each model. This means that make_grid gets called multiple times, but
    (hopefully) makes it a lot easier to debug when only a subset of tests
    fails.  It also makes sure that all models are run, even if some fail,
    which did not happen when they were all stuck in the same test function.

    Convolution models are skipped (easier to filter out here given the
    current design).
    """

    from sherpa.astro import xspec

    # use an identifier in case there is an error
    mdl = modelcls()
    if isinstance(mdl, xspec.XSConvolutionKernel):
        return

    elo, ehi, wlo, whi = make_grid()

    # The model checks that the values are all finite,
    # so there is no need to check that the output of
    # mdl does not contain non-finite values.
    #
    evals = mdl(elo, ehi)
    wvals = mdl(wlo, whi)

    assert_is_finite(evals, modelcls, "energy")
    assert_is_finite(wvals, modelcls, "wavelength")

    assert wvals == pytest.approx(evals)


@requires_xspec
@pytest.mark.parametrize("modelcls", get_xspec_models())
def test_evaluate_xspec_model_noncontiguous2(modelcls):
    """Can we evaluate an XSPEC model with a non-contiguous grid?

    Note that there is no test that the non-contiguous grid
    results are similar to the result from using a contiguous
    grid and then filtering out the missing bins. The way that
    some models are implemented make this a tricky test to
    write (due to numerical tolerances), as bins at the
    edges may not match well.
    """

    from sherpa.astro import xspec

    # use an identifier in case there is an error
    mdl = modelcls()
    if isinstance(mdl, xspec.XSConvolutionKernel):
        return

    elo, ehi, wlo, whi = make_grid_noncontig2()

    evals2 = mdl(elo, ehi)
    wvals2 = mdl(wlo, whi)

    assert_is_finite(evals2, modelcls, "energy")
    assert_is_finite(wvals2, modelcls, "wavelength")

    assert wvals2 == pytest.approx(evals2)


@requires_xspec
def test_apec_redshift_parameter_is_case_insensitive():
    """This should be a given, but just check."""

    import sherpa.astro.xspec as xs

    mdl = xs.XSapec()
    assert mdl.redshift.val == 0
    assert mdl.Redshift.val == 0

    mdl.redshift = 0.4
    assert mdl.redshift.val == 0.4
    assert mdl.Redshift.val == 0.4

    mdl.REDshift = 0.2
    assert mdl.redshift.val == 0.2
    assert mdl.Redshift.val == 0.2
