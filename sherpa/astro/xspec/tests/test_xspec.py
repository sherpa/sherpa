#
#  Copyright (C) 2007, 2015, 2016  Smithsonian Astrophysical Observatory
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
from distutils.version import LooseVersion
from numpy.testing import assert_allclose, assert_array_equal
from sherpa.astro import ui
from sherpa.utils import SherpaTestCase
from sherpa.utils import requires_data, requires_fits, requires_xspec

try:
    from sherpa.astro.xspec import get_xsversion
    GREATER_THAN_120901 = LooseVersion(get_xsversion()) > LooseVersion('12.9.1')
except:
    GREATER_THAN_120901 = True # on error, skip it anyway


# How many models should there be?
XSPEC_MODELS_COUNT = 166

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
def get_xspec_models(xs):
    """What are the XSpec model names to test."""

    # The alternate approach is to use that taken by
    # test_create_model_instances
    #
    model_names = [model_name for model_name in dir(xs) if model_name.startswith('XS')]

    # Could just exclude any names that end in 'Model', but this
    # could remove valid model names, so be explicit.
    remove_item(model_names, 'XSModel')
    remove_item(model_names, 'XSMultiplicativeModel')
    remove_item(model_names, 'XSAdditiveModel')
    remove_item(model_names, 'XSTableModel')

    # The sirf model - in 12.8.2 and up to 12.9.0d at least - includes
    # a read outside of an array. This has been seen to cause occasional
    # errors in the Sherpa test case, so it is removed from the test
    # for now. This problem has been reported to the XSPEC developers,
    # so it will hopefully be fixed in one of ther 12.9.0 patches.
    remove_item(model_names, 'XSsirf')

    models = list(map(lambda model_name: getattr(xs, model_name), model_names))
    models = [model for model in models if model.version_enabled]

    return models


def make_grid():
    """Return the 'standard' contiguous grid used in these tests.

    Returns egrid, elo, ehi, wgrid, wlo, whi where the e-xxx values
    are in keV and the w-xxx values in Angstrom, *grid are a
    single grid of values (i.e. one array) and the *lo/*hi values
    are this grid separated out into bin edges.

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

    egrid = numpy.arange(0.1, 5.01, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]

    wgrid = _hc / egrid
    whi = wgrid[:-1]
    wlo = wgrid[1:]

    return egrid, elo, ehi, wgrid, wlo, whi


def make_grid_noncontig2():
    """Return the 'standard' non-contiguous grid used in these tests,
    with 2 gaps.

    The grid from make_grid is modified to make one gap.
    """

    egrid, elo, ehi, wgrid, wlo, whi = make_grid()

    # remove two sections
    idx1 = (elo <= 1.1) | (elo >= 1.4)
    idx2 = (elo <= 2.3) | (elo >= 2.9)
    idx = idx1 & idx2

    # at the moment do not return the original grid and filter
    return elo[idx], ehi[idx], wlo[idx], whi[idx]


@requires_xspec
class test_xspec(SherpaTestCase):

    def setUp(self):
        from sherpa.astro import xspec
        ui.clean()
        self._defaults = xspec.get_xsstate()

    def tearDown(self):
        from sherpa.astro import xspec
        xspec.set_xsstate(self._defaults)

    # In Sherpa 4.8.0 development, the XSPEC model class included
    # an explicit check that all the bins were finite. This has
    # been removed as testing has shown there are complications when
    # using this with some real-world responses. So add in an
    # explicit - hopefully temporary - test here.
    #
    def assertFinite(self, vals, model, label):
        emsg = "model {} is finite [{}]".format(model, label)
        self.assertTrue(numpy.isfinite(vals).all(),
                        msg=emsg)

    def test_create_model_instances(self):
        import sherpa.astro.xspec as xs
        count = 0

        for cls in dir(xs):
            if not cls.startswith('XS'):
                continue

            cls = getattr(xs, cls)

            if is_proper_subclass(cls, (xs.XSAdditiveModel,
                                        xs.XSMultiplicativeModel)):
                # Ensure that we can create an instance, but do
                # nothing with it.
                cls()
                count += 1

        self.assertEqual(count, XSPEC_MODELS_COUNT)

    def test_norm_works(self):
        # Check that the norm parameter for additive models
        # works, as it is handled separately from the other
        # parameters.
        import sherpa.astro.xspec as xs

        # need an additive model
        mdl = xs.XSpowerlaw()
        mdl.PhoIndex = 2
        egrid = [0.1, 0.2, 0.3, 0.4]

        mdl.norm = 1.2
        y1 = mdl(egrid)

        mfactor = 2.1
        mdl.norm = mdl.norm.val * mfactor
        y2 = mdl(egrid)

        # check that the sum is not 0 and that it
        # scales as expected.
        s1 = y1.sum()
        s2 = y2.sum()
        self.assertGreater(s1, 0.0, msg='powerlaw is positive')
        self.assertAlmostEqual(s2, mfactor * s1,
                               msg='powerlaw norm scaling')

    def test_evaluate_model(self):
        import sherpa.astro.xspec as xs
        mdl = xs.XSbbody()
        out = mdl([1, 2, 3, 4])
        if mdl.calc.__name__.startswith('C_'):
            otype = numpy.float64
        else:
            otype = numpy.float32
        self.assertEqual(out.dtype.type, otype)
        self.assertEqual(int(numpy.flatnonzero(out == 0.0)), 3)

    def test_checks_input_length(self):
        import sherpa.astro.xspec as xs
        mdl = xs.XSpowerlaw()

        # Check when input array is too small (< 2 elements)
        self.assertRaises(TypeError, mdl, [0.1])

        # Check when input arrays are not the same size (when the
        # low and high bin edges are given)
        self.assertRaises(TypeError, mdl, [0.1, 0.2, 0.3], [0.2, 0.3])
        self.assertRaises(TypeError, mdl, [0.1, 0.2], [0.2, 0.3, 0.4])

    def test_xspec_models(self):
        import sherpa.astro.xspec as xs
        models = get_xspec_models(xs)

        egrid, elo, ehi, wgrid, wlo, whi = make_grid()
        for model in models:
            # use an identifier in case there is an error
            mdl = model()

            # The model checks that the values are all finite,
            # so there is no need to check that the output of
            # mdl does not contain non-finite values.
            # NOTE: this is no-longer the case, so include an
            #       explicit check for the single-argument forms
            #
            evals1 = mdl(egrid)
            evals2 = mdl(elo, ehi)

            wvals1 = mdl(wgrid)
            wvals2 = mdl(wlo, whi)

            self.assertFinite(evals1, model, "energy")
            self.assertFinite(wvals1, model, "wavelength")

            emsg = "{} model evaluation failed: ".format(model)

            # It might be expected that the test should be
            #   assert_allclose(evals1[:-1], evals2)
            # to ensure there's no floating-point issues,
            # but in this case the grid and parameter
            # values *should* be exactly the same, so the
            # results *should* be exactly equal, hence
            # the use of assert_array_equal
            assert_array_equal(evals1[:-1], evals2,
                               err_msg=emsg + "energy comparison")
            assert_array_equal(wvals1[:-1], wvals2,
                               err_msg=emsg + "wavelength comparison")

            # When comparing wavelength to energy values, have
            # to use allclose since the bins are not identically
            # equal.
            assert_allclose(evals1, wvals1,
                            err_msg=emsg + "energy to wavelength")

    def test_xspec_models_noncontiguous2(self):
        """Note that there is no test that the non-contiguous grid
        results are similar to the result from using a contiguous
        grid and then filtering out the missing bins. The way that
        some models are implemented make this a tricky test to
        write (due to numerical tolerances), as bins at the
        edges may not match well.
        """

        import sherpa.astro.xspec as xs
        models = get_xspec_models(xs)

        elo, ehi, wlo, whi = make_grid_noncontig2()
        for model in models:
            # use an identifier in case there is an error
            mdl = model()

            evals2 = mdl(elo, ehi)
            wvals2 = mdl(wlo, whi)

            self.assertFinite(evals2, model, "energy")
            self.assertFinite(wvals2, model, "wavelength")

            emsg = "{} non-contiguous model evaluation " + \
                "failed: ".format(model)
            assert_allclose(evals2, wvals2,
                            err_msg=emsg + "energy to wavelength")

    # Support for XSPEC support in load_table_model has been deprecated
    # in Sherpa 4.9.0; use load_xstable_model instead. For now the
    # tests are run with both functions, which means making the loading
    # function a parameter of the tests.
    #
    def _test_xspec_tablemodel_checks_input_length(self, loadfunc):

        loadfunc('mdl', self.make_path('xspec-tablemodel-RCS.mod'))
        mdl = ui.get_model_component('mdl')

        # Check when input array is too small (< 2 elements)
        self.assertRaises(TypeError, mdl, [0.1])

        # Check when input arrays are not the same size (when the
        # low and high bin edges are given)
        self.assertRaises(TypeError, mdl, [0.1, 0.2, 0.3], [0.2, 0.3])
        self.assertRaises(TypeError, mdl, [0.1, 0.2], [0.2, 0.3, 0.4])

    def _test_xspec_tablemodel(self, loadfunc):
        # Just test one table model; use the same scheme as
        # test_xspec_models_noncontiguous().
        #
        # The table model is from
        # https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/rcs.html
        # retrieved July 9 2015. The exact model is irrelevant for this
        # test, so this was chosen as it's relatively small.
        loadfunc('tmod', self.make_path('xspec-tablemodel-RCS.mod'))

        # when used in the test suite it appears that the tmod
        # global symbol is not created, so need to access the component
        tmod = ui.get_model_component('tmod')

        self.assertEqual(tmod.name, 'xstablemodel.tmod')

        egrid, elo, ehi, wgrid, wlo, whi = make_grid()

        evals1 = tmod(egrid)
        evals2 = tmod(elo, ehi)

        wvals1 = tmod(wgrid)
        wvals2 = tmod(wlo, whi)

        self.assertFinite(evals1, tmod, "energy")
        self.assertFinite(wvals1, tmod, "wavelength")

        emsg = "table model evaluation failed: "
        assert_array_equal(evals1[:-1], evals2,
                           err_msg=emsg + "energy comparison")

        assert_allclose(evals1, wvals1,
                        err_msg=emsg + "single arg")
        assert_allclose(evals2, wvals2,
                        err_msg=emsg + "two args")

    def _test_xspec_tablemodel_noncontiguous2(self, loadfunc):

        loadfunc('tmod', self.make_path('xspec-tablemodel-RCS.mod'))
        tmod = ui.get_model_component('tmod')

        elo, ehi, wlo, whi = make_grid_noncontig2()

        evals2 = tmod(elo, ehi)
        wvals2 = tmod(wlo, whi)

        self.assertFinite(evals2, tmod, "energy")
        self.assertFinite(wvals2, tmod, "wavelength")

        emsg = "table model non-contiguous evaluation failed: "
        rtol = 1e-3
        assert_allclose(evals2, wvals2, rtol=rtol,
                        err_msg=emsg + "energy to wavelength")

    @requires_data
    @requires_fits
    def test_xstablemodel_checks_input_length(self):
        self._test_xspec_tablemodel_checks_input_length(
            ui.load_xstable_model)

    @requires_data
    @requires_fits
    def test_tablemodel_checks_input_length(self):
        self._test_xspec_tablemodel_checks_input_length(
            ui.load_table_model)

    @requires_data
    @requires_fits
    def test_xspec_xstablemodel(self):
        self._test_xspec_tablemodel(ui.load_xstable_model)

    @requires_data
    @requires_fits
    def test_xspec_tablemodel(self):
        self._test_xspec_tablemodel(ui.load_table_model)

    @requires_data
    @requires_fits
    def test_xspec_xstablemodel_noncontiguous2(self):
        self._test_xspec_tablemodel_noncontiguous2(ui.load_xstable_model)

    @requires_data
    @requires_fits
    def test_xspec_tablemodel_noncontiguous2(self):
        self._test_xspec_tablemodel_noncontiguous2(ui.load_table_model)

    def test_convolution_model_cflux(self):
        # Use the cflux convolution model, since this gives
        # an easily-checked result. At present the only
        # interface to these models is via direct access
        # to the functions (i.e. there are no model classes
        # providing access to this functionality).
        #
        import sherpa.astro.xspec as xs

        if not hasattr(xs._xspec, 'C_cflux'):
            self.skipTest('cflux convolution model is missing')

        # The energy grid should extend beyond the energy grid
        # used to evaluate the model, to avoid any edge effects.
        # It also makes things easier if the elo/ehi values align
        # with the egrid bins.
        elo = 0.55
        ehi = 1.45
        egrid = numpy.linspace(0.5, 1.5, 101)

        mdl1 = xs.XSpowerlaw()
        mdl1.PhoIndex = 2

        # flux of mdl1 over the energy range of interest; converting
        # from a flux in photon/cm^2/s to erg/cm^2/s, when the
        # energy grid is in keV.
        y1 = mdl1(egrid)
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
        y2 = xs._xspec.C_cflux(pars, y1, egrid)

        self.assertFinite(y2, "cflux", "energy")

        elo = egrid[:-1]
        ehi = egrid[1:]
        wgrid = _hc / egrid
        whi = wgrid[:-1]
        wlo = wgrid[1:]

        expected = y1 * 10**lflux / f1
        numpy.testing.assert_allclose(expected, y2,
                                      err_msg='energy, single grid')

        y1 = mdl1(wgrid)
        y2 = xs._xspec.C_cflux(pars, y1, wgrid)
        self.assertFinite(y2, "cflux", "wavelength")
        numpy.testing.assert_allclose(expected, y2,
                                      err_msg='wavelength, single grid')

        expected = expected[:-1]

        y1 = mdl1(elo, ehi)
        y2 = xs._xspec.C_cflux(pars, y1, elo, ehi)
        numpy.testing.assert_allclose(expected, y2,
                                      err_msg='energy, two arrays')

        y1 = mdl1(wlo, whi)
        y2 = xs._xspec.C_cflux(pars, y1, wlo, whi)
        numpy.testing.assert_allclose(expected, y2,
                                      err_msg='wavelength, two arrays')

    def test_convolution_model_cpflux_noncontiguous(self):
        # The models should raise an error if given a non-contiguous
        # grid.
        import sherpa.astro.xspec as xs

        if not hasattr(xs._xspec, 'C_cpflux'):
            self.skipTest('cpflux convolution model is missing')

        elo, ehi, wlo, whi = make_grid_noncontig2()

        lflux = -5.0
        pars = [0.2, 0.8, lflux]
        y1 = numpy.zeros(elo.size)

        self.assertRaises(ValueError, xs._xspec.C_cpflux, pars, y1, elo, ehi)
        self.assertRaises(ValueError, xs._xspec.C_cpflux, pars, y1, wlo, whi)

    @requires_data
    @requires_fits
    def test_set_analysis_wave_fabrizio(self):
        rmf = self.make_path('3c273.rmf')
        arf = self.make_path('3c273.arf')

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

        self.assertAlmostEqual(y_m, y2_m)

    def test_xsxset_get(self):
        import sherpa.astro.xspec as xs
        # TEST CASE #1 Case insentitive keys
        xs.set_xsxset('fooBar', 'somevalue')
        self.assertEqual('somevalue', xs.get_xsxset('Foobar'))


@requires_xspec
@pytest.mark.skipif(GREATER_THAN_120901, reason="This test only makes sense with earlier versions of XSPEC")
def test_nonexistent_model():
    from sherpa.models import Parameter
    from sherpa.astro.xspec.utils import version_at_least
    from sherpa.astro.xspec import XSAdditiveModel

    @version_at_least("12.9.1")
    class XSbtapec(XSAdditiveModel):
        __function__ = "C_btapec"

        def __init__(self, name='btapec'):
            self.kT = Parameter(name, 'kT', 1.0)
            XSAdditiveModel.__init__(self, name, (self.kT))

    m = XSbtapec()

    with pytest.raises(AttributeError) as exc:
        m([])

    assert version_at_least.DISABLED_MODEL_MESSAGE.format("XSbtapec") == str(exc.value)


@requires_xspec
def test_existent_model():
    from sherpa.astro.xspec import XSzbabs

    m = XSzbabs()

    m([1, 2, 3])


@requires_xspec
def test_not_compiled_model():
    """
    Test the error handling case where a model is included according to the conditional decorator, but it wraps a
    function that had not been compiled.
    """
    from sherpa.models import Parameter
    from sherpa.astro.xspec.utils import include_if, ModelMeta
    from sherpa.astro.xspec import get_xsversion, XSAdditiveModel

    @include_if(get_xsversion() >= "0.0.0")  # Use the lowest version possible, so this is always true.
    class XSfoo(XSAdditiveModel):
        __function__ = "C_foo"

        def __init__(self, name='foo'):
            self.kT = Parameter(name, 'kT', 1.0)
            XSAdditiveModel.__init__(self, name, (self.kT))

    m = XSfoo()

    with pytest.raises(AttributeError) as exc:
        m([])

    assert ModelMeta.NOT_COMPILED_FUNCTION_MESSAGE == str(exc.value)


@requires_xspec
def test_old_style_xspec_class():
    """
    We changed the way xspec models are declared, but just in case let's make sure old-style declarations still work.
    """
    from sherpa.models import Parameter
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

    actual = m([1, 2, 3])

    expected = XSzbabs()([1, 2, 3])

    assert_array_equal(expected, actual)
