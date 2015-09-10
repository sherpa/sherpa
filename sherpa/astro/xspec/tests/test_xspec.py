#
#  Copyright (C) 2007, 2015  Smithsonian Astrophysical Observatory
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

import unittest
import numpy
from numpy.testing import assert_allclose, assert_array_equal
from sherpa.astro import ui
from sherpa.utils import SherpaTestCase, test_data_missing
from sherpa.utils import has_package_from_list, has_fits_support


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
    models = [model for model in dir(xs) if model.startswith('XS')]

    # Could just exclude any names that end in 'Model', but this
    # could remove valid model names, so be explicit.
    remove_item(models, 'XSModel')
    remove_item(models, 'XSMultiplicativeModel')
    remove_item(models, 'XSAdditiveModel')
    remove_item(models, 'XSTableModel')

    # In XSPEC 12.8.2, the nteea model is written in such a way that
    # it fails in Sherpa (but not from within XSPEC). This has
    # been fixed in 12.9.0, but can cause problems for Sherpa tests
    # when the model is evaluated multiple times.
    #
    version = [int(v) for v in xs.get_xsversion().split('.')[0:2]]
    if tuple(version) < (12, 9):
        remove_item(models, 'XSnteea')

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


@unittest.skipIf(not has_package_from_list('sherpa.astro.xspec'),
                 "required sherpa.astro.xspec module missing")
class test_xspec(SherpaTestCase):

    def setUp(self):
        from sherpa.astro import xspec
        self._defaults = xspec.get_xsstate()

    def tearDown(self):
        from sherpa.astro import xspec
        xspec.set_xsstate(self._defaults)

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

        self.assertEqual(count, 164)

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
        self.assertAlmostEqual(s2, mfactor * s1, msg='powerlaw norm scaling')

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

        self.assertRaises(TypeError, mdl, [0.1])

    def test_xspec_models(self):
        import sherpa.astro.xspec as xs
        models = get_xspec_models(xs)

        egrid, elo, ehi, wgrid, wlo, whi = make_grid()
        for model in models:
            cls = getattr(xs, model)
            mdl = cls(model)  # use an identifier in case there is an error

            # The model checks that the values are all finite,
            # so there is no need to check that the output of
            # mdl does not contain non-finite values.
            #
            evals1 = mdl(egrid)
            evals2 = mdl(elo, ehi)

            wvals1 = mdl(wgrid)
            wvals2 = mdl(wlo, whi)

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

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_xspec_tablemodel(self):
        # Just test one table model; use the same scheme as
        # test_xspec_models_noncontiguous().
        #
        # The table model is from
        # https://heasarc.gsfc.nasa.gov/xanadu/xspec/models/rcs.html
        # retrieved July 9 2015. The exact model is irrelevant for this
        # test, so this was chosen as it's relatively small.
        ui.load_table_model('tmod',
                            self.make_path('xspec/tablemodel/RCS.mod'))

        # when used in the test suite it appears that the tmod
        # global symbol is not created, so need to access the component
        tmod = ui.get_model_component('tmod')

        egrid, elo, ehi, wgrid, wlo, whi = make_grid()

        evals1 = tmod(egrid)
        evals2 = tmod(elo, ehi)

        # wvals1 = tmod(wgrid)
        # wvals2 = tmod(wlo, whi)

        # Direct comparison of wavelength and energy grid.
        #
        emsg = "table model evaluation failed: "
        rtol = 1e-3

        assert_allclose(evals1[:-1], evals2, rtol=rtol,
                        err_msg=emsg + "energy comparison")

        # assert_allclose(evals1, wvals1, rtol=rtol,
        #                 err_msg=emsg + "single arg")
        # assert_allclose(evals2, wvals2, rtol=rtol,
        #                 err_msg=emsg + "two args")

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, astropy.io.fits, or pyfits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_set_analysis_wave_fabrizio(self):
        rmf = self.make_path('ciao4.3/fabrizio/Data/3c273.rmf')
        arf = self.make_path('ciao4.3/fabrizio/Data/3c273.arf')

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


if __name__ == '__main__':
    import os
    import sys
    import sherpa.astro.xspec as xs
    from sherpa.utils import SherpaTest

    if len(sys.argv) > 1:
        datadir = sys.argv[1]
        if not os.path.exists(datadir):
            datadir = None
    else:
        datadir = None

    SherpaTest(xs).test(datadir=datadir)
