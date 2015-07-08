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
import numpy.testing # this was added numpy 1.5.0
from sherpa.astro import ui
from sherpa.utils import SherpaTestCase, test_data_missing
from sherpa.utils import has_package_from_list, has_fits_support

def is_proper_subclass(obj, cls):
    if type(cls) is not tuple:
        cls = (cls,)
    if obj in cls:
        return False
    return issubclass(obj, cls)

def _get_xspec_models(xs):
    """What are the XSpec model names to test."""

    models = [model for model in dir(xs) if model.startswith('XS')]
    models.remove('XSModel')
    models.remove('XSMultiplicativeModel')
    models.remove('XSAdditiveModel')
    models.remove('XSTableModel')

    # The nteea model is known to be broken in XSpec 12.8.2q
    # (it appears to work in XSpec, but does not in Sherpa due to
    # differences in how the models are evaluated). The broken
    # behavior is that a second call to evaluate the model without
    # changing the parameter values - as done in the test below -
    # will fail, returning 0's for all bin values.
    #
    # The model removal should be made contingent on the XSpec
    # library version - see xs.get_xsversion() - once a fix has
    # been made.
    #
    # Doug Burke, July 8, 2015
    models.remove('XSnteea')

    return models

def _make_noncontiguous_grid():
    """Return the grids for test_xspec_models_noncontiguous()."""

    egrid = numpy.arange(0.1, 11.01, 0.01, dtype=float)

    def grid(idx):
        """Return lo,hi values for these indexes."""
        return (egrid[idx][:-1], egrid[idx][1:])

    # gaps are somewhat-randomly chosen.
    idx1 = numpy.arange(0, 40)
    idx2 = numpy.arange(41, 458)
    idx3 = numpy.arange(499, 802)
    idx4 = numpy.arange(917, egrid.size)

    xlos = []
    xhis = []
    idxs = []
    for idx in [idx1, idx2, idx3, idx4]:
        (xlo, xhi) = grid(idx)
        xlos.append(xlo)
        xhis.append(xhi)
        idxs.append(idx[:-1])

    xlo = numpy.concatenate(xlos)
    xhi = numpy.concatenate(xhis)
    idx = numpy.concatenate(idxs)

    # indexes where differences have been seen.
    # TODO: this should be generated from the idxs array
    # rather than hard coded.
    #
    ends = [38, 454, 756]

    # What are the indexes for those points where we expect the
    # comparison to be "exact" and the end-points, which are
    # known to be problematic (i.e. have larger variation).
    #
    gidx = numpy.asarray([r for r in range(0, idx.size) if r not in ends])
    bidx = numpy.asarray(ends)

    return (egrid, xlo, xhi, idx, gidx, bidx)

@unittest.skipIf(not has_package_from_list('sherpa.astro.xspec'),
                 "required sherpa.astro.xspec module missing")
class test_xspec(SherpaTestCase):

    def test_create_model_instances(self):
        import sherpa.astro.xspec as xs
        count = 0

        for cls in dir(xs):
            if not cls.startswith('XS'):
                continue

            cls = getattr(xs, cls)

            if is_proper_subclass(cls, (xs.XSAdditiveModel,
                                        xs.XSMultiplicativeModel)):
                m = cls()
                count += 1

        self.assertEqual(count, 164)

    def test_evaluate_model(self):
        import sherpa.astro.xspec as xs
        m = xs.XSbbody()
        out = m([1,2,3,4])
        if m.calc.__name__.startswith('C_'):
            otype = numpy.float64
        else:
            otype = numpy.float32
        self.assert_(out.dtype.type is otype)
        self.assertEqual(int(numpy.flatnonzero(out == 0.0)), 3)


    def test_xspec_models(self):
        import sherpa.astro.xspec as xs
        models = _get_xspec_models(xs)

        # There are two ways to call an XSpec model; single
        # argument taken as a continuous grid or with two arguments
        # giving the low and high values of each bin. Check that they
        # give the same results.
        #
        egrid = numpy.arange(0.1, 11.01, 0.01, dtype=float)
        xlo = egrid[:-1]
        xhi = egrid[1:]
        for model in models:
            cls = getattr(xs, model)
            mdl = cls('foo')
            vals1 = mdl(egrid)
            vals2 = mdl(xlo,xhi)
            emsg = "{0} model evaluation failed [".format(model)
            self.assertTrue(numpy.isfinite(vals1).all(), msg=emsg + "1]")
            self.assertTrue(numpy.isfinite(vals2).all(), msg=emsg + "2]")

            # Ideally the two should be exactly the same.
            #numpy.testing.assert_allclose(vals1[:-1], vals2,
            #                              err_msg=emsg + "comparison]")
            self.assertTrue((vals1[:-1] == vals2).all(),
                            msg=emsg + "comparison]")

    # It would make sense to just use this test, rather than do both
    # the contiguous test (test_xspec_models) as well as this one.
    # However, due to the number of models that are skipped in this
    # test it was decided to keep both tests.
    #
    def test_xspec_models_noncontiguous(self):
        import sherpa.astro.xspec as xs
        models = _get_xspec_models(xs)

        # What extra models should be skipped?
        # With XSpec 12.8.2q, models that call sumdem -
        # so thermal-plasma related - will crash when
        # called using the current Sherpa code - see
        # https://github.com/sherpa/sherpa/issues/62 -
        # so these are removed.
        #
        models.remove("XSapec")
        models.remove("XSbapec")
        models.remove("XSbvapec")
        models.remove("XSbvvapec")
        models.remove("XSequil")
        models.remove("XSgadem")
        models.remove("XSgnei")
        models.remove("XSnei")
        models.remove("XSrnei")
        models.remove("XSnpshock")
        models.remove("XSpshock")
        models.remove("XSsedov")
        models.remove("XSvapec")
        models.remove("XSvequil")
        models.remove("XSvgadem")
        models.remove("XSvgnei")
        models.remove("XSvnei")
        models.remove("XSvnpshock")
        models.remove("XSvpshock")
        models.remove("XSvrnei")
        models.remove("XSvsedov")
        models.remove("XSvvapec")
        models.remove("XSvvgnei")
        models.remove("XSvvnei")
        models.remove("XSvvnpshock")
        models.remove("XSvvpshock")
        models.remove("XSvvrnei")
        models.remove("XSvvsedov")

        # The kerrdisk model is known to show large differences
        # here (XSpec 12.8.2q; but it's likely down to how the model
        # is created) and so it is skipped here. See
        # https://gist.github.com/DougBurke/d9a0074489b6d1de108e
        # for an example.
        #
        models.remove("XSkerrdisk")

        (egrid, elo, ehi, idx, gidx, bidx) = _make_noncontiguous_grid()

        # Ideally the same tolerances could be used for all the models,
        # but for now apply some model-specific values.
        #
        # The default rtol value is 1e-7 and atol is 0; these are
        # used to allow a difference of
        #     atol + rtol * abs(desired value)
        #
        rtols = { 'XSlaor': 1e-5, 'XSlaor2': 1e-5, 'XSnsa': 1e-5 }
        atols = { 'XSnsa': 1e-11 }

        # For the "edge" bins, use a default rtol of 1e-4
        rtols_edges = { }
        atols_edges = { 'XSpexmon': 1e-5,
                        'XSposm': 1e-8,
                        'XSswind1': 1e-3
        }

        # Skip models for which we know (using 12.8.2q)
        # that the check is not "worth it" (i.e. that the
        # tolerance has to be so large as to make it
        # non-informative).
        #
        skip_ends = [ 'XScompbb', 'XSlaor', 'XSlaor2' ]

        for model in models:
            cls = getattr(xs, model)
            mdl = cls('foo')

            vals1 = mdl(egrid)
            vals2 = mdl(elo,ehi)
            emsg = "{0} model evaluation failed [noncontig; ".format(model)
            self.assertTrue(numpy.isfinite(vals1).all(), msg=emsg + "1]")
            self.assertTrue(numpy.isfinite(vals2).all(), msg=emsg + "2]")

            # filter down the "contiguous" version
            vals1 = vals1[idx]

            # Check the "non-edge" bins.
            #
            kwargs = {}
            try:
                kwargs['rtol'] = rtols[model]
            except KeyError:
                pass

            try:
                kwargs['atol'] = atols[model]
            except KeyError:
                pass

            numpy.testing.assert_allclose(vals1[gidx], vals2[gidx],
                                          err_msg=emsg + "comparison]",
                                          **kwargs)

            if model in skip_ends:
                continue

            # Check the "edge" bins.
            #
            kwargs = {}
            try:
                kwargs['rtol'] = rtols_edges[model]
            except KeyError:
                kwargs['rtol'] = 1e-4

            try:
                kwargs['atol'] = atols_edges[model]
            except KeyError:
                pass

            numpy.testing.assert_allclose(vals1[bidx], vals2[bidx],
                                          err_msg=emsg + "comparison, edges]",
                                          **kwargs)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits')
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_set_analysis_wave_fabrizio(self):
        rmf = self.datadir + '/ciao4.3/fabrizio/Data/3c273.rmf'
        arf = self.datadir + '/ciao4.3/fabrizio/Data/3c273.arf'

        ui.set_model("fabrizio", "xspowerlaw.p1")
        ui.fake_pha("fabrizio", arf, rmf, 10000)

        model = ui.get_model("fabrizio")
        bare_model, _ = ui._session._get_model_status("fabrizio")
        y = bare_model.calc([1,1], model.xlo, model.xhi)
        y_m = numpy.mean(y)

        ui.set_analysis("fabrizio","wave")

        model2 = ui.get_model("fabrizio")
        bare_model2, _ = ui._session._get_model_status("fabrizio")
        y2 = bare_model2.calc([1,1], model2.xlo, model2.xhi)
        y2_m = numpy.mean(y2)

        self.assertAlmostEqual(y_m, y2_m)

    def test_xsxset_get(self):
        import sherpa.astro.xspec as xs
        # TEST CASE #1 Case insentitive keys
        xs.set_xsxset('fooBar', 'somevalue')
        self.assertEqual('somevalue', xs.get_xsxset('Foobar'))


if __name__ == '__main__':
    import sherpa.astro.xspec as xs
    from sherpa.utils import SherpaTest
    SherpaTest(xs).test()
