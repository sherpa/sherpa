#
#  Copyright (C) 2016  Smithsonian Astrophysical Observatory
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

"""Test the TableModel class.
"""

import sys

from sherpa.utils import SherpaTest, SherpaTestCase
from sherpa.utils import linear_interp, nearest_interp, neville
from sherpa.utils.err import ModelErr

from sherpa.data import Data1D
from sherpa.models.basic import TableModel

import numpy as np
from numpy.testing import assert_allclose


class test_table_model(SherpaTestCase):

    # What interpolation functions should be used. Note that
    # the code is somewhat hard-coded to these values (e.g. for
    # selecting tolerances), so perhaps more, or less, abstraction
    # is needed.
    #
    interp1d = [linear_interp, nearest_interp, neville]

    def setUp(self):
        # These values match sherpa/ui/tests/data/gauss1d-error.dat
        dtype = np.float32
        self.x = np.asarray([80, 95, 110, 125, 140, 155, 170,
                             185, 200], dtype=dtype)
        self.y = np.asarray([0, 0, 9, 35, 93, 96, 49, 15, 0],
                            dtype=dtype)
        self.dy = np.asarray([1.86603, 1.86603, 4.1225, 6.97913,
                              10.6825, 10.8362, 8.05337, 4.96863,
                              1.86603], dtype=dtype)

    def _make(self, x, y, name='tbl', filename='dummy', method=None):
        """Create an instance"""

        mdl = TableModel(name)
        mdl.filename = filename
        if method is not None:
            mdl.method = method
        mdl.load(x, y)
        return mdl

    def test_basic(self):
        """Check that a table can be created and has expected properties"""

        mdl = self._make(self.x, self.y)

        # Just check that data in == data out, but allow for
        # data type changes, so allow for a tolerance.
        assert_allclose(self.x, mdl.get_x())
        assert_allclose(self.y, mdl.get_y())

        pars = mdl.pars
        self.assertEqual(1, len(pars), msg='One parameter')

        par = pars[0]
        self.assertEqual("ampl", par.name)
        self.assertEqual(False, par.frozen, msg='Is ampl frozen')

        self.assertAlmostEqual(1.0, par.val)

    def test_unsorted(self):
        """Check what happens when x array is not sorted."""

        mdl = self._make(self.y, self.dy, filename='unsorted')

        idx = np.argsort(self.y)
        assert_allclose(self.y[idx], mdl.get_x())
        assert_allclose(self.dy[idx], mdl.get_y())

    def test_fail_if_rows_do_not_match(self):
        """Model evaluation fails if #rows does not match (no X grid)."""

        mdl = self._make(None, self.y, name='oned')
        self.assertRaises(ModelErr, mdl, self.x[1:], self.y[1:])

    def test_fail_if_filter_mismatch(self):
        """Model evaluation fails if #rows does not match (dataset)."""

        mdl = self._make(None, self.y, name='oned')
        dset = Data1D('oned', self.x[1:], self.y[1:])
        # ignore the first point (this is to make sure that there is
        # a mask, so as to trigger the check in fold).
        dset.notice(self.x[1], self.x[-1])
        self.assertRaises(ModelErr, mdl.fold, dset)

    def test_fail_if_not_a_function(self):
        """What happens if the interpolator is not a function?"""

        # The NoNewAttributesAfterInit base class should raise an
        # AttributeError if a callable attribute is replaced by
        # a non-callable. This is left here in case there are changes
        # to that logic.
        self.assertRaises(AttributeError, self._make, self.x, self.y,
                          name='fail', method=True)

    def test_eval1(self):
        """Evaluate the model when no X column is given."""

        nval = 2.1e6
        m = self._make(None, self.y, filename='eval1')
        m.ampl = nval
        for interp in self.interp1d:
            m.method = interp
            ymdl = m(self.x)
            assert_allclose(self.y * nval, ymdl)

    def test_eval2_exact(self):
        """Evaluate the model when the X values of the data set
        match those of the table."""

        nval = 2.1e6
        m = self._make(self.x, self.y, filename='eval2_exact')
        m.ampl = nval
        for interp in self.interp1d:
            m.method = interp
            ymdl = m(self.x)
            assert_allclose(self.y * nval, ymdl)

    def test_eval2_interp(self):
        """Evaluate the model when the X values of the data set
        do not match those of the table."""

        # Pick some values outside the data range; the expected values
        # are based on a visual check of the interpolation results,
        # and so are a regression test. The first/last points of
        # nevile are exclued from the final comparison, since the
        # tolerance on the comparison on these values is not really
        # well defined.
        #
        xvals = [75, 100, 130, 150, 190, 210]
        yexp = {}
        yexp[linear_interp] = [0, 3, 54.3, 95, 10, -10]
        yexp[nearest_interp] = [0, 0, 35, 96, 15, 0]
        yexp[neville] = [164.6, 9.5, 55.3, 103.6, 5.3, 199.6]

        m = self._make(self.x, self.y, filename='eval2_interp')
        for interp in self.interp1d:
            m.method = interp
            ydiff = m(xvals) - yexp[interp]
            if interp == neville:
                ydiff = ydiff[1:-1]

            ymax = np.abs(ydiff).max()
            self.assertLess(ymax, 0.1)

    def test_eval2_filter(self):
        """Can we filter the data?

        It is not entirely clear how useful this test is here,
        since it tests out the interplay of several concepts
        revolving around Sherpa objects and state management,
        but it fits in with the other tests so leave here for
        now.
        """

        # 120 and 140 will be filtered out
        xvals = [75, 100, 120, 130, 140, 150, 190, 210]
        yexp = {}

        # These were calculated using the table model, and so
        # are regression tests (i.e. not calculated manually).
        yexp[linear_interp] = [0, 3, 26.3, 54.3, 93, 95, 10, -10]
        yexp[nearest_interp] = [0, 0, 35, 35, 93, 96, 15, 0]
        yexp[neville] = [164.6, 9.5, 19.3, 55.3, 93, 103.6, 5.3, 199.6]

        m = self._make(self.x, self.y, filename='eval2_filter')
        for interp in self.interp1d:
            m.method = interp

            dset = Data1D('filt', xvals, yexp[interp])
            dset.ignore(110, 125)
            dset.ignore(135, 149)

            m.fold(dset)

            yfull = dset.eval_model(m)
            yfilt = dset.eval_model_to_fit(m)

            y = np.asarray(yexp[interp])
            assert_allclose(y, yfull, atol=0.1)
            assert_allclose(y[dset.mask], yfilt, atol=0.1)

    def test_eval1_integrated(self):
        """Test handling with integrated data sets."""

        # This is based on test_eval2_integrated, and included
        # for completeness.
        #
        xgrid = np.asarray([10, 17, 20, 22, 31, 35, 40])
        xlo = xgrid[:-1]
        xhi = xgrid[1:]

        ytbl = np.asarray([-10, 10, 15, 20, 20, 10])
        ymdl = ytbl / 2.0

        tmdl = self._make(None, ytbl, filename='tmp')
        tmdl.ampl = 0.5
        for interp in self.interp1d:
            tmdl.method = interp
            ycalc = tmdl(xlo, xhi)
            assert_allclose(ymdl, ycalc)

    def test_eval2_integrated(self):
        """Test handling with integrated data sets."""

        # make sure the grid is not evenly spaced, since this
        # should show up any issues with trying to integrate
        # the model over the grid.
        #
        xgrid = np.asarray([10, 17, 20, 22, 31, 35, 40])
        xlo = xgrid[:-1]
        xhi = xgrid[1:]

        ytbl = np.asarray([-10, 10, 15, 20, 20, 10])
        ymdl = ytbl / 2.0

        tmdl = self._make(xlo, ytbl, filename='tmp')
        tmdl.ampl = 0.5
        for interp in self.interp1d:
            tmdl.method = interp
            ycalc = tmdl(xlo, xhi)
            assert_allclose(ymdl, ycalc)

if __name__ == '__main__':

    from sherpa import models

    if len(sys.argv) > 1:
        datadir = sys.argv[1]
    else:
        datadir = None

    SherpaTest(models).test(datadir=datadir)
