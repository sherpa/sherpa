#
#  Copyright (C) 2015  Smithsonian Astrophysical Observatory
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
Provide tests for the sherpa.astro.io routines.

These are aimed at low-level tests of the API, covering parts that
are not handled by existing higher-level tests (i.e. of those parts
of the API that provide wrappers around the IO backend code).

"""

import unittest
import tempfile

import numpy
from numpy.testing import assert_allclose

from sherpa.utils import SherpaTestCase, SherpaTest, has_fits_support
from sherpa.astro import io
from sherpa.data import Data1D


class WriteArrays(SherpaTestCase):
    """io.write_arrays is not exposed to the user, but is used by
    multiple routines, so test it out.
    """

    # Technically do not need _fits and _read_func since
    # they are a paired together,but separate out for now.
    #
    # _colnames:  are column names given to write_arrays or is it
    #             to calculate them?
    # _fits:      use FITS format for output?

    _colnames = False
    _fits = None

    def _read_func(self, filename, ncols):
        "Only expose the ncols argument"
        raise NotImplementedError

    def write_arrays(self):
        """Write out a small set of data using io.write_arrays
        and then read it back in, to check it was written out
        correctly (or, at least, in a way that can be read back
        in).
        """

        # It looks like the input arrays to `write_arrays` should be numpy
        # arrays, so enforce that invariant.
        a = numpy.asarray([1, 3, 9])
        b = numpy.sqrt(numpy.asarray(a))
        c = b * 0.1

        if self._colnames:
            fields = ["x", "yy", "z"]
        else:
            fields = None

        ofh = tempfile.NamedTemporaryFile(suffix='sherpa_test')
        io.write_arrays(ofh.name, [a, b, c], fields=fields,
                        ascii=not self._fits, clobber=True)

        out = self._read_func(ofh.name, ncols=3)

        rtol = 0
        atol = 1e-5
        self.assertIsInstance(out, Data1D)
        self.assertEqual(out.name, ofh.name, msg="file name")
        assert_allclose(out.x, a, rtol=rtol, atol=atol, err_msg="x column")
        assert_allclose(out.y, b, rtol=rtol, atol=atol, err_msg="y column")
        assert_allclose(out.staterror, c, rtol=rtol, atol=atol,
                        err_msg="staterror")
        self.assertIsNone(out.syserror, msg="syserror")


@unittest.skipIf(not has_fits_support(),
                 'need pycrates, astropy.io.fits, or pyfits')
class TestWriteArraysNoColsFITS(WriteArrays):

    _fits = True

    def _read_func(self, filename, ncols):
        return io.read_table(filename, ncols=ncols)

    def test_write_arrays(self):
        self.write_arrays()


@unittest.skipIf(not has_fits_support(),
                 'need pycrates, astropy.io.fits, or pyfits')
class TestWriteArraysColsFITS(TestWriteArraysNoColsFITS):

    _colnames = True

    def test_write_arrays(self):
        self.write_arrays()


@unittest.skipIf(not has_fits_support(),
                 'need pycrates, astropy.io.fits, or pyfits')
class TestWriteArraysNoColsASCII(WriteArrays):

    _fits = False

    def _read_func(self, filename, ncols):
        return io.read_ascii(filename, ncols=ncols)

    def test_write_arrays(self):
        self.write_arrays()


@unittest.skipIf(not has_fits_support(),
                 'need pycrates, astropy.io.fits, or pyfits')
class TestWriteArraysColsASCII(TestWriteArraysNoColsASCII):

    _colnames = True

    def test_write_arrays(self):
        self.write_arrays()


if __name__ == '__main__':
    SherpaTest(io).test()
