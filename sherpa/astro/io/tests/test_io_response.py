#
#  Copyright (C) 2021, 2023
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

"""Check that we can read in and write out PHA response files.

"""

import numpy as np

import pytest

from sherpa.astro import io
from sherpa.astro.data import DataARF, DataRMF
from sherpa.data import Data1DInt
from sherpa.utils.testing import requires_data, requires_fits


def backend_is(name):
    """Are we using the specified backend?"""
    return io.backend.__name__ == f"sherpa.astro.io.{name}_backend"


@requires_data
@requires_fits
def test_read_arf(make_data_path):
    """Read in a Chandra grating ARF"""

    infile = make_data_path("3c120_meg_1.arf.gz")
    arf = io.read_arf(infile)
    assert isinstance(arf, DataARF)
    assert arf.name.endswith("/3c120_meg_1.arf.gz")
    assert arf.exposure == pytest.approx(77715.64389196)
    assert np.log10(arf.ethresh) == pytest.approx(-10)

    hdr = arf.header
    assert "HDUNAME" not in hdr
    assert hdr["CONTENT"] == "SPECRESP"
    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "SPECRESP"
    assert hdr["HDUVERS"] == "1.1.0"
    assert hdr["HDUVERS1"] == "1.0.0"
    assert hdr["HDUVERS2"] == "1.1.0"
    assert hdr["ARFVERSN"] == "1992a"

    assert hdr["MISSION"] == "AXAF"
    assert hdr["TELESCOP"] == "CHANDRA"
    assert hdr["INSTRUME"] == "ACIS"
    assert hdr["GRATING"] == "HETG"
    assert hdr["FILTER"] == ""

    assert hdr["OBJECT"] == "3C 120"
    assert hdr["RESPFILE"] == "grid(meg_1.rmf)"

    assert hdr["OBS_ID"] == "16221"
    assert hdr["OBI_NUM"] == 0
    assert hdr["CTI_CORR"]
    assert hdr["DTCOR"] == pytest.approx(0.98318749385508,)

    for field in ["energ_lo", "energ_hi", "specresp", "bin_lo", "bin_hi"]:
        attr = getattr(arf, field)
        assert len(attr) == 8192
        assert attr.dtype == np.float64

    assert (arf.energ_lo[1:] == arf.energ_hi[:-1]).all()
    assert arf.energ_lo[0] == pytest.approx(0.29548186064)
    assert arf.energ_hi[-1] == pytest.approx(12.398418427)

    assert arf.specresp.sum() == pytest.approx(67517.7296)
    assert np.argmax(arf.specresp) == 7041

    assert (arf.bin_lo[:-1] == arf.bin_hi[1:]).all()
    assert arf.bin_lo[-1] == pytest.approx(1.0)
    assert arf.bin_hi[0] == pytest.approx(41.959999084)


@requires_data
@requires_fits
def test_read_rmf(make_data_path):
    """Read in an ASCA RMF with slightly different structure to Chandra

    This has F_CHAN[2]/N_CHAN[2] rather than scalar values.
    """

    infile = make_data_path("s0_mar24.rmf")
    rmf = io.read_rmf(infile)
    print(rmf)
    assert isinstance(rmf, DataRMF)
    assert rmf.name.endswith("/s0_mar24.rmf")
    assert rmf.detchans == 512
    assert rmf.offset == 0
    assert np.log10(rmf.ethresh) == pytest.approx(-10)

    hdr = rmf.header
    assert "HDUNAME" not in hdr
    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "RSP_MATRIX"
    assert hdr["HDUCLAS3"] == "DETECTOR"
    assert hdr["HDUVERS1"] == "1.0.0"
    assert hdr["HDUVERS2"] == "1.1.0"
    assert hdr["RMFVERSN"] == "1992a"

    assert hdr["TELESCOP"] == "ASCA"
    assert hdr["INSTRUME"] == "SIS0"
    assert hdr["FILTER"] == "NONE"

    assert hdr["CHANTYPE"] == "PI"
    assert np.log10(hdr["LO_THRES"]) == pytest.approx(-7)

    assert len(rmf.energ_lo) == 1180
    assert len(rmf.energ_hi) == 1180
    assert len(rmf.n_grp) == 1180

    assert len(rmf.f_chan) == 1189
    assert len(rmf.n_chan) == 1189

    assert len(rmf.matrix) == 260775

    assert len(rmf.e_min) == 512
    assert len(rmf.e_max) == 512

    assert rmf.n_grp.dtype == np.uint64

    for field in ["f_chan", "n_chan"]:
        attr = getattr(rmf, field)
        if backend_is("crates"):
            assert attr.dtype == np.uint32
        elif backend_is("pyfits"):
            assert attr.dtype == np.uint64
        else:
            raise RuntimeError(f"unsupported I/O backend: {io.backend}")

    for field in ["energ_lo", "energ_hi", "matrix", "e_min", "e_max"]:
        attr = getattr(rmf, field)
        assert attr.dtype == np.float64

    assert (rmf.energ_lo[1:] == rmf.energ_hi[:-1]).all()
    assert rmf.energ_lo[0] == pytest.approx(0.20000000298)
    assert rmf.energ_hi[-1] == pytest.approx(12)

    assert rmf.n_grp.sum() == 1189
    assert np.argmax(rmf.n_grp) == 12

    # It is not obvious how the RMF is flattened here, so treat this
    # as a regression test rather than "from first principles"
    #
    assert rmf.f_chan.sum() == 16711
    assert rmf.n_chan.sum() == 260775

    assert rmf.matrix.sum() == pytest.approx(544.8837793875416)

    # e_min/max come from the EBOUNDS block
    #
    assert (rmf.e_min[1:] == rmf.e_max[:-1]).all()
    assert rmf.e_min[0] == pytest.approx(0.01329296827316)
    assert rmf.e_max[-1] == pytest.approx(14.678317070)


TEST_ARF_ELO_SUM = 6042.19
TEST_ARF_EHI_SUM = 6052.97
TEST_ARF_SPECRESP_SUM = 251968.38017

@requires_data
@requires_fits
def test_write_arf_fits(make_data_path, tmp_path):
    """Check we can write out an ARF as a FITS file.

    We check we can read-back the ARF we have written.
    We use a different ARF to that used in test_read_arf.
    """

    infile = make_data_path("9774.arf")
    orig = io.read_arf(infile)

    outpath = tmp_path / "out.arf"
    outfile = str(outpath)
    io.write_arf(outfile, orig, ascii=False)

    new = io.read_arf(outfile)
    assert isinstance(new, DataARF)
    assert new.exposure == pytest.approx(75141.23109910)
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "SPECRESP"
    assert hdr["HDUVERS"] == "1.1.0"
    assert hdr["TELESCOP"] == "CHANDRA"
    assert hdr["INSTRUME"] == "ACIS"
    assert hdr["FILTER"] == "NONE"

    assert new.energ_lo.sum() == pytest.approx(TEST_ARF_ELO_SUM)
    assert new.energ_hi.sum() == pytest.approx(TEST_ARF_EHI_SUM)
    assert new.specresp.sum() == pytest.approx(TEST_ARF_SPECRESP_SUM)


@requires_data
@requires_fits
def test_write_arf_ascii(make_data_path, tmp_path):
    """Check we can write out an ARF as an ASCII file.

    We check we can read-back the ARF we have written.
    We use a different ARF to that used in test_read_arf.
    """

    infile = make_data_path("9774.arf")
    orig = io.read_arf(infile)

    outpath = tmp_path / "out.arf"
    outfile = str(outpath)
    io.write_arf(outfile, orig, ascii=True)

    new = io.read_ascii(outfile, ncols=3, dstype=Data1DInt,
                        colkeys=["ENERG_LO", "ENERG_HI", "SPECRESP"])

    assert new.xlo.sum() == pytest.approx(TEST_ARF_ELO_SUM)
    assert new.xhi.sum() == pytest.approx(TEST_ARF_EHI_SUM)
    assert new.y.sum() == pytest.approx(TEST_ARF_SPECRESP_SUM)
