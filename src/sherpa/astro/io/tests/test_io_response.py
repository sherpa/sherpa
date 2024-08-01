#
#  Copyright (C) 2021, 2023, 2024
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

import logging
import warnings

import numpy as np

import pytest

from sherpa.astro.data import DataARF, DataPHA, DataRMF
from sherpa.astro.instrument import RMF1D, create_arf, create_delta_rmf
from sherpa.astro.instrument import create_delta_rmf
from sherpa.astro import io
from sherpa.astro.io.types import HeaderItem, Header, Column, TableBlock, BlockList
from sherpa.data import Data1DInt
from sherpa.models.basic import Gauss1D
from sherpa.utils.err import ArgumentErr, IOErr
from sherpa.utils.testing import requires_data, requires_fits


def backend_is(name: str) -> bool:
    """Are we using the specified backend?"""
    return io.backend.name == name


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
        assert attr.dtype == np.uint64

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

    assert new.bin_lo is None
    assert new.bin_hi is None


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


@requires_data
@requires_fits
def test_write_arf_fits_chandra_acis_hetg(make_data_path, tmp_path):
    """Check we can write out an ARF with BIN_LO/HI as a FITS file.
    """

    infile = make_data_path("3c120_meg_-1.arf.gz")
    orig = io.read_arf(infile)
    assert orig.exposure > 10
    assert orig.bin_lo is not None
    assert orig.bin_hi.size == orig.specresp.size

    outpath = tmp_path / "out.arf"
    outfile = str(outpath)
    io.write_arf(outfile, orig, ascii=False)

    new = io.read_arf(outfile)
    assert isinstance(new, DataARF)
    assert new.exposure == pytest.approx(orig.exposure)
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "SPECRESP"
    assert hdr["HDUVERS"] == "1.1.0"
    assert hdr["TELESCOP"] == "CHANDRA"
    assert hdr["INSTRUME"] == "ACIS"
    assert hdr["GRATING"] == "HETG"
    assert hdr["FILTER"] == ""
    assert hdr["TG_SRCID"] == 0
    assert hdr["TG_M"] == -1
    assert hdr["TG_PART"] == 2
    assert hdr["OBJECT"] == "3C 120"
    assert hdr["OBS_ID"] == "16221"  # note a string

    assert new.energ_lo == pytest.approx(orig.energ_lo)
    assert new.energ_hi == pytest.approx(orig.energ_hi)
    assert new.specresp == pytest.approx(orig.specresp)

    assert new.bin_lo == pytest.approx(orig.bin_lo)
    assert new.bin_hi == pytest.approx(orig.bin_hi)


@requires_data
@requires_fits
def test_write_arf_ascii_chandra_acis_hetg(make_data_path, tmp_path):
    """Check we can write out an ARF with BIN_LO/HI as an ASCII file.
    """

    infile = make_data_path("3c120_meg_-1.arf.gz")
    orig = io.read_arf(infile)

    outpath = tmp_path / "out.arf"
    outfile = str(outpath)
    io.write_arf(outfile, orig, ascii=True)

    new = io.read_ascii(outfile, ncols=3, dstype=Data1DInt,
                        colkeys=["BIN_LO", "BIN_HI", "SPECRESP"])

    assert new.xlo == pytest.approx(orig.bin_lo)
    assert new.xhi == pytest.approx(orig.bin_hi)
    assert new.y == pytest.approx(orig.specresp)


@requires_data
@requires_fits
def test_write_rmf_fits_chandra_acis(make_data_path, tmp_path):
    """Check we can write out an RMF as a FITS file.

    We check we can read-back the RMF we have written.
    We use a different RMF to that used in test_read_rmf.
    """

    NCHAN = 1024
    NENERGY = 1078
    NUMGRP = 1481
    NUMELT = 438482

    infile = make_data_path("9774.rmf")
    orig = io.read_rmf(infile)

    # Delete some keywords we want to check
    del orig.header["NUMELT"]
    del orig.header["NUMGRP"]
    store = [key for key in orig.header if key.startswith("HDU")]
    for key in store:
        del orig.header[key]

    outpath = tmp_path / "out.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, orig)

    new = io.read_rmf(outfile)
    assert isinstance(new, DataRMF)
    assert new.detchans == NCHAN
    assert new.offset == 1
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert "DETCHANS" not in hdr

    assert hdr["NUMGRP"] == NUMGRP
    assert hdr["NUMELT"] == NUMELT

    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "RSP_MATRIX"
    assert "HDUCLAS3" not in hdr
    assert "HDUVERS1" not in hdr
    assert "HDUVERS2" not in hdr
    assert hdr["HDUVERS"] == "1.3.0"
    assert hdr["RMFVERSN"] == "1992a"

    assert hdr["TELESCOP"] == "CHANDRA"
    assert hdr["INSTRUME"] == "ACIS"
    assert hdr["FILTER"] == "NONE"

    assert hdr["CHANTYPE"] == "PI"
    assert np.log10(hdr["LO_THRES"]) == pytest.approx(-5)

    assert len(new.energ_lo) == NENERGY
    assert len(new.energ_hi) == NENERGY
    assert len(new.n_grp) == NENERGY

    assert len(new.f_chan) == NUMGRP
    assert len(new.n_chan) == NUMGRP

    assert len(new.matrix) == NUMELT

    assert len(new.e_min) == NCHAN
    assert len(new.e_max) == NCHAN

    assert new.n_grp.dtype == np.uint64

    for field in ["f_chan", "n_chan"]:
        attr = getattr(new, field)
        assert attr.dtype == np.uint64

    for field in ["energ_lo", "energ_hi", "matrix", "e_min", "e_max"]:
        attr = getattr(new, field)
        assert attr.dtype == np.float64

    assert (new.energ_lo[1:] == new.energ_hi[:-1]).all()
    assert new.energ_lo[0] == pytest.approx(0.21999999880791)
    assert new.energ_hi[-1] == pytest.approx(11)

    assert new.n_grp.sum() == NUMGRP
    assert new.n_grp.max() == 3
    assert new.n_grp.min() == 1

    assert new.f_chan.sum() == 266767
    assert new.n_chan.sum() == NUMELT

    # Technically this should be NENERGY, since each row sums to 1 and
    # there's no missing rows, but we check the actual value.
    #
    assert new.matrix.sum() == pytest.approx(1078.000087318186)

    # e_min/max come from the EBOUNDS block
    #
    assert (new.e_min[1:] == new.e_max[:-1]).all()
    assert new.e_min[0] == pytest.approx(0.0014600000577047467)
    assert new.e_max[-1] == pytest.approx(14.950400352478027)


@requires_data
@requires_fits
def test_write_rmf_fits_asca_sis(make_data_path, tmp_path, caplog):
    """Check we can write out an RMF as a FITS file.

    Use the ASCA SIS as it has a different structure to
    the ACIS RMF - for instance, empty rows and each row only
    has one group (but it is not stored as an array, but a VLF).

    """

    NCHAN = 1024
    NENERGY = 1180
    NUMGRP = 1171
    NUMELT = 505483

    infile = make_data_path("sis0.rmf")
    assert len(caplog.record_tuples) == 0
    orig = io.read_rmf(infile)
    assert len(caplog.record_tuples) == 1

    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.astro.io"
    assert lvl == logging.ERROR
    assert msg.startswith("Failed to locate TLMIN keyword for F_CHAN column in RMF file '")
    assert msg.endswith("sis0.rmf'; Update the offset value in the RMF data set to the appropriate TLMIN value prior to fitting")

    # Unlike the ACIS case, this does not have NUMELT/GRP and
    # leave in the HDU keys to check they don't get over-written.
    #
    assert "NUMELT" not in orig.header
    assert "NUMGRP" not in orig.header
    assert "HDUVERS" not in orig.header

    outpath = tmp_path / "out.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, orig)

    new = io.read_rmf(outfile)
    assert isinstance(new, DataRMF)
    assert new.detchans == NCHAN
    assert new.offset == 1
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert "DETCHANS" not in hdr

    assert hdr["NUMGRP"] == NUMGRP
    assert hdr["NUMELT"] == NUMELT

    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "RSP_MATRIX"
    assert hdr["HDUCLAS3"] == "SPECRESP MATRIX"
    assert hdr["HDUVERS1"] == "1.0.0"
    assert hdr["HDUVERS2"] == "1.1.0"
    assert hdr["HDUVERS"] == "1.3.0"
    assert hdr["RMFVERSN"] == "1992a"

    assert hdr["TELESCOP"] == "ASCA"
    assert hdr["INSTRUME"] == "SIS0"
    assert hdr["FILTER"] == "NONE"

    assert hdr["CHANTYPE"] == "PI"
    assert np.log10(hdr["LO_THRES"]) == pytest.approx(-7)

    assert len(new.energ_lo) == NENERGY
    assert len(new.energ_hi) == NENERGY
    assert len(new.n_grp) == NENERGY

    assert len(new.f_chan) == NUMGRP
    assert len(new.n_chan) == NUMGRP

    assert len(new.matrix) == NUMELT

    assert len(new.e_min) == NCHAN
    assert len(new.e_max) == NCHAN

    assert new.n_grp.dtype == np.uint64

    for field in ["f_chan", "n_chan"]:
        attr = getattr(new, field)
        assert attr.dtype == np.uint64

    for field in ["energ_lo", "energ_hi", "matrix", "e_min", "e_max"]:
        attr = getattr(new, field)
        assert attr.dtype == np.float64

    assert (new.energ_lo[1:] == new.energ_hi[:-1]).all()
    assert new.energ_lo[0] == pytest.approx(0.20000000298023224)
    assert new.energ_hi[-1] == pytest.approx(12)

    assert new.n_grp.sum() == NUMGRP
    assert new.n_grp.max() == 1
    assert new.n_grp.min() == 0

    assert new.f_chan.sum() == 30446
    assert new.n_chan.sum() == NUMELT

    # Note that this is significantly less than NENERGY.
    #
    assert new.matrix.sum() == pytest.approx(552.932040504181)

    # e_min/max come from the EBOUNDS block
    #
    assert (new.e_min[1:] == new.e_max[:-1]).all()
    assert new.e_min[0] == pytest.approx(0.0)
    assert new.e_max[-1] == pytest.approx(14.240230560302734)


@requires_data
@requires_fits
def test_write_rmf_fits_swift_xrt(make_data_path, tmp_path):
    """Check we can write out an RMF as a FITS file.

    Use the Swift XRT as it has F/N_CHAN as fixed-length, not VLF,
    data.

    """

    NCHAN = 1024
    NENERGY = 2400
    NUMGRP = 2400
    NUMELT = NUMGRP * NCHAN

    infile = make_data_path("swxpc0to12s6_20130101v014.rmf.gz")
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter("always")
        orig = io.read_rmf(infile)

    assert len(ws) == 1
    msg = ws[0].message
    assert isinstance(msg, UserWarning)
    msg = str(msg)
    assert msg.startswith("The minimum ENERG_LO in the RMF '")
    assert msg.endswith(".rmf.gz' was 0 and has been replaced by 1e-10")

    assert "NUMELT" not in orig.header
    assert "NUMGRP" not in orig.header
    assert orig.header["HDUVERS"] == "1.3.0"
    assert orig.offset == 0

    outpath = tmp_path / "out.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, orig)

    # Note that we get no warning from this read as the min value has
    # been replaced by 1e-10.
    #
    with warnings.catch_warnings(record=True) as ws2:
        warnings.simplefilter("always")
        new = io.read_rmf(outfile)

    assert len(ws2) == 0

    assert isinstance(new, DataRMF)
    assert new.detchans == NCHAN
    assert new.offset == 0
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert "DETCHANS" not in hdr

    assert hdr["NUMGRP"] == NUMGRP
    assert hdr["NUMELT"] == NUMELT

    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "RSP_MATRIX"
    assert "HDUVERS1" not in hdr
    assert "HDUVERS2" not in hdr
    assert "HDUVERS3" not in hdr
    assert hdr["HDUVERS"] == "1.3.0"
    assert "RMFVERSN" not in hdr

    assert hdr["TELESCOP"] == "SWIFT"
    assert hdr["INSTRUME"] == "XRT"
    assert hdr["FILTER"] == "none"

    assert hdr["CHANTYPE"] == "PI"
    assert hdr["LO_THRES"] == 0

    assert len(new.energ_lo) == NENERGY
    assert len(new.energ_hi) == NENERGY
    assert len(new.n_grp) == NENERGY

    assert len(new.f_chan) == NUMGRP
    assert len(new.n_chan) == NUMGRP

    assert len(new.matrix) == NUMELT

    assert len(new.e_min) == NCHAN
    assert len(new.e_max) == NCHAN

    assert new.n_grp.dtype == np.uint64

    for field in ["f_chan", "n_chan"]:
        attr = getattr(new, field)
        assert attr.dtype == np.uint64

    for field in ["energ_lo", "energ_hi", "matrix", "e_min", "e_max"]:
        attr = getattr(new, field)
        assert attr.dtype == np.float64

    assert (new.energ_lo[1:] == new.energ_hi[:-1]).all()
    assert new.energ_lo[0] == pytest.approx(1e-10)  # 0 has been replaced
    assert new.energ_hi[-1] == pytest.approx(12)

    assert new.n_grp.sum() == NUMGRP
    assert new.n_grp.max() == 1
    assert new.n_grp.min() == 1

    assert new.f_chan.sum() == 0
    assert new.n_chan.sum() == NUMELT

    # Note that this is significantly less than NENERGY.
    #
    assert new.matrix.sum() == pytest.approx(1133.8128197300257)

    # e_min/max come from the EBOUNDS block
    #
    assert (new.e_min[1:] == new.e_max[:-1]).all()
    assert new.e_min[0] == pytest.approx(0.0)
    assert new.e_max[-1] == pytest.approx(10.239999771118164)


@requires_data
@requires_fits
def test_write_rmf_fits_rosat_pspc(make_data_path, tmp_path):
    """Check we can write out an RMF as a FITS file.

    Use the ROSAT PSPC-C example dataset as it has its own
    DataRosatRMF class.  It's also a RSP rather than a RMF (ie
    includes the ARF, I believe).

    """

    NCHAN = 256
    NENERGY = 729
    NUMGRP = 719
    NUMELT = 115269

    infile = make_data_path("pspcc_gain1_256.rsp")
    orig = io.read_rmf(infile)

    assert "NUMELT" not in orig.header
    assert "NUMGRP" not in orig.header
    assert "HDUVERS" not in orig.header
    assert orig.header["HDUVERS2"] == "1.2.0"
    assert orig.header["HDUCLAS3"] == "FULL"

    outpath = tmp_path / "out.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, orig)

    new = io.read_rmf(outfile)
    assert isinstance(new, DataRMF)
    assert new.detchans == NCHAN
    assert new.offset == 1
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert "DETCHANS" not in hdr

    assert hdr["NUMGRP"] == NUMGRP
    assert hdr["NUMELT"] == NUMELT

    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "RSP_MATRIX"
    assert hdr["HDUCLAS3"] == "FULL"
    assert "HDUVERS1" not in hdr
    assert hdr["HDUVERS2"] == "1.2.0"
    assert "HDUVERS3" not in hdr
    assert hdr["HDUVERS"] == "1.3.0"
    assert hdr["RMFVERSN"] == "1992a"

    assert hdr["TELESCOP"] == "ROSAT"
    assert hdr["INSTRUME"] == "PSPCC"
    assert hdr["FILTER"] == "NONE"

    assert hdr["CHANTYPE"] == "PI"
    assert hdr["LO_THRES"] == 0

    assert len(new.energ_lo) == NENERGY
    assert len(new.energ_hi) == NENERGY
    assert len(new.n_grp) == NENERGY

    assert len(new.f_chan) == NUMGRP
    assert len(new.n_chan) == NUMGRP

    assert len(new.matrix) == NUMELT

    assert len(new.e_min) == NCHAN
    assert len(new.e_max) == NCHAN

    assert new.n_grp.dtype == np.uint64

    for field in ["f_chan", "n_chan"]:
        attr = getattr(new, field)
        assert attr.dtype == np.uint64

    for field in ["energ_lo", "energ_hi", "matrix", "e_min", "e_max"]:
        attr = getattr(new, field)
        assert attr.dtype == np.float64

    assert (new.energ_lo[1:] == new.energ_hi[:-1]).all()
    assert new.energ_lo[0] == pytest.approx(0.054607998579740524)
    assert new.energ_hi[-1] == pytest.approx(3.009999990463257)

    assert new.n_grp.sum() == NUMGRP
    assert new.n_grp.max() == 1
    assert new.n_grp.min() == 0

    assert new.f_chan.sum() == 7732
    assert new.n_chan.sum() == NUMELT

    assert new.matrix.sum() == pytest.approx(42994.35767966656)

    # e_min/max come from the EBOUNDS block
    #
    assert (new.e_min[1:] == new.e_max[:-1]).all()
    assert new.e_min[0] == pytest.approx(0.01495600026100874)
    assert new.e_max[-1] == pytest.approx(2.552633047103882)


@requires_data
@requires_fits
def test_write_rmf_fits_xmm_rgs(make_data_path, tmp_path):
    """Check we can write out an RMF as a FITS file.

    Try out a grating RMF. Use XMM since it is smaller than Chandra.

    """

    NCHAN = 3600
    NENERGY = 4000
    NUMGRP = 82878
    NUMELT = 8475822

    infile = make_data_path("xmmrgs/P0112880201R1S004RSPMAT1003.FTZ")
    orig = io.read_rmf(infile)

    assert "NUMELT" not in orig.header
    assert "NUMGRP" not in orig.header
    assert orig.header["HDUVERS"] == "1.3.0"
    assert orig.header["HDUVERS1"] == "1.1.0"
    assert orig.header["HDUVERS2"] == "1.2.0"
    assert "HDUVERS3" not in orig.header
    assert orig.header["HDUCLAS3"] == "FULL"

    # Just to check, as this is different to other RMNFs
    assert (orig.e_min[:-1] == orig.e_max[1:]).all()
    assert orig.e_min[-1] == pytest.approx(0.30996060371398926)
    assert orig.e_max[0] == pytest.approx(3.0996062755584717)

    outpath = tmp_path / "out.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, orig)

    new = io.read_rmf(outfile)
    assert isinstance(new, DataRMF)
    assert new.detchans == NCHAN
    assert new.offset == 1
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert "DETCHANS" not in hdr

    assert hdr["NUMGRP"] == NUMGRP
    assert hdr["NUMELT"] == NUMELT

    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "RSP_MATRIX"
    assert hdr["HDUCLAS3"] == "FULL"
    assert hdr["HDUVERS1"] == "1.1.0"
    assert hdr["HDUVERS2"] == "1.2.0"
    assert "HDUVERS3" not in hdr
    assert hdr["HDUVERS"] == "1.3.0"
    assert hdr["RMFVERSN"] == "1992a"

    assert hdr["TELESCOP"] == "XMM"
    assert hdr["INSTRUME"] == "RGS1"
    assert hdr["FILTER"] == "NONE"

    assert hdr["CHANTYPE"] == "PI"
    assert np.log10(hdr["LO_THRES"]) == pytest.approx(-7.0)

    assert len(new.energ_lo) == NENERGY
    assert len(new.energ_hi) == NENERGY
    assert len(new.n_grp) == NENERGY

    assert len(new.f_chan) == NUMGRP
    assert len(new.n_chan) == NUMGRP

    assert len(new.matrix) == NUMELT

    assert len(new.e_min) == NCHAN
    assert len(new.e_max) == NCHAN

    assert new.n_grp.dtype == np.uint64

    for field in ["f_chan", "n_chan"]:
        attr = getattr(new, field)
        assert attr.dtype == np.uint64

    for field in ["energ_lo", "energ_hi", "matrix", "e_min", "e_max"]:
        attr = getattr(new, field)
        assert attr.dtype == np.float64

    assert (new.energ_lo[1:] == new.energ_hi[:-1]).all()
    assert new.energ_lo[0] == pytest.approx(0.30000001192092896)
    assert new.energ_hi[-1] == pytest.approx(2.7998998165130615)

    assert new.n_grp.sum() == NUMGRP
    assert new.n_grp.max() == 48
    assert new.n_grp.min() == 12

    assert new.f_chan.sum() == 173869713
    assert new.n_chan.sum() == NUMELT

    assert new.matrix.sum() == pytest.approx(102430.35533461263)

    # e_min/max come from the EBOUNDS block; note the ordering
    # which is unlike the other cases (i.e. descending)
    #
    assert (new.e_min[:-1] == new.e_max[1:]).all()
    assert new.e_min[-1] == pytest.approx(0.30996060371398926)
    assert new.e_max[0] == pytest.approx(3.0996062755584717)


@pytest.mark.parametrize("rsptype", ["arf", "rmf"])
def test_write_response_sent_invalid_data(rsptype, tmp_path):
    """Check we error out. It is not a complete check."""

    pha = DataPHA("ex", [1, 2], [1, 2])
    savefunc = getattr(io, f"write_{rsptype}")
    outpath = tmp_path / "no.fits"

    emsg = "data set is not a"
    if rsptype == "arf":
        emsg += "n ARF"
    else:
        emsg += " RMF"

    with pytest.raises(IOErr, match=f"^{emsg}$"):
        savefunc(str(outpath), pha)

    assert not outpath.exists()


# At the moment there's no equivalent check for the RMF.
def test_write_arf_missing_cols(tmp_path):
    """Limited check."""

    empty = DataARF("empty", np.asarray([1, 2]),
                    np.asarray([2, 3]), None)
    outpath = tmp_path / "no.arf"

    with pytest.raises(ArgumentErr,
                       match="^Invalid ARF: 'SPECRESP column is missing'$"):
        io.write_arf(str(outpath), empty)

    assert not outpath.exists()


@requires_fits
def test_write_fake_arf(tmp_path):
    """Check we can write out an ARF we create.

    Just check the FITS path.
    """

    ebins = np.arange(0.1, 1.2, 0.1)
    elo = ebins[:-1]
    ehi = ebins[1:]
    y = np.arange(10, 30, 2)
    arf = create_arf(elo, ehi, y)

    outpath = tmp_path / "arf.faked"
    outfile = str(outpath)
    io.write_arf(outfile, arf, ascii=False)

    new = io.read_arf(outfile)
    assert isinstance(new, DataARF)
    assert new.exposure is None
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "SPECRESP"
    assert hdr["HDUVERS"] == "1.1.0"
    assert hdr["TELESCOP"] == "none"
    assert hdr["INSTRUME"] == "none"
    assert hdr["FILTER"] == "none"
    assert hdr["CREATOR"].startswith("sherpa ")

    # We may add keywords, which will mean this needs updating.  It
    # also relies on the removal of "structural" keywords to make this
    # match between different backends.
    #
    assert len(hdr) == 8

    assert new.energ_lo == pytest.approx(elo)
    assert new.energ_hi == pytest.approx(ehi)
    assert new.specresp == pytest.approx(y)

    assert new.bin_lo is None
    assert new.bin_hi is None


# Is offset=10 valid? As far as I can see there's no reason to not
# allow this - from reading the OGIP documentation - even if no-one
# uses it.
#
@requires_fits
@pytest.mark.parametrize("offset", [0, 1, 10])
def test_write_fake_perfect_rmf(offset, tmp_path):
    """Check we can write out a RMF we create."""

    ebins = np.arange(0.1, 1.2, 0.1)
    elo = ebins[:-1]
    ehi = ebins[1:]
    rmf = create_delta_rmf(elo, ehi, offset=offset, e_min=elo, e_max=ehi)

    outpath = tmp_path / "rmf.faked"
    outfile = str(outpath)
    io.write_rmf(outfile, rmf)

    new = io.read_rmf(outfile)
    assert isinstance(new, DataRMF)
    assert new.detchans == 10
    assert new.offset == offset
    assert np.log10(new.ethresh) == pytest.approx(-10)

    hdr = new.header
    assert "HDUNAME" not in hdr
    assert hdr["HDUCLASS"] == "OGIP"
    assert hdr["HDUCLAS1"] == "RESPONSE"
    assert hdr["HDUCLAS2"] == "RSP_MATRIX"
    assert hdr["HDUVERS"] == "1.3.0"
    assert hdr["TELESCOP"] == "none"
    assert hdr["INSTRUME"] == "none"
    assert hdr["FILTER"] == "none"
    assert hdr["CHANTYPE"] == "PI"
    assert hdr["NUMGRP"] == 10
    assert hdr["NUMELT"] == 10
    assert hdr["CREATOR"].startswith("sherpa ")

    # We may add keywords, which will mean this needs updating.  It
    # also relies on the removal of "structural" keywords to make this
    # match between different backends.
    #
    assert len(hdr) == 11

    assert new.energ_lo == pytest.approx(elo)
    assert new.energ_hi == pytest.approx(ehi)
    assert new.n_grp == pytest.approx([1] * 10)
    assert new.f_chan == pytest.approx(np.arange(1, 11) + offset - 1)
    assert new.n_chan == pytest.approx([1] * 10)

    assert new.matrix == pytest.approx([1] * 10)

    assert new.e_min == pytest.approx(elo)
    assert new.e_max == pytest.approx(ehi)


@requires_data
@requires_fits
def test_write_rmf_fits_xmm_epn(make_data_path, tmp_path, caplog):
    """This XMM dataset has been useful for edge cases. So reuse it.

    See sherpa-test-data/sherpatest/xmm-fit.py
    """

    NCHAN = 4096
    NENERGY = 1319
    NUMELT = 487198

    def check(rmf):
        assert rmf.detchans == NCHAN
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
        assert "RMFVERSN" not in hdr

        assert hdr["TELESCOP"] == "XMM"
        assert hdr["INSTRUME"] == "EPN"
        assert hdr["FILTER"] == "NONE"

        assert hdr["CHANTYPE"] == "PI"
        assert hdr["LO_THRES"] == pytest.approx(2e-6)

        assert len(rmf.energ_lo) == NENERGY
        assert len(rmf.energ_hi) == NENERGY
        assert len(rmf.n_grp) == NENERGY

        assert len(rmf.f_chan) == NENERGY
        assert len(rmf.n_chan) == NENERGY

        assert len(rmf.matrix) == NUMELT

        assert len(rmf.e_min) == NCHAN
        assert len(rmf.e_max) == NCHAN

        assert rmf.n_grp.dtype == np.uint64

        for field in ["f_chan", "n_chan"]:
            attr = getattr(rmf, field)
            assert attr.dtype == np.uint64

        for field in ["energ_lo", "energ_hi", "matrix", "e_min", "e_max"]:
            attr = getattr(rmf, field)
            assert attr.dtype == np.float64

        assert (rmf.energ_lo[1:] == rmf.energ_hi[:-1]).all()
        assert rmf.energ_lo[0] == pytest.approx(0.05000000074505806)
        assert rmf.energ_hi[-1] == pytest.approx(18.049999237060547)

        assert rmf.n_grp.sum() == NENERGY
        assert np.argmax(rmf.n_grp) == 0

        # It is not obvious how the RMF is flattened here, so treat this
        # as a regression test rather than "from first principles"
        #
        assert rmf.f_chan.sum() == 544580
        assert rmf.n_chan.sum() == NUMELT

        assert rmf.matrix.sum() == pytest.approx(1200.610821738967)

        # e_min/max come from the EBOUNDS block
        #
        assert (rmf.e_min[1:] == rmf.e_max[:-1]).all()
        assert rmf.e_min[0] == pytest.approx(0)
        assert rmf.e_max[-1] == pytest.approx(20.48000144958496)


    rmfname = "epn_ff20_dY9.rmf"
    infile = make_data_path(rmfname)

    # In 4.16.0 and earlier this gave a TLMIN warning, even though we
    # managed to read the data in from the EBOUNDS/CHANNEL column. The
    # code now doesn't warn in this case.
    #
    assert len(caplog.record_tuples) == 0
    orig = io.read_rmf(infile)
    assert len(caplog.record_tuples) == 0

    assert isinstance(orig, DataRMF)
    assert orig.name.endswith(f"/{rmfname}")
    check(orig)

    outpath = tmp_path / "out.rmf"
    outfile = str(outpath)
    io.write_rmf(outfile, orig)
    assert len(caplog.record_tuples) == 0

    new = io.read_rmf(outfile)
    assert isinstance(new, DataRMF)
    check(new)


@requires_fits
def test_read_multi_matrix_rmf(tmp_path, caplog):
    """Check reading in a multi-MATRIX RMF block.

    As we don't currently have one "in the wild" we manually create
    one.

    This is a regression test since we currently only read in one of
    the matrices.

    """

    # First block is "perfect" and the second adds a "blob"
    # component.
    #
    egrid = np.linspace(0.1, 2.1, 21)
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf1 = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi,
                            name="perfect")
    rmf1.header["TESTBLCK"] = 1

    # reduce the "energy" resolution of the second block.
    #
    e4grid = np.linspace(0.1, 2.1, 6)
    e4lo = e4grid[:-1]
    e4hi = e4grid[1:]
    matrix = np.zeros((5, 4), dtype=np.float32)
    blur = [0.1, 0.3, 0.4, 0.2]
    for i in range(5):
        matrix[i] = blur

    rmf2 = DataRMF("blob", detchans=rmf1.detchans,
                   energ_lo=e4lo, energ_hi=e4hi,
                   n_grp=np.ones(5, dtype=np.int32),
                   f_chan=np.asarray([2, 4, 6, 8, 10], dtype=np.int32),
                   n_chan=4 * np.ones(5, dtype=np.int32),
                   matrix=matrix.flatten(),
                   e_min=elo, e_max=ehi,
                   header={"TESTBLCK": 2})

    # Extract the data and then reconstruct a "multi" RMF.
    #
    blocks1 = io._pack_rmf(rmf1)
    blocks2 = io._pack_rmf(rmf2)

    header = Header([HeaderItem(name="TESTKEY", value=12)])
    blocks = [blocks1.blocks[0], blocks2.blocks[0], blocks1.blocks[1]]
    blist = BlockList(blocks=blocks, header=header)

    outpath = tmp_path / "multi.rmf"
    outfile = str(outpath)
    io.backend.set_hdus(outfile, blist)

    # What happens if we try to read this in?
    #
    assert len(caplog.record_tuples) == 0
    rmf = io.read_rmf(outfile)
    assert len(caplog.record_tuples) == 1

    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.astro.io"
    assert lvl == logging.ERROR
    assert msg.startswith("RMF in ")
    assert msg.endswith("/multi.rmf contains 2 MATRIX blocks; Sherpa only uses the first block!")

    # What happens if we apply the RMF to a model?
    #
    mdl = Gauss1D()
    mdl.pos = 1.1
    mdl.fwhm = 0.8
    mdl.ampl = 1e4

    resp = RMF1D(rmf)
    conv = resp(mdl)
    chans = np.arange(1, 21, dtype=np.int16)
    y = conv(chans)

    # We create the responses from both matrices as eventually we will
    # need both (i.e. when multi-matrix RMF are supported).
    #
    expected_perfect = mdl(elo, ehi)

    # Create a 2D array to represent the blurry matrix
    #
    blurry_matrix = np.zeros((5, 20))
    for idx, fchan in enumerate([2, 4, 6, 8, 10]):
        blurry_matrix[idx, fchan - 1:fchan + 3] = blur

    expected_blurry = mdl(e4lo, e4hi) @ blurry_matrix

    assert y == pytest.approx(expected_perfect)


@requires_fits
@requires_data
def test_read_arf_object(make_data_path):
    """Check we can send in a 'object'.

    This support is not well advertised so it could perhaps be
    removed, but let's add a basic test at least.

    """

    # Could create the "object" manually, but let's just read one in
    #
    infile = make_data_path("MNLup_2138_0670580101_EMOS1_S001_spec.arf")
    close = False

    if backend_is("crates"):
        import pycrates  # type: ignore
        arg = pycrates.read_file(infile)

    elif backend_is("pyfits"):
        from astropy.io import fits  # type: ignore
        arg = fits.open(infile)
        close = True

    else:
        assert False, f"unknown backend: {io.backend}"

    try:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            arf = io.read_arf(arg)

    finally:
        if close:
            arg.close()

    assert isinstance(arf, DataARF)
    assert arf.energ_lo[0:4] == pytest.approx([1e-10, 0.005, 0.01, 0.015])
    assert arf.energ_hi[0:4] == pytest.approx([0.005, 0.01, 0.015, 0.02])
    assert arf.specresp.max() == pytest.approx(449.2618408203125)

    assert len(ws) == 1
    msg = ws[0].message
    assert isinstance(msg, UserWarning)
    msg = str(msg)
    assert msg.startswith("The minimum ENERG_LO in the ARF '")
    assert msg.endswith("_spec.arf' was 0 and has been replaced by 1e-10")


@requires_fits
@requires_data
def test_read_rmf_object(make_data_path):
    """Check we can send in a 'object'.

    This support is not well advertised so it could perhaps be
    removed, but let's add a basic test at least.

    """

    # Could create the "object" manually, but let's just read one in
    #
    infile = make_data_path("source.rmf")
    close = False

    if backend_is("crates"):
        import pycrates  # type: ignore
        arg = pycrates.RMFCrateDataset(infile, mode='r')

    elif backend_is("pyfits"):
        from astropy.io import fits  # type: ignore
        arg = fits.open(infile)
        close = True

    else:
        assert False, f"unknown backend: {io.backend}"

    try:
        rmf = io.read_rmf(arg)

    finally:
        if close:
            arg.close()

    assert isinstance(rmf, DataRMF)
    assert rmf.energ_lo[0:4] == pytest.approx([0.1, 0.11, 0.12, 0.13])
    assert rmf.energ_hi[0:4] == pytest.approx([0.11, 0.12, 0.13, 0.14])
    assert rmf.matrix.max() == pytest.approx(0.12998242676258087)


def test_rmf_fails_if_no_ebounds(tmp_path):
    """Check we error out if e_min/max has not been set up.

    This is a regression test.
    """

    egrid = np.asarray([0.1, 0.2, 0.3])
    rmf = create_delta_rmf(egrid[:-1], egrid[1:])
    rmf.e_min = None
    rmf.e_max = None

    outfile = tmp_path / 'do-not-create'
    with pytest.raises(IOErr,
                       match="^RMF delta-rmf has no E_MIN or E_MAX data$"):
        io.write_rmf(str(outfile), rmf, clobber=True)
