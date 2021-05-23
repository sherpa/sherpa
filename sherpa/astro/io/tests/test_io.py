#
#  Copyright (C) 2016, 2017, 2018, 2020, 2021  Smithsonian Astrophysical Observatory
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

from tempfile import NamedTemporaryFile
import struct

import numpy as np

import pytest

from sherpa.utils.testing import requires_data, requires_fits, requires_xspec
from sherpa.astro import ui
from sherpa.models.basic import Box1D, Const1D
from sherpa.astro.io import print_meta


def test_str_singleton():
    """stringification: single value"""

    store = {}

    # Could loop over this, but a bit easier to check the string
    # output this way.
    #
    store['key'] = ""
    assert print_meta(store) == '\n key           = '

    store['key'] = "  "
    assert print_meta(store) == '\n key           =   '

    store['key'] = " a string "
    assert print_meta(store) == '\n key           =  a string '

    store['key'] = 23
    assert print_meta(store) == '\n key           = 23'

    store['key'] = 23.0
    assert print_meta(store) == '\n key           = 23.0'

    store['key'] = False
    assert print_meta(store) == '\n key           = False'

    # Now some special cases
    for val in [None, "None", "NONE", "none"]:
        store['key'] = val
        assert print_meta(store) == ''


@requires_fits
def test_str_multi():
    """Multiple keys are displayed as expected"""

    store = {}
    store['Xkey'] = 'X X'
    store['xkey'] = ' y  y'
    store['a'] = 23
    store['INFILE'] = 'none'
    store['outfile'] = '/tmp/b.fits'

    lines = str(store).split('\n')
    assert len(lines) == 5
    assert lines[0] == ''
    assert lines[1] == ' Xkey          = X X'
    assert lines[2] == ' a             = 23'
    assert lines[3] == ' outfile       = /tmp/b.fits'
    assert lines[4] == ' xkey          =  y  y'


@requires_data
@requires_fits
@requires_xspec
def test_mod_fits(make_data_path):
    tablemodelfile = make_data_path("xspec-tablemodel-RCS.mod")
    ui.load_table_model("tmod", tablemodelfile)
    tmod = ui.get_model_component("tmod")
    assert tmod.name == "xstablemodel.tmod"


@requires_fits
def test_warnings_are_gone_arrays():
    ui.load_arrays(1, [1, 2, 3], [4, 5, 6])
    #  We now have logic in conftest.py to catch white-listed warnings and fail on unexpected ones.
    #  We just need to make any warnings bubble up, here and in the following test.
    with NamedTemporaryFile() as f:
        ui.save_data(1, f.name, ascii=True, clobber=True)
    with NamedTemporaryFile() as f:
        ui.save_data(1, f.name, ascii=False, clobber=True)


@requires_fits
@requires_data
def test_warnings_are_gone_pha(make_data_path):
    pha = make_data_path("3c273.pi")
    ui.load_pha(pha)
    with NamedTemporaryFile() as f:
        ui.save_data(1, f.name, ascii=False, clobber=True)


def assert_staterr(use_errors):
    assert np.all(ui.get_data("phacounts").counts ==
                  pytest.approx(ui.get_data("pharate").counts))
    if use_errors is True:
        assert np.all(ui.get_data("phacounts").staterror ==
                      pytest.approx(ui.get_data("pharate").staterror))
    else:
        assert ui.get_data("phacounts").staterror is None
        assert ui.get_data("pharate").staterror is None

    for n in ['phacounts', 'pharate']:
        ui.group_bins(n, 16)
    ui.set_analysis('energy')
    ui.ignore(None, 3.)
    if use_errors is True:
        assert np.all(ui.get_data("phacounts").get_error() ==
                      pytest.approx(ui.get_data("pharate").get_error()))
    else:
        assert ui.get_data("phacounts").get_error() is None
        assert ui.get_data("pharate").get_error() is None


@requires_fits
@requires_data
@pytest.mark.parametrize("use_errors", [True, False])
def test_scaling_staterr(make_data_path, use_errors):
    '''Regression test for https://github.com/sherpa/sherpa/issues/800

    Notes
    -----
    Test files are made with (in sherpa-test-data/sherpatest)
    dmtcalc source.pi source_edit.pi expression="stat_err=stat_err/exposure"
    dmcopy "source_edit.pi[cols channel,pi,RATE=count_rate,stat_err]" source_rate.pi clobber=yes
    # Then fix units to counts/s manually in the header
    '''
    ui.load_pha("phacounts", make_data_path("source.pi"),
                use_errors=use_errors)
    ui.load_pha("pharate", make_data_path("source_rate.pi"),
                use_errors=use_errors)
    for n in ['phacounts', 'pharate']:
        ui.load_arf(n, make_data_path("source.arf"))
        ui.load_rmf(n, make_data_path("source.rmf"))

    assert_staterr(use_errors)


@requires_fits
@requires_data
@pytest.mark.parametrize("use_errors", [True, False])
def test_scaling_staterr_pha2(make_data_path, use_errors):
    '''Regression test for https://github.com/sherpa/sherpa/issues/800

    Notes
    -----
    Test files are made with (in sherpa-test-data/sherpatest)
    dmtcalc 3c120_pha2.gz 3c120_pha2_edit.gz expression="RATE=(float)counts/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit.gz 3c120_pha2_edit2.gz expression="stat_err=stat_err/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit2.gz 3c120_pha2_edit3.gz expression="backgroun_up=(float)background_up/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit3.gz 3c120_pha2_edit4.gz expression="backgroun_down=(float)background_down/exposure" clobber=yes

    dmcopy "3c120_pha2_edit4.gz[cols -COUNTS,-BACKGROUND_UP,-BACKGROUND_DOWN][cols *,BACKGROUND_UP=backgroun_up,BACKGROUND_DOWN=backgroun_down][col -COUNTS]" 3c120_pha2_rate.gz clobber=yes
    Then fix units manually.
    '''
    ui.load_pha(make_data_path("3c120_pha2_rate.gz"),
                use_errors=use_errors)
    ui.copy_data(9, "phacounts")
    ui.load_pha(make_data_path("3c120_pha2_rate.gz"),
                use_errors=use_errors)
    ui.copy_data(9, "pharate")
    for n in ['phacounts', 'pharate']:
        ui.load_arf(n, make_data_path("3c120_meg_-1.arf.gz"))
        ui.load_rmf(n, make_data_path("3c120_meg_-1.rmf.gz"))
    assert_staterr(use_errors)


def fake_rmf(outfile):
    """Create a "perfect" RMF with a scalar MATRIX column

    We do this rather than add another file to the test
    data directory as it is easier.
    """

    ebins = np.arange(0.15, 0.2, 0.01, dtype=np.float32)
    elo = ebins[:-1]
    ehi = ebins[1:]
    nchan = elo.size
    chans = np.arange(1, nchan + 1, dtype=np.int16)

    def hdr(key, value):
        if isinstance(value, str):
            value = "'{}'".format(value)
        elif isinstance(value, bool):
            value = 'T' if value else 'F'

        out = "{:8s}= {}".format(key, value)
        return out.ljust(80)

    hdr_img = [hdr('SIMPLE', True),
               hdr('BITPIX', 8),
               hdr('NAXIS', 0),
               hdr('EXTEND', True)]

    hdr_ebounds = [hdr('XTENSION', 'BINTABLE'),
                   hdr('BITPIX', 8),
                   hdr('NAXIS', 2),
                   hdr('NAXIS1', 10),
                   hdr('NAXIS2', nchan),
                   hdr('PCOUNT', 0),
                   hdr('GCOUNT', 1),
                   hdr('TFIELDS', 3),
                   hdr('TTYPE1', 'CHANNEL'),
                   hdr('TFORM1', 'I'),
                   hdr('TTYPE2', 'E_MIN'),
                   hdr('TFORM2', 'E'),
                   hdr('TTYPE3', 'E_MAX'),
                   hdr('TFORM3', 'E'),
                   hdr('EXTNAME', 'EBOUNDS'),
                   hdr('HDUCLASS', 'OGIP'),
                   hdr('HDUCLAS1', 'RESPONSE'),
                   hdr('HDUCLAS2', 'EBOUNDS'),
                   hdr('HDUVERS', '1.3.0'),
                   hdr('CHANTYPE', 'PI'),
                   hdr('DETCHANS', nchan)]

    hdr_matrix = [hdr('XTENSION', 'BINTABLE'),
                  hdr('BITPIX', 8),
                  hdr('NAXIS', 2),
                  hdr('NAXIS1', 18),
                  hdr('NAXIS2', nchan),
                  hdr('PCOUNT', 0),
                  hdr('GCOUNT', 1),
                  hdr('TFIELDS', 6),
                  hdr('TTYPE1', 'ENERG_LO'),
                  hdr('TFORM1', 'E'),
                  hdr('TTYPE2', 'ENERG_HI'),
                  hdr('TFORM2', 'E'),
                  hdr('TTYPE3', 'N_GRP'),
                  hdr('TFORM3', 'I'),
                  hdr('TTYPE4', 'F_CHAN'),
                  hdr('TFORM4', 'I'),
                  hdr('TTYPE5', 'N_CHAN'),
                  hdr('TFORM5', 'I'),
                  hdr('TTYPE6', 'MATRIX'),
                  hdr('TFORM6', 'E'),
                  hdr('EXTNAME', 'MATRIX'),
                  hdr('HDUCLASS', 'OGIP'),
                  hdr('HDUCLAS1', 'RESPONSE'),
                  hdr('HDUCLAS2', 'RSP_MATRIX'),
                  hdr('HDUCLAS3', 'REDIST'),
                  hdr('HDUVERS', '1.3.0'),
                  hdr('CHANTYPE', 'PI'),
                  hdr('TLMIN4', 1),
                  hdr('DETCHANS', nchan),
                  hdr('NUMGRP', nchan),
                  hdr('NUMELT', nchan),
                  hdr('LO_THRES', '1E-06')]

    ngrps = np.ones(nchan, dtype=np.int16)
    fchans = np.arange(1, nchan + 1, dtype=np.int16)
    nchans = ngrps
    matrix = np.ones(nchan, dtype=np.float32)

    def row1(chan, e1, e2):
        return struct.pack('>hff', chan, e1, e2)

    def row2(e1, e2, ngrp, fchan, nchan, mat):
        return struct.pack('>ffhhhf', e1, e2, ngrp, fchan, nchan, mat)

    with open(outfile, 'wb') as fh:

        def extend():
            """extend to next 2880 block"""
            p = fh.tell()
            n = p % 2880
            if n != 0:
                fh.write(b' ' * (2880 - n))

        def mkhdr(cards):
            for card in cards:
                fh.write(card.encode('ascii'))

            fh.write(b'END' + b' ' * 77)
            extend()

        mkhdr(hdr_img)

        mkhdr(hdr_ebounds)
        for row in zip(chans, elo, ehi):
            fh.write(row1(*row))

        extend()

        mkhdr(hdr_matrix)
        for row in zip(elo, ehi, ngrps, fchans, nchans, matrix):
            fh.write(row2(*row))

        extend()


@requires_fits
def test_read_ideal_rmf():
    """Can a RMF similar to issue #862 be read in?

    The MATRIX column in this file is a scalar rather than array,
    and let's do EBOUNDS then MATRIX blocks.
    """

    from sherpa.astro.io import read_rmf

    ebins = np.arange(0.15, 0.2, 0.01)
    elo = ebins[:-1]
    ehi = ebins[1:]

    with NamedTemporaryFile() as f:
        fake_rmf(f.name)
        r = read_rmf(f.name)

    # Can we read in the data
    #
    assert r.detchans == 5
    assert r.energ_lo == pytest.approx(elo)
    assert r.energ_hi == pytest.approx(ehi)
    assert (r.n_grp == [1, 1, 1, 1, 1]).all()
    assert (r.f_chan == [1, 2, 3, 4, 5]).all()
    assert (r.n_chan == [1, 1, 1, 1, 1]).all()
    assert r.offset == 1
    assert r.e_min == pytest.approx(elo)
    assert r.e_max == pytest.approx(ehi)
    assert r.ethresh == 1e-10

    # Can we apply it?
    #
    # The cmdl evalutes to a value of 2 * bin width
    # The bmdl evalates to the bin width * x
    # where x = [0, 1, 0.5, 0, 0]
    #
    cmdl = Const1D()
    cmdl.c0 = 2

    bmdl = Box1D()
    bmdl.xlow = 0.16
    bmdl.xhi = 0.175

    mdl = bmdl + cmdl

    # Multiply by 100 so numbers are close to unity
    expected = 100 * 0.01 * np.asarray([2, 3, 2.5, 2, 2])
    y = 100 * r.eval_model(mdl)
    assert y == pytest.approx(expected, rel=2e-6)
