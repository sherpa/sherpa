#
#  Copyright (C) 2015, 2016, 2018, 2021 - 2024
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

# break out some of the I/O tests in the UI layer into a separate set
# of tests, in part to reduce file size, but also because some of these
# may be better placed in tests of the sherpa.astro.io module, once that
# becomes possible

import pytest

from sherpa.astro import ui
from sherpa.astro.data import DataPHA
from sherpa.astro.instrument import ARF1D, RMF1D
from sherpa.utils.err import ArgumentErr
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_group, requires_xspec


FILE_NAME = 'acisf01575_001N001_r0085_pha3.fits'


def validate_pha(idval):
    """Check that the PHA dataset in id=idval is
    as expected.
    """

    assert ui.list_data_ids() == [idval]

    pha = ui.get_data(idval)
    assert isinstance(pha, DataPHA)

    arf = ui.get_arf(idval)
    assert isinstance(arf, ARF1D)

    rmf = ui.get_rmf(idval)
    assert isinstance(rmf, RMF1D)

    bpha = ui.get_bkg(idval, bkg_id=1)
    assert isinstance(bpha, DataPHA)

    barf = ui.get_arf(idval, bkg_id=1)
    assert isinstance(barf, ARF1D)

    brmf = ui.get_rmf(idval, bkg_id=1)
    assert isinstance(brmf, RMF1D)

    # normally the background data set would have a different name,
    # but this is a  PHA Type 3 file.
    # assert pha.name == bpha.name
    assert arf.name == barf.name
    assert rmf.name == brmf.name


@requires_fits
@requires_data
def test_pha3_read_explicit(make_data_path, clean_astro_ui):
    """Include .gz in the file name"""

    fname = make_data_path(FILE_NAME + '.gz')
    idval = 12
    ui.load_pha(idval, fname)

    validate_pha(idval)

    # TODO: does this indicate that the file name, as read in,
    #       should have the .gz added to it to match the data
    #       read in, or left as is?
    pha = ui.get_data(idval)
    bpha = ui.get_bkg(idval, bkg_id=1)
    assert pha.name == bpha.name + '.gz'
    assert pha.name == fname


@requires_fits
@requires_data
def test_pha3_read_implicit(make_data_path, clean_astro_ui):
    """Exclude .gz from the file name"""

    idval = "13"
    fname = make_data_path(FILE_NAME)
    ui.load_pha(idval, fname)

    validate_pha(idval)

    pha = ui.get_data(idval)
    bpha = ui.get_bkg(idval, bkg_id=1)
    assert pha.name == bpha.name
    assert pha.name == fname


@requires_fits
@requires_data
def test_hrci_imaging_mode_spectrum(make_data_path, clean_astro_ui):
    """This is a follow-on test based on issue #1830"""

    infile = make_data_path("chandra_hrci/hrcf24564_000N030_r0001"
                            "_pha3.fits.gz")
    ui.load_pha(infile)

    # Just check we can evaluate the model folding through the
    # response, and compare to the data. Getting a "nice" fit
    # is surprisingly hard so do a non-physical fit.
    #
    ui.notice(4, 6)
    pl = ui.create_model_component("powlaw1d", "pl")
    pl.gamma = 2.24
    pl.ampl = 4.6e-3
    ui.set_source(ui.powlaw1d.pl)

    ui.set_stat("chi2gehrels")

    # regression tests (the expected values were calculated using CIAO
    # 4.15/crates).
    #
    assert ui.calc_stat() == pytest.approx(39.58460766890158)

    ui.subtract()
    assert ui.calc_stat() == pytest.approx(39.611346193696164)


def check_fit_stats():
    """Do we get the expected results?

    Use of this means @requires_xspec, @requires_group
    """

    ui.set_method("levmar")
    ui.set_stat("chi2xspecvar")
    ui.set_xsabund("angr")
    ui.set_xsxsect("vern")

    gal = ui.create_model_component("xsphabs", "gal")
    pl = ui.create_model_component("powlaw1d", "pl")
    ui.set_source(gal * pl)

    ui.notice(0.5, 6)
    ui.group_counts(15, tabStops=~ui.get_data().get_mask())

    ui.subtract()
    ui.fit()

    # This is a regression test, as all we care about is that the
    # results match for the direct route and after saving the files
    # and loading them in. Therefore it's okay for these values to
    # change.
    #
    # The test could be changed to just fit a power law, which would
    # mean no need for the requires_xspec decorator, but given we need
    # I/O and the data files anyway it doesn't really improve things
    # if we did that.
    #
    # This test produces different results on macOS to Linux, hence
    # the relaxed tolerances.
    #
    assert pl.gamma.val == pytest.approx(1.6870563856336693,
                                         rel=3.5e-5)
    assert pl.ampl.val == pytest.approx(3.831373278007354e-05,
                                        rel=6.5e-5)
    assert gal.nH.val == pytest.approx(0.24914805790330877,
                                       rel=1.1e-4)
    assert ui.calc_stat() == pytest.approx(41.454087606314054,
                                           rel=1e-5)


@requires_xspec
@requires_group
@requires_fits
@requires_data
def test_pha3_original_check(make_data_path, clean_astro_ui):
    """Check we get the expected fit statistic with the PHA3 file.

    See also test_pha3_roundtrip_check.
    """

    fname = make_data_path(FILE_NAME)
    ui.load_pha(fname)

    check_fit_stats()


@requires_xspec
@requires_group
@requires_fits
@requires_data
def test_pha3_roundtrip_check(make_data_path, clean_astro_ui, tmp_path):
    """Check we get the expected fit statistic with the PHA3 file.

    See also test_pha3_original_check. This version writes out the
    PHA, background PHA, ARF, and RMF and then reads those back in and
    fits with those.

    """

    fname = make_data_path(FILE_NAME)
    ui.load_pha(fname)

    src = str(tmp_path / "src.pi")
    bkg = str(tmp_path / "bkg.pi")
    arf = str(tmp_path / "arf")
    rmf = str(tmp_path / "rmf")

    # This will write out invalid names for the ARF and RMF files (as
    # it refers to the versions from the original file which do not
    # exist in the temporary directory) but that's okay and we do not
    # want to check the screen output to validate this behavior as
    # this is not the point of this test.
    #
    ui.save_pha(1, src)
    ui.save_pha(1, bkg, bkg_id=1)
    ui.save_arf(1, arf)
    ui.save_rmf(1, rmf)

    ui.clean()

    ui.load_pha(src)
    ui.load_bkg(bkg)
    ui.load_rmf(rmf)
    ui.load_arf(arf)

    check_fit_stats()


@requires_fits
@pytest.mark.parametrize("rsptype", ["arf", "rmf"])
@pytest.mark.parametrize("bkg", [None, 1])
def test_save_response_no_arf(rsptype, bkg, clean_astro_ui, tmp_path):
    """Check error handling"""

    # Create PHA with no responses
    ui.load_arrays("pha", [1, 2, 3], [4, 0, 2], DataPHA)

    # Create a background component
    ui.set_bkg("pha", DataPHA("x", [1, 2, 3], [1, 0, 0]))

    if bkg is None:
        emsg = ""
    else:
        emsg = "background '1' of "

    emsg += r"data set 'pha' \(response 1\) does not contain a"
    if rsptype == "arf":
        emsg += "n "
    else:
        emsg += " "

    emsg += rsptype.upper()

    savefunc = getattr(ui, f"save_{rsptype}")
    outpath = tmp_path / "should-not-exist"
    with pytest.raises(ArgumentErr, match=f"^{emsg}$"):
        savefunc("pha", str(outpath), bkg_id=bkg)

    assert not outpath.exists()
