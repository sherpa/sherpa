#
#  Copyright (C) 2019 - 2024
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

"""This is labelled ui but really tests astro/ui

"""

import logging

import numpy as np

import pytest

from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, \
    requires_xspec, requires_fits, requires_group

from sherpa.astro import ui
from sherpa.models.basic import Gauss2D
from sherpa.utils.err import StatErr


def startswith(got, expected):
    """Provide useful information if the strings do not match"""

    # We can use got.startswith(expected) but if there is an error
    # there is no information on where the difference is.
    #
    n = len(expected)
    assert got[:n] == expected


@requires_xspec
@requires_group
@requires_fits
@requires_data
def test_plot_pvalue(make_data_path, clean_astro_ui, caplog):
    """Check plot_pvalue with PHA data."""

    # Set the seed. However, there is a limited amount of
    # repeatability here as LikelihoodRatioTest.run uses
    # parallel_map_rng which means we can not rely on using the
    # RandomState class [as internally it ends up creating other
    # objects which are not guaranteed to be long-term repeatable],
    # which limits the tests we can make.
    #
    ui.set_rng(np.random.RandomState(123))

    ui.set_stat('cstat')
    ui.set_method("neldermead")

    fname = make_data_path("qso.pi")
    with SherpaVerbosity('WARN'):
        ui.load_pha(fname)
        # Why are we grouping to so high-a-value with CStat?
        ui.group_counts(10)
        ui.notice(0.3, 8)

    ui.set_model(ui.xsphabs.abs1 * (ui.xspowerlaw.p1 + ui.gauss1d.g1))

    # move the fit close to the best fit to save a small amount
    # of time.
    abs1.nh = 0.05
    p1.phoindex = 1.28
    p1.norm = 2e-4
    g1.ampl = 1.8e-5

    g1.pos = 3.
    g1.fwhm = 0.1
    ui.freeze(g1.pos, g1.fwhm)

    # Pick a small number of bins to try to reduce the runtime of the
    # test while still exercising the code.
    #
    NUM = 20
    NBINS = 8

    ui.fit()
    assert len(caplog.records) == 2

    ui.plot_pvalue(p1, p1 + g1, num=NUM, bins=NBINS)

    LR = 2.679487496941789

    tmp = ui.get_pvalue_results()

    assert tmp.null == pytest.approx(210.34566845619273)
    assert tmp.alt == pytest.approx(207.66618095925094)
    assert tmp.lr == pytest.approx(LR)

    # Have we returned the correct info?
    #
    assert tmp.samples.shape == (NUM, 2)
    assert tmp.stats.shape == (NUM, 2)
    assert tmp.ratios.shape == (NUM, )

    # Check the screen output
    #
    assert len(caplog.records) == 4

    r = caplog.records[2]
    assert r.name == "sherpa.astro.ui.utils"
    assert r.levelname == "WARNING"
    msg = "data set 1 has associated backgrounds, but they " + \
        "have not been subtracted, nor have background " + \
        "models been set"
    assert r.getMessage() == msg

    # For this, we do not do an exact comparison due to numerical
    # differences.
    #
    r = caplog.records[3]
    assert r.name == "sherpa.ui.utils"
    assert r.levelname == "INFO"
    msgs = r.getMessage().split("\n")
    assert msgs[0] == "Likelihood Ratio Test"
    startswith(msgs[1], "null statistic   =  210.3456")
    startswith(msgs[2], "alt statistic    =  207.6661")
    startswith(msgs[3], "likelihood ratio =  2.67948")
    startswith(msgs[4], "p-value       ")  # may be '= ..." or '< ...'
    assert len(msgs) == 5

    # Check the plot
    #
    tmp = ui.get_pvalue_plot()

    assert tmp.lr == pytest.approx(LR)

    assert tmp.xlabel == 'Likelihood Ratio'
    assert tmp.ylabel == 'Frequency'
    assert tmp.title == 'Likelihood Ratio Distribution'

    assert tmp.ratios.shape == (NUM, )
    assert tmp.xlo.shape == (NBINS + 1, )
    assert tmp.xhi.shape == (NBINS + 1, )
    assert tmp.y.shape == (NBINS + 1, )

    assert tmp.ratios.shape == (20, )
    assert tmp.xlo.shape == (9, )
    assert tmp.xhi.shape == (9, )
    assert tmp.y.shape == (9, )

    # The code may use parallel_map_rng which would mean this test is
    # not guaranteed to be repeatable.
    #
    # assert tmp.y == pytest.approx([0.35, 0.35, 0, 0, 0.2, 0.05, 0, 0, 0.05])


@requires_xspec
@requires_group
@requires_fits
@requires_data
def test_plot_pvalue_with_wstat(make_data_path, clean_astro_ui):
    """Check plot_pvalue with PHA data and wstat

    These values being tested were generated by Sherpa, so they only
    check that the code is doing the same thing. It is not clear what
    the correct values should be here.

    """

    ui.set_rng(np.random.RandomState(123))

    ui.set_stat('wstat')
    ui.set_method("neldermead")

    fname = make_data_path("3c273.pi")
    with SherpaVerbosity("ERROR"):
        ui.load_pha(fname)
        ui.notice(0.5, 6)

    ui.set_source("xsphabs.abs1 * xspowerlaw.p1")

    # move the fit close to the best fit to save a small amount
    # of time.
    abs1.nh = 0.05
    p1.phoindex = 2.05
    p1.norm = 2.1e-4

    # The number of iterations has been reduced so it's less than the
    # fit time. This leads to using a small value for the number of
    # bins.
    #
    ui.fit()
    ui.plot_pvalue(p1, abs1 * p1, num=20, bins=8)

    LR = 1.353942679847087

    tmp = ui.get_pvalue_results()

    assert tmp.null == pytest.approx(37.216613841917265)
    assert tmp.alt == pytest.approx(35.86267116207018)
    assert tmp.lr == pytest.approx(LR)

    assert tmp.samples.shape == (20, 2)
    assert tmp.stats.shape == (20, 2)
    assert tmp.ratios.shape == (20, )

    tmp = ui.get_pvalue_plot()

    assert tmp.lr == pytest.approx(LR)

    assert tmp.xlabel == 'Likelihood Ratio'
    assert tmp.ylabel == 'Frequency'
    assert tmp.title == 'Likelihood Ratio Distribution'

    assert tmp.ratios.shape == (20, )
    assert tmp.xlo.shape == (9, )
    assert tmp.xhi.shape == (9, )
    assert tmp.y.shape == (9, )

    # The code may use parallel_map_rng which would mean this test is
    # not guaranteed to be repeatable.
    #
    # assert tmp.y == pytest.approx([0.5, 0.1, 0, 0.1, 0.15, 0.05, 0.05, 0, 0.05])


@requires_xspec
@requires_group
@requires_fits
@requires_data
def test_plot_pvalue_with_bkg(make_data_path, clean_astro_ui):
    """Check plot_pvalue with PHA data.

    This is known to fail: #1377

    Similar to test_plot_value but uses a different dataset.
    As we don't know what the results should be we treat the failure
    as the "correct" value so that when it gets fixed the test will
    point out it needs to be updated.

    """

    ui.set_rng(np.random.RandomState(123))
    np.random.seed(123)

    ui.set_stat('cstat')
    ui.set_method("neldermead")

    fname = make_data_path("3c273.pi")
    with SherpaVerbosity("ERROR"):
        ui.load_pha(fname)
        ui.notice(0.5, 6)

    ui.set_source("xsphabs.abs1 * xspowerlaw.p1")
    ui.set_bkg_source("const1d.bmdl")

    # move the fit close to the best fit to save a small amount
    # of time.
    abs1.nh = 0.05
    p1.phoindex = 2.05
    p1.norm = 2.1e-4

    bmdl.c0.min = 0
    bmdl.c0 = 3.8e-6

    # The number of iterations has been reduced so it's less than the
    # fit time. This leads to using a small value for the number of
    # bins.
    #
    ui.fit()
    with pytest.raises(StatErr,
                       match=r"size mismatch between number of data sets \(2\) and model expressions \(1\)"):
        ui.plot_pvalue(p1, abs1 * p1, num=20, bins=8)


@pytest.fixture
def setup_imgdata_model():
    """Use a model for the PSF"""

    # Fake an image
    #
    x1, x0 = np.mgrid[0:8, 0:10]
    ymod = 10 + 100 / ((x0 - 5.5)**2 + (x1 - 3.5)**2)

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = ymod.flatten()

    y = y.astype(int)  # convert to integers

    ui.load_arrays(1, x0, x1, y, (8, 10), ui.DataIMG)

    pmodel = ui.create_model_component('gauss2d', 'pmodel')
    pmodel.xpos = 0
    pmodel.ypos = 0
    pmodel.fwhm = 3

    ui.load_psf('psf', pmodel)
    psf.size = [10, 10]  # not sure if this is useful but leave in
    psf.center = [5, 4]
    ui.set_psf(psf)

    setup_imgdata_source()


def setup_imgdata_source():
    """Set the model and fit setup"""
    ui.set_source(ui.gauss2d.g1 + ui.const2d.c1)

    c1.c0 = 10
    g1.fwhm = 4
    g1.xpos = 4
    g1.ypos = 3
    g1.ampl = 100

    ui.set_stat('cstat')
    ui.set_method("simplex")

    # Fit and check we are in the expected place (so this is a
    # pre-condition for the tests).
    #
    ui.fit()
    assert ui.calc_stat() == pytest.approx(362.9264257040166)


def fake_psf_image(outfile):
    """Create a PSF image

    This is based on the output of image_psf() saved to a file
    for the test_plot_pvalue_imgpsf_convolved case. As this is
    a very-small file and it could change (to make the test more
    meaningful) we do not add the data to the repository but
    create it directly.

    """

    hdr = ['SIMPLE  =                    T /',
           'BITPIX  =                  -64 / Bits per pixel',
           'NAXIS   =                    2 /',
           'NAXIS1  =                   10 /',
           'NAXIS2  =                    8 /']

    x1, x0 = np.mgrid[0:8, 0:10]
    dx0 = (x0 - 5).flatten()
    dx1 = (x1 - 4).flatten()

    mdl = Gauss2D()
    mdl.fwhm = 3
    y = mdl(dx0, dx1)
    y = y / y.max()

    # Create the FITS structure
    #
    with open(outfile, "wb") as fh:
        for h in hdr:
            fh.write(f'{h:80s}'.encode('ascii'))

        n = 2880 - (len(hdr) + 1) * 80

        # assume we don't have enough lines to pass 1 block!
        assert n > 0

        fmt = f'{{:{n}s}}'
        fh.write(fmt.format(' ').encode('ascii'))

        fh.write('{:80s}'.format('END').encode('ascii'))

        # FITS is big-endian
        yb = y.astype('>f8').tobytes()
        fh.write(yb)

        n = 2880 - len(yb)
        assert n > 0

        fh.write(b'\0' * n)


@pytest.fixture
def setup_imgdata_file(tmp_path, recwarn):
    """Use a model for the PSF"""

    # Fake an image
    #
    x1, x0 = np.mgrid[0:8, 0:10]
    ymod = 10 + 100 / ((x0 - 5.5)**2 + (x1 - 3.5)**2)

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = ymod.flatten()

    y = y.astype(int)  # convert to integers

    ui.load_arrays(1, x0, x1, y, (8, 10), ui.DataIMG)

    # Sherpa can confuse astropy if you send in pathlib objects so
    # convert to the raw filesystem path.
    #
    fname = tmp_path / "psf.fits"
    infile = str(fname)
    fake_psf_image(infile)

    ui.load_psf('psf', infile)
    psf.size = [10, 8]
    psf.center = [5, 4]

    # Just check where the warning comes from (pytest.warns didn't
    # seem to remove the warning so it still triggered the general
    # "check no warning" fixture, so use recwarn like this).
    #
    assert len(recwarn) == 0
    ui.set_psf(psf)
    assert len(recwarn) == 1

    setup_imgdata_source()

    assert len(recwarn) == 1
    w = recwarn.pop()
    assert issubclass(w.category, UserWarning)

    wmsg = 'PSF Image does not have a pixel size. Sherpa ' + \
        'will assume the pixel size is the same as the data'
    assert str(w.message) == wmsg


def check_imgdata_unconvolved(caplog):
    """What is the behavior when we do not add the PSF to plot_pvalue?

    Note we add a check of the screen output here.
    """

    # include a check of the screen output
    #
    with caplog.at_level('INFO', logger='sherpa'):
        with SherpaVerbosity('INFO'):
            ui.plot_pvalue(c1, c1 + g1, num=40, bins=5)

    assert len(caplog.records) == 1

    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == 'sherpa.ui.utils'
    assert lvl == logging.INFO

    # Do not use equality tests for the numeric values in case
    # there are numpy-version differences in the number of
    # significant figures.
    #
    toks = msg.split('\n')
    assert len(toks) == 5
    assert toks[0] == 'Likelihood Ratio Test'
    assert toks[1].startswith('null statistic   =  2391.2696')
    assert toks[2].startswith('alt statistic    =  353.82')
    assert toks[3].startswith('likelihood ratio =  2037.446')
    assert toks[4] == 'p-value          <  0.025'

    tmp = ui.get_pvalue_results()

    assert tmp.null == pytest.approx(2391.2696310023503)
    assert tmp.alt == pytest.approx(353.8235336370698)
    assert tmp.lr == pytest.approx(2037.4460973652804)

    assert tmp.samples.shape == (40, 1)
    assert tmp.stats.shape == (40, 2)
    assert tmp.ratios.shape == (40, )

    tmp = ui.get_pvalue_plot(2037.4460973652804)

    assert tmp.lr == pytest.approx(2037.4460973652804)

    assert tmp.xlabel == 'Likelihood Ratio'
    assert tmp.ylabel == 'Frequency'
    assert tmp.title == 'Likelihood Ratio Distribution'

    assert tmp.ratios.shape == (40, )
    assert tmp.xlo.shape == (6, )
    assert tmp.xhi.shape == (6, )
    assert tmp.y.shape == (6, )


def check_imgdata_convolved():
    """What is the behavior when we add the PSF to plot_pvalue?"""

    r1 = ui.get_psf()
    ui.plot_pvalue(c1, c1 + g1, conv_model=r1, num=40, bins=5)

    tmp = ui.get_pvalue_results()

    # these values are different to test_plot_pvalue_imgpsf_unconvolved
    #
    assert tmp.null == pytest.approx(2391.26963100235)
    assert tmp.alt == pytest.approx(563.3992697080881)
    assert tmp.lr == pytest.approx(1827.8703612942618)

    assert tmp.samples.shape == (40, 1)
    assert tmp.stats.shape == (40, 2)
    assert tmp.ratios.shape == (40, )

    tmp = ui.get_pvalue_plot()

    assert tmp.lr == pytest.approx(1827.8703612942618)

    assert tmp.xlabel == 'Likelihood Ratio'
    assert tmp.ylabel == 'Frequency'
    assert tmp.title == 'Likelihood Ratio Distribution'

    assert tmp.ratios.shape == (40, )
    assert tmp.xlo.shape == (6, )
    assert tmp.xhi.shape == (6, )
    assert tmp.y.shape == (6, )


@requires_fits
@requires_data
def test_plot_pvalue_imgpsf_model_unconvolved(clean_astro_ui, hide_logging,
                                              setup_imgdata_model, caplog,
                                              plot_backends):
    """Test of issue #1214 but with no explicit convolution

    This is an extra check added while working on 1214:
    this does not add an explicit convolution model.

    I've taken the liberty to add an explicit check of
    the screen output with this routine which we don't
    do in the convolved case.

    Since this test checks for a warning that is only triggered in the
    plotting process itself, it only makes sense for functional plotting
    backends, not for dummies that skip the plotting process.
    """

    check_imgdata_unconvolved(caplog)


@requires_fits
@requires_data
def test_plot_pvalue_imgpsf_model_convolved(clean_astro_ui, hide_logging, setup_imgdata_model):
    """Test of issue #1214

    It is not obvious if this is a meaningful test. It also does not
    trigger #1214.

    """

    check_imgdata_convolved()


@requires_fits
@requires_data
def test_plot_pvalue_imgpsf_file_unconvolved(clean_astro_ui, hide_logging,
                                             setup_imgdata_file, caplog,
                                             plot_backends):
    """Added as part of #1214 but just to exercise the system.

    As with the imgpsf_model version we check the screen
    output, which should be the same.

    Since this test checks for an exception that is only triggered in the
    plotting process itself, it only makes sense for functional plotting
    backends, not for dummies that skip the plotting process.
    """

    check_imgdata_unconvolved(caplog)


@requires_fits
@requires_data
def test_plot_pvalue_imgpsf_file_convolved(clean_astro_ui, hide_logging, setup_imgdata_file):
    """Test of issue #1214. This does trigger the problem."""

    check_imgdata_convolved()
