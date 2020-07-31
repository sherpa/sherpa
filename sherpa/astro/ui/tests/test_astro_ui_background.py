#
#  Copyright (C) 2020  Smithsonian Astrophysical Observatory
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

"""Test handling of fitting a background to PHA data.

There are a few tests scattered around that do this, but this
attempts to provide a solid test base. There is likely to be
overlap in the WSTAT tests, for one.
"""

import numpy as np

import pytest

from sherpa.astro import ui
from sherpa.astro.data import DataARF, DataPHA
from sherpa.utils.err import DataErr, ModelErr
from sherpa.utils.testing import requires_data, requires_fits, requires_group


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_filter(id, make_data_path, clean_astro_ui, hide_logging):
    """Can I change filtering of the background?"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    sdata = ui.get_data(id)
    bdata = ui.get_bkg(id)

    # Need an argument for notice/ignore_id
    iid = 1 if id is None else id

    # The idea here is to check that we can have
    # different filters for the source and background
    # regions.
    #
    assert sdata.get_dep(filter=True).size == 46
    assert bdata.get_dep(filter=True).size == 46

    ui.ignore_id(iid, None, 0.5)
    ui.ignore_id(iid, 7, None)

    assert sdata.get_dep(filter=True).size == 40
    assert bdata.get_dep(filter=True).size == 40

    # need to reset the background filter
    ui.notice_id(iid, None, None, bkg_id=1)

    assert sdata.get_dep(filter=True).size == 40
    assert bdata.get_dep(filter=True).size == 46

    ui.ignore_id(iid, None, 0.3, bkg_id=1)
    ui.ignore_id(iid, 8, None, bkg_id=1)

    assert sdata.get_dep(filter=True).size == 40
    assert bdata.get_dep(filter=True).size == 42


@requires_data
@requires_fits
@requires_group
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_group(id, make_data_path, clean_astro_ui, hide_logging):
    """Can I change grouping of the background?"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    sdata = ui.get_data(id)
    bdata = ui.get_bkg(id)

    dtype = np.dtype('>i2')  # why not np.dtype(np.int16)?
    assert sdata.grouping.dtype == dtype
    assert bdata.grouping.dtype == dtype

    assert sdata.quality.dtype == dtype
    assert bdata.quality.dtype == dtype

    assert sdata.grouped
    assert bdata.grouped

    # Check default grouping (taken from source)
    sgrouping = sdata.grouping.copy()
    assert sgrouping.sum() == -932
    assert (bdata.grouping == sgrouping).all()

    # Change the grouping of the background daga
    #
    ui.group_counts(id, 5, bkg_id=1)
    assert (sdata.grouping == sgrouping).all()
    assert bdata.grouping.sum() == -952

    # check we can ungroup just the background
    ui.ungroup(id, bkg_id=1)
    assert sdata.grouped
    assert not bdata.grouped

    assert bdata.grouping.sum() == -952


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_models(id, make_data_path, clean_astro_ui, hide_logging):
    """Is the model set correctly for the background"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    # Check the source models
    #
    ui.set_source(id, ui.powlaw1d.pl)
    ui.set_bkg_source(id, ui.powlaw1d.bpl)

    smdl = ui.get_source(id)
    assert smdl.name == 'powlaw1d.pl'

    bmdl = ui.get_bkg_source(id)
    assert bmdl.name == 'powlaw1d.bpl'

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf((38564.608926889 * (powlaw1d.pl + 0.134921 * (powlaw1d.bpl)))))'

    bmdl = ui.get_bkg_model(id)
    assert bmdl.name == 'apply_rmf(apply_arf((38564.608926889 * powlaw1d.bpl)))'


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_models_two(id, make_data_path, clean_astro_ui, hide_logging):
    """Is the model set correctly for two background components"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    # Tweak scaling values to make output easier to check.
    # The ARF is changed as well, to a different value, to
    # check what value is used in the model expression.
    #
    sdata = ui.get_data(id)
    sdata.exposure = 100
    sdata.get_arf().exposure = 150
    sdata.backscal = 0.1
    sdata.areascal = 0.8

    # scale factor to correct to source is
    #   = (100 / 1000) * (0.1 / 0.4) * (0.8 / 0.4)
    #     exp = 0.1   back = 0.25  area = 2
    #   = 0.025 if exclude the area scaling
    #   = 0.05  if include the area scaling
    #
    bdata1 = ui.get_bkg(id, bkg_id=1)
    bdata1.exposure = 1000
    # bdata1.get_arf().exposure = 1500  no ARF
    bdata1.backscal = 0.4
    bdata1.areascal = 0.4

    # add a second dataset (copy of the background)
    #
    # scale factor
    #   = (100 / 2000) * (0.1 / 0.8) * (0.8 / 0.5)
    #     exp = 0.05   back = 0.125  area = 1.6
    #   = 0.00625 if exclude the area scaling
    #   = 0.01    if include the area scaling
    #
    bdata2 = ui.unpack_pha(make_data_path('3c273_bg.pi'))
    bdata2.exposure = 2000
    # bdata2.get_arf().exposure = 2500  no ARF
    bdata2.backscal = 0.8
    bdata2.areascal = 0.5
    ui.set_bkg(id, bdata2, bkg_id=2)

    # Check the source models:
    #  - first with one
    #  - then with two
    #
    ui.set_source(id, ui.powlaw1d.pl)
    ui.set_bkg_source(id, ui.powlaw1d.bpl)

    smdl = ui.get_source(id)
    assert smdl.name == 'powlaw1d.pl'

    bmdl = ui.get_bkg_source(id)
    assert bmdl.name == 'powlaw1d.bpl'

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf((100 * (powlaw1d.pl + 0.03 * (powlaw1d.bpl)))))'

    bmdl = ui.get_bkg_model(id)
    assert bmdl.name == 'apply_rmf(apply_arf((1000 * powlaw1d.bpl)))'

    ui.set_bkg_source(id, ui.polynom1d.bpl2, bkg_id=2)

    bmdl = ui.get_bkg_model(id, bkg_id=1)
    assert bmdl.name == 'apply_rmf(apply_arf((1000 * powlaw1d.bpl)))'

    bmdl = ui.get_bkg_model(id, bkg_id=2)
    assert bmdl.name == 'apply_rmf(apply_arf((2000 * polynom1d.bpl2)))'

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf((100 * (powlaw1d.pl + 0.03 * (powlaw1d.bpl + polynom1d.bpl2)))))'


def setup_pha1(exps, bscals, ascals):
    """Create multiple PHA files (source + backgrounds).

    There's an ARF (for each dataset) but no RMF.

    """

    channels = np.arange(1, 20, dtype=np.int16)
    counts = np.ones(19, dtype=np.int16)

    egrid = 0.5 + np.arange(1, 21) / 10
    elo = egrid[:-1]
    ehi = egrid[1:]

    dset = None
    for i, (exp, bscal, ascal) in enumerate(zip(exps, bscals, ascals)):
        d = DataPHA('tst{}'.format(i),
                    channel=channels, counts=counts,
                    exposure=exp, backscal=bscal, areascal=ascal)

        a = DataARF('arf{}'.format(i), elo, ehi,
                    np.ones(19) * 100 - i * 10)
        d.set_response(arf=a)

        if dset is None:
            dset = d
        else:
            dset.set_background(d, id=i)

    dset.set_analysis('energy')
    return dset


# Try and pick routines with different pathways through the code
# to catch all cases.
#
@pytest.mark.parametrize("func,etype,emsg",
                         [(ui.calc_stat, ModelErr, "background model 2 for data set 1 has not been set"),
                          (ui.get_model_plot, DataErr, "background 2 has no associated model")])
def test_evaluation_requires_models(func, etype, emsg, clean_astro_ui):
    """We need a background model for each dataset"""

    exps = (100.0, 1000.0, 1500.0)
    bscales = (0.01, 0.05, 0.02)
    ascales = (0.8, 0.4, 0.5)
    ui.set_data(setup_pha1(exps, bscales, ascales))

    ui.set_source(ui.box1d.smdl)
    ui.set_bkg_source(ui.box1d.bmdl)

    with pytest.raises(etype) as exc:
        func()

    assert str(exc.value) == emsg


SCALING = np.ones(19)
SCALING[2:5] = 0.8
SCALING[15] = 0.7

@pytest.mark.parametrize("exps,bscales,ascales,result",
                         [((100, 1000), (0.01, 0.05), (0.4, 0.5),
                           (100 * 0.01 * 0.4) / (1000 * 0.05 * 0.5)),
                          ((100, 1000, 1500), (0.01, 0.05, 0.04), (0.2, 0.5, 0.4),
                           0.5 * (100 * 0.01 * 0.2) *
                           (1 / (1000 * 0.05 * 0.5) + 1 / (1500 * 0.04 * 0.4))),
                          ((100, 1000), (0.01, 0.05 * SCALING), (0.4, 0.5),
                           (100 * 0.01 * 0.4) / (1000 * 0.05 * SCALING * 0.5)),
                          ((100, 1000, 1500), (0.01, 0.05 * SCALING, 0.04), (0.4, 0.5, 0.6 * SCALING),
                           0.5 * (100 * 0.01 * 0.4) *
                           (1 / (1000 * 0.05 * SCALING * 0.5) +
                            1 / (1500 * 0.04 * 0.6 * SCALING))),
                         ])
@pytest.mark.parametrize("id", [None, "x"])
def test_pha1_get_bkg_scale(id, exps, bscales, ascales, result, clean_astro_ui):
    """Check we can calculate the scaling factor"""

    ui.set_data(id, setup_pha1(exps, bscales, ascales))
    bscal = ui.get_bkg_scale(id)
    assert bscal == pytest.approx(result)


def test_pha1_eval(clean_astro_ui):
    """Check we can evaluate the source/bgnd values.

    A faked-up PHA file is created, with only ARFs,
    and two background components.

    """
    exps = (100.0, 1000.0, 1500.0)
    bscales = (0.01, 0.05, 0.02)
    ascales = (0.8, 0.4, 0.5)
    ui.set_data(setup_pha1(exps, bscales, ascales))

    ui.set_source(ui.box1d.smdl)
    ui.set_bkg_source(ui.box1d.bmdl1)
    ui.set_bkg_source(ui.box1d.bmdl2, bkg_id=2)

    smdl.ampl.max = 10
    smdl.ampl = 10
    smdl.xlow = 0.95
    smdl.xhi = 1.59

    bmdl1.ampl.max = 2
    bmdl1.ampl = 2
    bmdl1.xlow = 0.74
    bmdl1.xhi = 1.71

    bmdl2.ampl = 1
    bmdl2.xlow = 1.85
    bmdl2.xhi = 2.21

    # Check the evaluation of the source (no instrument)
    #
    splot = ui.get_source_plot()
    bplot1 = ui.get_bkg_source_plot()

    assert splot.title == 'Source Model of tst0'
    assert bplot1.title == 'Source Model of tst1'

    # check the model evaluates correctly
    #
    # source: bins 3-9
    # bgnd 1:      1-11
    # bgnd 2:      12-16
    #
    sy = np.zeros(19)
    sy[3] = 5
    sy[4:9] = 10
    sy[9] = 9

    by1 = np.zeros(19)
    by1[1] = 1.2
    by1[2:11] = 2
    by1[11] = 0.2

    by2 = np.zeros(19)
    by2[12] = 0.5
    by2[13:16] = 1
    by2[16] = 0.1

    assert splot.y == pytest.approx(sy)
    assert bplot1.y == pytest.approx(by1)

    # only check bplot2 after finishing with bplot1 as the
    # return value actually refers to the same object.
    #
    bplot2 = ui.get_bkg_source_plot(bkg_id=2)
    assert bplot2.title == 'Source Model of tst2'
    assert bplot2.y == pytest.approx(by2)

    # Check the evaluation of the source (no instrument)
    #
    splot = ui.get_model_plot()
    bplot1 = ui.get_bkg_model_plot()

    assert splot.title == 'Model'
    assert bplot1.title == 'Model'

    # check the model evaluates correctly
    # - backgrond is just multiplied by the arf
    # - source needs to include the scaled background
    #
    sy *= 100
    by1 *= 90
    by2 *= 80

    # need to correct by exposure time, backscal,
    # area scaling, and to correct for the ARF used to
    # calculate by, and then divide by the number of
    # backgrounds.
    #
    def r1(sval, bval1, bval2):
        return sval / bval1

    def r2(sval, bval1, bval2):
        return sval / bval2

    sy += r1(*exps) * r1(*bscales) * r1(*ascales) * (100 / 90) * by1 / 2
    sy += r2(*exps) * r2(*bscales) * r2(*ascales) * (100 / 80) * by2 / 2

    assert splot.y == pytest.approx(sy)
    assert bplot1.y == pytest.approx(by1)

    bplot2 = ui.get_bkg_model_plot(bkg_id=2)
    assert bplot2.title == 'Model'
    assert bplot2.y == pytest.approx(by2)


# The calc_stat call fails because of
# TypeError: only size-1 arrays can be converted to Python scalars
#
# The current expected stats are "best guess" (based on developoing
# the code) and may need to be adjusted once the code is fixed.
#
@pytest.mark.xfail
@pytest.mark.parametrize('dofilter,expected',
                         [(False, 15.059210673609382),
                          (True, 2.6057323845799996)])
def test_pha1_eval_vector_stat(dofilter, expected, clean_astro_ui):
    """Compare statistic, with and without filtering.

    The statistic values are calculated from the code, and
    so only act as a regression test.
    """

    scale = np.ones(19)
    scale[9:15] = 0.8

    exps = (100.0, 1000.0)
    bscales = (0.01, 0.05 * scale)
    ascales = (0.8, 0.4)
    ui.set_data(setup_pha1(exps, bscales, ascales))

    # change counts to reduce the statistic difference
    #
    y = np.zeros(19, dtype=np.int16)
    y[1] = 4
    y[2] = 10
    y[3] = 507
    y[4:9] = 1007
    y[9] = 910
    y[10] = 11
    y[11] = 2
    ui.set_dep(y * 8)

    y = np.zeros(19, dtype=np.int16)
    y[1] = 110
    y[2:11] = 181
    y[11] = 17
    ui.set_dep(y * 40, bkg_id=1)

    ui.set_source(ui.box1d.smdl)
    ui.set_bkg_source(ui.box1d.bmdl1)

    smdl.ampl.max = 10
    smdl.ampl = 10
    smdl.xlow = 0.95
    smdl.xhi = 1.59

    bmdl1.ampl.max = 2
    bmdl1.ampl = 2
    bmdl1.xlow = 0.74
    bmdl1.xhi = 1.71

    ui.set_stat('chi2datavar')

    if dofilter:
        ui.ignore(None, 0.79)
        ui.ignore(1.38, None)

    s = ui.calc_stat()
    assert s == pytest.approx(expected)
