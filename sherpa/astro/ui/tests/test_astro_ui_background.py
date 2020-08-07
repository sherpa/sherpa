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

from io import StringIO

import numpy as np

import pytest

from sherpa.astro import ui
from sherpa.astro.data import DataARF, DataPHA
from sherpa.astro.instrument import ARFModelPHA
from sherpa.models.model import ArithmeticConstantModel
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

    # check the UI get_grouping method agrees
    #
    assert (ui.get_grouping(id) == ui.get_grouping(id, bkg_id=1)).all()

    # The quality array is all 0's, so not an illuminating test
    assert (ui.get_quality(id) == ui.get_quality(id, bkg_id=1)).all()

    # Change the grouping of the background daga
    #
    ui.group_counts(id, 5, bkg_id=1)
    assert (sdata.grouping == sgrouping).all()
    assert bdata.grouping.sum() == -952

    assert (ui.get_grouping(id) != ui.get_grouping(id, bkg_id=1)).any()

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


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_models_two_single(id, make_data_path, clean_astro_ui, hide_logging):
    """The common use case is for the same background model

    See test_setup_pha1_file_models_two
    """

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    sdata = ui.get_data(id)
    sdata.exposure = 100
    sdata.get_arf().exposure = 150
    sdata.backscal = 0.1
    sdata.areascal = 0.8

    bdata1 = ui.get_bkg(id, bkg_id=1)
    bdata1.exposure = 1000
    bdata1.backscal = 0.4
    bdata1.areascal = 0.4

    bdata2 = ui.unpack_pha(make_data_path('3c273_bg.pi'))
    bdata2.exposure = 2000
    bdata2.backscal = 0.8
    bdata2.areascal = 0.5
    ui.set_bkg(id, bdata2, bkg_id=2)

    ui.set_source(id, ui.powlaw1d.pl)
    ui.set_bkg_source(id, ui.powlaw1d.bpl)
    ui.set_bkg_source(id, ui.powlaw1d.bpl, bkg_id=2)

    iid = 1 if id is None else id

    bmdl = ui.get_bkg_model(id, bkg_id=1)
    assert bmdl.name == 'apply_rmf(apply_arf((1000 * powlaw1d.bpl)))'

    bmdl = ui.get_bkg_model(id, bkg_id=2)
    assert bmdl.name == 'apply_rmf(apply_arf((2000 * powlaw1d.bpl)))'

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf((100 * (powlaw1d.pl + 0.03 * (powlaw1d.bpl + powlaw1d.bpl)))))'

    assert ui.list_model_components() == ['bpl', 'pl']


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


@pytest.mark.parametrize("exps,bscales,ascales,results",
                         [((100, 200), (2, 4), (0.5, 0.1),
                           [(100 * 2 * 1) / (200 * 4 * 1)]),
                          ((100, 200, 400), (2, 4, 0.5), (0.5, 0.1, 0.8),
                           [0.5 * (100 * 2 * 1) / (200 * 4 * 1),
                            0.5 * (100 * 2 * 1) / (400 * 0.5 * 1)]),
                          ((100, 200, 400),
                           (2, 4 * np.ones(19) * 0.8, 0.5 * np.ones(19) * 1.1),
                           (0.5, 0.1, 0.8),
                           [0.5 * (100 * 2 * 1) / (200 * 4 * np.ones(19) * 0.8 * 1),
                            0.5 * (100 * 2 * 1) / (400 * 0.5 * np.ones(19) * 1.1 * 1)]),
                          ])
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_pha1_subtract(id, exps, bscales, ascales, results, clean_astro_ui, hide_logging):
    """Check we can subtract the background.

    NOTE: the correction factor does NOT include the
          areascal correction - is this a bug?

    """

    ui.set_data(id, setup_pha1(exps, bscales, ascales))

    sdata = ui.get_dep(id)
    expected = np.zeros(19)
    for i, r in enumerate(results, 1):
        expected += r * ui.get_dep(id, bkg_id=i)

    expected = sdata - expected

    ui.subtract(id)
    data = ui.get_dep(id)

    assert (data < sdata).all()  # this just checks the subtraction does subtract
    assert data == pytest.approx(expected)


@pytest.mark.parametrize("exps,bscales,ascales,results",
                         [((100, 200), (2, 4), (0.5, 0.1),
                           [(100 * 2 * 0.5) / (200 * 4 * 0.1)]),
                          ((100, 200, 400), (2, 4, 0.5), (0.5, 0.1, 0.8),
                           [0.5 * (100 * 2 * 0.5) / (200 * 4 * 0.1),
                            0.5 * (100 * 2 * 0.5) / (400 * 0.5 * 0.8)]),
                          ((100, 200, 1.0), (2, 4, 1 * np.ones(19)), (0.5, 0.1, 1.0),
                           [0.5 * (100 * 2 * 0.5) / (200 * 4 * 0.1),
                            None]),
                         ])
@pytest.mark.parametrize("id", [None, 1, "x"])
def test_pha1_show_data(id, exps, bscales, ascales, results, clean_astro_ui):
    """Check we can show the data and get the scaling.

    This *only* checks the background scaling, not the
    full output.
    """

    iid = 1 if id is None else id

    ui.set_data(id, setup_pha1(exps, bscales, ascales))

    out = StringIO()
    ui.show_data(outfile=out)
    msg = out.getvalue().split('\n')

    assert msg[0] == 'Data Set: {}'.format(iid)
    assert msg[1] == 'Filter: 0.6500-2.4500 Energy (keV)'

    if None not in results:
        bscale = sum([r for r in results if r is not None])
        assert msg[2] == 'Bkg Scale: {}'.format(bscale)
        notice = 3
    else:
        notice = 2

    assert msg[notice] == 'Noticed Channels: 1-19'
    assert msg[notice + 1] == 'name           = tst0'


def test_pha1_instruments(clean_astro_ui):
    """Check we get the correct responses for various options.
    """

    # We don't check the scaling here
    scales = (1, 1, 1)
    ui.set_data(1, setup_pha1(scales, scales, scales))
    ui.set_source(1, ui.stephi1d.smdl)
    ui.set_bkg_source(1, ui.steplo1d.bmdl)
    ui.set_bkg_source(1, bmdl, bkg_id=2)

    smdl = ui.get_model()
    bmdl1 = ui.get_bkg_model(bkg_id=1)
    bmdl2 = ui.get_bkg_model(bkg_id=2)

    # check the response contains the correct names
    for i, mdl in enumerate([smdl, bmdl1, bmdl2]):
        assert isinstance(mdl, ARFModelPHA)
        assert mdl.pha.name == 'tst{}'.format(i)
        assert mdl.arf.name == 'arf{}'.format(i)


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
                          (True, 4.344584636218002)])
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


# Fails from ui.get_model raising
# TypeError: only size-1 arrays can be converted to Python scalars
#
@pytest.mark.xfail
def test_pha1_eval_vector_show(clean_astro_ui):
    """Check we can show the source/bgnd models with vector scaling
    test_pha1_eval does most of the work; this test is
    to check when there's a vector, not scalar, for
    scaling the background to match the source.
    """

    scale = np.ones(19)
    scale[9:15] = 0.8

    exps = (100.0, 1000.0)
    bscales = (0.01, 0.05 * scale)
    ascales = (0.8, 0.4)
    ui.set_data(setup_pha1(exps, bscales, ascales))

    ui.set_source(ui.box1d.smdl)
    ui.set_bkg_source(ui.box1d.bmdl1)

    def r(sval, bval):
        return sval / bval

    bmdl = ui.get_bkg_model()
    assert bmdl.name == 'apply_arf((1000.0 * box1d.bmdl1))'

    smdl = ui.get_model()

    array = r(*exps) * r(*bscales) * r(*ascales)
    src = '(apply_arf((100.0 * box1d.smdl))'
    src += ' + (apply_arf((100.0 * box1d.bmdl1))'
    src += ' * {}))'.format(array)

    assert smdl.name == src


# Fails from ui.get_model_plot raising
# TypeError: only size-1 arrays can be converted to Python scalars
#
@pytest.mark.xfail
def test_pha1_eval_vector(clean_astro_ui):
    """Check we can evaluate the source/bgnd values and vector scaling
    test_pha1_eval does most of the work; this test is
    to check when there's a vector, not scalar, for
    scaling the background to match the source.
    """

    scale = np.ones(19)
    scale[9:15] = 0.8

    exps = (100.0, 1000.0)
    bscales = (0.01, 0.05 * scale)
    ascales = (0.8, 0.4)
    ui.set_data(setup_pha1(exps, bscales, ascales))

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
    #
    sy = np.zeros(19)
    sy[3] = 5
    sy[4:9] = 10
    sy[9] = 9

    by1 = np.zeros(19)
    by1[1] = 1.2
    by1[2:11] = 2
    by1[11] = 0.2

    assert splot.y == pytest.approx(sy)
    assert bplot1.y == pytest.approx(by1)

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

    # need to correct by exposure time, backscal,
    # area scaling, and to correct for the ARF used to
    # calculate by, and then divide by the number of
    # backgrounds.
    #
    def r(sval, bval):
        return sval / bval

    sy += r(*exps) * r(*bscales) * r(*ascales) * (100 / 90) * by1

    assert splot.y == pytest.approx(sy)
    assert bplot1.y == pytest.approx(by1)


@requires_data
@requires_fits
def test_use_source_data(make_data_path, clean_astro_ui, hide_logging):
    """Check a fit uses the sousrce data to define its fit

    The source region binning and filtering is used here
    as nothing has been set for the background.

    It does not check the statistic value

    See also test_use_background_data
    """

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    # Use the source to define grouping and filtering
    ui.notice(0.5, 7)

    ui.set_source(ui.powlaw1d.smdl)
    ui.set_bkg_source(ui.const1d.bmdl)

    stats = ui.get_stat_info()
    assert len(stats) == 3

    src = stats[0]
    bgnd = stats[1]
    comb = stats[2]

    assert src.name == 'Dataset 1'
    assert bgnd.name == 'Background 1 for Dataset 1'
    assert comb.name == 'Dataset [1]'

    assert src.ids == (1,)
    assert bgnd.ids == (1,)
    assert comb.ids == [1]

    assert src.bkg_ids is None
    assert bgnd.bkg_ids == (1,)
    assert comb.bkg_ids is None

    assert src.numpoints == 42
    assert bgnd.numpoints == 42
    assert comb.numpoints == 84

    # Note: comb.dof > src.dof + bgnd.dof here
    #
    # source dataset is fit with powlaw1d + const1d, so has
    # three free parameters.
    # background dataset is fit with const1d, so has 1 free
    # parameter.
    # the combined dataset is fit with powlae1d+const1d, const1d
    # so has three free parameters
    #
    assert src.dof == 39
    assert bgnd.dof == 41
    assert comb.dof == 81


@requires_data
@requires_fits
def test_use_source_data_manual(make_data_path, clean_astro_ui, hide_logging):
    """Check a fit uses the source data to define its fit (load in bkg)

    The source region binning and filtering is used here
    as nothing has been set for the background. Unlike test_use_source_data
    we load in the background separately (ie not via BACKFILE keyword)
    to see if this makes any difference.

    It does not check the statistic value

    See also test_use_background_data
    """

    infile = make_data_path('12845.pi')
    ui.load_pha(infile)

    assert ui.list_bkg_ids() == []

    # Add in a background
    bpha = ui.unpack_pha(infile)
    bpha.delete_response()

    ui.set_bkg(bpha)

    assert ui.list_bkg_ids() == [1]

    # Use the source to define grouping and filtering
    ui.notice(0.5, 7)

    ui.set_source(ui.powlaw1d.smdl)
    ui.set_bkg_source(ui.const1d.bmdl)

    stats = ui.get_stat_info()
    assert len(stats) == 3

    src = stats[0]
    bgnd = stats[1]
    comb = stats[2]

    assert src.name == 'Dataset 1'
    assert bgnd.name == 'Background 1 for Dataset 1'
    assert comb.name == 'Dataset [1]'

    assert src.ids == (1,)
    assert bgnd.ids == (1,)
    assert comb.ids == [1]

    assert src.bkg_ids is None
    assert bgnd.bkg_ids == (1,)
    assert comb.bkg_ids is None

    assert src.numpoints == 446
    assert bgnd.numpoints == 446
    assert comb.numpoints == 892

    assert src.dof == 443
    assert bgnd.dof == 445
    assert comb.dof == 889


@requires_group
@requires_data
@requires_fits
def test_use_background_data(make_data_path, clean_astro_ui, hide_logging):
    """Check a fit uses the background data to define its fit

    See also test_use_source_data and test_use_backgroud_data_two
    """

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    bexpval = ui.get_exposure(bkg_id=1)

    # Use different grouping and filtering for source and
    # background apertures.
    #
    ui.group_width(5, bkg_id=1)

    ui.notice(0.5, 7)

    ui.notice(None, None, bkg_id=1)
    ui.notice(0.3, 8, bkg_id=1)

    ui.set_source(ui.powlaw1d.smdl)
    ui.set_bkg_source(ui.const1d.bmdl)

    stats = ui.get_stat_info()
    assert len(stats) == 3

    src = stats[0]
    bgnd = stats[1]
    comb = stats[2]

    assert src.name == 'Dataset 1'
    assert bgnd.name == 'Background 1 for Dataset 1'
    assert comb.name == 'Dataset [1]'

    assert src.ids == (1,)
    assert bgnd.ids == (1,)
    assert comb.ids == [1]

    assert src.bkg_ids is None
    assert bgnd.bkg_ids == (1,)
    assert comb.bkg_ids is None

    assert src.numpoints == 42
    assert bgnd.numpoints == 106
    assert comb.numpoints == 148

    assert src.dof == 39
    assert bgnd.dof == 105
    assert comb.dof == 145


@requires_group
@requires_data
@requires_fits
def test_use_background_data_two(make_data_path, clean_astro_ui, hide_logging):
    """Check a fit uses the background data to define its fit (2 backgrounds)

    See also test_use_source_data and  and test_use_backgroud_data
    """

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    # Fake a second background, and use one with a response
    # (this one doesn't have a background)
    bgfile = make_data_path('12845.pi')
    ui.load_bkg(bgfile, bkg_id=2)

    # Use different grouping and filtering for source and
    # background apertures.
    #
    ui.group_width(5, bkg_id=1)
    ui.group_width(10, bkg_id=2)

    ui.notice(0.5, 7)

    ui.notice(None, None, bkg_id=1)
    ui.notice(0.3, 8, bkg_id=1)

    ui.notice(None, None, bkg_id=2)
    ui.notice(1, 6, bkg_id=2)

    ui.set_source(ui.powlaw1d.smdl)
    ui.set_bkg_source(ui.const1d.bmdl)
    ui.set_bkg_source(bmdl, bkg_id=2)

    stats = ui.get_stat_info()
    assert len(stats) == 4

    src = stats[0]
    bgnd1 = stats[1]
    bgnd2 = stats[2]
    comb = stats[3]

    assert src.name == 'Dataset 1'
    assert bgnd1.name == 'Background 1 for Dataset 1'
    assert bgnd2.name == 'Background 2 for Dataset 1'
    assert comb.name == 'Dataset [1]'

    assert src.ids == (1,)
    assert bgnd1.ids == (1,)
    assert bgnd2.ids == (1,)
    assert comb.ids == [1]

    assert src.bkg_ids is None
    assert bgnd1.bkg_ids == (1,)
    assert bgnd2.bkg_ids == (2,)
    assert comb.bkg_ids is None

    assert src.numpoints == 42
    assert bgnd1.numpoints == 106
    assert bgnd2.numpoints == 36
    assert comb.numpoints == 184

    assert src.dof == 39
    assert bgnd1.dof == 105
    assert bgnd2.dof == 35
    assert comb.dof == 181


@requires_data
@requires_fits
def test_response_source_data(make_data_path, clean_astro_ui, hide_logging):
    """Check the full model expressions"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    bexpval = ui.get_exposure(bkg_id=1)

    ui.set_source(ui.powlaw1d.smdl)
    ui.set_bkg_source(ui.const1d.bmdl)

    # What is used as the background response?
    #
    mdl = ui.get_bkg_model()

    pname = mdl.pha.name.split('/')[-1]
    assert pname == '3c273_bg.pi'

    aname = mdl.arf.name.split('/')[-1]
    assert aname == '3c273.arf'

    rname = mdl.rmf.name.split('/')[-1]
    assert rname == '3c273.rmf'

    # Check the exposure time (this is sensitive to the
    # structure of the model expression).
    #
    eterm = mdl.parts[0].parts[0]
    assert isinstance(eterm, ArithmeticConstantModel)
    assert eterm.val == pytest.approx(bexpval)

    # Source response
    #
    mdl = ui.get_model()

    pname = mdl.pha.name.split('/')[-1]
    assert pname == '3c273.pi'

    aname = mdl.arf.name.split('/')[-1]
    assert aname == '3c273.arf'

    rname = mdl.rmf.name.split('/')[-1]
    assert rname == '3c273.rmf'

    # The exposure time is the same as the background
    eterm = mdl.parts[0].parts[0]
    assert isinstance(eterm, ArithmeticConstantModel)
    assert eterm.val == pytest.approx(bexpval)


def validate_response_source_data_manual(direct=True):
    """Tests for test_response_source_data_manual[_datapha]

    It checks out issue #880.
    """

    expval = ui.get_exposure()
    bexpval = ui.get_exposure(bkg_id=1)
    assert expval == pytest.approx(37845.662207644)
    assert bexpval == pytest.approx(37845.662207644 * 5)

    # What is used as the response?
    #
    mdl = ui.get_bkg_model()

    # The value here should be 'fake' but thanks to #880 it
    # depends. Rather than make this an XFAIL for one case,
    # adapt to the current behavior.
    #
    pname = mdl.pha.name.split('/')[-1]
    if direct:
        assert pname == 'fake'
    else:
        assert pname == '12845.pi'

    aname = mdl.arf.name.split('/')[-1]
    assert aname == '12845.warf'

    rname = mdl.rmf.name.split('/')[-1]
    assert rname == '12845.wrmf'

    eterm = mdl.parts[0].parts[0]
    assert isinstance(eterm, ArithmeticConstantModel)
    if direct:
        assert eterm.val == pytest.approx(bexpval)
    else:
        assert eterm.val == pytest.approx(expval)

    # Source response
    #
    mdl = ui.get_model()

    pname = mdl.pha.name.split('/')[-1]
    assert pname == '12845.pi'

    aname = mdl.arf.name.split('/')[-1]
    assert aname == '12845.warf'

    rname = mdl.rmf.name.split('/')[-1]
    assert rname == '12845.wrmf'

    # The exposure time is the same as not the same as the background
    eterm = mdl.parts[0].parts[0]
    assert isinstance(eterm, ArithmeticConstantModel)
    assert eterm.val == pytest.approx(expval)


@requires_data
@requires_fits
def test_response_source_data_manual(make_data_path, clean_astro_ui, hide_logging):
    """Check the full model expressions: explicitly load background"""

    infile = make_data_path('12845.pi')
    ui.load_pha(infile)

    assert ui.list_bkg_ids() == []

    expval = ui.get_exposure()

    # Change scaling and remove response from the fake background
    bpha = ui.unpack_pha(infile)
    bpha.name = 'fake'
    bpha.exposure *= 5
    bpha.backscal *= 2
    bpha.delete_response()

    ui.set_bkg(bpha)

    bexpval = ui.get_exposure(bkg_id=1)
    assert bexpval == pytest.approx(expval * 5)

    ui.set_source(ui.powlaw1d.smdl)
    ui.set_bkg_source(ui.const1d.bmdl)

    validate_response_source_data_manual()


@requires_data
@requires_fits
def test_response_source_data_manual_datapha(make_data_path, clean_astro_ui, hide_logging):
    """Check the full model expressions: explicitly load background: DataPHA

    This is test_response_source_data_manual but using the
    DataPHA method to add the background to the source, since it
    leads to different behavior with the background instrument
    response. This is a combination of issue #879 and #880
    """

    infile = make_data_path('12845.pi')
    ui.load_pha(infile)

    assert ui.list_bkg_ids() == []

    expval = ui.get_exposure()

    # Change scaling and remove response from the fake background
    bpha = ui.unpack_pha(infile)
    bpha.name = 'fake'
    bpha.exposure *= 5
    bpha.backscal *= 2
    bpha.delete_response()

    ui.get_data().set_background(bpha)

    bexpval = ui.get_exposure(bkg_id=1)
    assert bexpval == pytest.approx(expval * 5)

    ui.set_source(ui.powlaw1d.smdl)
    ui.set_bkg_source(ui.const1d.bmdl)

    validate_response_source_data_manual(direct=False)


@requires_data
@requires_fits
def test_response_two_backgrounds(make_data_path, clean_astro_ui, hide_logging):
    """Check the full model expressions: two backgrounds"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    # Fake a second background, and use one with a response
    # (this one doesn't have a background)
    bgfile = make_data_path('12845.pi')
    ui.load_bkg(bgfile, bkg_id=2)

    bexpval1 = ui.get_exposure(bkg_id=1)
    bexpval2 = ui.get_exposure(bkg_id=2)

    assert bexpval1 == pytest.approx(38564.608926889)
    assert bexpval2 == pytest.approx(37845.662207644)

    ui.set_source(ui.powlaw1d.smdl)
    ui.set_bkg_source(ui.const1d.bmdl)
    ui.set_bkg_source(bmdl, bkg_id=2)

    # What is used as the response?
    #
    mdl1 = ui.get_bkg_model()

    pname = mdl1.pha.name.split('/')[-1]
    assert pname == '3c273_bg.pi'

    aname = mdl1.arf.name.split('/')[-1]
    assert aname == '3c273.arf'

    rname = mdl1.rmf.name.split('/')[-1]
    assert rname == '3c273.rmf'

    eterm = mdl1.parts[0].parts[0]
    assert isinstance(eterm, ArithmeticConstantModel)
    assert eterm.val == pytest.approx(bexpval1)

    mdl2 = ui.get_bkg_model(bkg_id=2)

    pname = mdl2.pha.name.split('/')[-1]
    assert pname == '12845.pi'

    aname = mdl2.arf.name.split('/')[-1]
    assert aname == '12845.warf'

    rname = mdl2.rmf.name.split('/')[-1]
    assert rname == '12845.wrmf'

    eterm = mdl2.parts[0].parts[0]
    assert isinstance(eterm, ArithmeticConstantModel)
    assert eterm.val == pytest.approx(bexpval2)

    # Source response
    #
    mdl = ui.get_model()

    pname = mdl.pha.name.split('/')[-1]
    assert pname == '3c273.pi'

    aname = mdl.arf.name.split('/')[-1]
    assert aname == '3c273.arf'

    rname = mdl.rmf.name.split('/')[-1]
    assert rname == '3c273.rmf'

    # The exposure time is the same as the first background
    # component
    #
    eterm = mdl.parts[0].parts[0]
    assert isinstance(eterm, ArithmeticConstantModel)
    assert eterm.val == pytest.approx(bexpval1)


@requires_data
@requires_fits
def test_source_overwrites_background(make_data_path, clean_astro_ui, hide_logging):
    """Changing the source clears out the background filters"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)

    ui.notice(2, 4, bkg_id=1)

    assert ui.get_dep(filter=True).size == 46
    assert ui.get_dep(filter=True, bkg_id=1).size == 13

    # This over-writes the background filter
    ui.notice(0.5, 7)

    assert ui.get_dep(filter=True).size == 42
    assert ui.get_dep(filter=True, bkg_id=1).size == 42


def fake_pha(idval, direct, response=True):
    """Create a fake PHA file with a response and background"""

    chans = np.arange(3, dtype=np.int16)
    d = ui.DataPHA('ex', chans, chans * 0, exposure=20, backscal=0.1)
    b = ui.DataPHA('bg', chans, chans * 0, exposure=200, backscal=0.2)

    # add in a RMF
    if response:
        ebins = np.asarray([0.1, 0.2, 0.4, 0.5])
        elo = ebins[:-1]
        ehi = ebins[1:]
        r = ui.create_rmf(elo, ehi, e_min=elo, e_max=ehi)
        d.set_rmf(r)

    if direct:
        d.set_background(b)

    if idval is None:
        ui.set_data(d)
    else:
        ui.set_data(idval, d)

    if not direct:
        if idval is None:
            ui.set_bkg(b)
        else:
            ui.set_bkg(idval, b)


@pytest.mark.parametrize("idval", [None, 1, "one"])
@pytest.mark.parametrize("direct", [True, False])
def test_bkg_analysis_setting_no_response(idval, direct, clean_astro_ui):
    """There's no response, so setting the analysis will error out.
    """

    fake_pha(idval, direct, response=False)

    if idval is None:
        with pytest.raises(DataErr) as exc:
            ui.set_analysis('energy')
    else:
        with pytest.raises(DataErr) as exc:
            ui.set_analysis(idval, 'energy')

    assert str(exc.value) == 'No instrument model found for dataset ex'


@pytest.mark.parametrize("idval", [None, 1, "one"])
@pytest.mark.parametrize("direct", [True, False])
def test_bkg_analysis_setting_default(idval, direct, clean_astro_ui):
    """Check the analysis setting of the background.

    The default setting is channel

    """

    fake_pha(idval, direct)

    src = ui.get_data(idval)
    bkg = ui.get_bkg(idval)

    assert src.units == 'channel'
    assert bkg.units == 'channel'


@pytest.mark.parametrize("idval", [None, 1, "one"])
@pytest.mark.parametrize("direct", [True, False])
@pytest.mark.parametrize("analysis", ['channel', 'energy', 'wavelength'])
def test_bkg_analysis_setting_changed(idval, direct, analysis, clean_astro_ui):
    """Change the analysis setting of the background.

    This compares datapha.set_background to ui.set_bkg to
    check issue #879

    """

    fake_pha(idval, direct)

    if idval is None:
        ui.set_analysis(analysis)
    else:
        ui.set_analysis(idval, analysis)

    src = ui.get_data(idval)
    bkg = ui.get_bkg(idval)

    # Want to be able to mark some tests as XFAIL, but as it
    # involves a parameter setting, how best to do this and
    # get to see a XPASS when the code is fixed?
    #
    if direct and analysis != 'channel':
        if bkg.units == analysis:
            assert False, "Test was expected to fail, but can not mark XPASS"

        pytest.xfail("Using set_background does not set units")

    assert src.units == analysis
    assert bkg.units == analysis
