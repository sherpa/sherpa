#
#  Copyright (C) 2020 - 2024
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

"""Test handling of fitting a background to PHA data.

There are a few tests scattered around that do this, but this
attempts to provide a solid test base. There is likely to be
overlap in the WSTAT tests, for one.
"""

from io import StringIO
import logging

import numpy as np

import pytest

from sherpa.astro import ui
from sherpa.astro.data import DataARF, DataPHA
from sherpa.astro.instrument import ARFModelPHA, create_arf, create_delta_rmf
from sherpa.models.model import ArithmeticConstantModel
from sherpa.utils.err import DataErr, IdentifierErr, ModelErr
from sherpa.utils.logging import SherpaVerbosity
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

    assert sdata.grouping.dtype.type == np.int16
    assert bdata.grouping.dtype.type == np.int16

    assert sdata.quality.dtype.type == np.int16
    assert bdata.quality.dtype.type == np.int16

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

    # Using the actual exposure time can make the check below
    # sensitive to how the value is represented, in particular
    # with various Python versions. Since the model code doesn't
    # make any guarantees about the precision, let's use as
    # (hopefully) easy-to-display value. It's also done for
    # the BACKSCAL values to pick a nice ratio.
    #
    ui.set_exposure(id, 1234.5)
    ui.set_exposure(id, 5432.1, bkg_id=1)

    ui.set_backscal(id, 0.1)
    ui.set_backscal(id, 0.2, bkg_id=1)

    assert ui.list_model_components() == []

    # Check the source models
    #
    ui.set_source(id, ui.powlaw1d.pl)
    assert ui.list_model_components() == ['pl']

    ui.set_bkg_source(id, ui.powlaw1d.bpl)
    assert ui.list_model_components() == ['bpl', 'pl']

    smdl = ui.get_source(id)
    assert smdl.name == 'powlaw1d.pl'
    assert ui.list_model_components() == ['bpl', 'pl']

    bmdl = ui.get_bkg_source(id)
    assert bmdl.name == 'powlaw1d.bpl'
    assert ui.list_model_components() == ['bpl', 'pl']

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf(1234.5 * (powlaw1d.pl + 0.5 * powlaw1d.bpl)))'
    assert ui.list_model_components() == ['bpl', 'pl']

    bmdl = ui.get_bkg_model(id)
    assert bmdl.name == 'apply_rmf(apply_arf(5432.1 * powlaw1d.bpl))'
    assert ui.list_model_components() == ['bpl', 'pl']


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
    sdata.backscal = 0.1
    sdata.areascal = 0.8

    # Check that the ARF exposure isn't used
    sdata.get_arf().exposure = 15000

    # scale factor to correct to source is
    #   = (100 / 1000) * (0.1 / 0.4) * (0.8 / 0.4)
    #     exp = 0.1   back = 0.25  area = 2
    #   = 0.25  if only have BACKSCAL
    #   = 0.05  if include all scaling
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
    #   = 0.125   if only have BACKSCAL
    #   = 0.01    if include all scaling
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

    # Now try the "model" model:
    # - fail when only 1 background dataset has a model
    #
    iid = 1 if id is None else id
    estr = f"^background model 2 for data set {iid} has not been set$"
    with pytest.raises(ModelErr, match=estr):
        ui.get_model(id)

    bmdl = ui.get_bkg_model(id)
    assert bmdl.name == 'apply_rmf(apply_arf(1000.0 * powlaw1d.bpl))'

    ui.set_bkg_source(id, ui.polynom1d.bpl2, bkg_id=2)

    bmdl = ui.get_bkg_model(id, bkg_id=1)
    assert bmdl.name == 'apply_rmf(apply_arf(1000.0 * powlaw1d.bpl))'

    bmdl = ui.get_bkg_model(id, bkg_id=2)
    assert bmdl.name == 'apply_rmf(apply_arf(2000.0 * polynom1d.bpl2))'

    smdl = ui.get_model(id)

    toks = smdl.name.split('(')
    assert toks[0] == 'apply_rmf'
    assert toks[1] == 'apply_arf'
    assert toks[2] == '100.0 * '

    # The scale factors here are half of the values above, as there
    # are now two background components. The ordering used to depend
    # on the Python version but this should no-longer be an issue.
    #
    x1 = 'powlaw1d.pl + 0.125 * powlaw1d.bpl + 0.0625 * polynom1d.bpl2)))'
    assert toks[3] == x1

    assert len(toks) == 4

    assert ui.list_model_components() == ['bpl', 'bpl2', 'pl']


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

    bmdl = ui.get_bkg_model(id, bkg_id=1)
    assert bmdl.name == 'apply_rmf(apply_arf(1000.0 * powlaw1d.bpl))'

    bmdl = ui.get_bkg_model(id, bkg_id=2)
    assert bmdl.name == 'apply_rmf(apply_arf(2000.0 * powlaw1d.bpl))'

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf(100.0 * (powlaw1d.pl + 0.1875 * powlaw1d.bpl)))'

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
                            None])])
@pytest.mark.parametrize("id", [None, 1, "x"])
def test_pha1_show_data(id, exps, bscales, ascales, results, clean_astro_ui):
    """Check we can show the data and get the scaling.

    This *only* checks the background scaling, not the
    full output. Since the scaling value isn't shown
    when it's a vector we don't need to care about the values
    """

    iid = 1 if id is None else id

    ui.set_data(id, setup_pha1(exps, bscales, ascales))

    out = StringIO()
    ui.show_data(outfile=out)
    msg = out.getvalue().split('\n')

    assert msg[0] == 'Data Set: {}'.format(iid)
    assert msg[1] == 'Filter: 0.6000-2.5000 Energy (keV)'

    def check(line, result, bkg_id=None):
        out = 'Bkg Scale'
        if bkg_id is not None:
            out += ' {}'.format(bkg_id)

        out += ': '
        if result is None:
            out += 'float64[19]'
        else:
            out += '{}'.format(result)

        assert line == out

    nbkg = len(results)
    notice = 2 + nbkg
    if len(results) == 1:
        check(msg[2], results[0])
    else:
        for i, r in enumerate(results, 1):
            check(msg[1 + i], r, bkg_id=i)

    # just check we have the correct following lines
    assert msg[notice] == 'Noticed Channels: 1-19'
    assert msg[notice + 1] == 'name           = tst0'


def test_pha1_show_data_missing(clean_astro_ui):
    """Check we can show the data and get the scaling (with a missing component).

    This does not check all the cases that
    test_pha1_show_data does.
    """

    exps = (10, 20, 100)
    bscales = (0.1, 0.2, 0.4)
    ascales = (1, 1, 1)
    ui.set_data('x', setup_pha1(exps, bscales, ascales))

    d = ui.get_data('x')
    d.delete_background(id=1)

    out = StringIO()
    ui.show_data('x', outfile=out)
    msg = out.getvalue().split('\n')

    assert msg[0] == 'Data Set: x'
    assert msg[1] == 'Filter: 0.6000-2.5000 Energy (keV)'
    assert msg[2] == 'Bkg Scale 2: 0.025'
    assert msg[3] == 'Noticed Channels: 1-19'
    assert msg[4] == 'name           = tst0'


def test_pha1_show_data_no_bkg(clean_astro_ui):
    """Check we don't get a scaling value if there is no background,
    """

    ui.load_arrays(2, np.arange(1, 4), np.ones(3), ui.DataPHA)

    out = StringIO()
    ui.show_data(2, outfile=out)
    msg = out.getvalue().split('\n')

    assert msg[0] == 'Data Set: 2'
    assert msg[1] == 'Filter: 1-3 Channel'
    assert msg[2] == 'Noticed Channels: 1-3'
    assert msg[3] == 'name           = '


def test_pha1_instruments(clean_astro_ui):
    """Check we get the correct responses for various options.
    """

    # We don't check the scaling here
    scales = (1, 1, 1)
    ui.set_data(1, setup_pha1(scales, scales, scales))

    s = ui.stephi1d.s
    b = ui.steplo1d.b
    ui.set_source(1, s)
    ui.set_bkg_source(1, b)
    ui.set_bkg_source(1, b, bkg_id=2)

    ssrc = ui.get_source()
    bsrc1 = ui.get_bkg_source(bkg_id=1)
    bsrc2 = ui.get_bkg_source(bkg_id=2)

    assert ssrc == s
    assert bsrc1 == b
    assert bsrc2 == b

    smdl = ui.get_model()
    bmdl1 = ui.get_bkg_model(bkg_id=1)
    bmdl2 = ui.get_bkg_model(bkg_id=2)

    # check the response contains the correct names
    for i, mdl in enumerate([smdl, bmdl1, bmdl2]):
        assert isinstance(mdl, ARFModelPHA)
        assert mdl.pha.name == 'tst{}'.format(i)
        assert mdl.arf.name == 'arf{}'.format(i)


@pytest.mark.parametrize('bid', [1, 2])
def test_pha1_instruments_missing(bid, clean_astro_ui):
    """Check we get the correct responses for various options.
    """

    # We don't check the scaling here
    scales = (1, 1, 1)
    ui.set_data(1, setup_pha1(scales, scales, scales))

    # Why can we not say ui.stephi1d.bmdl here?
    #
    bmdl = ui.create_model_component('stephi1d', 'bmdl')
    ui.set_bkg_source(1, bmdl)
    ui.set_bkg_source(1, bmdl, bkg_id=2)

    bkg = ui.get_bkg(1, bkg_id=bid)
    bkg.delete_response()

    # We test the pha name field as a simple way to check what
    # PHA is being used here.
    #
    # Uses the source dataset for the response.
    bmdl = ui.get_bkg_model(bkg_id=bid)
    assert isinstance(bmdl, ARFModelPHA)
    assert bmdl.pha.name == "tst0"

    # Uses the background dataset for the response.
    oid = 1 if bid == 2 else 2
    bmdl = ui.get_bkg_model(bkg_id=oid)
    assert isinstance(bmdl, ARFModelPHA)
    assert bmdl.pha.name == f"tst{oid}"


# Try and pick routines with different pathways through the code
# to catch all cases.
#
@pytest.mark.parametrize("func", [ui.calc_stat, ui.get_model_plot])
def test_evaluation_requires_models(func, clean_astro_ui):
    """We need a background model for each dataset"""

    exps = (100.0, 1000.0, 1500.0)
    bscales = (0.01, 0.05, 0.02)
    ascales = (0.8, 0.4, 0.5)
    ui.set_data(setup_pha1(exps, bscales, ascales))

    ui.set_source(ui.box1d.smdl)
    ui.set_bkg_source(ui.box1d.bmdl)

    with pytest.raises(ModelErr,
                       match='^background model 2 for data set 1 has not been set$'):
        func()


SCALING = np.ones(19)
SCALING[2:5] = 0.8
SCALING[15] = 0.7


@pytest.mark.parametrize("exps,bscales,ascales,bkg_id,result",
                         [((100, 1000), (0.01, 0.05), (0.4, 0.5),
                           1, (100 * 0.01 * 0.4) / (1000 * 0.05 * 0.5)),
                          ((100, 1000, 1500), (0.01, 0.05, 0.04), (0.2, 0.5, 0.4),
                           1, 0.5 * (100 * 0.01 * 0.2) / (1000 * 0.05 * 0.5)),
                          ((100, 1000, 1500), (0.01, 0.05, 0.04), (0.2, 0.5, 0.4),
                           2, 0.5 * (100 * 0.01 * 0.2) / (1500 * 0.04 * 0.4)),
                          ((100, 1000, 1500), (0.01, 0.05 * SCALING, 0.04), (0.4, 0.5, 0.6 * SCALING),
                           1, 0.5 * (100 * 0.01 * 0.4) / (1000 * 0.05 * SCALING * 0.5)),
                          ((100, 1000, 1500), (0.01, 0.05 * SCALING, 0.04), (0.4, 0.5, 0.6 * SCALING),
                           2, 0.5 * (100 * 0.01 * 0.4) / (1500 * 0.04 * 0.6 * SCALING))])
@pytest.mark.parametrize("id", [None, "x"])
def test_pha1_get_bkg_scale(id, exps, bscales, ascales, bkg_id, result, clean_astro_ui):
    """Check we can calculate the scaling factor"""

    ui.set_data(id, setup_pha1(exps, bscales, ascales))
    bscal = ui.get_bkg_scale(id, bkg_id)
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
    assert bplot1.title == 'Background Model Contribution'

    # check the model evaluates correctly
    # - background is just multiplied by the arf
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

    # sy += r1(*exps) * r1(*bscales) * r1(*ascales) * (100 / 90) * by1 / 2
    # sy += r2(*exps) * r2(*bscales) * r2(*ascales) * (100 / 80) * by2 / 2

    sy += r1(*bscales) * (100 / 90) * by1 / 2
    sy += r2(*bscales) * (100 / 80) * by2 / 2

    assert splot.y == pytest.approx(sy)
    assert bplot1.y == pytest.approx(by1)

    bplot2 = ui.get_bkg_model_plot(bkg_id=2)
    assert bplot2.title == 'Background Model Contribution'
    assert bplot2.y == pytest.approx(by2)


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

    assert ui.list_model_components() == ['bmdl1', 'smdl']

    bmdl = ui.get_bkg_model()
    assert bmdl.name == 'apply_arf(1000.0 * box1d.bmdl1)'

    assert ui.list_model_components() == ['bmdl1', 'smdl']

    smdl = ui.get_model()

    src = 'apply_arf(100.0 * box1d.smdl)' + \
        ' + scale1 * apply_arf(100.0 * box1d.bmdl1)'

    assert smdl.name == src

    assert ui.list_model_components() == ['bmdl1', 'smdl']

    # deconstruct the model to extract the components of the scale factor
    #
    scale = smdl.parts[1].parts[0]
    assert isinstance(scale, ArithmeticConstantModel)

    def r(sval, bval):
        return sval / bval

    # expected = r(*exps) * r(*bscales) * r(*ascales)
    expected = r(*bscales)
    assert scale.val == pytest.approx(expected)


def test_multi_show_data(clean_astro_ui):
    """What happens if have PHA and non-PHA data?

    The behavior is somewhat different if we use show_bkg: see
    test_multi_show_bkg
    """

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
    ui.load_arrays("pha", [1, 2, 3], [1, 2, 3], ui.DataPHA)
    ui.load_arrays(7, [1, 2, 3], [1, 2, 3], ui.DataPHA)

    bkg = ui.DataPHA("x", [1, 2, 3], [12, 13, 0])
    ui.set_bkg(7, bkg)

    out = StringIO()
    ui.show_data(outfile=out)
    msg = out.getvalue().split('\n\n')

    # There are four useful chunks. We do not bother checking the
    # individual data, just that we get the data in an 'expected'
    # order.
    #
    assert len(msg) == 5

    assert msg[0].startswith("Data Set: 1\n")
    assert msg[1].startswith("Data Set: 7\n")
    assert msg[2].startswith("Background Data Set: 7:1\n")
    assert msg[3].startswith("Data Set: pha\n")
    assert msg[4] == "\n"


def test_multi_show_bkg(clean_astro_ui):
    """This used to fail because of #1645
    """

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
    ui.load_arrays("pha", [1, 2, 3], [1, 2, 3], ui.DataPHA)
    ui.load_arrays(7, [1, 2, 3], [1, 2, 3], ui.DataPHA)

    bkg = ui.DataPHA("x", [1, 2, 3], [12, 13, 0])
    ui.set_bkg(7, bkg)

    out = StringIO()
    ui.show_bkg(outfile=out)

    msg = out.getvalue().split('\n\n')

    # We do not do a full check, just enough to check the data is
    # what we expect.
    #
    assert len(msg) == 2

    assert msg[0].startswith("Background Data Set: 7:1\nFilter: 1-3 Channel\n")
    assert msg[1] == "\n"


def test_multi_show_bkg_explicit(clean_astro_ui):
    """Ensure coverage
    """

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
    ui.load_arrays("pha", [1, 2, 3], [1, 2, 3], ui.DataPHA)
    ui.load_arrays(7, [1, 2, 3], [1, 2, 3], ui.DataPHA)

    bkg = ui.DataPHA("x", [1, 2, 3], [12, 13, 0])
    ui.set_bkg(7, bkg)

    out = StringIO()
    ui.show_bkg(7, bkg_id=1, outfile=out)

    # We just want to check the correct data is being shown and do not
    # do a full check here.
    #
    msg = out.getvalue().split('\n\n')
    assert len(msg) == 2
    assert msg[0].startswith("Background Data Set: 7:1\n")
    assert msg[1] == "\n"


def test_pha1_eval_vector_show_two(clean_astro_ui):
    """Check we can show the source/bgnd models with vector scaling

    test_pha1_eval_vector_show but with two background
    components with different scaling (but the same model).

    """

    scale1 = np.ones(19)
    scale2 = np.ones(19)
    scale1[9:15] = 0.8
    scale2[3:12] = 0.7

    exps = (100.0, 1000.0, 750.0)
    bscales = (0.01, 0.05 * scale1, 0.02 * scale2)
    ascales = (0.8, 0.4, 0.2)
    ui.set_data(setup_pha1(exps, bscales, ascales))

    ui.set_source(ui.box1d.smdl)
    ui.set_bkg_source(ui.box1d.bmdl1)
    ui.set_bkg_source(bmdl1, bkg_id=2)

    assert ui.list_model_components() == ['bmdl1', 'smdl']

    bmdl = ui.get_bkg_model()
    assert bmdl.name == 'apply_arf(1000.0 * box1d.bmdl1)'

    bmdl = ui.get_bkg_model(bkg_id=2)
    assert bmdl.name == 'apply_arf(750.0 * box1d.bmdl1)'

    assert ui.list_model_components() == ['bmdl1', 'smdl']

    smdl = ui.get_model()

    src = 'apply_arf(100.0 * box1d.smdl)' + \
        ' + scale1 * apply_arf(100.0 * box1d.bmdl1)'

    assert smdl.name == src

    assert ui.list_model_components() == ['bmdl1', 'smdl']

    # deconstruct the model to extract the components of the scale factor
    #
    scale = smdl.parts[1].parts[0]
    assert isinstance(scale, ArithmeticConstantModel)

    def r1(sval, bval1, bval2):
        return sval / bval1

    def r2(sval, bval1, bval2):
        return sval / bval2

    # expected = (r1(*exps) * r1(*bscales) * r1(*ascales) + \
    #             r2(*exps) * r2(*bscales) * r2(*ascales)) / 2
    expected = (r1(*bscales) + r2(*bscales)) / 2
    assert scale.val == pytest.approx(expected)


def test_pha1_eval_vector_show_two_separate(clean_astro_ui):
    """Check we can show the source/bgnd models with vector scaling

    test_pha1_eval_vector_show_two with different background
    models.

    """

    scale1 = np.ones(19)
    scale2 = np.ones(19)
    scale1[9:15] = 0.8
    scale2[3:12] = 0.7

    exps = (100.0, 1000.0, 750.0)
    bscales = (0.01, 0.05 * scale1, 0.02 * scale2)
    ascales = (0.8, 0.4, 0.2)
    ui.set_data('foo', setup_pha1(exps, bscales, ascales))

    ui.set_source('foo', ui.box1d.smdl)
    ui.set_bkg_source('foo', ui.box1d.bmdl1)
    ui.set_bkg_source('foo', ui.const1d.bmdl2, bkg_id=2)

    assert ui.list_model_components() == ['bmdl1', 'bmdl2', 'smdl']

    bmdl = ui.get_bkg_model('foo')
    assert bmdl.name == 'apply_arf(1000.0 * box1d.bmdl1)'

    bmdl = ui.get_bkg_model('foo', bkg_id=2)
    assert bmdl.name == 'apply_arf(750.0 * const1d.bmdl2)'

    assert ui.list_model_components() == ['bmdl1', 'bmdl2', 'smdl']

    smdl = ui.get_model('foo')

    # The ordering of this test can depend on the Python version
    # (before Python 3.7).
    #
    src = 'apply_arf(100.0 * box1d.smdl)'

    src1 = src
    src1 += ' + scalefoo_1 * apply_arf(100.0 * box1d.bmdl1)'
    src1 += ' + scalefoo_2 * apply_arf(100.0 * const1d.bmdl2)'

    src2 = src
    src2 += ' + scalefoo_1 * apply_arf(100.0 * const1d.bmdl2)'
    src2 += ' + scalefoo_2 * apply_arf(100.0 * box1d.bmdl1)'

    assert smdl.name in [src1, src2]

    assert ui.list_model_components() == ['bmdl1', 'bmdl2', 'smdl']

    # deconstruct the model to extract the components of the scale factor
    # (and because of the ordering it requires some work to get
    # the path pre Python 3.7)
    #
    lterm, rterm = smdl.parts

    scale1 = lterm.parts[1].parts[0]
    scale2 = rterm.parts[0]
    assert isinstance(scale1, ArithmeticConstantModel)
    assert isinstance(scale2, ArithmeticConstantModel)

    def r1(sval, bval1, bval2):
        return sval / bval1

    def r2(sval, bval1, bval2):
        return sval / bval2

    # expected1 = r1(*exps) * r1(*bscales) * r1(*ascales) / 2
    # expected2 = r2(*exps) * r2(*bscales) * r2(*ascales) / 2

    # assert scale1.val == pytest.approx(expected1)
    # assert scale2.val == pytest.approx(expected2)

    expected1 = r1(*bscales) / 2
    expected2 = r2(*bscales) / 2

    check1_a = scale1.val == pytest.approx(expected1)
    check2_a = scale2.val == pytest.approx(expected2)

    check1_b = scale1.val == pytest.approx(expected2)
    check2_b = scale2.val == pytest.approx(expected1)

    assert (check1_a and check2_a) or (check1_b and check2_b)


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

    # Check the evaluation of the source (with instrument)
    #
    splot = ui.get_model_plot()
    bplot1 = ui.get_bkg_model_plot()

    assert splot.title == 'Model'
    assert bplot1.title == 'Background Model Contribution'

    # check the model evaluates correctly
    # - background is just multiplied by the arf
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

    # sy += r(*exps) * r(*bscales) * r(*ascales) * (100 / 90) * by1
    sy += r(*bscales) * (100 / 90) * by1

    # Just check a problem I encountered during development
    ndim = splot.y.ndim
    assert ndim == 1

    assert splot.y == pytest.approx(sy)
    assert bplot1.y == pytest.approx(by1)


@pytest.mark.parametrize('dofilter,expected',
                         [(False, 15.059210673609382),
                          (True, 4.344584636218002)])
def test_pha1_eval_vector_two(dofilter, expected, clean_astro_ui):
    """Compare statistic, with and without filtering.

    test_pha1_eval_vector but with two background components.

    The statistic values are calculated from the code, and
    so only act as a regression test.
    """

    scale1 = np.ones(19)
    scale2 = np.ones(19)
    scale1[9:15] = 0.8
    scale2[3:12] = 0.7

    exps = (100.0, 1000.0, 750.0)
    bscales = (0.01, 0.05 * scale1, 0.02 * scale2)
    ascales = (0.8, 0.4, 0.2)
    ui.set_data(setup_pha1(exps, bscales, ascales))

    ui.set_source(ui.box1d.smdl)
    ui.set_bkg_source(ui.box1d.bmdl1)
    ui.set_bkg_source(bmdl1, bkg_id=2)

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

    # check the second background component
    # (have to do after the first since the same object is returned)
    #
    bplot2 = ui.get_bkg_source_plot(bkg_id=2)

    assert bplot2.title == 'Source Model of tst2'

    # As the model is plotted as a rate it has the same values as
    # bkg_id=1.
    #
    assert bplot2.y == pytest.approx(by1)

    # Check the evaluation of the source (with instrument)
    #
    splot = ui.get_model_plot()
    bplot1 = ui.get_bkg_model_plot()

    assert splot.title == 'Model'
    assert bplot1.title == 'Background Model Contribution'

    # check the model evaluates correctly
    # - background is just multiplied by the arf
    # - source needs to include the scaled background
    #
    sy *= 100

    # need to correct by exposure time, backscal,
    # area scaling, and to correct for the ARF used to
    # calculate by, and then divide by the number of
    # backgrounds.
    #
    def r1(sval, bval1, bval2):
        return sval / bval1

    def r2(sval, bval1, bval2):
        return sval / bval2

    # sy += r1(*exps) * r1(*bscales) * r1(*ascales) * 100 * by1 / 2
    # sy += r2(*exps) * r2(*bscales) * r2(*ascales) * 100 * by1 / 2

    sy += r1(*bscales) * 100 * by1 / 2
    sy += r2(*bscales) * 100 * by1 / 2

    # Just check a problem I encountered during development
    assert splot.y.ndim == 1

    assert splot.y == pytest.approx(sy)
    assert bplot1.y == pytest.approx(by1 * 90)

    bplot2 = ui.get_bkg_model_plot(bkg_id=2)

    assert bplot2.title == 'Background Model Contribution'

    assert bplot2.y == pytest.approx(by1 * 80)


SCALING2 = np.ones(19)
SCALING2[9:15] = 0.8


@pytest.mark.parametrize('exps,bscales,ascales,dofilter,expected',
                         [((100.0, 1000.0), (0.01, 0.05 * SCALING2), (0.8, 0.4), False, 2742.484820710452),
                          ((100.0, 1000.0), (0.01, 0.05 * SCALING2), (0.8, 0.4), True, 764.242724090022),
                          ((100.0, 1000.0, 750.0), (0.01, 0.05 * SCALING2, 0.02 * SCALING), (0.8, 0.4, 0.2), False, 53935156.203287244),
                          ((100.0, 1000.0, 750.0), (0.01, 0.05 * SCALING2, 0.02 * SCALING), (0.8, 0.4, 0.2), True, 28780409.397312585)])
def test_pha1_eval_vector_stat(exps, bscales, ascales,
                               dofilter, expected, clean_astro_ui):
    """Compare statistic, with and without filtering.

    The statistic values are calculated from the code, and
    so only act as a regression test.
    """

    ui.set_data(setup_pha1(exps, bscales, ascales))

    # change counts to reduce the statistic difference (but doesn't
    # make much difference with some cases, thanks to the made-up
    # nature of the test data).
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

    bmdl1 = ui.create_model_component('box1d', 'bmdl1')
    for i, _ in enumerate(exps[1:], 1):
        ui.set_bkg_source(bmdl1, bkg_id=i)

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


def test_jdpileup_no_warning(caplog, clean_astro_ui):
    """jdpileup model has no warning when scalar scaling"""

    exps = (100.0, 1000.0)
    bscales = (0.01, 0.05)
    ascales = (0.8, 0.4)
    ui.set_data('x', setup_pha1(exps, bscales, ascales))

    ui.set_source('x', ui.box1d.smdl)
    ui.set_bkg_source('x', ui.box1d.bmdl1)

    ui.set_pileup_model('x', ui.jdpileup.jdp)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.get_model('x')

    assert len(caplog.records) == 0


def test_jdpileup_warning(caplog, clean_astro_ui):
    """jdpileup model has warning when vector scaling"""

    # The vector needs to be added to BACKSCAL
    exps = (100.0, 1000.0, 200)
    bscales = (0.01, 0.02, 0.05 * np.ones(19))
    ascales = (0.8, 0.8, 0.4)
    ui.set_data('x', setup_pha1(exps, bscales, ascales))

    ui.set_source('x', ui.box1d.smdl)
    ui.set_bkg_source('x', ui.box1d.bmdl1)
    ui.set_bkg_source('x', bmdl1, bkg_id=2)

    ui.set_pileup_model('x', ui.jdpileup.jdp)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.get_model('x')

    assert len(caplog.records) == 1
    name, level, msg = caplog.record_tuples[0]
    assert name == 'sherpa.astro.background'
    assert level == logging.WARNING
    assert msg == 'model results for dataset x likely wrong: use of pileup model and array scaling for the background'


def test_get_bkg_scale_nodata(clean_astro_ui):
    """Can we call get_bkg_scale with no background?"""

    ui.set_data(ui.DataPHA('foo', np.arange(3), np.arange(3)))

    with pytest.raises(DataErr,
                       match="^data set '1' does not have any associated backgrounds$"):
        ui.get_bkg_scale()


def test_get_bkg_scale_invalid(clean_astro_ui):
    """Can we call get_bkg_scale with invalid units?"""

    exps = (100.0, 1000.0)
    bscales = (0.01, 0.02)
    ascales = (0.8, 0.8)
    ui.set_data(setup_pha1(exps, bscales, ascales))

    with pytest.raises(ValueError,
                       match='^Invalid units argument: subtract$'):
        ui.get_bkg_scale(units='subtract')


def test_get_bkg_scale(clean_astro_ui):
    """Can we call get_bkg_scale with units=counts and rate?"""

    # Due to the choice of scalar vs vector values the scale
    # factor depends on whether units=counts or rate.
    #
    exps = (100.0, 1000.0, 200)
    bscales = (0.01, 0.02, 0.05 * np.ones(19))
    ascales = (0.8, 0.8 * np.ones(19), 0.4)
    ui.set_data('x', setup_pha1(exps, bscales, ascales))

    def r(vals, bkg):
        return vals[0] / vals[bkg]

    bscale1 = ui.get_bkg_scale('x')
    assert bscale1.size == 19
    assert bscale1 == pytest.approx(0.5 * r(exps, 1) * r(bscales, 1) * r(ascales, 1))

    bscale2 = ui.get_bkg_scale('x', bkg_id=2)
    assert bscale2.size == 19
    assert bscale2 == pytest.approx(0.5 * r(exps, 2) * r(bscales, 2) * r(ascales, 2))

    bscale1 = ui.get_bkg_scale('x', units='rate')
    assert np.isscalar(bscale1)
    assert bscale1 == pytest.approx(0.5 * r(bscales, 1))

    bscale2 = ui.get_bkg_scale('x', bkg_id=2, units='rate')
    assert bscale2.size == 19
    assert bscale2 == pytest.approx(0.5 * r(bscales, 2))


def check_stat_info(sinfo, ids, bkg_ids, numpoints, dof):
    """Check the stat-info structure."""

    assert sinfo.ids == ids
    if bkg_ids is None:
        assert sinfo.bkg_ids is None
    else:
        assert sinfo.bkg_ids == bkg_ids

    assert sinfo.numpoints == numpoints
    assert sinfo.dof == dof


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

    # Note: stats[2].dof > stats[0].dof + stats[1].dof here
    #
    # source dataset is fit with powlaw1d + const1d, so has
    # three free parameters.
    # background dataset is fit with const1d, so has 1 free
    # parameter.
    # the combined dataset is fit with powlaw1d+const1d, const1d
    # so has three free parameters
    #
    assert stats[0].name == 'Dataset 1'
    assert stats[1].name == 'Background 1 for Dataset 1'
    assert stats[2].name == 'Dataset [1]'

    check_stat_info(stats[0], (1, ), None, 42, 39)
    check_stat_info(stats[1], (1, ), (1, ), 42, 41)
    check_stat_info(stats[2], [1], None, 84, 81)


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

    assert stats[0].name == 'Dataset 1'
    assert stats[1].name == 'Background 1 for Dataset 1'
    assert stats[2].name == 'Dataset [1]'

    check_stat_info(stats[0], (1, ), None, 446, 443)
    check_stat_info(stats[1], (1, ), (1, ), 446, 445)
    check_stat_info(stats[2], [1], None, 892, 889)


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

    assert stats[0].name == 'Dataset 1'
    assert stats[1].name == 'Background 1 for Dataset 1'
    assert stats[2].name == 'Dataset [1]'

    check_stat_info(stats[0], (1, ), None, 42, 39)
    check_stat_info(stats[1], (1, ), (1, ), 106, 105)
    check_stat_info(stats[2], [1], None, 148, 145)


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

    assert stats[0].name == 'Dataset 1'
    assert stats[1].name == 'Background 1 for Dataset 1'
    assert stats[2].name == 'Background 2 for Dataset 1'
    assert stats[3].name == 'Dataset [1]'

    check_stat_info(stats[0], (1, ), None, 42, 39)
    check_stat_info(stats[1], (1, ), (1, ), 106, 105)
    check_stat_info(stats[2], (1, ), (2, ), 36, 35)
    check_stat_info(stats[3], [1], None, 184, 181)


def test_get_stat_info_multi_backgrounds(clean_astro_ui):
    """A complex case to check the ordering. Testing #1721

    Three PHA datasets:
      - first with no background
      - second with one background
      - third with two backgrounds

    """

    # Create the data sets
    #   a - 2 channels, no background
    #   b - 3 channels, one background
    #   c - 4 channels, two backgrounds
    #
    ui.dataspace1d(1, 2, id="a", dstype=ui.DataPHA)
    ui.dataspace1d(1, 3, id="b", dstype=ui.DataPHA)
    ui.dataspace1d(1, 3, id="b", bkg_id=1, dstype=ui.DataPHA)
    ui.dataspace1d(1, 4, id="c", dstype=ui.DataPHA)
    ui.dataspace1d(1, 4, id="c", bkg_id="up", dstype=ui.DataPHA)
    ui.dataspace1d(1, 4, id="c", bkg_id="down", dstype=ui.DataPHA)

    # set some data
    #
    ui.set_counts("a", [1, 2])
    ui.set_counts("b", [1, 2, 1])
    ui.set_counts("b", [1, 0, 1], bkg_id=1)
    ui.set_counts("c", [2, 1, 0, 2])
    ui.set_counts("c", [0, 0, 1, 1], bkg_id="up")
    ui.set_counts("c", [1, 1, 0, 0], bkg_id="down")

    # Create the responses. Use the same grid and just subset
    # for the different datasets.
    #
    egrid = np.asarray([0.1, 0.2, 0.4, 0.5, 0.8])

    elo = egrid[0:2]
    ehi = egrid[1:3]
    ui.set_rmf("a", create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi))

    elo = egrid[0:3]
    ehi = egrid[1:4]
    ui.set_rmf("b", create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi))
    ui.set_rmf("b", create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi), bkg_id=1)

    elo = egrid[0:4]
    ehi = egrid[1:5]
    ui.set_rmf("c", create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi))
    ui.set_rmf("c", create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi), bkg_id="up")
    ui.set_rmf("c", create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi), bkg_id="down")

    # Very simple model
    #
    ui.set_source("a", ui.const1d.smdl)
    ui.set_source("b", smdl)
    ui.set_source("c", smdl)

    ui.set_bkg_source("b", ui.const1d.bmdl, bkg_id=1)
    ui.set_bkg_source("c", bmdl, bkg_id="up")
    ui.set_bkg_source("c", bmdl, bkg_id="down")

    smdl.c0 = 0.2
    bmdl.c0 = 0.1
    stats = ui.get_stat_info()
    assert len(stats) == 7

    assert stats[0].name == 'Dataset a'
    assert stats[1].name == 'Dataset b'
    assert stats[2].name == 'Background 1 for Dataset b'
    assert stats[3].name == 'Dataset c'
    assert stats[4].name == 'Background up for Dataset c'
    assert stats[5].name == 'Background down for Dataset c'
    assert stats[6].name == "Datasets ['a', 'b', 'c']"

    check_stat_info(stats[0], ("a", ), None, 2, 1)
    check_stat_info(stats[1], ("b", ), None, 3, 1)
    check_stat_info(stats[2], ("b", ), (1, ), 3, 2)
    check_stat_info(stats[3], ("c", ), None, 4, 2)
    check_stat_info(stats[4], ("c", ), ("up", ), 4, 3)
    check_stat_info(stats[5], ("c", ), ("down", ), 4, 3)
    check_stat_info(stats[6], ["a", "b", "c"], None, 20, 18)


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


def validate_response_source_data_manual():
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

    pname = mdl.pha.name.split('/')[-1]
    assert pname == 'fake'

    aname = mdl.arf.name.split('/')[-1]
    assert aname == '12845.warf'

    rname = mdl.rmf.name.split('/')[-1]
    assert rname == '12845.wrmf'

    eterm = mdl.parts[0].parts[0]
    assert isinstance(eterm, ArithmeticConstantModel)
    assert eterm.val == pytest.approx(bexpval)

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

    validate_response_source_data_manual()


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

    emsg = '^No instrument response found for dataset ex$'
    if idval is None:
        with pytest.raises(DataErr, match=emsg):
            ui.set_analysis('energy')
    else:
        with pytest.raises(DataErr, match=emsg):
            ui.set_analysis(idval, 'energy')


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

    assert src.units == analysis
    assert bkg.units == analysis


@pytest.mark.parametrize("idval", [None, 1, "one"])
def test_partially_set_bkg_models(idval, clean_astro_ui):
    """We do not set a background model for all the backgrounds.

    Regression test what happens.
    """

    ui.set_stat("leastsq")
    ui.set_method("simplex")

    idopt = 1 if idval is None else idval
    ui.load_arrays(idopt, [1, 2, 3], [1, 1, 1], ui.DataPHA)
    ui.get_data(idopt).name = f"fake-{idval}"

    ui.set_bkg(idopt, DataPHA("up", [1, 2, 3], [0, 1, 1]), bkg_id="up")
    ui.set_bkg(idopt, DataPHA("down", [1, 2, 3], [1, 0, 0]), bkg_id="down")

    # Fake a response
    #
    egrid = np.asarray([0.2, 0.3, 0.4, 0.5])

    ui.set_arf(idopt, create_arf(egrid[:-1], egrid[1:]))

    # Set the source model
    #
    ui.set_source(idopt, ui.const1d.mdl)

    # Only set one of the background models
    ui.set_bkg_source(idopt, ui.const1d.bmdldn, bkg_id="down")

    emsg = f"^background model up for data set {idopt} has not been set$"
    with pytest.raises(ModelErr, match=emsg):
        ui.get_stat_info()

    with pytest.raises(ModelErr, match=emsg):
        ui.fit()

    # Note that fit_bkg will skip the bkg_id=up data set.
    #
    ui.fit_bkg()
    fres = ui.get_fit_results()
    assert fres.datasets == (idopt, )
    assert fres.succeeded
    assert fres.parnames == ("bmdldn.c0", )
    assert fres.numpoints == 3
    assert fres.dof == 2
    assert fres.istatval == pytest.approx(0.83)
    assert fres.statval == pytest.approx(2 / 3)
    assert fres.parvals == pytest.approx([3 + 1/3])


def test_fit_with_bkg_models_missing(clean_astro_ui):
    """Check what happens if a bkg model is missing.

    This is a regression test and should be compared to
    test_fit_bkg_with_bkg_models_missing.
    """

    ui.load_arrays(1, [1, 2], [2, 4], ui.DataPHA)
    ui.set_bkg(1, DataPHA("up", [1, 2], [2, 1]), bkg_id=1)
    ui.set_bkg(1, DataPHA("down", [1, 2], [2, 1]), bkg_id=2)

    egrid = np.asarray([0.2, 0.3, 0.4])
    arf = create_arf(egrid[:-1], egrid[1:])
    ui.set_arf(1, arf)
    ui.set_arf(1, arf, bkg_id=1)
    ui.set_arf(1, arf, bkg_id=2)

    ui.set_source(1, ui.const1d.smdl)
    ui.set_bkg_source(1, ui.const1d.smdl, bkg_id=1)

    with pytest.raises(ModelErr,
                       match="^background model 2 for data set 1 has not been set"):
        ui.fit()


def test_fit_bkg_with_bkg_models_missing(clean_astro_ui):
    """Check what happens if a bkg model is missing.

    This is a regression test and should be compared to
    test_fit_with_bkg_models_missing.
    """

    ui.set_stat("leastsq")

    ui.load_arrays(1, [1, 2], [2, 4], ui.DataPHA)
    ui.set_bkg(1, DataPHA("up", [1, 2], [2, 1]), bkg_id=1)
    ui.set_bkg(1, DataPHA("down", [1, 2], [2, 2]), bkg_id=2)

    egrid = np.asarray([0.2, 0.3, 0.4])
    arf = create_arf(egrid[:-1], egrid[1:])
    ui.set_arf(1, arf)
    ui.set_arf(1, arf, bkg_id=1)
    ui.set_arf(1, arf, bkg_id=2)

    # Note that we do not set a source model as it is not required
    # by fit_bkg:
    #    ui.set_source(1, ui.const1d.smdl)
    #
    ui.set_bkg_source(1, ui.const1d.bmdl, bkg_id=1)

    # Unlike ui.fit, this does not complain about a missing model for
    # bkg_id=2.
    #
    ui.fit_bkg()
    fres = ui.get_fit_results()
    assert fres.parnames == ("bmdl.c0", )
    assert fres.parvals == pytest.approx([15])


@requires_group
@requires_data
@requires_fits
def test_get_bkg_stat_info(make_data_path, clean_astro_ui, caplog):
    """Check out the "background" stats.

    """

    ui.set_stat("cstat")

    # We gave 4 observations we can use, but obs2 looks like the
    # spectral shape is different, so we load it in but do not fit it
    # (which adds a check for how we handle un-used data).
    #
    pl = ui.create_model_component("powlaw1d", "pl")
    bpl = ui.create_model_component("powlaw1d", "bpl")
    ngrp = 2

    with SherpaVerbosity("ERROR"):
        for idx in range(1, 5):
            ui.load_pha(idx, make_data_path(f"obs{idx}.pi"))
            if idx == 2:
                continue

            ui.set_source(idx, pl)
            ui.set_bkg_source(idx, bpl)

            ui.notice_id(idx, 1, 6)
            data = ui.get_data(idx)
            ui.group_counts(idx, ngrp, tabStops=~data.get_mask())

            ui.notice_id(idx, 1, 8, bkg_id=1)
            bkg = ui.get_bkg(idx, bkg_id=1)
            ui.group_counts(idx, ngrp, bkg_id=1, tabStops=~bkg.get_mask())

        ui.fit_bkg()
        ui.fit()

    expected_stats = [18.227080107057652,
                      34.552754127308,
                      17.105502808259327]
    expected_sum = 69.88533704262498
    assert ui.calc_bkg_stat() == pytest.approx(expected_sum)

    bstats = ui.get_bkg_stat_info()
    assert len(bstats) == 4
    for bstat, estat in zip(bstats, expected_stats):
        assert bstat.statval == pytest.approx(estat)

    assert bstats[3].statval == pytest.approx(expected_sum)

    for bstat, vs in zip(bstats,
                         [{"idx": 1, "num": 13, "dof": 11},
                          {"idx": 3, "num": 23, "dof": 21},
                          {"idx": 4, "num": 14, "dof": 12}
                          ]):
        assert bstat.name == f"Background 1 for Dataset {vs['idx']}"
        assert bstat.ids == (vs['idx'], )
        assert bstat.bkg_ids == (1,)
        assert bstat.numpoints == vs['num']
        assert bstat.dof == vs['dof']
        assert bstat.statname == "cstat"

    bstat = bstats[3]
    assert bstat.name == "Backgrounds for Datasets [1, 3, 4]"
    assert bstat.ids == [1, 3, 4]
    assert bstat.bkg_ids == (1,)
    assert bstat.numpoints == 50
    assert bstat.dof == 48
    assert bstat.statname == "cstat"

    # Just check we get some output from calc_bkg_stat_info()
    # but do not do too many checks.
    #
    assert len(caplog.records) == 0
    ui.calc_bkg_stat_info()
    assert len(caplog.records) == 1

    msg = caplog.records[0].getMessage()
    toks = msg.split("\n\n")
    assert len(toks) == 4
    assert toks[0].startswith("Background 1 in Dataset = 1\n")
    assert toks[1].startswith("Background 1 in Dataset = 3\n")
    assert toks[2].startswith("Background 1 in Dataset = 4\n")
    assert toks[3].startswith("Backgrounds in Datasets = [1, 3, 4]\n")
