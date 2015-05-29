# 
#  Copyright (C) 2014  Smithsonian Astrophysical Observatory
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


from sherpa.utils import SherpaTestCase
import os
import unittest
from sherpa.utils import has_fits_support
from sherpa.astro.ui import *
from sherpa.astro.datastack import *
from sherpa.astro import datastack
from acis_bkg_model import acis_bkg_model
import numpy as np
import tempfile

logger = logging.getLogger('sherpa')


class test_design(SherpaTestCase):
    def setUp(self):
        clear_stack()
        ui.clean()
        set_template_id("ID")
        logger.setLevel(logging.ERROR)
        set_stack_verbosity(logging.ERROR)
        self._this_dir = os.path.dirname(sys.modules[self.__module__].__file__)

    def tearDown(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits')
    def test_case_1(self):
        datadir = '/'.join((self._this_dir, 'data'))
        ls = '@'+'/'.join((datadir, '3c273.lis'))
        load_pha(ls, use_errors=True)

        assert 2 == len(DATASTACK.datasets)
        assert 2 == len(ui._session._data)

        load_pha("myid", '/'.join((datadir, "3c273.pi")))

        assert 2 == len(DATASTACK.datasets)
        assert 3 == len(ui._session._data)

        load_pha('/'.join((datadir, "3c273.pi")))

        assert 3 == len(DATASTACK.datasets)
        assert 4 == len(ui._session._data)

        load_pha([], '/'.join((datadir, "3c273.pi")))

        assert 4 == len(DATASTACK.datasets)
        assert 5 == len(ui._session._data)

        ds = DataStack()

        load_pha(ds, ls)

        assert 4 == len(DATASTACK.datasets)
        assert 7 == len(ui._session._data)
        assert 2 == len(ds.datasets)

        load_pha(ds, '/'.join((datadir, "3c273.pi")))

        assert 4 == len(DATASTACK.datasets)
        assert 8 == len(ui._session._data)
        assert 3 == len(ds.datasets)

        dids = DATASTACK.get_stack_ids()
        assert dids == [1,2,3,4]

        sids = ui._session._data.keys()
        assert sids == [1,2,3,4,5,6,7, "myid"]

        set_source([1,2], "powlaw1d.pID")
        set_source([3,4], "brokenpowerlaw.bpID")

        dsids = ds.get_stack_ids()
        assert dsids == [5,6,7]

        p1 = ui._session._model_components['p1']
        p2 = ui._session._model_components['p2']
        bp3 = ui._session._model_components['bp3']
        bp4 = ui._session._model_components['bp4']

        assert p1 is not None
        assert p2 is not None
        assert bp3 is not None
        assert bp4 is not None

        set_source(1, "polynom1d.poly1")
        set_source([2,3,4], "atten.attID")

        poly1 = ui._session._model_components['poly1']
        a2 = ui._session._model_components['att2']
        a3 = ui._session._model_components['att3']
        a4 = ui._session._model_components['att4']

        assert poly1 is not None
        assert a2 is not None
        assert a3 is not None
        assert a4 is not None

        clean()

        assert 0 == len(DATASTACK.datasets)
        assert 0 == len(ui._session._data)
        assert 3 == len(ds.datasets)

class test_global(SherpaTestCase):
    def setUp(self):
        clear_stack()
        ui.clean()
        logger.setLevel(logging.ERROR)
        set_stack_verbosity(logging.ERROR)
        set_template_id("__ID")

    def tearDown(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")



    def test_case_2(self):
        x1 = np.arange(50)+100
        y1 = 2*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        load_arrays([[x1, y1], [x2, y2], [x3, y3]])

        set_source([], 'const1d.const * polynom1d.poly__ID')

        poly1 = ui._session._get_model_component('poly1')
        poly2 = ui._session._get_model_component('poly2')
        poly3 = ui._session._get_model_component('poly3')
        const = ui._session._get_model_component('const')

        freeze([], 'poly')
        thaw([], 'poly.c0')
        thaw([], 'poly.c1')
        thaw([], 'poly.c2')
        thaw([], 'const')

        assert poly1.c0.frozen is False
        assert poly1.c1.frozen is False
        assert poly1.c2.frozen is False
        assert poly1.c3.frozen is True

        assert poly2.c0.frozen is False
        assert poly2.c1.frozen is False
        assert poly2.c2.frozen is False
        assert poly2.c3.frozen is True

        assert poly3.c0.frozen is False
        assert poly3.c1.frozen is False
        assert poly3.c2.frozen is False
        assert poly3.c3.frozen is True

        set_par([], 'poly.c1', 0.45)

        assert poly1.c1.val == 0.45
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        set_par([1], 'poly.c1', 0.1)

        assert poly1.c1.val == 0.1
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        set_par([], 'const.c0', 2)

        assert const.c0.val == 2

        set_par([], 'const.integrate', False)
        freeze([], 'const.c0')

        vals = get_par([], 'poly.c1.val')
        assert ([0.1, 0.45, 0.45] == vals).all()

        pars = get_par([], 'const.c0')

        fit([])

        assert round(poly1.c1.val) == 1
        assert round(poly1.c2.val) == 3
        assert round(poly2.c1.val) == 3
        assert round(poly2.c2.val) == 1
        assert round(poly3.c1.val) == 3
        assert round(poly3.c2.val) == 1

        clear_stack()

        x1 = np.arange(50)+100
        y1 = 7*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        load_arrays([[x1, y1], [x2, y2], [x3, y3]])

        set_template_id("foo")

        set_source([], 'const1d.constfoo * polynom1d.polyfoo')

        const1 = ui._session._get_model_component('const1')
        const2 = ui._session._get_model_component('const2')
        const3 = ui._session._get_model_component('const3')

        link([2,3], 'const.c0')

        set_par([2], 'const.c0', 3)
        set_par([1], 'const.c0', 7)

        freeze([1], 'const.c0')

        assert const2.c0.frozen is False

        fit([])

        assert const2.c0.val == const3.c0.val
        assert const3.c0._link is const2.c0

        unlink([], "const.c0")

        assert const3.c0._link is not const2.c0


class test_load(SherpaTestCase):
    def setUp(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")
        self._this_dir = os.path.dirname(sys.modules[self.__module__].__file__)
        self.create_files()

    def tearDown(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")
        os.remove(self.lisname)
        os.remove(self.name1)
        os.remove(self.name2)
        set_stack_verbose(False)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits')
    def test_case_3(self):
        load_ascii("@{}".format(self.lisname))
        assert len(ui._session._data) == 2
        assert len(DATASTACK.datasets) == 2

        load_data("@{}".format(self.lisname))
        assert len(ui._session._data) == 4
        assert len(DATASTACK.datasets) == 4

        load_data("@"+"/".join((self._this_dir, 'data', 'pha.lis')))
        assert len(ui._session._data) == 6
        assert len(DATASTACK.datasets) == 6


    def create_files(self):
        fd1, self.name1 = tempfile.mkstemp()
        fd2, self.name2 = tempfile.mkstemp()

        fdlis, self.lisname = tempfile.mkstemp()

        ascii1 = open(self.name1, 'w')
        ascii2 = open(self.name2, 'w')
        lis = open(self.lisname, 'w')

        x1 = np.arange(50)+100
        y1 = 2*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)

        for x, y in zip(x1, y1):
            ascii1.write("{} {}\n".format(x, y))
        ascii1.close()
        os.close(fd1)

        for x, y in zip(x2, y2):
            ascii2.write("{} {}\n".format(x, y))
        ascii2.close()
        os.close(fd2)

        lis.write("{}\n".format(self.name1))
        lis.write("{}\n".format(self.name2))
        lis.close()
        os.close(fdlis)


class test_partial_oo(SherpaTestCase):

    def setUp(self):
        self.ds = DataStack()
        set_template_id("__ID")
        clear_stack()
        ui.clean()

    def tearDown(self):
        self.ds.clear_stack()
        clear_stack()
        set_template_id("__ID")
        ui.clean()

    def test_case_4(self):
        x1 = np.arange(50)+100
        y1 = 2*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        ds = self.ds

        load_arrays(ds, [[x1, y1], [x2, y2], [x3, y3]])

        set_source(ds, 'const1d.const * polynom1d.poly__ID')

        poly1 = ui._session._get_model_component('poly1')
        poly2 = ui._session._get_model_component('poly2')
        poly3 = ui._session._get_model_component('poly3')
        const = ui._session._get_model_component('const')

        freeze(ds, 'poly')
        thaw(ds, 'poly.c0')
        thaw(ds, 'poly.c1')
        thaw(ds, 'poly.c2')
        thaw(ds, 'const')

        assert poly1.c0.frozen is False
        assert poly1.c1.frozen is False
        assert poly1.c2.frozen is False
        assert poly1.c3.frozen is True

        assert poly2.c0.frozen is False
        assert poly2.c1.frozen is False
        assert poly2.c2.frozen is False
        assert poly2.c3.frozen is True

        assert poly3.c0.frozen is False
        assert poly3.c1.frozen is False
        assert poly3.c2.frozen is False
        assert poly3.c3.frozen is True

        set_par(ds, 'poly.c1', 0.45)

        assert poly1.c1.val == 0.45
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        set_par(ds[1], 'poly.c1', 0.1)

        assert poly1.c1.val == 0.1
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        set_par(ds, 'const.c0', 2)

        assert const.c0.val == 2

        set_par(ds, 'const.integrate', False)
        freeze(ds, 'const.c0')

        vals = get_par(ds, 'poly.c1.val')
        assert ([0.1, 0.45, 0.45] == vals).all()

        pars = get_par(ds, 'const.c0')

        fit(ds)

        assert round(poly1.c1.val) == 1
        assert round(poly1.c2.val) == 3
        assert round(poly2.c1.val) == 3
        assert round(poly2.c2.val) == 1
        assert round(poly3.c1.val) == 3
        assert round(poly3.c2.val) == 1

        ds.clear_stack()

        x1 = np.arange(50)+100
        y1 = 7*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        load_arrays(ds, [[x1, y1], [x2, y2], [x3, y3]])

        set_template_id("foo")

        set_source(ds, 'const1d.constfoo * polynom1d.polyfoo')

        const1 = ui._session._get_model_component('const1')
        const2 = ui._session._get_model_component('const2')
        const3 = ui._session._get_model_component('const3')

        link(ds[2,3], 'const.c0')

        set_par(ds[2], 'const.c0', 3)
        set_par(ds[1], 'const.c0', 7)

        freeze(ds[1], 'const.c0')

        assert const2.c0.frozen is False

        fit(ds)

        assert const2.c0.val == const3.c0.val
        assert const3.c0._link is const2.c0

        unlink(ds, "const.c0")

        assert const3.c0._link is not const2.c0


class test_oo(SherpaTestCase):

    def setUp(self):
        datastack.clear_stack()
        datastack.set_template_id("__ID")
        ui.clean()
        self.ds = datastack.DataStack()

    def tearDown(self):
        self.ds.clear_stack()
        datastack.clear_stack()
        datastack.set_template_id("__ID")
        ui.clean()

    def test_case_5(self):
        x1 = np.arange(50)+100
        y1 = 2*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        ds = self.ds

        ds.load_arrays([[x1, y1], [x2, y2], [x3, y3]])

        ds.set_source('const1d.const * polynom1d.poly__ID')

        poly1 = ui._session._get_model_component('poly1')
        poly2 = ui._session._get_model_component('poly2')
        poly3 = ui._session._get_model_component('poly3')
        const = ui._session._get_model_component('const')

        ds.freeze('poly')
        ds.thaw('poly.c0')
        ds.thaw('poly.c1')
        ds.thaw('poly.c2')
        ds.thaw('const')

        assert poly1.c0.frozen is False
        assert poly1.c1.frozen is False
        assert poly1.c2.frozen is False
        assert poly1.c3.frozen is True

        assert poly2.c0.frozen is False
        assert poly2.c1.frozen is False
        assert poly2.c2.frozen is False
        assert poly2.c3.frozen is True

        assert poly3.c0.frozen is False
        assert poly3.c1.frozen is False
        assert poly3.c2.frozen is False
        assert poly3.c3.frozen is True

        ds.set_par('poly.c1', 0.45)

        assert poly1.c1.val == 0.45
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        ds[1].set_par('poly.c1', 0.1)

        assert poly1.c1.val == 0.1
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        ds.set_par('const.c0', 2)

        assert const.c0.val == 2

        ds.set_par('const.integrate', False)
        ds.freeze('const.c0')

        vals = ds.get_par('poly.c1.val')
        assert ([0.1, 0.45, 0.45] == vals).all()

        pars = ds.get_par('const.c0')

        ds.fit()

        assert round(poly1.c1.val) == 1
        assert round(poly1.c2.val) == 3
        assert round(poly2.c1.val) == 3
        assert round(poly2.c2.val) == 1
        assert round(poly3.c1.val) == 3
        assert round(poly3.c2.val) == 1

        ds.clear_stack()

        x1 = np.arange(50)+100
        y1 = 7*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        ds.load_arrays([[x1, y1], [x2, y2], [x3, y3]])

        datastack.set_template_id("foo")

        ds.set_source('const1d.constfoo * polynom1d.polyfoo')

        const1 = ui._session._get_model_component('const1')
        const2 = ui._session._get_model_component('const2')
        const3 = ui._session._get_model_component('const3')

        ds[2,3].link('const.c0')

        ds[2].set_par('const.c0', 3)
        ds[1].set_par('const.c0', 7)

        ds[1].freeze('const.c0')

        assert const2.c0.frozen is False

        ds.fit()

        assert const2.c0.val == const3.c0.val
        assert const3.c0._link is const2.c0

        ds.unlink("const.c0")

        assert const3.c0._link is not const2.c0

class test_pha(SherpaTestCase):
    def setUp(self):
        clear_stack()
        ui.clean()
        set_template_id("ID")
        self._this_dir = os.path.dirname(sys.modules[self.__module__].__file__)

    def tearDown(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits')
    def test_case_6(self):
        datadir = '/'.join((self._this_dir, 'data'))
        ls = '@'+'/'.join((datadir, 'pha.lis'))
        load_pha(ls)

        load_bkg_rmf([], '/'.join((datadir, "acisf04938_000N002_r0043_rmf3.fits")))
        load_bkg_rmf([], '/'.join((datadir, "acisf07867_000N001_r0002_rmf3.fits")))

        load_bkg_arf([], '/'.join((datadir, "acisf04938_000N002_r0043_arf3.fits")))
        load_bkg_arf([], '/'.join((datadir, "acisf07867_000N001_r0002_arf3.fits")))


        # Define background models
        bkg_arfs = get_bkg_arf([])
        bkg_scales = get_bkg_scale([])
        bkg_models = [const1d.c1 * acis_bkg_model('acis7s'),
                      const1d.c2 * acis_bkg_model('acis7s')]
        bkg_rsps = get_response([], bkg_id=1)
        for i in range(2):
            id_ = i + 1
            # Make the ARF spectral response flat.  This is required for using
            # the acis_bkg_model.
            bkg_arfs[i].specresp = bkg_arfs[i].specresp * 0 + 1.
            set_bkg_full_model(id_, bkg_rsps[i](bkg_models[i]))

        # Fit background
        notice(0.5, 8.)
        set_method("neldermead")
        set_stat("cash")

        thaw(c1.c0)
        thaw(c2.c0)
        fit_bkg()
        freeze(c1.c0)
        freeze(c2.c0)

        # Define source models
        rsps = get_response([])
        src_model = powlaw1d.pow1
        src_models = [src_model,
                      src_model * const1d.ratio_12]
        for i in range(2):
            id_ = i + 1
            set_full_model(id_, (rsps[i](src_models[i])
                                 + bkg_scales[i] * bkg_rsps[i](bkg_models[i])))

        fit()

class test_query(SherpaTestCase):
    def setUp(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")
        self._this_dir = os.path.dirname(sys.modules[self.__module__].__file__)

    def tearDown(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")
        set_stack_verbose(False)

    @unittest.skipIf(not has_fits_support(),
                     'need pycrates, pyfits')
    def test_case_7(self):
        load_pha('@'+'/'.join((self._this_dir, 'data', 'pha.lis')))

        f = query_by_header_keyword('INSTRUME', 'ACIS')

        assert f == [1,2]

        f = query_by_obsid('7867')

        assert f == [2]

        ds = DataStack()

        ds.load_pha('@'+'/'.join((self._this_dir, 'data', 'pha.lis')))

        f = ds.query_by_obsid('4938')

        assert f == [3]

        f = query_by_obsid(ds, '7867')

        assert f == [4]
