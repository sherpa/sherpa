#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.utils import SherpaTestCase
import os
import sys
from acis_bkg_model import acis_bkg_model
from sherpa.astro.ui import *
from sherpa.astro.datastack import *
logger = logging.getLogger('sherpa')
logger.setLevel(logging.ERROR)


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



    def test_case_1(self):
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



