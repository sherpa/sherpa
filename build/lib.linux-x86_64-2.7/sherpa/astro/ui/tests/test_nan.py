from sherpa.utils import SherpaTest, SherpaTestCase, needs_data
import sherpa.astro.ui as ui
import logging
import os
import numpy
logger = logging.getLogger("sherpa")

class test_more_ui(SherpaTestCase):
    def assign_model(self, name, obj):
        self.locals[name] = obj

    def run_thread(self, name, scriptname='fit.py'):
        ui.clean()
        ui.set_model_autoassign_func(self.assign_model)
        self.locals = {}
        os.chdir(os.path.join(self.datadir, 'ciao4.3', name))
        execfile(scriptname, {}, self.locals)

    @needs_data
    def setUp(self):
        self.img = self.datadir + '/img.fits'
        self.pha = self.datadir + '/threads/simultaneous/pi2286.fits'
        self.rmf = self.datadir + '/threads/simultaneous/rmf2286.fits'
        self.nan = self.datadir + '/ciao4.3/filternan/with_nan.fits'
        logger.setLevel(logging.ERROR)

    #bug 12784
    @needs_data
    def test_filter_nan(self):
        self.run_thread('filternan')
        self.assertFalse(numpy.isnan(ui.get_fit_results().statval))

if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        SherpaTest(ui).test(datadir=sys.argv[1])