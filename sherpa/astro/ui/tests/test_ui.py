# 
#  Copyright (C) 2012  Smithsonian Astrophysical Observatory
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


from sherpa.utils import SherpaTest, SherpaTestCase, needs_data
import sherpa.astro.ui as ui
import numpy
import logginglogger = logging.getLogger("sherpa")

class test_ui(SherpaTestCase):

    @needs_data
    def setUp(self):
        self.ascii = self.datadir + '/threads/ascii_table/sim.poisson.1.dat'
        self.fits = self.datadir + '/1838_rprofile_rmid.fits'
        self.singledat = self.datadir + '/single.dat'
        self.singletbl = self.datadir + '/single.fits'
        self.doubledat = self.datadir + '/double.dat'
        self.doubletbl = self.datadir + '/double.fits'
        self.img = self.datadir + '/img.fits'
        self.filter_single_int_ascii = self.datadir + '/filter_single_integer.dat'
        self.filter_single_int_table = self.datadir + '/filter_single_integer.fits'
        self.filter_single_log_table = self.datadir + '/filter_single_logical.fits'

        self.func = lambda x: x
        ui.dataspace1d(1,1000,dstype=ui.Data1D)


    @needs_data
    def test_ascii(self):
        ui.load_ascii(1, self.ascii)
        ui.load_ascii(1, self.ascii, 2)
        ui.load_ascii(1, self.ascii, 2, ("col2", "col1"))


    @needs_data
    def test_table(self):
        ui.load_table(1, self.fits)
        ui.load_table(1, self.fits, 3)
        ui.load_table(1, self.fits, 3, ["RMID","SUR_BRI","SUR_BRI_ERR"])
        ui.load_table(1, self.fits, 4, ('R',"SUR_BRI",'SUR_BRI_ERR'),
                      ui.Data1DInt)


    # Test table model
    @needs_data
    def test_table_model_ascii_table(self):
        ui.load_table_model('tbl', self.singledat)
        ui.load_table_model('tbl', self.doubledat)


    @needs_data
    def test_table_model_fits_table(self):
        ui.load_table_model('tbl', self.singletbl)
        ui.load_table_model('tbl', self.doubletbl)


    @needs_data
    def test_table_model_fits_image(self):
        ui.load_table_model('tbl', self.img)


    # Test user model
    @needs_data
    def test_user_model_ascii_table(self):
        ui.load_user_model(self.func, 'mdl', self.singledat)
        ui.load_user_model(self.func, 'mdl', self.doubledat)


    @needs_data
    def test_user_model_fits_table(self):
        ui.load_user_model(self.func, 'mdl', self.singletbl)
        ui.load_user_model(self.func, 'mdl', self.doubletbl)


    @needs_data
    def test_filter_ascii(self):
        ui.load_filter(self.filter_single_int_ascii)
        ui.load_filter(self.filter_single_int_ascii, ignore=True)


    # Test load_filter
    @needs_data
    def test_filter_table(self):
        ui.load_filter(self.filter_single_int_table)
        ui.load_filter(self.filter_single_int_table, ignore=True)

        ui.load_filter(self.filter_single_log_table)
        ui.load_filter(self.filter_single_log_table, ignore=True)

class test_more_ui(SherpaTestCase):
    @needs_data
    def setUp(self):
        self.img = self.datadir + '/img.fits'
	self.pha = self.datadir + '/threads/simultaneous/pi2286.fits'
	self.rmf = self.datadir + '/threads/simultaneous/rmf2286.fits'
	logger.setLevel(logging.ERROR)

    #bug #12732
    @needs_data
    def test_string_model_with_rmf(self):
	ui.load_pha("foo", self.pha)
	ui.load_rmf("foo", self.rmf)
	# Check that get_rmf(id)('modelexpression') works
	caught = False
	try:
		m=ui.get_rmf("foo")("powlaw1d.pl1")
	except:
		caught = True
	if caught:
		self.fail("Exception caught when it shouldn't")
	from sherpa.astro.instrument import RMFModelPHA
	self.assertTrue(isinstance(m, RMFModelPHA))
    


class test_image_12578(SherpaTestCase):
    @needs_data
    def setUp(self):
        self.img = self.datadir + '/img.fits'
	logger.setLevel(logging.ERROR)
	ui.clean()
    
    #bug #12578
    @needs_data
    def test_set_coord_bad_coord(self):
	from sherpa.utils.err import IdentifierErr, DataErr

	# Test Case #1: if the list of ids is empty, raise a IdentifierErr['nodatasets']
	caught = False
	try:
		ui.set_coord('image')
	except IdentifierErr:
		caught = True
	if not caught:
		self.fail("Test Case #1: IdentifierErr Exception not caught")

 	# Test Case #2: check the user expected behavior. The call set_coord("sky")
	# will result in the error message DataErr: unknown coordinates: 'sky'\n \
	# Valid coordinates: logical, image, physical, world, wcs
	ui.load_image(self.img)

	caught = False
	try:
		ui.set_coord("sky")
	except DataErr as e:
		okmsg = "unknown coordinates: 'sky'\nValid options: logical, image, physical, world, wcs"
		self.assertEqual(okmsg, e.message)
		caught = True
	if not caught:
		self.fail("Test Case #2: DataErr Exception not caught")

class test_psf_ui(SherpaTestCase):

    models1d = ['beta1d', 'lorentz1d', 'normbeta1d']
    models2d = ['beta2d', 'devaucouleurs2d', 'hubblereynolds', 'lorentz2d']

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_psf_model2d(self):
        ui.dataspace1d(1, 10)
        for model in self.models1d:
            try:
                ui.load_psf('psf1d', model+'.mdl')
                ui.set_psf('psf1d')
                mdl = ui.get_model_component('mdl')
                self.assert_( (numpy.array(mdl.get_center()) ==
                               numpy.array([4])).all() )
            except:
                print model
                raise


    def test_psf_model2d(self):
        ui.dataspace2d([216,261])
        for model in self.models2d:
            try:
                ui.load_psf('psf2d', model+'.mdl')
                ui.set_psf('psf2d')
                mdl = ui.get_model_component('mdl')
                self.assert_( (numpy.array(mdl.get_center()) ==
                               numpy.array([108,130])).all() )
            except:
                print model
                raise
    
    #bug #12503
    def test_psf_pars_are_frozen(self):
	ui.load_psf('psf', ui.beta2d.p1)
	self.assertEqual([], p1.thawedpars)
    
class test_stats_ui(SherpaTestCase):

    @needs_data
    def setUp(self):
        self.data = self.datadir + '/threads/chi2/SWIFTJ0840.1+2946.pha.gz'
	ui.clean()

    #bugs #11400, #13297, #12365    
    @needs_data
    def test_chi2(self):

	#Case 1: first ds has no error, second has, chi2-derived (chi2gehrels) statistic
	#I expect stat.name to be chi2gehrels for ds1, chi2 for ds2, chi2gehrels for ds1,2
	ui.load_data(1, self.data)
	ui.load_data(2, self.data, use_errors=True)

	ui.set_source(1, "gauss1d.g1")
	ui.set_source(2, "gauss1d.g1")

	ui.set_stat("chi2gehrels")

	si=ui.get_stat_info()

	stat1 = si[0].statname
	stat2 = si[1].statname
	stat12 = si[2].statname

	self.assertEqual('chi2gehrels', stat1)
	self.assertEqual('chi2', stat2)
	self.assertEqual('chi2gehrels', stat12)

	#Case 2: first ds has errors, second has not, chi2-derived (chi2gehrels) statistic
	#I expect stat.name to be chi2 for ds1, chi2gehrels for ds2, chi2gehrels for ds1,2
	ui.load_data(2, self.data)
	ui.load_data(1, self.data, use_errors=True)
	
	si=ui.get_stat_info()

	stat1 = si[0].statname
	stat2 = si[1].statname
	stat12 = si[2].statname	

	self.assertEqual('chi2gehrels', stat2)
	self.assertEqual('chi2', stat1)
	self.assertEqual('chi2gehrels', stat12)

	#Case 3: both datasets have errors, chi2-derived (chi2gehrels) statistic
	#I expect stat.name to be chi2 for ds1, chi2 for ds2, chi2 for ds1,2
	ui.load_data(2, self.data, use_errors=True)
	ui.load_data(1, self.data, use_errors=True)
	
	si=ui.get_stat_info()

	stat1 = si[0].statname
	stat2 = si[1].statname
	stat12 = si[2].statname	

	self.assertEqual('chi2', stat2)
	self.assertEqual('chi2', stat1)
	self.assertEqual('chi2', stat12)

	#Case 4: first ds has errors, second has not, LeastSq statistic
	#I expect stat.name to be leastsq for ds1, leastsq for ds2, leastsq for ds1,2
	ui.load_data(2, self.data)
	ui.load_data(1, self.data, use_errors=True)
	
	ui.set_stat("leastsq")

	si=ui.get_stat_info()

	stat1 = si[0].statname
	stat2 = si[1].statname
	stat12 = si[2].statname	

	self.assertEqual('leastsq', stat2)
	self.assertEqual('leastsq', stat1)
	self.assertEqual('leastsq', stat12)

	#Case 5: both ds have errors, LeastSq statistic
	#I expect stat.name to be leastsq for ds1, leastsq for ds2, leastsq for ds1,2
	ui.load_data(2, self.data, use_errors=True)
	ui.load_data(1, self.data, use_errors=True)
	
	ui.set_stat("leastsq")

	si=ui.get_stat_info()

	stat1 = si[0].statname
	stat2 = si[1].statname
	stat12 = si[2].statname	

	self.assertEqual('leastsq', stat2)
	self.assertEqual('leastsq', stat1)
	self.assertEqual('leastsq', stat12)

	#Case 6: first ds has errors, second has not, CStat statistic
	#I expect stat.name to be cstat for ds1, cstat for ds2, cstat for ds1,2
	ui.load_data(2, self.data)
	ui.load_data(1, self.data, use_errors=True)
	
	ui.set_stat("cstat")

	si=ui.get_stat_info()

	stat1 = si[0].statname
	stat2 = si[1].statname
	stat12 = si[2].statname	

	self.assertEqual('cstat', stat2)
	self.assertEqual('cstat', stat1)
	self.assertEqual('cstat', stat12)

	#Case7: select chi2 as statistic. One of the ds does not provide errors
	#I expect sherpa to raise a StatErr exception.
	ui.set_stat('chi2')

	caught=False

	from sherpa.utils.err import StatErr
	try:
		ui.get_stat_info()
	except StatErr:
		caught=True
	
	self.assertTrue(caught)

	#Case8: select chi2 as statistic. Both datasets provide errors
	#I expect stat to be 'chi2'
	ui.load_data(2, self.data, use_errors=True)
	si=ui.get_stat_info()

	stat1 = si[0].statname
	stat2 = si[1].statname
	stat12 = si[2].statname

	self.assertEqual('chi2', stat2)
	self.assertEqual('chi2', stat1)
	self.assertEqual('chi2', stat12)
	
	

if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        SherpaTest(ui).test(datadir=sys.argv[1])
