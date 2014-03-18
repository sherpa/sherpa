#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import os.path
import sherpa
from sherpa.utils import SherpaTest, SherpaTestCase, needs_data
from sherpa import ui

class test_sherpa(SherpaTestCase):

    def test_include_dir(self):
        incdir = os.path.join(sherpa.get_include(), 'sherpa')
        self.assert_(os.path.isdir(incdir))

    @needs_data
    def setUp(self):
	self.agn2 = self.datadir + '/ciao4.3/faulty_load_data/agn2'
	self.agn2_fixed = self.datadir + '/ciao4.3/faulty_load_data/agn2_fixed'
	self.template_idx = self.datadir + '/ciao4.3/faulty_load_data/table.txt'

    @needs_data
    def test_not_reading_header_without_comment(self):
	self.assertRaises(ValueError, ui.load_data, self.agn2)

    @needs_data
    def test_reading_floats(self):
	ui.load_data(self.agn2_fixed)

    @needs_data
    def test_reading_strings(self):
	ui.load_data(self.template_idx, require_floats=False)

    @needs_data
    def test_require_float(self):
	self.assertRaises(ValueError, ui.load_data, self.agn2)

