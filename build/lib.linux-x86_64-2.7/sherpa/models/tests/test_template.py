#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2011)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.models import TableModel, Gauss1D
from sherpa.models.template import TemplateModel, create_template_model
from sherpa.utils import SherpaTest, SherpaTestCase, needs_data
from sherpa.utils.err import ModelErr
from sherpa import ui
import numpy, logging, os

logger = logging.getLogger("sherpa")

class test_new_templates_ui(SherpaTestCase):
    def assign_model(self, name, obj):
        self.locals[name] = obj

    def run_thread(self, name, scriptname='fit.py'):
        ui.clean()
        ui.set_model_autoassign_func(self.assign_model)
        self.locals = {}
        os.chdir(os.path.join(self.datadir, 'ciao4.3', name))
        execfile(scriptname, {}, self.locals)

    def setUp(self):
        logger.setLevel(logging.ERROR)


    # Regression reported by Aneta (original report by Fabrizio Nicastro)
    # When restoring a file saved with an older version of sherpa,
    # we need to make sure models are assigned an is_discrete field.
    # For model that do not have the is_discrete field, fallback to False.
    @needs_data
    def test_restore_is_discrete(self):
	self.run_thread('template_restore', 'test.py')
	

    # TestCase 1 load_template_model enables interpolation by default
    @needs_data
    def test_load_template_with_interpolation(self):
        self.run_thread('load_template_with_interpolation')
        try:   
                self.assertEqualWithinTol(2023.46, ui.get_fit_results().parvals[0], 0.001)
                self.assertEqualWithinTol(2743.47, ui.get_fit_results().parvals[1], 0.001)
        except:
                self.assertEqualWithinTol(2743.47, ui.get_fit_results().parvals[0], 0.001)
                self.assertEqualWithinTol(2023.46, ui.get_fit_results().parvals[1], 0.001)

    @needs_data
    def test_load_template_interpolator(self):
        self.run_thread('load_template_interpolator')
        self.assertEqualWithinTol(2743.91, ui.get_fit_results().parvals[0], 0.001)

    # TestCase 2 load_template_model with template_interpolator_name=None disables interpolation
    # TestCase 3.1 discrete templates fail when probed for values they do not represent (gridsearch with wrong values) 
    @needs_data
    def test_load_template_model_without_interpolation(self):
        try:   
                self.run_thread('load_template_without_interpolation', scriptname='test_case_2_and_3.1.py')
        except ModelErr:
                return
        self.fail('Fit should have failed: using gridsearch with wrong parvals')

    # TestCase 3.2 discrete templates fail when probed for values they do not represent (continuous method with discrete template)
    @needs_data
    def test_load_template_model_without_interpolation(self):
        try:   
                self.run_thread('load_template_without_interpolation', scriptname='test_case_3.2.py')
        except ModelErr:
                return
        self.fail('Fit should have failed: using gridsearch with wrong parvals')


    # TestCase 4 gridsearch with right values succeeds
    @needs_data
    def test_grid_search_with_discrete_template(self):
	self.run_thread('load_template_without_interpolation', scriptname='test_case_4.py')

    # TestCase 5 user can access interpolators' parvals
    @needs_data
    def test_grid_search_with_discrete_template(self):
	self.run_thread('load_template_with_interpolation', scriptname='test_case_5.py')



class test_template(SherpaTestCase):

    def setUp(self):
        self.num = 4
        self.ncoords = 100
        self.ntemplates = 2**self.num
        self.x = numpy.linspace(0.1, 5, 50)
        g1 = Gauss1D('g1')

        # create a 4-dimensional grid from 0 to 1 inclusive, shape = (16,4)
        grid = numpy.mgrid[ [slice(0,2,1) for ii in range(self.num)] ]
        grid = numpy.asarray(map(numpy.ravel, grid)).T
        coords = numpy.linspace(0.01, 6, 100)
        names = ["p%i" % i for i in range(self.num)]
        templates = []
        for ii in range(self.ntemplates):
            t = TableModel()
            g1.fwhm = numpy.random.uniform(0.5, 2.0)
            g1.pos  = numpy.random.uniform(1.0, 4.5)
            g1.ampl = numpy.random.uniform(1.0, 50.)
            t.load(coords, g1(coords))
            templates.append(t)

        self.model = create_template_model("mdl", names, grid, templates)

    def tearDown(self):
        self.model = None


    def test_template_model_evaluation(self):
        self.model.thawedpars = [0,1,0,1]
        vals = self.model(self.x)

#    def test_template_query_index(self):
#        expected = 5
#        result = self.model.query_index([0,1,0,1])
#        self.assertEqual(expected, result,
#                         "Expected %s instead of %s" % (str(expected), str(result)))

#    def test_template_query(self):
#        result = self.model.query([0,1,0,1])



if __name__ == '__main__':

    import sherpa.models as models
    SherpaTest(models).test()
