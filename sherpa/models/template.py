#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2011)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from parameter import Parameter, tinyval
from model import ArithmeticModel, modelCacher1d, CompositeModel, \
    ArithmeticFunctionModel
from basic import TableModel
import numpy, operator
from sherpa.utils.err import ModelErr

__all__ = ('create_template_model', 'TemplateModel', 'KNNInterpolator', 'Template')



def create_template_model(modelname, names, parvals, templates, template_interpolator_name='default'):
    """
    Create a TemplateModel model class from template input


    `modelname`  - name of the template model.

    `names`      - list of strings that define the order of the 
                   named parameters.

    `parvals`    - 2-D ndarray of parameter vectors, index corresponds
                   to the spectrum in `templates`. The parameter grid.

    `templates`  - list of TableModel objects that contain a spectrum
                   at a specific parameter vector (corresponds to a row
                   in `parvals`).

    `template_interpolator_name` - name of the template interpolator, or None
		   for disabling interpolation *between* templates.
                   See load_template_model for more information.

    """
    # Create a list of parameters from input
    pars = []   
    for ii, name in enumerate(names):
        minimum = min(parvals[:,ii])
        maximum = max(parvals[:,ii])
        initial = parvals[:,ii][0]
        # Initial parameter value is always first parameter value listed
        par = Parameter(modelname, name, initial,
                        minimum, maximum,
                        minimum, maximum)
        pars.append(par)

    # Create the templates table from input
    tm = TemplateModel(modelname, pars, parvals, templates)
    if template_interpolator_name is not None:
	if interpolators.has_key(template_interpolator_name):
		interp = interpolators[template_interpolator_name]
		args = interp[1]
		args['template_model'] = tm
		args['name'] = modelname
		return interp[0](**args)
    else:
	return tm


class InterpolatingTemplateModel(ArithmeticModel):
    def __init__(self, name, template_model):
        self.template_model = template_model
        for par in template_model.pars:
            self.__dict__[par.name] = par
	    self.parvals = template_model.parvals
	ArithmeticModel.__init__(self, name, template_model.pars)

    def fold(self, data):
        for template in self.template_model.templates:
            template.fold(data)

    @modelCacher1d
    def calc(self, p, x0, x1=None, *args, **kwargs):
        interpolated_template = self.interpolate(p, x0)
        return interpolated_template(x0, x1, *args, **kwargs)

class KNNInterpolator(InterpolatingTemplateModel):
    def __init__(self, name, template_model, k=None, order=2):
        self._distances = {}
        if k is None:
            self.k = 2*template_model.parvals[0].size
        else:
            self.k = k
        self.order = order
	InterpolatingTemplateModel.__init__(self, name, template_model)

    def _calc_distances(self, point):
        self._distances = {}
        for i, t_point in enumerate(self.template_model.parvals):
            self._distances[i] = numpy.linalg.norm(point - t_point, self.order)
        self._distances = sorted(self._distances.iteritems(), key=operator.itemgetter(1))
    
    def interpolate(self, point, x_out):
        self._calc_distances(point)
        if self._distances[0][1]==0:
            return self.template_model.templates[self._distances[0][0]]
        k_distances = self._distances[:self.k]
        weights = [(idx, 1/numpy.array(distance)) for idx, distance in k_distances]
        sum_weights = sum([1/weight for idx, weight in k_distances])
        y_out = numpy.zeros(len(x_out))
        for idx, weight in weights:
            y_out += self.template_model.templates[idx].calc((weight,), x_out)
        y_out /= sum_weights
        tm = TableModel('interpolated')
        tm.load(x_out, y_out)
        return tm

class Template(KNNInterpolator):
	def __init__(self, *args, **kwargs):
		KNNInterpolator.__init__(self, *args, **kwargs)


class TemplateModel(ArithmeticModel):

    def __init__(self, name='templatemodel', pars=(), parvals=[], templates=[]):
        self.parvals = parvals
        self.templates = templates
	self.index = {}

	for par in pars:
            self.__dict__[par.name] = par

        for ii, parval in enumerate(parvals):
            self.index[tuple(parval)] = templates[ii]        

        ArithmeticModel.__init__(self, name, pars)
	self.is_discrete = True

    def fold(self, data):
        for template in self.templates:
            template.fold(data)


    def get_x(self):
	p = tuple(par.val for par in self.pars)
	template = self.query(p)
        return template.get_x()

    def get_y(self):
	p = tuple(par.val for par in self.pars)
	template = self.query(p)
        return template.get_y()


    def query(self, p):
	try:
        	return self.index[tuple(p)]
	except:
		raise ModelErr("Interpolation of template parameters was disabled for this model, but parameter values not in the template library have been requested. Please use gridsearch method and make sure the sequence option is consistent with the template library")


    @modelCacher1d
    def calc(self, p, x0, x1=None, *args, **kwargs):
        table_model = self.query(p)

        # return interpolated the spectrum according to the input grid (x0, [x1])
        return table_model(x0, x1, *args, **kwargs)



interpolators = {
	'default' : (Template, {'k':2, 'order':2})
}
