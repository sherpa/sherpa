#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import sherpa.all
import sherpa.ui.utils
from sherpa.utils import calc_mlr, calc_ftest, rebin, histogram1d, \
    histogram2d, gamma, lgam, erf, igamc, igam, incbet, multinormal_pdf, \
    multit_pdf
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.logposterior import Prior

# We build up __all__ as we go along
__all__ = ['calc_mlr', 'calc_ftest','Data1D', 'Data1DInt',
           'Data2D', 'Data2DInt', 'rebin', 'histogram1d',
           'histogram2d','gamma', 'lgam', 'erf', 'igamc',
           'igam', 'incbet', 'Prior', 'multinormal_pdf', 'multit_pdf']

_session = utils.Session()
_session._add_model_types(sherpa.models.basic)
# To get PSFModel in list of models -- doesn't inherit from ArithmeticModel
_session._add_model_types(sherpa.instrument, baselist=(sherpa.models.Model,))
__all__.extend(_session._export_names(globals()))


__all__ = tuple(__all__)
