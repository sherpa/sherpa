#
#  Copyright (C) 2007, 2018  Smithsonian Astrophysical Observatory
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

import sherpa.all
import sherpa.ui.utils
from sherpa.utils import calc_mlr, calc_ftest, rebin, histogram1d, \
    histogram2d, gamma, lgam, igamc, igam, incbet, multinormal_pdf, \
    multit_pdf
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.logposterior import Prior

# We build up __all__ as we go along
__all__ = ['calc_mlr', 'calc_ftest', 'Data1D', 'Data1DInt',
           'Data2D', 'Data2DInt', 'rebin', 'histogram1d',
           'histogram2d', 'gamma', 'lgam', 'igamc',
           'igam', 'incbet', 'Prior', 'multinormal_pdf', 'multit_pdf']

_session = utils.Session()
_session._add_model_types(sherpa.models.basic)
_session._add_model_types(sherpa.models.template)
# To get PSFModel in list of models -- doesn't inherit from ArithmeticModel
_session._add_model_types(sherpa.instrument, baselist=(sherpa.models.Model,))
__all__.extend(_session._export_names(globals()))


__all__ = tuple(__all__)
