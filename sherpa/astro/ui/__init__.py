#
#  Copyright (C) 2007, 2018, 2020, 2021, 2024
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

import sherpa.all
import sherpa.astro.all
import sherpa.astro.ui.utils
from sherpa.utils import calc_mlr, calc_ftest, rebin, histogram1d, \
    histogram2d, gamma, lgam, igamc, igam, incbet, multinormal_pdf, \
    multit_pdf
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.astro.data import DataARF, DataRMF, DataPHA, DataIMG, DataIMGInt
from sherpa.logposterior import Prior


# We build up __all__ as we go along
__all__ = ['DataARF', 'DataRMF', 'DataPHA', 'DataIMG', 'DataIMGInt', 'Data1D',
           'Data1DInt', 'Data2D', 'Data2DInt', 'calc_mlr', 'calc_ftest', 'rebin',
           'histogram1d', 'histogram2d', 'gamma', 'lgam', 'igamc',
           'igam', 'incbet', 'Prior', 'multinormal_pdf', 'multit_pdf']

_session = utils.Session()
_session._add_model_types(sherpa.models.basic)
_session._add_model_types(sherpa.astro.models)
_session._add_model_types(sherpa.astro.optical)
_session._add_model_types(sherpa.models.template)

# To add PSFModel to list -- doesn't inherit from ArithmeticModel
_session._add_model_types(sherpa.instrument, baselist=(sherpa.models.Model,))

# Get RMFModel, ARFModel in list of models
_session._add_model_types(sherpa.astro.instrument)

if hasattr(sherpa.astro, 'xspec'):
    _session._add_model_types(sherpa.astro.xspec,
                              (sherpa.astro.xspec.XSAdditiveModel,
                               sherpa.astro.xspec.XSMultiplicativeModel,
                               sherpa.astro.xspec.XSConvolutionKernel))

    from sherpa.astro.xspec import get_xsabund, get_xscosmo, get_xsxsect, \
        set_xsabund, set_xscosmo, set_xsxsect, set_xsxset, get_xsxset, \
        get_xschatter, set_xschatter, get_xsabundances, set_xsabundances
    __all__.extend(('get_xsabund', 'get_xschatter', 'get_xscosmo',
                    'get_xsxsect', 'set_xsabund', 'set_xschatter',
                    'set_xscosmo', 'set_xsxsect', 'set_xsxset',
                    'get_xsxset', 'get_xsabundances', 'set_xsabundances'))

__all__.extend(_session._export_names(globals()))

__all__.append('_session')


__all__ = tuple(__all__)
