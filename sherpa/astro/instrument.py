# 
#  Copyright (C) 2010  Smithsonian Astrophysical Observatory
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

import numpy
import sherpa
from sherpa.utils.err import InstrumentErr, DataErr, PSFErr, ArgumentTypeErr
from sherpa.models.model import ArithmeticFunctionModel, NestedModel, \
    ArithmeticModel, CompositeModel, Model
WCS = None
try:
    from sherpa.astro.io.wcs import WCS
except:
    WCS = None

from sherpa.instrument import PSFModel as _PSFModel
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.data import BaseData, Data1D
from sherpa.astro.data import DataARF, DataRMF, DataPHA, _notice_resp, \
    DataIMG
from sherpa.utils import sao_fcmp, sum_intervals, sao_arange
from sherpa.astro.utils import compile_energy_grid
from itertools import izip

_tol = numpy.finfo(numpy.float32).eps

__all__ = ('RMFModel', 'ARFModel', 'RSPModel',
           'RMFModelPHA', 'RMFModelNoPHA',
           'ARFModelPHA', 'ARFModelNoPHA',
           'RSPModelPHA', 'RSPModelNoPHA',
           'MultiResponseSumModel', 'PileupRMFModel', 'RMF1D', 'ARF1D',
           'Response1D', 'MultipleResponse1D','PileupResponse1D',
           'PSFModel')


class RMFModel(CompositeModel, ArithmeticModel):
    """
    Base class for expressing RMF convolution in model expressions
    """
    def __init__(self, rmf, model):
        self.rmf = rmf
        self.model = model

        # Logic for ArithmeticModel.__init__
        self.pars = ()

        # FIXME: group pairs of coordinates with one attribute

        self.elo = None; self.ehi = None  # Energy space
        self.lo = None;  self.hi = None   # Wavelength space
        self.xlo = None; self.xhi = None  # Current Spectral coordinates

        # Used to rebin against finer or coarser energy grids
        self.rmfargs = ()

        CompositeModel.__init__(self, 'apply_rmf(%s)' % model.name, (model,))
        self.filter()


    def filter(self):
        # Energy grid (keV)
        self.elo, self.ehi = self.rmf.get_indep()

        # Wavelength grid (angstroms)
        self.lo, self.hi = DataPHA._hc/self.ehi, DataPHA._hc/self.elo

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi

        # Used to rebin against finer or coarser energy grids
        self.rmfargs = ()


    def startup(self):
        self.model.startup()
        CompositeModel.startup(self)


    def teardown(self):
        self.model.teardown()
        CompositeModel.teardown(self)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        raise NotImplementedError


class ARFModel(CompositeModel, ArithmeticModel):
    """
    Base class for expressing ARF convolution in model expressions
    """
    def __init__(self, arf, model):
        self.arf = arf
        self.model = model

        self.elo = None; self.ehi = None  # Energy space
        self.lo = None;  self.hi = None   # Wavelength space
        self.xlo = None; self.xhi = None  # Current Spectral coordinates

        # Used to rebin against finer or coarser energy grids
        self.arfargs = ()

        # Logic for ArithmeticModel.__init__
        self.pars = ()

        CompositeModel.__init__(self, 'apply_arf(%s)' % model.name, (model,))
        self.filter()


    def filter(self):
        # Energy grid (keV)
        self.elo, self.ehi = self.arf.get_indep()

        # Wavelength grid (angstroms)
        self.lo, self.hi = DataPHA._hc/self.ehi, DataPHA._hc/self.elo

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi

        # Used to rebin against finer or coarser energy grids
        self.arfargs = ()


    def startup(self):
        self.model.startup()
        CompositeModel.startup(self)


    def teardown(self):
        self.model.teardown()
        CompositeModel.teardown(self)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        raise NotImplementedError


class RSPModel(CompositeModel, ArithmeticModel):
    """
    Base class for expressing RMF + ARF convolution in model expressions
    """
    def __init__(self, arf, rmf, model):
        self.arf = arf
        self.rmf = rmf
        self.model = model

        self.elo = None; self.ehi = None  # Energy space
        self.lo = None; self.hi = None    # Wavelength space
        self.xlo = None; self.xhi = None  # Current Spectral coordinates

        # Used to rebin against finer or coarser energy grids
        self.rmfargs = ()
        self.arfargs = ()

        # Logic for ArithmeticModel.__init__
        self.pars = ()

        CompositeModel.__init__(self, 'apply_rmf(apply_arf(%s))' % model.name,
                                (model,))
        self.filter()


    def filter(self):
        # Energy grid (keV), ARF grid breaks tie
        self.elo, self.ehi = self.arf.get_indep()

        # Wavelength grid (angstroms)
        self.lo, self.hi = DataPHA._hc/self.ehi, DataPHA._hc/self.elo

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi

        # Used to rebin against finer or coarser energy grids
        self.rmfargs = ()
        self.arfargs = ()


    def startup(self):
        self.model.startup()
        CompositeModel.startup(self)


    def teardown(self):
        self.model.teardown()
        CompositeModel.teardown(self)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        raise NotImplementedError



class RMFModelPHA(RMFModel):
    """
    RMF convolution model with associated PHA
    """
    def __init__(self, rmf, pha, model):
        self.pha = pha
        self._rmf = rmf # store a reference to original
        RMFModel.__init__(self, rmf, model)


    def filter(self):

        RMFModel.filter(self)

        pha = self.pha
        # If PHA is a finer grid than RMF, evaluate model on PHA and
        # rebin down to the granularity that the RMF expects.
        if pha.bin_lo is not None and pha.bin_hi is not None:
            bin_lo, bin_hi = pha.bin_lo, pha.bin_hi

            # If PHA grid is in angstroms then convert to keV for
            # consistency
            if (bin_lo[0] > bin_lo[-1]) and (bin_hi[0] > bin_hi[-1]):
                bin_lo, bin_hi = DataPHA._hc/pha.bin_hi, DataPHA._hc/pha.bin_lo

            # FIXME: What about filtered option?? bin_lo, bin_hi are unfiltered??

            # Compare disparate grids in energy space
            self.rmfargs = ((self.elo, self.ehi), (bin_lo, bin_hi))

            # FIXME: Compute on finer energy grid?  Assumes that PHA has
            # finer grid than RMF
            self.elo, self.ehi = bin_lo, bin_hi

            # Wavelength grid (angstroms)
            self.lo, self.hi = DataPHA._hc/self.ehi, DataPHA._hc/self.elo

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi
        if self.pha.units == 'wavelength':
            self.xlo, self.xhi = self.lo, self.hi


    def startup(self):
        rmf = self._rmf # original

        # Create a view of original RMF
        self.rmf = DataRMF(rmf.name, rmf.detchans, rmf.energ_lo, rmf.energ_hi,
                            rmf.n_grp, rmf.f_chan, rmf.n_chan, rmf.matrix,
                            rmf.offset, rmf.e_min, rmf.e_max, rmf.header)

        # Filter the view for current fitting session
        _notice_resp(self.pha.get_noticed_channels(), None, self.rmf)

        self.filter()

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi
        if self.pha.units == 'wavelength':
            self.xlo, self.xhi = self.lo, self.hi

        RMFModel.startup(self)


    def teardown(self):
        self.rmf = self._rmf

        self.filter()
        RMFModel.teardown(self)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        # x is noticed/full channels here

        src = self.model.calc(p, self.xlo, self.xhi)
        return self.rmf.apply_rmf(src, *self.rmfargs)


class RMFModelNoPHA(RMFModel):
    """
    RMF convolution model without associated PHA
    """
    def __init__(self, rmf, model):
        RMFModel.__init__(self, rmf, model)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        # x is noticed/full channels here

        # Always evaluates source model in keV!
        src = self.model.calc(p, self.xlo, self.xhi)
        return self.rmf.apply_rmf(src)


class ARFModelPHA(ARFModel):
    """
    ARF convolution model with associated PHA
    """
    def __init__(self, arf, pha, model):
        self.pha = pha
        self._arf = arf # store a reference to original
        ARFModel.__init__(self, arf, model)


    def filter(self):

        ARFModel.filter(self)

        pha = self.pha
        # If PHA is a finer grid than ARF, evaluate model on PHA and
        # rebin down to the granularity that the ARF expects.
        if pha.bin_lo is not None and pha.bin_hi is not None:
            bin_lo, bin_hi = pha.bin_lo, pha.bin_hi

            # If PHA grid is in angstroms then convert to keV for
            # consistency
            if (bin_lo[0] > bin_lo[-1]) and (bin_hi[0] > bin_hi[-1]):
                bin_lo, bin_hi = DataPHA._hc/pha.bin_hi, DataPHA._hc/pha.bin_lo

            # FIXME: What about filtered option?? bin_lo, bin_hi are unfiltered??

            # Compare disparate grids in energy space
            self.arfargs = ((self.elo, self.ehi), (bin_lo, bin_hi))

            # FIXME: Assumes ARF grid is finest

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi
        if self.pha.units == 'wavelength':
            self.xlo, self.xhi = self.lo, self.hi


    def startup(self):
        arf = self._arf # original
        pha = self.pha

        # Create a view of original ARF
        self.arf = DataARF(arf.name, arf.energ_lo, arf.energ_hi, arf.specresp,
                           arf.bin_lo, arf.bin_hi, arf.exposure, arf.header)

        # Filter the view for current fitting session
        if numpy.iterable(pha.mask):
            mask = pha.get_mask()
            if len(mask) == len(self.arf.specresp):
                self.arf.notice(mask)

        self.filter()

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi
        if pha.units == 'wavelength':
            self.xlo, self.xhi = self.lo, self.hi

        ARFModel.startup(self)


    def teardown(self):
        self.arf = self._arf # restore original

        self.filter()
        ARFModel.teardown(self)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        # x could be channels or x, xhi could be energy|wave

        src = self.model.calc(p, self.xlo, self.xhi)
        return self.arf.apply_arf(src, *self.arfargs)


class ARFModelNoPHA(ARFModel):
    """
    ARF convolution model without associated PHA
    """
    def __init__(self, arf, model):
        ARFModel.__init__(self, arf, model)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        # x could be channels or x, xhi could be energy|wave

        #if (xhi is not None and
        #    x[0] > x[-1] and xhi[0] > xhi[-1]):
        #    xlo, xhi = self.lo, self.hi
        #else:

        # Always evaluates source model in keV!
        src = self.model.calc(p, self.xlo, self.xhi)
        return self.arf.apply_arf(src)



class RSPModelPHA(RSPModel):
    """
    RMF + ARF convolution model with associated PHA
    """
    def __init__(self, arf, rmf, pha, model):
        self.pha = pha
        self._arf = arf
        self._rmf = rmf
        RSPModel.__init__(self, arf, rmf, model)


    def filter(self):

        RSPModel.filter(self)

        pha = self.pha
        # If PHA is a finer grid than RMF, evaluate model on PHA and
        # rebin down to the granularity that the RMF expects.
        if pha.bin_lo is not None and pha.bin_hi is not None:
            bin_lo, bin_hi = pha.bin_lo, pha.bin_hi

            # If PHA grid is in angstroms then convert to keV for
            # consistency
            if (bin_lo[0] > bin_lo[-1]) and (bin_hi[0] > bin_hi[-1]):
                bin_lo, bin_hi = DataPHA._hc/pha.bin_hi, DataPHA._hc/pha.bin_lo

            # FIXME: What about filtered option?? bin_lo, bin_hi are unfiltered??

            # Compare disparate grids in energy space
            self.arfargs = ((self.elo, self.ehi), (bin_lo, bin_hi))

            # FIXME: Assumes ARF grid is finest

        elo, ehi = self.rmf.get_indep()
        # self.elo, self.ehi are from ARF
        if len(elo) != len(self.elo) and len(ehi) != len(self.ehi):

            self.rmfargs = ((elo, ehi), (self.elo, self.ehi))

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi
        if self.pha.units == 'wavelength':
            self.xlo, self.xhi = self.lo, self.hi


    def startup(self):
        arf = self._arf
        rmf = self._rmf

        # Create a view of original RMF
        self.rmf = DataRMF(rmf.name, rmf.detchans, rmf.energ_lo, rmf.energ_hi,
                           rmf.n_grp, rmf.f_chan, rmf.n_chan, rmf.matrix,
                           rmf.offset, rmf.e_min, rmf.e_max, rmf.header)

        # Create a view of original ARF
        self.arf = DataARF(arf.name, arf.energ_lo, arf.energ_hi, arf.specresp,
                           arf.bin_lo, arf.bin_hi, arf.exposure, arf.header)

        # Filter the view for current fitting session
        _notice_resp(self.pha.get_noticed_channels(), self.arf, self.rmf)

        self.filter()

        # Assume energy as default spectral coordinates
        self.xlo, self.xhi = self.elo, self.ehi
        if self.pha.units == 'wavelength':
            self.xlo, self.xhi = self.lo, self.hi

        RSPModel.startup(self)


    def teardown(self):
        self.arf = self._arf  # restore originals
        self.rmf = self._rmf

        self.filter()
        RSPModel.teardown(self)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        # x could be channels or x, xhi could be energy|wave

        src = self.model.calc(p, self.xlo, self.xhi)
        src = self.arf.apply_arf(src, *self.arfargs)
        return self.rmf.apply_rmf(src, *self.rmfargs)


class RSPModelNoPHA(RSPModel):
    """
    RMF + ARF convolution model without associated PHA
    """
    def __init__(self, arf, rmf, model):
        RSPModel.__init__(self, arf, rmf, model)


    def calc(self, p, x, xhi=None, *args, **kwargs):
        # x could be channels or x, xhi could be energy|wave

        # Always evaluates source model in keV!
        src = self.model.calc(p, self.xlo, self.xhi)
        src = self.arf.apply_arf(src, *self.arfargs)
        return self.rmf.apply_rmf(src, *self.rmfargs)


class ARF1D(NoNewAttributesAfterInit):

    def __init__(self, arf, pha=None, rmf=None):
        self._arf = arf
        self._pha = pha
        NoNewAttributesAfterInit.__init__(self)

    def __getattr__(self, name):
        arf = None
        try:
            arf = ARF1D.__getattribute__(self, '_arf')
        except:
            pass

        if name in ('_arf', '_pha'):
            return self.__dict__[name]

        if arf is not None:
            return DataARF.__getattribute__(arf, name)

        return ARF1D.__getattribute__(self, name)

    def __setattr__(self, name, val):
        arf = None
        try:
            arf = ARF1D.__getattribute__(self, '_arf')
        except:
            pass

        if arf is not None and hasattr(arf, name):
            DataARF.__setattr__(arf, name, val)
        else:
            NoNewAttributesAfterInit.__setattr__(self, name, val)

    def __dir__(self):
        return dir(self._arf)

    def __str__(self):
        return str(self._arf)

    def __repr__(self):
        return repr(self._arf)

    def __call__(self, model, session=None):
        arf = self._arf
        pha = self._pha

	if isinstance(model, basestring):
		if session is None:
			model = sherpa.astro.ui._session._eval_model_expression(model)
		else:
			model = session._eval_model_expression(model)
		

        # Automatically add exposure time to source model
        if pha is not None and pha.exposure is not None:
            model = pha.exposure * model
        elif arf.exposure is not None:
            model = arf.exposure * model
        # FIXME: display a warning if exposure is None?

        if pha is not None:
            return ARFModelPHA(arf, pha, model)

        return ARFModelNoPHA(arf, model)


class RMF1D(NoNewAttributesAfterInit):

    def __init__(self, rmf, pha=None, arf=None):
        self._rmf = rmf
        self._arf = arf
        self._pha = pha
        NoNewAttributesAfterInit.__init__(self)

    def __getattr__(self, name):
        rmf = None
        try:
            rmf = RMF1D.__getattribute__(self, '_rmf')
        except:
            pass

        if name in ('_rmf', '_pha'):
            return self.__dict__[name]

        if rmf is not None:
            return DataRMF.__getattribute__(rmf, name)

        return RMF1D.__getattribute__(self, name)


    def __setattr__(self, name, val):
        rmf = None
        try:
            rmf = RMF1D.__getattribute__(self, '_rmf')
        except:
            pass

        if rmf is not None and hasattr(rmf, name):
            DataRMF.__setattr__(rmf, name, val)
        else:
            NoNewAttributesAfterInit.__setattr__(self, name, val)

    def __dir__(self):
        return dir(self._rmf)

    def __str__(self):
        return str(self._rmf)

    def __repr__(self):
        return repr(self._rmf)

    def __call__(self, model, session=None):
        arf = self._arf
        rmf = self._rmf
        pha = self._pha

	if isinstance(model, basestring):
		if session is None:
			model = sherpa.astro.ui._session._eval_model_expression(model)
		else:
			model = session._eval_model_expression(model)

        # Automatically add exposure time to source model for RMF-only analysis
        if type(model) not in (ARFModel,ARFModelPHA,ARFModelNoPHA):

            if pha is not None and pha.exposure is not None:
                model = pha.exposure * model
            elif arf is not None and arf.exposure is not None:
                model = arf.exposure * model
        elif pha is not None and arf is not None:
            # If model is an ARF?
            # Replace RMF(ARF(SRC)) with RSP(SRC) for efficiency
            return RSPModelPHA(arf, rmf, pha, model.model)

        if pha is not None:
            return RMFModelPHA(rmf, pha, model)

        return RMFModelNoPHA(rmf, model)


class Response1D(NoNewAttributesAfterInit):

    def __init__(self, pha):
        self.pha = pha
        arf, rmf = pha.get_response()
        if arf is None and rmf is None:
            raise DataErr('norsp', pha.name)

        NoNewAttributesAfterInit.__init__(self)

    def __call__(self, model, session=None):
        pha = self.pha
        arf, rmf = pha.get_response()

	if isinstance(model, basestring):
		if session is None:
			model = sherpa.astro.ui._session._eval_model_expression(model)
		else:
			model = session._eval_model_expression(model)

        # Automatically add exposure time to source model
        if pha.exposure is not None:
            model = pha.exposure * model
        elif arf is not None and arf.exposure is not None:
            model = arf.exposure * model

        if arf is not None and rmf is not None:
            return RSPModelPHA(arf, rmf, pha, model)

        if arf is not None:
            return ARFModelPHA(arf, pha, model)

        if rmf is not None:
            return RMFModelPHA(rmf, pha, model)

        raise DataErr('norsp', pha.name)


class ResponseNestedModel(Model):

    def __init__(self, arf=None, rmf=None):
        self.arf = arf
        self.rmf = rmf

        name = ''
        if arf is not None and rmf is not None:
            name = 'apply_rmf(apply_arf('
        elif arf is not None:
            name = 'apply_arf('
        elif rmf is not None:
            name = 'apply_rmf('
        Model.__init__(self, name)


    def calc(self, p, *args, **kwargs):
        arf = self.arf
        rmf = self.rmf

        if arf is not None and rmf is not None:
            return rmf.apply_rmf(arf.apply_arf(*args, **kwargs))
        elif self.arf is not None:
            return arf.apply_arf(*args, **kwargs)

        return rmf.apply_rmf(*args, **kwargs)


class MultiResponseSumModel(CompositeModel, ArithmeticModel):

    def __init__(self, source, pha):
        self.channel = pha.channel
        self.mask = numpy.ones(len(pha.channel), dtype=bool)
        self.pha = pha
        self.source = source
        self.elo = None; self.ehi = None
        self.lo = None; self.hi = None
        self.table = None
        self.orders = None

        models = []
        grid = []

        for id in pha.response_ids:
            arf, rmf = pha.get_response(id)

            if arf is None and rmf is None:
                raise DataErr('norsp', pha.name)

            m = ResponseNestedModel(arf, rmf)
            indep = None

            if arf is not None:
                indep = arf.get_indep()

            if rmf is not None:
                indep = rmf.get_indep()

            models.append(m)
            grid.append(indep)

        self.models = models
        self.grid = grid

        name = '%s(%s)' % (type(self).__name__,
                           ','.join(['%s(%s)' % (m.name, source.name)
                                     for m in models]))
        CompositeModel.__init__(self, name, (source,))


    def _get_noticed_energy_list(self):
        grid = []
        for id in self.pha.response_ids:
            arf, rmf = self.pha.get_response(id)
            indep = None
            if arf is not None:
                indep = arf.get_indep()
            elif rmf is not None:
                indep = rmf.get_indep()
            grid.append(indep)

        self.elo, self.ehi, self.table = compile_energy_grid(grid)
        self.lo, self.hi = DataPHA._hc/self.ehi, DataPHA._hc/self.elo


    def startup(self):
        pha = self.pha
        if numpy.iterable(pha.mask):
            pha.notice_response(True)
        self.channel = pha.get_noticed_channels()
        self.mask = pha.get_mask()
        self._get_noticed_energy_list()
        CompositeModel.startup(self)


    def teardown(self):
        pha = self.pha
        if numpy.iterable(pha.mask):
            pha.notice_response(False)
        self.channel = pha.channel
        self.mask = numpy.ones(len(pha.channel), dtype=bool)
        self.elo = None; self.ehi = None; self.table = None
        self.lo = None; self.hi = None
        CompositeModel.teardown(self)


    def _check_for_user_grid(self, x, xhi=None):
        return (len(self.channel) != len(x) or
                not (sao_fcmp(self.channel, x, _tol)==0).all())


    def _startup_user_grid(self, x, xhi=None):
        # fit() never comes in here b/c it calls startup()
        pha = self.pha
        self.mask = numpy.zeros(len(pha.channel), dtype=bool)
        self.mask[numpy.searchsorted(pha.channel, x)]=True
        pha.notice_response(True, x)
        self._get_noticed_energy_list()


    def _teardown_user_grid(self):
        # fit() never comes in here b/c it calls startup()
        pha = self.pha
        self.mask = numpy.ones(len(pha.channel), dtype=bool)
        pha.notice_response(False)
        self.elo = None; self.ehi = None; self.table = None
        self.lo = None; self.hi = None


    def calc(self, p, x, xhi=None, *args, **kwargs):
        pha = self.pha

        user_grid = False
        try:

            if self._check_for_user_grid(x, xhi):
                user_grid = True
                self._startup_user_grid(x, xhi)

            # Slow
            if self.table is None:
                # again, fit() never comes in here b/c it calls startup()
                src = self.source
                vals = []
                for model, args in izip(self.models, self.grid):
                    elo,ehi = lo,hi = args
                    if pha.units == 'wavelength':
                        lo = DataPHA._hc / ehi
                        hi = DataPHA._hc / elo
                    vals.append(model(src(lo, hi)))
                self.orders = vals
            # Fast
            else:
                xlo,xhi = self.elo, self.ehi
                if pha.units == 'wavelength':
                    xlo, xhi = self.lo, self.hi

                src = self.source(xlo, xhi)  # hi-res grid of all ARF grids

                # Fold summed intervals through the associated response.
                self.orders = \
                    [model(sum_intervals(src, interval[0], interval[1]))
                     for model, interval in izip(self.models, self.table)]

            vals = sum(self.orders)
            if self.mask is not None:
                vals = vals[self.mask]

        finally:
            if user_grid:
                self._teardown_user_grid()


        return vals


class MultipleResponse1D(Response1D):

    def __call__(self, model, session=None):
        pha = self.pha

	if isinstance(model, basestring):
		if session is None:
			model = sherpa.astro.ui._session._eval_model_expression(model)
		else:
			model = session._eval_model_expression(model)

        pha.notice_response(False)

        model = MultiResponseSumModel(model, pha)

        if pha.exposure:
            model = pha.exposure * model

        return model


class PileupRMFModel(CompositeModel, ArithmeticModel):

    def __init__(self, rmf, model, pha=None):
        self.pha = pha
        self.channel = sao_arange(1, rmf.detchans)  # sao_arange is inclusive
        self.mask = numpy.ones(rmf.detchans, dtype=bool)
        self.rmf = rmf

        self.elo, self.ehi = rmf.get_indep()
        self.lo, self.hi = DataPHA._hc/self.ehi, DataPHA._hc/self.elo
        self.model = model
        self.otherargs = None
        self.otherkwargs = None
        self.pars = ()
        CompositeModel.__init__(self,
                                ('%s(%s)' % ('apply_rmf', self.model.name)),
                                (self.model,))

    def startup(self):
        pha = self.pha
        pha.notice_response(False)
        self.channel = pha.get_noticed_channels()
        self.mask = pha.get_mask()
        self.model.startup()
        CompositeModel.startup(self)


    def teardown(self):
        pha = self.pha
        rmf = self.rmf
        self.channel = sao_arange(1, rmf.detchans)
        self.mask = numpy.ones(rmf.detchans, dtype=bool)
        self.model.teardown()
        CompositeModel.teardown(self)


    def _check_for_user_grid(self, x):
        return (len(self.channel) != len(x) or
                not (sao_fcmp(self.channel, x, _tol)==0).all())


    def _startup_user_grid(self, x):
        # fit() never comes in here b/c it calls startup()
        self.mask = numpy.zeros(self.rmf.detchans, dtype=bool)
        self.mask[numpy.searchsorted(self.pha.channel, x)]=True


    def _calc(self, p, xlo, xhi):
        # Evaluate source model on RMF energy/wave grid OR
        # model.calc --> pileup_model
        src = self.model.calc(p, xlo, xhi)

        # rmf_fold
        return self.rmf.apply_rmf(src)


    def calc(self, p, x, xhi=None, **kwargs):
        pha = self.pha
        # x is noticed/full channels here

        user_grid = False
        try:
            if self._check_for_user_grid(x):
                user_grid = True
                self._startup_user_grid(x)

            xlo, xhi = self.elo,self.ehi
            if pha is not None and pha.units == 'wavelength':
                xlo, xhi = self.lo, self.hi

            vals = self._calc(p, xlo, xhi)
            if self.mask is not None:
                vals = vals[self.mask]

        finally:
            if user_grid:
                self.mask = numpy.ones(self.rmf.detchans, dtype=bool)

        return vals


class PileupResponse1D(NoNewAttributesAfterInit):

    def __init__(self, pha, pileup_model):
        self.pha = pha
        self.pileup_model = pileup_model
        NoNewAttributesAfterInit.__init__(self)

    def __call__(self, model, session=None):
        pha = self.pha
        # clear out any previous response filter
        pha.notice_response(False)

	if isinstance(model, basestring):
		if session is None:
			model = sherpa.astro.ui._session._eval_model_expression(model)
		else:
			model = session._eval_model_expression(model)

        arf, rmf = pha.get_response()
        err_msg = None

        if arf is None and rmf is None:
            raise DataErr('norsp', pha.name)

        if arf is None:
            err_msg = 'does not have an associated ARF'
        elif pha.exposure is None:
            err_msg = 'does not specify an exposure time'

        if err_msg:
            raise InstrumentErr('baddata', pha.name, err_msg)

        # Currently, the response is NOT noticed using pileup

        # ARF convolution done inside ISIS pileup module
        # on finite grid scale
        model = model.apply(self.pileup_model, pha.exposure, arf.energ_lo,
                            arf.energ_hi, arf.specresp, model)

        if rmf is not None:
            model = PileupRMFModel(rmf, model, pha)
        return model


class PSFModel(_PSFModel):

    def fold(self, data):
        _PSFModel.fold(self, data)

        # Set WCS coordinates of kernel data set to match source data set.
        if (isinstance(data, DataIMG) and
            isinstance(self.kernel, DataIMG)):
            self.kernel.set_coord(data.coord)


    def get_kernel(self, data, subkernel=True):

        indep, dep, kshape, lo, hi = self._get_kernel_data(data, subkernel)

        # Use kernel data set WCS if available
        eqpos = getattr(self.kernel, 'eqpos', None)
        sky   = getattr(self.kernel, 'sky', None)

        # If kernel is a model, use WCS from data if available
        if callable(self.kernel):
            eqpos = getattr(data, 'eqpos', None)
            sky   = getattr(data, 'sky', None)

        dataset = None
        ndim = len(kshape)
        if ndim == 1:
            dataset = Data1D('kernel', indep[0], dep)

        elif ndim == 2:

            # Edit WCS to reflect the subkernel extraction in
            # physical coordinates.
            if (subkernel and sky is not None and
                lo is not None and hi is not None):

                if (WCS != None):
                    sky = WCS(sky.name, sky.type, sky.crval,
                              sky.crpix - lo, sky.cdelt, sky.crota,
                              sky.epoch, sky.equinox)

                # FIXME: Support for WCS only (non-Chandra) coordinate
                # transformations?

            dataset = DataIMG('kernel', indep[0], indep[1], dep,
                              kshape[::-1], sky=sky, eqpos=eqpos)
        else:
            raise PSFErr('ndim')

        return dataset
