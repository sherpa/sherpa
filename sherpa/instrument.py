from __future__ import division
#
#  Copyright (C) 2008, 2016, 2018, 2019, 2020
#        Smithsonian Astrophysical Observatory
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
import math
import warnings

import numpy
import logging

from sherpa.astro.data import DataIMG
from sherpa.data import Data, Data1D, Data2D
from sherpa.models import ArithmeticModel, ArithmeticConstantModel, \
    ArithmeticFunctionModel, CompositeModel, Model
from sherpa.models.parameter import Parameter
from sherpa.models.regrid import EvaluationSpace2D, rebin_2d
from sherpa.utils import bool_cast, NoNewAttributesAfterInit
from sherpa.utils.err import PSFErr
from sherpa.utils._psf import extract_kernel, get_padsize, normalize, \
    pad_data, set_origin, tcdData, unpad_data

import sherpa
info = logging.getLogger(__name__).info

string_types = (str, )


__all__ = ('Kernel', 'PSFKernel', 'RadialProfileKernel', 'PSFModel',
           'ConvolutionModel', 'PSFSpace2D')


class ConvolutionModel(CompositeModel, ArithmeticModel):

    @staticmethod
    def wrapobj(obj):
        if isinstance(obj, ArithmeticModel):
            return obj
        return ArithmeticFunctionModel(obj)

    @staticmethod
    def wrapkern(obj):
        if isinstance(obj, ArithmeticModel):
            return obj
        elif callable(obj):
            return ArithmeticFunctionModel(obj)
        return ArithmeticConstantModel(obj, 'kernel')

    def __init__(self, lhs, rhs, psf):
        self.lhs = self.wrapkern(lhs)
        self.rhs = self.wrapobj(rhs)
        self.psf = psf
        CompositeModel.__init__(self,
                                ('%s(%s)' %
                                 (self.psf.name, self.rhs.name)),
                                (self.psf, self.lhs, self.rhs))

    def calc(self, p, *args, **kwargs):
        nlhs = len(self.lhs.pars)
        return self.psf.calc(p[:nlhs], p[nlhs:],
                             self.lhs.calc, self.rhs.calc, *args, **kwargs)


class Kernel(NoNewAttributesAfterInit):
    "Base class for convolution kernels"

    def __init__(self, dshape, kshape, norm=False, frozen=True,
                 center=None, args=[], kwargs={},
                 do_pad=False, pad_mask=None, origin=None):

        if origin is None:
            origin = numpy.zeros(len(kshape))
        self.dshape = dshape
        self.kshape = kshape
        self.kernel = None
        self.skshape = None
        self.norm   = norm
        self.origin = origin
        self.frozen = frozen
        self.center = center
        self.args = args
        self.kwargs = kwargs
        self.renorm_shape = None
        self.renorm = None
        self.do_pad = do_pad
        self.pad_mask = pad_mask
        self.frac = None
        self._tcd = tcdData()
        NoNewAttributesAfterInit.__init__(self)

    def __setstate__(self, state):
        state['_tcd'] = tcdData()
        self.__dict__.update(state)

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop('_tcd')
        return state

    def __repr__(self):
        return "<%s kernel instance>" % type(self).__name__

    def __str__(self):
        ss = [
            'dshape   = %s' % str(self.dshape),
            'kshape   = %s' % str(self.kshape),
            #            'kernel   = %s' % type(self.kernel).__name__,
            'skshape  = %s' % str(self.skshape),
            'norm     = %s' % str(self.norm),
            'origin   = %s' % str(self.origin),
            'frozen   = %s' % str(self.frozen),
            'center   = %s' % str(self.center),
            'args     = %s' % str(self.args),
            'kwargs   = %s' % str(self.kwargs),
            'renorm_shape  = %s' % str(self.renorm_shape),
            'renorm   = %s' % str(self.renorm),
            'do_pad   = %s' % str(self.do_pad),
            'pad_mask = %s' % str(self.pad_mask),
            'frac     = %s' % str(self.frac)
        ]
        return '\n'.join(ss)

    def init_kernel(self, kernel):
        if not self.frozen:
            self._tcd.clear_kernel_fft()

        renorm_shape = []
        for axis in self.dshape:
            renorm_shape.append(get_padsize(2 * axis))
        self.renorm_shape = tuple(renorm_shape)

        kernpad = pad_data(kernel, self.dshape, self.renorm_shape)

        self.renorm = self._tcd.convolve(numpy.ones(len(kernel)), kernpad,
                                         self.dshape, renorm_shape,
                                         self.origin)
        self.renorm = unpad_data(self.renorm, renorm_shape, self.dshape)
        return (kernel, self.dshape)

    def init_data(self, data):
        if self.renorm_shape is None:
            renorm_shape = []
            for axis in self.dshape:
                renorm_shape.append(get_padsize(2 * axis))
            self.renorm_shape = tuple(renorm_shape)

        # pad the data and convolve with unpadded kernel
        datapad = pad_data(data, self.dshape, self.renorm_shape)
        return (datapad, self.renorm_shape)

    def deinit(self, vals):
        if self.renorm is not None:
            vals = unpad_data(vals, self.renorm_shape, self.dshape)
            vals = vals / self.renorm
        if self.do_pad:
            vals = vals[self.pad_mask]
        return vals

    def convolve(self, data, dshape, kernel, kshape):
        return self._tcd.convolve(data, kernel, dshape, kshape, self.origin)

    def calc(self, pl, pr, lhs, rhs, *args, **kwargs):
        if self.do_pad and len(args[0]) == numpy.prod(self.dshape):
            self.do_pad = False

        data = rhs(pr, *self.args, **self.kwargs)
        (data, dshape) = self.init_data(data)

        if self.kernel is None or not self.frozen:
            kernel = lhs(pl, *self.args, **self.kwargs)
            (self.kernel, self.skshape) = self.init_kernel(kernel)

        vals = self.convolve(data, dshape, self.kernel, self.skshape)

        return self.deinit(vals)


class ConvolutionKernel(Model):

    def __init__(self, kernel, name='conv'):
        self.kernel = kernel
        self.name = name
        self._tcd = tcdData()
        Model.__init__(self, name)

    def __setstate__(self, state):
        state['_tcd'] = tcdData()
        self.__dict__.update(state)

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop('_tcd')
        return state

    def __repr__(self):
        return "<%s kernel instance>" % type(self).__name__

    def __str__(self):
        if self.kernel is None:
            raise PSFErr('notset')
        return "Convolution Kernel:\n" + self.kernel.__str__()

    def __call__(self, model, session=None):
        if self.kernel is None:
            raise PSFErr('notset')
        kernel = self.kernel
        if isinstance(kernel, Data):
            kernel = numpy.asarray(kernel.get_dep())

        if isinstance(model, string_types):
            if session is None:
                model = sherpa.astro.ui._session._eval_model_expression(model)
            else:
                model = session._eval_model_expression(model)

        return ConvolutionModel(kernel, model, self)

    def set_kernel(self, kernel):
        self.kernel = kernel

    def calc(self, pl, pr, lhs, rhs, *args, **kwargs):

        self._tcd.clear_kernel_fft()

        data = numpy.asarray(rhs(pr, *args, **kwargs))
        kern = numpy.asarray(lhs(pl, *args, **kwargs))

        size = data.size
        return self._tcd.convolve(data, kern, size, kern.size,
                                  int(size / 2))[:size]


class PSFKernel(Kernel):
    "class for PSF convolution kernels"

    def __init__(self, dshape, kshape, is_model=False, norm=True, frozen=True,
                 center=None, size=None, lo=None, hi=None, width=None,
                 args=[], kwargs={},
                 pad_mask=None, do_pad=False, origin=None):

        self.is_model = is_model
        self.size = size
        self.lo = lo
        self.hi = hi
        self.width = width
        self.radial = 0
        Kernel.__init__(self, dshape, kshape, norm, frozen,
                        center, args, kwargs,
                        do_pad, pad_mask, origin)
        self.origin = origin

    def __str__(self):
        ss = [
            'is_model = %s' % str(self.is_model),
            'size     = %s' % str(self.size),
            'lo       = %s' % str(self.lo),
            'hi       = %s' % str(self.hi),
            'width    = %s' % str(self.width),
            'radial   = %s' % str(self.radial)
        ]
        return Kernel.__str__(self) + '\n' + '\n'.join(ss)

    def init_kernel(self, kernel):
        # If PSF dataset, normalize before kernel extraction
        # if not self.is_model and self.norm:
        if self.norm:
            kernel = normalize(kernel)

        (kernel, kshape, self.frac,
         lo, hi) = extract_kernel(kernel, self.kshape, self.size, self.center,
                                  self.lo, self.hi, self.width, self.radial)

        # If PSF model, then normalize integrated volume to 1, after
        # kernel extraction
        # if self.is_model and self.norm:
        #    self.frac = 1.0
        #    kernel = normalize(kernel)

        # Find brightest pixel of PSF--assume that is the origin
        # Just assuming that the origin is half of szs1 can lead to
        # unwanted pixel shifts--but this assumes that origin should
        # be centered on brightest pixel.
        brightPixel = list(numpy.where(kernel == kernel.max())).pop()

        # if more than one pixel qualifies as brightest, such as const2D
        # use the middle of subkernel -- assumes the user provided center at
        # time of kernel extraction, so that should be middle of subkernel.
        origin = None
        if (not numpy.isscalar(brightPixel)) and len(brightPixel) != 1:
            origin = set_origin(kshape)
        else:
            # brightPixel is a NumPy index (int64) which - as of NumPy 1.18
            # and Python 3.8 - causes a TypeError with the message
            # "only integer scalar arrays can be converted to a scalar index"
            # to be thrown here if sent directly to set_origin. So
            # we convert to a Python integer type.
            #
            origin = set_origin(kshape, int(brightPixel))

        if self.origin is None:
            self.origin = origin

        if self.is_model and not self.frozen:
            # if the kernel model has thawed parameters, clear the old FFT and
            # recompute the kernel FFT at each model evaluation
            self._tcd.clear_kernel_fft()

        return (kernel, kshape)

    def init_data(self, data):
        return (data, self.dshape)


class RadialProfileKernel(PSFKernel):
    "class for 1D radial profile PSF convolution kernels"

    def __init__(self, dshape, kshape, is_model=False,
                 norm=True, frozen=True,
                 center=None, size=None, lo=None, hi=None, width=None,
                 args=[], kwargs={},
                 pad_mask=None, do_pad=False, origin=None):

        self.radialsize = None
        PSFKernel.__init__(self, dshape, kshape, is_model, norm,
                           frozen, center, size, lo, hi, width, args, kwargs,
                           pad_mask, do_pad, origin)
        self.radial = 1

    def __str__(self):
        return (PSFKernel.__str__(self) + '\n' +
                'radialsize = %s' % str(self.radialsize))

    def init_data(self, data):
        data, dshape = PSFKernel.init_data(self, data)
        # NOTICE: radial profile is 1D only!
        if self.radialsize is None:
            self.radialsize = self.dshape[0]
        return data, dshape

    def deinit(self, vals):
        # NOTICE: radial profile is 1D only!
        vals = vals[:self.radialsize]
        return PSFKernel.deinit(self, vals)

    def convolve(self, data, dshape, kernel, kshape):
        origin = self.origin
        if self.radialsize is not None:
            origin = self.origin + (numpy.asarray(dshape) -
                                    numpy.asarray(self.radialsize))
        return self._tcd.convolve(data, kernel, dshape, kshape, origin)

    def calc(self, pl, pr, lhs, rhs, *args, **kwargs):
        if self.do_pad and len(args[0]) == numpy.prod(self.dshape):
            self.do_pad = False

        data = rhs(pr, *self.args, **self.kwargs)
        (data, dshape) = self.init_data(data)

        # NOTICE: radial profile is 1D only!
        # old sherpa source model grid extension to zero
        # (radial profile add core)
        tail_grid = _create_tail_grid(self.args)
        if tail_grid is not None:
            tail = rhs(pr, *tail_grid, **self.kwargs)
            data = numpy.concatenate([tail, data])
            dshape = (len(data),)

        if self.kernel is None or not self.frozen:
            kernel = lhs(pl, *self.args, **self.kwargs)
            (self.kernel, self.skshape) = self.init_kernel(kernel)

        vals = self.convolve(data, dshape, self.kernel, self.skshape)

        return self.deinit(vals)


def _create_tail_grid(axis_list):
    if len(axis_list) == 1:
        # non-binned axis
        grid = axis_list[0]
        # origsize = len(grid)
        width = grid[1] - grid[0]
        tail = numpy.arange(grid[0] - width, 0., -width)[::-1]
        return (tail,)

    elif len(axis_list) == 2:
        # binned axis
        gridlo, gridhi = axis_list
        # origsize = len(gridlo)
        width = (gridhi[0] - gridlo[0])
        mid = (gridlo[0] + gridhi[0]) / 2.
        mids = numpy.arange(mid, 0., -width)[::-1]
        taillo = mids - width / 2.
        tailhi = mids + width / 2.
        return (taillo, tailhi)

    return None


class PSFModel(Model):
    def __init__(self, name='psfmodel', kernel=None):
        self._name = name
        self._size = None
        self._origin = None
        self._center = None
        self._must_rebin = False
        self.radial = Parameter(name, 'radial', 0, 0, 1, hard_min=0,
                                hard_max=1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1, 0, 1, hard_min=0, hard_max=1,
                              alwaysfrozen=True)
        self.kernel = kernel
        self.model = None
        self.data_space = None
        self.psf_space = None
        Model.__init__(self, name)

    def _get_center(self):
        if self._center is not None:
            if len(self._center) == 1:
                return self._center[0]
        return self._center

    def _set_center(self, vals):
        par = vals
        if type(vals) in (str, numpy.string_):
            raise PSFErr('nostr')
        elif type(vals) not in (list, tuple, numpy.ndarray):
            par = [vals]
        self._center = tuple(par)
        if par is None:
            self._center = None

    center = property(_get_center, _set_center, doc='array of size parameters')

    def _get_size(self):
        if self._size is not None:
            if len(self._size) == 1:
                return self._size[0]
        return self._size

    def _set_size(self, vals):
        par = vals
        if type(vals) in (str, numpy.string_):
            raise PSFErr('notstr')
        elif type(vals) not in (list, tuple, numpy.ndarray):
            par = [vals]
        self._size = tuple(par)
        if par is None:
            self._size = None

    size = property(_get_size, _set_size, doc='array of size parameters')

    def _get_origin(self):
        if self._origin is not None:
            if len(self._origin) == 1:
                return self._origin[0]
        return self._origin

    def _set_origin(self, vals):
        par = vals
        if type(vals) in (str, numpy.string_):
            raise PSFErr('notstr')
        elif type(vals) not in (list, tuple, numpy.ndarray):
            par = [vals]
        self._origin = tuple(par)
        if par is None:
            self._origin = None

    origin = property(_get_origin, _set_origin, doc='FFT origin')

    def _get_str(self):
        s = ''
        if self.kernel is not None:
            s += ('\n   %-12s %-6s %12s' %
                  ('%s.kernel' % self._name, 'frozen',
                   self.kernel.name))
        if self.size is not None:
            s += ('\n   %-12s %-6s %12s %12s %12s' %
                  ('%s.size' % self._name, 'frozen',
                   self.size, self.size, self.size))
        if self.center is not None:
            s += ('\n   %-12s %-6s %12s %12s %12s' %
                  ('%s.center' % self._name, 'frozen',
                   self.center, self.center, self.center))
        if self.origin is not None:
            s += ('\n   %-12s %-6s %12s %12s %12s' %
                  ('%s.origin' % self._name, 'frozen',
                   self.origin, self.origin, self.origin))
        for p in [self.radial, self.norm]:
            s += ('\n   %-12s %-6s %12g %12g %12g %10s' %
                  (p.fullname, 'frozen', p.val, p.min, p.max, p.units))
        return s

    def __str__(self):
        s = self.name
        hfmt = '\n   %-12s %-6s %12s %12s %12s %10s'
        s += hfmt % ('Param', 'Type', 'Value', 'Min', 'Max', 'Units')
        s += hfmt % ('-' * 5, '-' * 4, '-' * 5, '-' * 3, '-' * 3, '-' * 5)
        s += self._get_str()
        return s

    def __call__(self, model, session=None):
        if self.kernel is None:
            raise PSFErr('notset')
        kernel = self.kernel
        if isinstance(kernel, Data):
            kernel = numpy.asarray(kernel.get_dep())

        if isinstance(model, string_types):
            if session is None:
                model = sherpa.astro.ui._session._eval_model_expression(model)
            else:
                model = session._eval_model_expression(model)

        return ConvolutionModel(kernel, model, self)

    def calc(self, *args, **kwargs):
        if self.model is None:
            raise PSFErr('nofold')
        psf_space_evaluation = self.model.calc(*args, **kwargs)

        if self._must_rebin:
            return rebin_2d(psf_space_evaluation, self.psf_space, self.data_space).ravel()
        else:
            return psf_space_evaluation

    def fold(self, data):
        # FIXME how will we know the native dimensionality of the
        # raveled model without the values?
        kargs = {}

        kshape = None
        dshape = data.get_dims()

        (size, center, origin,
         kargs['norm'], radial) = (self.size, self.center, self.origin,
                                   bool_cast(self.norm.val),
                                   int(self.radial.val))

        kargs['size'] = size
        kargs['center'] = center
        kargs['origin'] = origin
        kargs['is_model'] = False
        kargs['do_pad'] = False

        kargs['args'] = data.get_indep()

        pixel_size_comparison = self._check_pixel_size(data)

        if pixel_size_comparison == self.SAME_RESOLUTION:  # Don't do anything special
            self.data_space = EvaluationSpace2D(*data.get_indep())
            self._must_rebin = False
        elif pixel_size_comparison == self.BETTER_RESOLUTION:  # Evaluate model in PSF space
            self.data_space = EvaluationSpace2D(*data.get_indep())
            self.psf_space = PSFSpace2D(self.data_space, self, data.sky.cdelt)
            kargs['args'] = self.psf_space.grid
            dshape = self.psf_space.shape
            self._must_rebin = True
        else:  # PSF has worse resolution, error out
            raise AttributeError("The PSF has a worse resolution than the data.")

        if isinstance(self.kernel, Data):

            kshape = self.kernel.get_dims()
            # (kargs['lo'], kargs['hi'],
            # kargs['width']) = _get_axis_info(self.kernel.get_indep(), kshape)

            kargs['lo'] = [1] * len(kshape)
            kargs['hi'] = kshape
            kargs['width'] = [1] * len(kshape)

            if center is None:
                kargs['center'] = [int(dim / 2.) for dim in kshape]
                # update center param to default
                self.center = kargs['center']

            if size is None:
                kargs['size'] = kshape
                # update size param to default
                self.size = kargs['size']

        else:
            if (self.kernel is None) or (not callable(self.kernel)):
                raise PSFErr('nopsf', self._name)
            kshape = data.get_dims()
            # (kargs['lo'], kargs['hi'],
            # kargs['width']) = _get_axis_info(kargs['args'], dshape)

            kargs['lo'] = [1] * len(kshape)
            kargs['hi'] = kshape
            kargs['width'] = [1] * len(kshape)

            if center is None:
                kargs['center'] = [int(dim / 2.) for dim in dshape]
                # update center param to default
                self.center = kargs['center']

            if size is None:
                kargs['size'] = dshape
                # update size param to default
                self.size = kargs['size']

            kargs['is_model'] = True
            if hasattr(self.kernel, 'pars'):
                # freeze all PSF model parameters if not already.
                for par in self.kernel.pars:
                    par.freeze()

            if hasattr(self.kernel, 'thawedpars'):
                kargs['frozen'] = (len(self.kernel.thawedpars) == 0)

        is_kernel = (kargs['is_model'] and not kargs['norm'] and
                     len(kshape) == 1)
        # Handle noticed regions for convolution
        if numpy.iterable(data.mask):
            kargs['do_pad'] = True
            kargs['pad_mask'] = data.mask

        if is_kernel:
            for id in ['is_model', 'lo', 'hi', 'width', 'size']:
                kargs.pop(id)
            self.model = Kernel(dshape, kshape, **kargs)
            return

        if radial:
            self.model = RadialProfileKernel(dshape, kshape, **kargs)
            return

        self.model = PSFKernel(dshape, kshape, **kargs)
        return

    def _get_kernel_data(self, data, subkernel=True):
        self.fold(data)
        if self.kernel is None:
            raise PSFErr('notset')
        kernel = self.kernel
        dep = None
        indep = None
        lo = None
        hi = None

        if isinstance(kernel, Data):
            dep = numpy.asarray(kernel.get_dep())
            indep = kernel.get_indep()

        elif callable(kernel):
            dep = kernel(*self.model.args, **self.model.kwargs)
            indep = self.model.args

        kshape = self.model.kshape
        if subkernel:
            (dep, newshape) = self.model.init_kernel(dep)

            if (numpy.array(kshape) != numpy.array(newshape)).any():
                newindep = []
                for axis in indep:
                    args = extract_kernel(axis,
                                          self.model.kshape,
                                          self.model.size,
                                          self.model.center,
                                          self.model.lo,
                                          self.model.hi,
                                          self.model.width,
                                          self.model.radial)
                    newaxis = args[0]
                    lo = args[3]  # subkernel offsets (lower bound)
                    hi = args[4]  # subkernel offsets (upper bound)
                    newindep.append(newaxis)
                indep = newindep

            kshape = newshape

        if self.model.frac is not None:
            info('PSF frac: %s' % self.model.frac)

        if numpy.isscalar(kshape):
            kshape = [kshape]

        return (indep, dep, kshape, lo, hi)

    def get_kernel(self, data, subkernel=True):

        indep, dep, kshape, lo, hi = self._get_kernel_data(data, subkernel)

        ndim = len(kshape)
        if ndim == 1:
            dataset = Data1D('kernel', indep[0], dep)
        elif ndim == 2:
            dataset = Data2D('kernel', indep[0], indep[1], dep,
                             kshape[::-1])
        else:
            raise PSFErr('ndim')

        return dataset

    def _check_pixel_size(self, data):
        """
        If the data and PSF dot not have WCS information, assume the PSF has the same resolution
        as the image.
        Otherwise check the WCS information to determine the relative resolutions.

        We only check the resolution in one dimention and assume they are the same
        """
        if hasattr(self.kernel, "sky"):
            # This corresponds to the case when the kernel is actually a psf image, not just a model.

            try:
                psf_pixel_size = self.kernel.sky.cdelt
            except AttributeError: # If the kernel does not have a pixel size, issue a warning and keep going
                warnings.warn("PSF Image does not have a pixel size. Sherpa will assume the pixel size is the same"
                              "as the data")
                return self.SAME_RESOLUTION

            try:
                data_pixel_size = data.sky.cdelt
            except AttributeError:
                warnings.warn("Data Image does not have a pixel size. Sherpa will assume the pixel size is the same"
                              "as the PSF")
                return self.SAME_RESOLUTION

            if numpy.allclose(psf_pixel_size, data_pixel_size):
                return self.SAME_RESOLUTION
            if psf_pixel_size[0] < data_pixel_size[0]:
                return self.BETTER_RESOLUTION
            if psf_pixel_size[0] > data_pixel_size[0]:
                return self.WORSE_RESOLUTION

        return self.SAME_RESOLUTION

    SAME_RESOLUTION = 0
    BETTER_RESOLUTION = 1
    WORSE_RESOLUTION = -1


class PSFSpace2D(EvaluationSpace2D):
    """
    This class defines a special evaluation space that has the same boundaries as the data space but a number of pixels
    consistent with the pixel size of the PSF. This space is used when the data image and the PSF have a different
    pixel size so the model is evaluated on the "PSF space" and then rebinned back to the data space for the calculation
    of the statistic during a fit.
    """
    def __init__(self, data_space, psf_model, data_pixel_size=None):
        if data_pixel_size is None:
            data_pixel_size = [1, 1]

        x_start, y_start = data_space.start
        x_end, y_end = data_space.end
        psf_pixel_size_axis0 = psf_model.kernel.sky.cdelt[0]
        psf_pixel_size_axis1 = psf_model.kernel.sky.cdelt[1]
        data_pixel_size_axis0 = data_pixel_size[0]
        data_pixel_size_axis1 = data_pixel_size[1]
        step_x = psf_pixel_size_axis0 / data_pixel_size_axis0
        step_y = psf_pixel_size_axis1 / data_pixel_size_axis1
        x_range_end, y_range_end = x_end + 1, y_end + 1
        x = numpy.arange(x_start, x_range_end, step_x)
        y = numpy.arange(y_start, y_range_end, step_y)
        self.data_2_psf_pixel_size_ratio = (step_x, step_y)
        super(PSFSpace2D, self).__init__(x, y)
