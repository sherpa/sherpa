#
#  Copyright (C) 2008, 2016, 2018, 2019, 2020, 2021, 2022, 2023
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

import logging
import warnings

import numpy

from sherpa.data import Data, Data1D, Data2D
from sherpa.models import ArithmeticModel, ArithmeticConstantModel, \
    ArithmeticFunctionModel, CompositeModel, Model
from sherpa.models.parameter import Parameter
from sherpa.models.regrid import EvaluationSpace1D, EvaluationSpace2D, rebin_2d
from sherpa.utils import bool_cast, NoNewAttributesAfterInit
from sherpa.utils.err import PSFErr
from sherpa.utils._psf import extract_kernel, get_padsize, normalize, \
    pad_data, set_origin, tcdData, unpad_data

import sherpa
info = logging.getLogger(__name__).info

string_types = (str, )


__all__ = ('Kernel', 'PSFKernel', 'RadialProfileKernel', 'PSFModel',
           'ConvolutionModel', 'PSFSpace2D')


def make_renorm_shape(shape):
    """Given a shape, calculate the appropriate renorm_shape."""

    out = []
    for axis in shape:
        out.append(get_padsize(2 * axis))

    return out


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

        if callable(obj):
            return ArithmeticFunctionModel(obj)

        return ArithmeticConstantModel(obj, 'kernel')

    def __init__(self, lhs, rhs, psf):
        self.lhs = self.wrapkern(lhs)
        self.rhs = self.wrapobj(rhs)
        self.psf = psf
        CompositeModel.__init__(self,
                                f"{self.psf.name}({self.rhs.name})",
                                (self.psf, self.lhs, self.rhs))

    def calc(self, p, *args, **kwargs):
        nlhs = len(self.lhs.pars)
        return self.psf.calc(p[:nlhs], p[nlhs:],
                             self.lhs.calc, self.rhs.calc, *args, **kwargs)


class Kernel(NoNewAttributesAfterInit):
    """Base class for convolution kernels

    There are some validation checks made when the object is created
    but not when fields are changed. The assumption is that concepts
    like the dimensionality of the kernel are not going to be changed.

    """

    def __init__(self, dshape, kshape, norm=False, frozen=True,
                 center=None, args=[], kwargs={},
                 do_pad=False, pad_mask=None, origin=None):

        # As these are low-level routines use Python exceptions
        # rather than the Sherpa-specific ones.
        #
        try:
            nd = len(dshape)
        except TypeError:
            raise TypeError("dshape must be a sequence")

        try:
            nk = len(kshape)
        except TypeError:
            raise TypeError("kshape must be a sequence")

        if nd != nk:
            raise ValueError(f"dshape and kshape must be the same size, not {nd} and {nk}")

        if nd == 0:
            raise ValueError("0D kernel is not supported")

        # There is a specific PSFErr here so use it.
        if nd > 2:
            raise PSFErr("ndim")

        self.ndim = nd

        if origin is None:
            origin = numpy.zeros(self.ndim)

        self.dshape = dshape
        self.kshape = kshape
        self.kernel = None
        self.skshape = None
        self.norm = norm
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
        super().__init__()

    def __setstate__(self, state):
        state['_tcd'] = tcdData()
        self.__dict__.update(state)

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop('_tcd')
        return state

    def __repr__(self):
        return f"<{type(self).__name__} kernel instance>"

    def __str__(self):
        ss = [
            f"dshape   = {self.dshape}",
            f"kshape   = {self.kshape}",
            #            f"kernel   = {type(self.kernel).__name__}",
            f"skshape  = {self.skshape}",
            f"norm     = {self.norm}",
            f"origin   = {self.origin}",
            f"frozen   = {self.frozen}",
            f"center   = {self.center}",
            f"args     = {self.args}",
            f"kwargs   = {self.kwargs}",
            f"renorm_shape  = {self.renorm_shape}",
            f"renorm   = {self.renorm}",
            f"do_pad   = {self.do_pad}",
            f"pad_mask = {self.pad_mask}",
            f"frac     = {self.frac}"
        ]
        return "\n".join(ss)

    # The kernel is a 1D array and does not know it's original
    # dimensions.
    #
    def init_kernel(self, kernel):
        if not self.frozen:
            self._tcd.clear_kernel_fft()

        renorm_shape = make_renorm_shape(self.dshape)
        self.renorm_shape = tuple(renorm_shape)

        kernpad = pad_data(kernel, self.dshape, self.renorm_shape)

        renorm = self._tcd.convolve(numpy.ones(len(kernel)), kernpad,
                                    self.dshape, renorm_shape,
                                    self.origin)
        self.renorm = unpad_data(renorm, renorm_shape, self.dshape)
        return (kernel, self.dshape)

    def init_data(self, data):
        if self.renorm_shape is None:
            renorm_shape = make_renorm_shape(self.dshape)
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
        super().__init__(name)

    def __setstate__(self, state):
        state['_tcd'] = tcdData()
        self.__dict__.update(state)

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop('_tcd')
        return state

    def __repr__(self):
        return f"<{type(self).__name__} kernel instance>"

    def __str__(self):
        if self.kernel is None:
            raise PSFErr('notset')

        return f"Convolution Kernel:\n{self.kernel}"

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
        super().__init__(dshape, kshape, norm, frozen,
                         center, args, kwargs,
                         do_pad, pad_mask, origin)

        # The super-class handles a missing origin in a different
        # manner to this class.
        #
        if origin is None:
            self.origin = origin

    def __str__(self):
        ss = [
            f"is_model = {self.is_model}",
            f"size     = {self.size}",
            f"lo       = {self.lo}",
            f"hi       = {self.hi}",
            f"width    = {self.width}",
            f"radial   = {self.radial}"
        ]
        return Kernel.__str__(self) + "\n" + "\n".join(ss)

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
            # we convert to a Python integer type. In NumPy 1.25 it became
            # a deprecation error to call int on an array with ndim > 0.
            #
            # assume there is only one element in brightPixel if not
            # a scalar
            #
            if not numpy.isscalar(brightPixel):
                loc = brightPixel[0]
            else:
                loc = brightPixel
            origin = set_origin(kshape, int(loc))

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
        super().__init__(dshape, kshape, is_model, norm,
                         frozen, center, size, lo, hi, width, args, kwargs,
                         pad_mask, do_pad, origin)
        self.radial = 1  # over-ride super-class

        if self.ndim != 1:
            raise PSFErr(f"Radial profile requires 1D data, not {self.ndim}D")

    def __str__(self):
        return (PSFKernel.__str__(self) + "\n" +
                f"radialsize = {self.radialsize}")

    def init_data(self, data):
        data, dshape = super().init_data(data)

        # NOTICE: radial profile is 1D only!
        if self.radialsize is None:
            self.radialsize = self.dshape[0]

        return data, dshape

    def deinit(self, vals):
        # NOTICE: radial profile is 1D only!
        vals = vals[:self.radialsize]
        return super().deinit(vals)

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

    if len(axis_list) == 2:
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
    """Convolve a model by another model or data set.

    At the moment the code does not distinguish between 1D and 2D data
    and models.

    Parameters
    ----------
    name : str
        The name for the model.
    kernel : sherpa.data.Data instance, callable, or None, optional
        The kernel used to convolve models. This can be changed.

    Notes
    -----
    A number of attributes are displayed as parameters, if set, but
    are not handled as parameters. The attributes are: kernel, size,
    origin, and center. The size, center, and origin values are only
    displayed when set (and will be set by the `fold` method if
    needed). There is an attempt to ensure that the size, origin, and
    center fields have the correct size - that is they match the
    dimensionality of the kernel - but it is possible for an invalid
    combination to be set.

    """

    def __init__(self, name='psfmodel', kernel=None):

        # store the name without the leading "psfmodel." term that Model adds.
        self._name = name

        self._size = None
        self._origin = None
        self._center = None
        self._must_rebin = False
        self._model = None

        self._kernel = None

        self.radial = Parameter(name, 'radial', 0, 0, 1, hard_min=0,
                                hard_max=1, alwaysfrozen=True)
        self.norm = Parameter(name, 'norm', 1, 0, 1, hard_min=0, hard_max=1,
                              alwaysfrozen=True)

        self.kernel = kernel
        self.data_space = None
        self.psf_space = None
        super().__init__(name)

    def _get_kernel(self):
        return self._kernel

    def _set_kernel(self, kernel):

        # Always clear the model
        self._model = None

        odim = self.ndim
        if odim is None:
            # avid having to check for None in code below
            odim = 0

        def clear_fields():
            "Clear the array-like fields"
            if self.ndim is not None and self.ndim == odim:
                return

            self.size = None
            self.origin = None
            self.center = None

        if kernel is None:
            self._kernel = None
            self.ndim = None
            clear_fields()
            return

        if isinstance(kernel, Data):
            self._kernel = kernel
            self.ndim = kernel.ndim
            clear_fields()
            return

        if not callable(kernel):
            raise PSFErr('nopsf', self._name)

        # We could only allow a sherpa.models.model.Model instance here,
        # but allowable any callable, but that means the dimensionality
        # may not be set.
        #
        self._kernel = kernel
        try:
            self.ndim = getattr(kernel, "ndim")
        except AttributeError:
            # It's hard to trigger this case so we have no test coverage.
            self.ndim = None

        clear_fields()

    kernel = property(_get_kernel, _set_kernel,
                      doc="""The kernel (sherpa.data.Data or sherpa.models.model.Model instance, callable, or None).

The kernel determines the dimensionality of the model (the
`ndim` attribute), although it can remain as `None` for
callable arguments. The size, origin, and center fields, if
set, must match the `ndim` field, and will be cleared if
they do not match.
""")

    def _get_field(self, name):
        """Return the field value"""
        return getattr(self, name)

    def _set_field(self, name, vals):
        """Set the field value

        The value is checked for the correct size when possible.
        However we have to support self.ndim being None because a
        callable (or some specialised model) has been used, which
        makes it hard to ensure that everything matches. Setting the
        kernel will clear these fields if ndim changes (or is None)
        but it is still possible to get into cases where the fields do
        not match. These fields are only really meaningful after the
        fold method has been called.

        Parameters
        ----------
        name : str
            The name of the field (the attribute name).
        vals
            The value to set the field. It can not be a string, but does
            not need to be a sequence (so a scalar can be used).

        Notes
        -----
        The stored value is converted to a tuple if is not a tuple,
        list, or ndarray. This is an attempt to make sure that the
        field can be accessed as if it is a sequence (it is unclear
        what the original design intended for these fields, but
        existing code seems to require a tuple-like interface).

        """
        if vals is None:
            setattr(self, name, None)
            return

        # TODO: do we still expect to get bytes here?
        if isinstance(vals, (str, numpy.bytes_)):
            raise PSFErr('nostr')

        if not isinstance(vals, (list, tuple, numpy.ndarray)):
            vals = [vals]

        nvals = len(vals)
        if self.ndim is not None and nvals != self.ndim:
            # remove leading underscore from the name when reporting an error
            raise PSFErr("mismatch_dims", self.name, name[1:], self.ndim, nvals)

        if self.ndim == 1:
            vals = vals[0]
        else:
            vals = tuple(vals)

        setattr(self, name, vals)

    def _get_center(self):
        return self._get_field("_center")

    def _set_center(self, vals):
        self._set_field("_center", vals)

    center = property(_get_center, _set_center, doc='array of size parameters')

    def _get_size(self):
        return self._get_field("_size")

    def _set_size(self, vals):
        self._set_field("_size", vals)

    size = property(_get_size, _set_size, doc='array of size parameters')

    def _get_origin(self):
        return self._get_field("_origin")

    def _set_origin(self, vals):
        self._set_field("_origin", vals)

    origin = property(_get_origin, _set_origin, doc='FFT origin')

    @property
    def model(self):
        """The model that applies the convolution.

        This is set by the `fold` method and can not be changed
        directly.
        """
        return self._model

    def _set_model(self, model):
        """Set the model, after checking the dimensions

        This is not made into the model.setter property as the idea is
        that external users do not change this setting.

        """

        # This is tricky to trigger (in fact, the only existing uses
        # of this code does not set the model to None), so we do not
        # add a test.
        #
        if model is None:
            self._model = None
            return

        # Is this worthwhile, e.g.:
        # It's hard to trigger this case so we have no test coverage.
        if self.ndim is not None and self.ndim != model.ndim:
            raise PSFErr(f"Dimension of model do not match the kernel: {model.ndim}D and {self.ndim}D")

        self._model = model

    def _get_array_str(self, name, value):
        """Display the 'array-like' fields"""

        name = f"{self._name}.{name}"
        flag = "frozen"

        # We could be a bit-more clever about this conversion but does
        # not seem worth it.
        value = str(value)

        # Do we need to return the min/max values?
        return f"\n   {name:12s} {flag:6s} {value:>12s} {value:>12s} {value:>12s}"

    def _get_str(self):
        s = ''
        if self.kernel is not None:
            s += ('\n   %-12s %-6s %12s' %
                  ('%s.kernel' % self._name, 'frozen',
                   self.kernel.name))

        if self.size is not None:
            s += self._get_array_str("size", self.size)

        if self.center is not None:
            s += self._get_array_str("center", self.center)

        if self.origin is not None:
            s += self._get_array_str("origin", self.origin)

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

    def calc(self, p, *args, **kwargs):
        if self.model is None:
            raise PSFErr('nofold')

        psf_space_evaluation = self.model.calc(p, *args, **kwargs)

        if self._must_rebin:
            return rebin_2d(psf_space_evaluation, self.psf_space, self.data_space).ravel()

        return psf_space_evaluation

    def fold(self, data):
        """The data to be convolved by the PSF.

        Parameters
        ----------
        data : sherpa.data.Data or sherpa.models.model.Model instance or a callable
            It must match the dimensionality of the kernel.

        Raises
        ------
        sherpa.utils.err.PSFErr
            The kernel has not been set.

        """

        if self.kernel is None:
            raise PSFErr('nopsf', self._name)

        # TODO: Should we treat origin as we do center and size?
        kwargs = {"norm": bool_cast(self.norm.val),
                  "origin": self.origin
                  }

        # This validates the dimensionality of data. We can probably
        # remove the ndim checks below because of this, but I am not
        # convinced yet, so leave them in. This means that several
        # error conditions do not have test coverage.
        #
        (args, dshape) = self._create_spaces(data)
        kwargs['args'] = args

        if isinstance(self.kernel, Data):
            kwargs['is_model'] = False

            kshape = self.kernel.get_dims()
            nkernel = self.ndim

            if nkernel != data.ndim:
                raise PSFErr("mismatch_dims", self.kernel.name,
                             data.name, nkernel, data.ndim)

            if self.center is None:
                self.center = [int(dim / 2.) for dim in kshape]

            if self.size is None:
                self.size = kshape

        else:
            kwargs['is_model'] = True

            kshape = data.get_dims()
            nkernel = len(kshape)

            # To support using any callable, not just a model, we need
            # to allow the kernel dimensions to be unknown. The
            # alternative is to require the kernel to have a ndim
            # attribute, but there are a number of places the code
            # allows the kernel to not be a model instance (e.g. the
            # checks for pars/thawedpars).
            #
            if self.ndim is not None and nkernel != self.ndim:
                raise PSFErr("mismatch_dims", self.kernel.name, data.name, self.ndim, nkernel)

            if self.center is None:
                self.center = [int(dim / 2.) for dim in dshape]

            if self.size is None:
                self.size = dshape

            if hasattr(self.kernel, 'pars'):
                # freeze all PSF model parameters if not already.
                for par in self.kernel.pars:
                    par.freeze()

            # TODO: shouldn't thawedpars always be True with the above?
            #
            if hasattr(self.kernel, 'thawedpars'):
                kwargs['frozen'] = (len(self.kernel.thawedpars) == 0)

        kwargs['center'] = self.center
        kwargs['size'] = self.size

        # TODO: Why is this restricted to 1D?
        is_kernel = (kwargs['is_model'] and not kwargs['norm'] and
                     nkernel == 1)

        # Handle noticed regions for convolution
        if numpy.iterable(data.mask):
            kwargs['do_pad'] = True
            kwargs['pad_mask'] = data.mask
        else:
            kwargs['do_pad'] = False

        if is_kernel:
            for kwarg in ['is_model', 'size']:
                kwargs.pop(kwarg)

            self._set_model(Kernel(dshape, kshape, **kwargs))
            return

        # TODO:
        # If these are not set some tests seem to go into an infinite loop
        # eg calling a convolved model in
        # sherpa/models/tests/test_regrid_unit.py::test_regrid1d_works_with_convolution_style
        # Does this indicate that there should be better argument checking
        # or defaults?
        #
        kwargs['lo'] = numpy.ones(nkernel)
        kwargs['hi'] = kshape
        kwargs['width'] = numpy.ones(nkernel)

        # TODO: why is this not just checking 'self.radial.val > 0'
        # instead of 'int(self.radial.val)'? Aren't they the same, and
        # the former is clearer?
        #
        if int(self.radial.val):
            self._set_model(RadialProfileKernel(dshape, kshape, **kwargs))
            return

        self._set_model(PSFKernel(dshape, kshape, **kwargs))

    def _get_kernel_data(self, data, subkernel=True):
        if self.kernel is None:
            raise PSFErr('notset')

        self.fold(data)

        kernel = self.kernel
        if isinstance(kernel, Data):
            dep = numpy.asarray(kernel.get_dep())
            indep = kernel.get_indep()

        else:
            dep = kernel(*self.model.args, **self.model.kwargs)
            indep = self.model.args

        kshape = self.model.kshape
        lo = None
        hi = None

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
                    newindep.append(args[0])

                    # TODO: shouldn't we store these values like we do
                    # newindep?  This looks to be a bug but we do not
                    # have tests to check the expected behavior (it
                    # requires DataIMG + WCS data).
                    #
                    lo = args[3]
                    hi = args[4]

                indep = newindep

            kshape = newshape

        if self.model.frac is not None:
            info('PSF frac: %s' % self.model.frac)

        if numpy.isscalar(kshape):
            kshape = [kshape]

        return (indep, dep, kshape, lo, hi)

    def get_kernel(self, data, subkernel=True):
        """Return a data object representing the kernel.

        Parameters
        ----------
        data : sherpa.data.Data or sherpa.models.model.Model instance
            The data to apply the kernel to. This routine will pass
            `data` to the `fold` method.
        subkernel : bool, optional

        Returns
        -------
        data : sherpa.data.Data1D or sherpa.data.Data2D instance

        """

        kdata = self._get_kernel_data(data, subkernel)
        indep = kdata[0]
        dep = kdata[1]
        kshape = kdata[2]

        # ndim should be the same as self.ndim
        ndim = len(kshape)

        # TODO: what happens with integrated datasets? This currently
        # returns the low edge of each bin. Should it use the center of
        # the bin or use the full ranges? Is this even an issue?
        #
        if ndim == 1:
            return Data1D('kernel', indep[0], dep)

        if ndim == 2:
            # Note that the shape order is reversed.
            return Data2D('kernel', indep[0], indep[1], dep, kshape[::-1])

        raise PSFErr('ndim')

    def _create_spaces(self, data):
        """Setup the data space based on the pixel size."""

        # This has been pulled out of fold so is currently lacking in documentation.
        #
        # To support using any callable, not just a model, we need
        # to allow the "kernel dimensionality" to be unknown. The
        # alternative is to require the kernel to nave a ndim
        # attribute, but there are a number of places the code
        # allows the kernel to not be a model instance (e.g. the
        # checks for pars/thawedpars).
        #
        if self.ndim is not None and hasattr(data, "ndim") and data.ndim != self.ndim:
            raise PSFErr("mismatch_dims", self.kernel.name, data.name,
                         self.ndim, data.ndim)

        pixel_size_comparison = self._check_pixel_size(data)

        indep = data.get_indep()

        # Don't do anything special
        if pixel_size_comparison == self.SAME_RESOLUTION:
            if self.ndim == 1:
                self.data_space = EvaluationSpace1D(*indep)
            elif self.ndim == 2:
                self.data_space = EvaluationSpace2D(*indep)
            else:
                # leave in in case we support higher dimensions or one
                # of the rare dimensionless models is in use.
                raise PSFErr("ndim")

            self._must_rebin = False
            return (indep, data.get_dims())

        # Evaluate model in PSF space. Note that if we get here then
        # we have to have 2D data.
        #
        if pixel_size_comparison == self.BETTER_RESOLUTION:
            self.data_space = EvaluationSpace2D(*indep)
            self.psf_space = PSFSpace2D(self.data_space, self, data.sky.cdelt)
            self._must_rebin = True
            return (self.psf_space.grid, self.psf_space.shape)

        # PSF has worse resolution, error out
        raise AttributeError("The PSF has a worse resolution than the data.")

    def _check_pixel_size(self, data):
        """
        If the data and PSF dot not have WCS information, assume the PSF has the same resolution
        as the image.
        Otherwise check the WCS information to determine the relative resolutions.

        We only check the resolution in one dimension and assume they are the same.
        """

        if not hasattr(self.kernel, "sky"):
            return self.SAME_RESOLUTION

        # This corresponds to the case when the kernel is actually a psf image, not just a model.
        try:
            psf_pixel_size = self.kernel.sky.cdelt
        except AttributeError:
            # If the kernel does not have a pixel size, issue a warning and keep going
            warnings.warn("PSF Image does not have a pixel size. Sherpa will assume "
                          "the pixel size is the same as the data")
            return self.SAME_RESOLUTION

        try:
            data_pixel_size = data.sky.cdelt
        except AttributeError:
            warnings.warn("Data Image does not have a pixel size. Sherpa will assume "
                          "the pixel size is the same as the PSF")
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
        super().__init__(x, y)
