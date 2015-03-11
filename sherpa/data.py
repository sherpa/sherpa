# 
#  Copyright (C) 2008  Smithsonian Astrophysical Observatory
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

"""
Tools for creating, storing, inspecting, and manipulating data sets
"""


import sys
import inspect
from itertools import izip
import numpy
from sherpa.utils.err import DataErr, NotImplementedErr
from sherpa.utils import SherpaFloat, NoNewAttributesAfterInit, \
     print_fields, create_expr, calc_total_error, bool_cast, \
     filter_bins


_all__ = ('Data', 'DataSimulFit', 'Data1D', 'Data1DInt', 'Data2D', 'Data2DInt')


class BaseData(NoNewAttributesAfterInit):
    "Base class for all data set types"

    def _get_filter(self):
        return self._filter
    def _set_filter(self, val):
        self._filter = val
        self._mask = True
    filter = property(_get_filter, _set_filter,
                      doc='Filter for dependent variable')

    def _get_mask(self):
        return self._mask
    def _set_mask(self, val):
        if (val is True) or (val is False):
            self._mask = val
        elif (val is None) or numpy.isscalar(val):
            raise DataErr('ismask')
        else:
            self._mask = numpy.asarray(val, numpy.bool_)
        self._filter = None
    mask = property(_get_mask, _set_mask,
                    doc='Mask array for dependent variable')

    def __init__(self):
        """

        Initialize a data object.  This method can only be called from
        a derived class constructor.  Attempts to create a BaseData
        instance will raise NotImplementedErr.

        Derived class constructors must call this method directly (and
        not indirectly through a superclass constructor).  When thus
        invoked, this method will extract the argument names and
        values from the derived class constructor invocation and set
        corresponding attributes on the instance (thereby eliminating
        the need for the derived class constructor to do its own
        attribute setting).  If the name of an argument matches the
        name of a DataProperty of the derived class, then the
        corresponding attribute name will have an underscore prepended
        (meaning the property will use the value directly instead of
        relying on _get_*/_set_* methods).

        """

        if type(self) is BaseData:
            raise NotImplementedErr('noinstanceallowed', 'BaseData')

        frame = sys._getframe().f_back
        cond = (frame.f_code is self.__init__.im_func.func_code)
        assert cond, (('%s constructor must call BaseData constructor ' +
                       'directly') % type(self).__name__)
        args = inspect.getargvalues(frame)
        
        self._fields = tuple(args[0][1:])
        for f in self._fields:
            cond = (f not in vars(self))
            assert cond, (("'%s' object already has attribute '%s'") %
                          (type(self).__name__, f))
            setattr(self, f, args[3][f])

        self.filter = None
        self.mask = True

        NoNewAttributesAfterInit.__init__(self)

    def __str__(self):
        """

        Return a listing of the attributes listed in self._fields and,
        if present, self._extra_fields.

        """

        fields = self._fields + getattr(self, '_extra_fields', ())
        fdict = dict(izip(fields, [getattr(self, f) for f in fields]))
        return print_fields(fields, fdict)

    def apply_filter(self, data):
        if data is not None:
            if self.filter is not None:
                if callable(self.filter):
                    data = self.filter(data)
                else:
                    data = data[self.filter]
            elif self.mask is not True:
                if self.mask is False:
                    raise DataErr('notmask')
                data = numpy.asarray(data)
                if data.shape != self.mask.shape:
                    raise DataErr('mismatch', 'mask', 'data array')
                data = data[self.mask]
        return data

    def ignore(self, *args, **kwargs):
        kwargs['ignore'] = True
        self.notice(*args, **kwargs)

    def notice(self, mins, maxes, axislist, ignore=False):

        ignore = bool_cast(ignore)
        if( str in [type(min) for min in mins] ):
            raise DataErr('typecheck', 'lower bound')
        elif( str in [type(max) for max in maxes] ):
            raise DataErr('typecheck', 'upper bound')
        elif( str in [type(axis) for axis in axislist] ):
            raise DataErr('typecheck', 'grid')

        mask = filter_bins(mins, maxes, axislist)

        if mask is None:
            self.mask = not ignore
        elif not ignore:
            if self.mask is True:
                self.mask = mask
            else:
                self.mask |= mask
        else:
            mask = ~mask
            if self.mask is False:
                self.mask = mask
            else:
                self.mask &= mask


class Data(BaseData):
    "Generic data set"

    def __init__(self, name, indep, dep, staterror=None, syserror=None):
        """

        Initialize a Data instance.  indep should be a tuple of
        independent axis arrays, dep should be an array of dependent
        variable values, and staterror and syserror should be arrays
        of statistical and systematic errors, respectively, in the
        dependent variable (or None).

        """

        BaseData.__init__(self)

    def __repr__(self):
        r = '<%s data set instance' % type(self).__name__
        if hasattr(self, 'name'):
            r += " '%s'" % self.name
        r += '>'
        return r

    def eval_model(self, modelfunc):
        return modelfunc(*self.get_indep())

    def eval_model_to_fit(self, modelfunc):
        return modelfunc(*self.get_indep(filter=True))

    #
    # Primary properties.  These can depend only on normal attributes (and not
    # other properties).
    #

    def get_indep(self, filter=False):
        "Return a tuple containing the independent variables/axes"
        indep = getattr(self, 'indep', None)
        filter=bool_cast(filter)
        if filter:
            indep = tuple([self.apply_filter(x) for x in indep])
        return indep

    def get_dep(self, filter=False):
        "Return an array of dependent variable values"
        dep = getattr(self, 'dep', None)
        filter=bool_cast(filter)
        if filter:
            dep = self.apply_filter(dep)
        return dep

    def get_staterror(self, filter=False, staterrfunc=None):
        "Return the statistical error array"

        staterror = getattr(self, 'staterror', None)
        filter=bool_cast(filter)
        if filter:
            staterror = self.apply_filter(staterror)
        
        if (staterror is None) and (staterrfunc is not None):
            dep = self.get_dep()
            if filter:
                dep = self.apply_filter(dep)
            staterror = staterrfunc(dep)
        return staterror

    def get_syserror(self, filter=False):
        "Return the systematic error array"
        syserr = getattr(self, 'syserror', None)
        filter=bool_cast(filter)
        if filter:
            syserr = self.apply_filter(syserr)
        return syserr

    #
    # Utility methods
    #

    def _wrong_dim_error(self, baddim):
	raise DataErr('wrongdim', self.name, baddim)

    def _no_image_error(self):
	raise DataErr('notimage', self.name)

    def _no_dim_error(self):
        raise DataErr('nodim', self.name)

    #
    # Secondary properties.  To best support subclasses, these should depend
    # only on the primary properties whenever possible, though there may be
    # instances when they depend on normal attributes.
    #

    def get_dims(self):
        self._no_dim_error()

    def get_error(self, filter=False, staterrfunc=None):
        "Return total error in dependent variable"
        return calc_total_error(self.get_staterror(filter, staterrfunc),
                                self.get_syserror(filter))

    def get_x(self, filter=False):
        "Return linear view of independent axis/axes"
	self._wrong_dim_error(1)

    def get_xerr(self, filter=False):
        "Return linear view of bin size in independent axis/axes"
	return None

    def get_xlabel(self):
        "Return label for linear view ofindependent axis/axes"
	return 'x'

    def get_y(self, filter=False, yfunc=None):
        "Return dependent axis in N-D view of dependent variable"
        y = self.get_dep(filter)

        if yfunc is not None:
            if filter:
                yfunc = self.eval_model_to_fit(yfunc)
            else:
                yfunc = self.eval_model(yfunc)
            y = (y, yfunc)

        return y

    def get_yerr(self, filter=False, staterrfunc=None):
        "Return errors in dependent axis in N-D view of dependent variable"
	return self.get_error(filter, staterrfunc) 

    def get_ylabel(self, yfunc=None):
        "Return label for dependent axis in N-D view of dependent variable"
	return 'y'

    def get_x0(self, filter=False):
        "Return first dimension in 2-D view of independent axis/axes"
	self._wrong_dim_error(2)

    def get_x0label(self):
        "Return label for first dimension in 2-D view of independent axis/axes"
	return 'x0'

    def get_x1(self, filter=False):
        "Return second dimension in 2-D view of independent axis/axes"
	self._wrong_dim_error(2)

    def get_x1label(self):
        """

        Return label for second dimension in 2-D view of independent axis/axes

        """
	return 'x1'

    # For images, only need y-array
    # Also, we do not filter, as imager needs M x N (or
    # L x M x N) array
    def get_img(self, yfunc=None):
        "Return dependent variable as an image"
	self._no_image_error()

    def get_imgerr(self, yfunc=None):
        "Return total error in dependent variable as an image"
	self._no_image_error()

    def to_guess(self):
        arrays = [self.get_y(True)]
        arrays.extend(self.get_indep(True))
        return tuple(arrays)

    def to_fit(self, staterrfunc=None):
        return (self.get_dep(True),
                self.get_staterror(True, staterrfunc),
                self.get_syserror(True))

    def to_plot(self, yfunc=None, staterrfunc=None):
        return (self.get_x(True),
                self.get_y(True, yfunc),
                self.get_yerr(True, staterrfunc),
                self.get_xerr(True),
                self.get_xlabel(),
                self.get_ylabel())

    def to_contour(self, yfunc=None):
        return (self.get_x0(True),
                self.get_x1(True),
                self.get_y(True, yfunc),
                self.get_x0label(),
                self.get_x1label())


class DataSimulFit(Data):

    def __init__(self, name, datasets):
        if len(datasets) == 0:
            raise DataErr('zerodatasimulfit', type(self).__name__)
        datasets = tuple(datasets)
        BaseData.__init__(self)

    def eval_model_to_fit(self, modelfuncs):
        total_model = []

        for func, data in izip(modelfuncs, self.datasets):
            total_model.append(data.eval_model_to_fit(func))
        
        return numpy.concatenate(total_model)

    def to_fit(self, staterrfunc=None):        
        total_dep = []
        total_staterror = []
        total_syserror = []

        no_staterror = True
        no_syserror  = True

        for data in self.datasets:
            dep, staterror, syserror = data.to_fit(staterrfunc)

            total_dep.append(dep)

            if staterror is not None:
                no_staterror = False
            total_staterror.append(staterror)

            if syserror is not None:
                no_syserror = False
            else:
                syserror = numpy.zeros_like(dep)
            total_syserror.append(syserror)

        total_dep = numpy.concatenate(total_dep)

        if no_staterror:
            total_staterror = None
        elif None in total_staterror:
            raise DataErr('staterrsimulfit')
        else:
            total_staterror = numpy.concatenate(total_staterror)

        if no_syserror:
            total_syserror = None
        else:
            total_syserror = numpy.concatenate(total_syserror)

        return (total_dep, total_staterror, total_syserror)

    def to_plot(self, yfunc=None, staterrfunc=None):
        return self.datasets[0].to_plot(yfunc.parts[0], staterrfunc)


class DataND(Data):
    "Base class for Data1D, Data2D, etc."

    def get_dep(self, filter=False):
        y = self.y
        filter=bool_cast(filter)
        if filter:
            y = self.apply_filter(y)
        return y

    def set_dep(self, val):
        "Set the dependent variable values"
        dep = None
        if numpy.iterable(val):
            dep = numpy.asarray(val, SherpaFloat)
        else:
            val = SherpaFloat(val)
            dep = numpy.array([val]*len(self.get_indep()[0]))
        setattr(self, 'y', dep)


class Data1D(DataND):
    "1-D data set"

    def _set_mask(self, val):
        DataND._set_mask(self, val)
        try:
            self._x = self.apply_filter(self.x)
        except DataErr:
            self._x = self.x

    mask = property(DataND._get_mask, _set_mask,
                    doc='Mask array for dependent variable')

    def __init__(self, name, x, y, staterror=None, syserror=None):
        self._x = x
	BaseData.__init__(self)

    def get_indep(self, filter=False):
        filter=bool_cast(filter)
        if filter:
            return (self._x,)
        return (self.x,)

    def get_x(self, filter=False):
        return self.get_indep(filter)[0]

    def get_dims(self, filter=False):
        return (len(self.get_x(filter)),)

    def get_filter(self, format='%.4f', delim=':'):
        # for derived intergrated classes, this will return values in center of
        # bin.
        x = self.get_x(filter=True)
        mask = numpy.ones(len(x), dtype=bool)
        if numpy.iterable(self.mask):
            mask = self.mask
        return create_expr(x, mask, format, delim)

    def get_filter_expr(self):
        return (self.get_filter(delim='-') + ' ' + self.get_xlabel())

    def get_bounding_mask(self):
        mask = self.mask
        size = None
        if numpy.iterable(self.mask):
            # create bounding box around noticed image regions
            mask = numpy.array(self.mask)
#            xi = numpy.where(mask == True)[0]
#            xlo = xi.min()
#            xhi = xi.max()
#            size = (mask[xlo:xhi+1].size,)
#            mask = mask[xlo:xhi+1]
            size = (mask.size,)
        return mask, size

    def get_img(self, yfunc=None):
        "Return 1D dependent variable as a 1 x N image"
        y_img = self.get_y(False, yfunc)
        if yfunc is not None:
            y_img = (y_img[0].reshape(1,y_img[0].size),
                     y_img[1].reshape(1,y_img[1].size))
        else:
            y_img = y_img.reshape(1,y_img.size)
        return y_img
    
    def get_imgerr(self):
	err = self.get_error()
	if err is not None:
	    err = err.reshape(1,err.size)
	return err

    def notice(self, xlo=None, xhi=None, ignore=False):
        BaseData.notice(self, (xlo,), (xhi,), self.get_indep(), ignore)


class Data1DInt(Data1D):
    "1-D integrated data set"

    def _set_mask(self, val):
        DataND._set_mask(self, val)
        try:
            self._lo = self.apply_filter(self.xlo)
            self._hi = self.apply_filter(self.xhi)
        except DataErr:
            self._lo = self.xlo
            self._hi = self.xhi

    mask = property(DataND._get_mask, _set_mask,
                    doc='Mask array for dependent variable')

    def __init__(self, name, xlo, xhi, y, staterror=None, syserror=None):
        self._lo = xlo
        self._hi = xhi
	BaseData.__init__(self)

    def get_indep(self, filter=False):
        filter=bool_cast(filter)
        if filter:
            return (self._lo, self._hi)
	return (self.xlo, self.xhi)

    def get_x(self, filter=False):
        indep = self.get_indep(filter)
        return (indep[0] + indep[1]) / 2.0

    def get_xerr(self, filter=False):
        xlo,xhi = self.get_indep(filter)
        return xhi-xlo

    def notice(self, xlo=None, xhi=None, ignore=False):
        BaseData.notice(self, (None, xlo), (xhi, None), self.get_indep(),
                        ignore)


class Data2D(DataND):
    "2-D data set"

    def _set_mask(self, val):
        DataND._set_mask(self, val)
        try:
            self._x0 = self.apply_filter(self.x0)
            self._x1 = self.apply_filter(self.x1)
        except DataErr:
            self._x0 = self.x0
            self._x1 = self.x1

    mask = property(DataND._get_mask, _set_mask,
                    doc='Mask array for dependent variable')

    def __init__(self, name, x0, x1, y, shape=None, staterror=None,
                 syserror=None):
        self._x0 = x0
        self._x1 = x1
	BaseData.__init__(self)

    def get_indep(self, filter=False):
        filter=bool_cast(filter)
        if filter:
            return (self._x0, self._x1)
	return (self.x0, self.x1)

    def get_x0(self, filter=False):
        return self.get_indep(filter)[0]


    def get_x1(self, filter=False):
        return self.get_indep(filter)[1]

    def get_axes(self):
        self._check_shape()
        # FIXME: how to filter an axis when self.mask is size of self.y?
        return (numpy.arange(self.shape[1])+1, numpy.arange(self.shape[0])+1)

    def get_dims(self, filter=False):
        #self._check_shape()
        if self.shape is not None:
            return self.shape[::-1]
        return (len(self.get_x0(filter)), len(self.get_x1(filter)))

    def get_filter_expr(self):
        return ''

    get_filter = get_filter_expr

    def _check_shape(self):
        if self.shape is None:
            raise DataErr('shape',self.name)

    def get_max_pos(self, dep=None):
        if dep is None:
            dep = self.get_dep(True)
        x0 = self.get_x0(True)
        x1 = self.get_x1(True)

        pos = numpy.asarray(numpy.where(dep == dep.max())).squeeze()
        if pos.ndim == 0:
            pos = int(pos)
            return (x0[pos], x1[pos])

        return [(x0[index], x1[index]) for index in pos]

    def get_img(self, yfunc=None):
        self._check_shape()
        y_img = self.get_y(False, yfunc)
        if yfunc is not None:
            y_img = (y_img[0].reshape(*self.shape),
                     y_img[1].reshape(*self.shape))
        else:
            y_img = y_img.reshape(*self.shape)
	return y_img

    def get_imgerr(self):
        self._check_shape()
	err = self.get_error()
	if err is not None:
	    err = err.reshape(*self.shape)
	return err

    def notice(self, x0lo=None, x0hi=None, x1lo=None, x1hi=None, ignore=False):
        BaseData.notice(self, (x0lo, x1lo), (x0hi, x1hi), self.get_indep(),
                        ignore)


class Data2DInt(Data2D):
    "2-D integrated data set"

    def _set_mask(self, val):
        DataND._set_mask(self, val)
        try:
            self._x0lo = self.apply_filter(self.x0lo)
            self._x0hi = self.apply_filter(self.x0hi)
            self._x1lo = self.apply_filter(self.x1lo)
            self._x1hi = self.apply_filter(self.x1hi)
        except DataErr:
            self._x0lo = self.x0lo
            self._x1lo = self.x1lo
            self._x0hi = self.x0hi
            self._x1hi = self.x1hi

    mask = property(DataND._get_mask, _set_mask,
                    doc='Mask array for dependent variable')

    def __init__(self, name, x0lo, x1lo, x0hi, x1hi, y, shape=None,
                 staterror=None, syserror=None):
        self._x0lo = x0lo
        self._x1lo = x1lo
        self._x0hi = x0hi
        self._x1hi = x1hi
	BaseData.__init__(self)

    def get_indep(self, filter=False):
        filter=bool_cast(filter)
        if filter:
            return (self._x0lo, self._x1lo, self._x0hi, self._x1hi)
	return (self.x0lo, self.x1lo, self.x0hi, self.x1hi)

    def get_x0(self, filter=False):
        indep = self.get_indep(filter)
        return (indep[0] + indep[2]) / 2.0

    def get_x1(self, filter=False):
        indep = self.get_indep(filter)
        return (indep[1] + indep[3]) / 2.0

    def notice(self, x0lo=None, x0hi=None, x1lo=None, x1hi=None, ignore=False):
        BaseData.notice(self, (None, None, x0lo, x1lo),
                        (x0hi, x1hi, None, None), self.get_indep(), ignore)
