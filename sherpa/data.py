#
#  Copyright (C) 2008, 2015, 2016, 2017, 2019, 2020
#     Smithsonian Astrophysical Observatory
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
import warnings
from abc import ABCMeta

import numpy

from sherpa.models.regrid import EvaluationSpace1D, EvaluationSpace2D
from sherpa.utils.err import DataErr
from sherpa.utils import SherpaFloat, NoNewAttributesAfterInit, \
    print_fields, create_expr, calc_total_error, bool_cast, \
    filter_bins, parallel_map_funcs


__all__ = ('Data', 'DataSimulFit', 'Data1D', 'Data1DInt',
           'Data1DAsymmetricErrs', 'Data2D', 'Data2DInt')


def _check(array):
    if array is None:
        # There may be valid reasons for the array to be None, e.g. that's what we do in fake_pha
        return array

    if hasattr(array, "shape"):
        if len(array.shape) != 1:
            raise TypeError("Data arrays should be 1-dimensional. Did you call 'flatten()' on {})".format(array))
    else:
        warnings.warn("Converting array {} to numpy array".format(array))
        array = numpy.asanyarray(array)
        return _check(array)
    return array


def _check_indep(array):
    if hasattr(array, 'mask'):
        warnings.warn('Dropping mask attribute for array {}. Masks are supported for dependent variables only.'.format(array))
    return array

def _check_dep(array):
    if not hasattr(array, 'mask'):
        return _check(array), True
    else:
        # We know the mask convention is opposite to sherpa
        if isinstance(array, numpy.ma.MaskedArray):
            return _check(array), ~array.mask
        # We don't know what the mask convention is
        else:
            warnings.warn('Dropping mask for array {}. Set .mask attribute manually or use "set_filter" function.'.format(array))
            return _check(array), True


class DataSpace1D(EvaluationSpace1D):
    """
    Class for representing 1-D Data Space. Data Spaces are spaces that describe the data domain. As models can be
    evaluated over data spaces, data spaces can be considered evaluation spaces themselves. However this "is-a"
    relationship is in the code mostly for convenience and could be removed in future versions.
    """
    def __init__(self, filter, x):
        """
        Parameters
        ----------
        filter : Filter
            a filter object that initialized this data space
        x : array_like
            the x axis of this data space
        """
        self.filter = filter
        EvaluationSpace1D.__init__(self, _check_indep(x))

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        DataSpace1D
        """
        filter = bool_cast(filter)

        if not filter:
            return self

        data = self.grid[0]

        data = self.filter.apply(data)
        return DataSpace1D(self.filter, data)

    def for_model(self, model):
        """
        Models can be defined over arbitrary evaluation spaces. However, 
        at evaluation time during a fit, the model's evaluation space shall
        be done at the user's request space only and set to 0 every where else.

        Parameters
        ----------
        model : The model whose evaluation space needs to be joined with the dataset's data space.

        Returns
        -------
        DataSpace1D
            A data space that joins this data space with the model's evaluation space. if the model does not have an
            evaluation space assigned to itself then `self` is returned.
        """
        evaluation_space = None

        if model is not None and hasattr(model, "evaluation_space"):
            if self not in model.evaluation_space:
                evaluation_space = self

        return self if evaluation_space is None else evaluation_space


class IntegratedDataSpace1D(EvaluationSpace1D):
    """
    Same as DataSpace1D, but for supporting integrated data sets.
    """
    def __init__(self, filter, xlo, xhi):
        """
        Parameters
        ----------
        filter : Filter
            a filter object that initialized this data space
        xlo : array_like
            the lower bounds array of this data space
        xhi : array_like
            the higher bounds array of this data space
        """
        self.filter = filter
        EvaluationSpace1D.__init__(self, _check_indep(xlo), _check_indep(xhi))

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        IntegratedDataSpace1D
        """
        filter = bool_cast(filter)

        if not filter:
            return self

        data = self.grid

        data = tuple(self.filter.apply(axis) for axis in data)
        return IntegratedDataSpace1D(self.filter, *data)

    def for_model(self, model):
        """
        Models can be defined over arbitrary evaluation spaces. However, at evaluation time during a fit, the model's
        evaluation space and the data space will be joined together and the model will be evaluated over the joined
        domain. This makes sure that when the models are rebinned back to the data space the evaluation does not have
        to be extrapolated from the model's evaluation space alone.

        Parameters
        ----------
        model : The model whose evaluation space needs to be joined with the dataset's data space.

        Returns
        -------
        IntegratedDataSpace1D
            A data space that joins this data space with the model's evaluation space. if the model does not have an
            evaluation space assigned to itself then `self` is returned.
        """
        evaluation_space = None

        if model is not None and hasattr(model, "evaluation_space"):
            if self not in model.evaluation_space:
                evaluation_space = self

        return self if evaluation_space is None else evaluation_space


class DataSpace2D():
    """
    Class for representing 2-D Data Spaces. Data Spaces are spaces that describe the data domain.
    """
    def __init__(self, filter, x0, x1):
        """
        Parameters
        ----------
        filter : Filter
            a filter object that initialized this data space
        x0 : array_like
            the first axis of this data space
        x1 : array_like
            the second axis of this data space
        """
        self.filter = filter
        self.x0 = _check(_check_indep(x0))
        self.x1 = _check(_check_indep(x1))

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        DataSpace2D
        """
        filter = bool_cast(filter)

        if not filter:
            return self

        data = self.grid

        data = tuple(self.filter.apply(axis) for axis in data)
        return DataSpace2D(self.filter, *data)

    @property
    def grid(self):
        """
        Return the grid representation of this dataset.

        The x0 and x1 arrays in the grid are one-dimensional representations of the meshgrid obtained
        from the x and y axis arrays, as in `numpy.meshgrid(x, y)[0].ravel()`

        Returns
        -------
        tuple
            A tuple representing the x0 and x1 axes. The tuple will contain two arrays.
        """
        return self.x0, self.x1


class IntegratedDataSpace2D():
    """
    Same as DataSpace2D, but for supporting integrated data sets.
    """
    def __init__(self, filter, x0lo, x1lo, x0hi, x1hi):
        """
        Parameters
        ----------
        filter : Filter
            a filter object that initialized this data space
        x0lo : array_like
            the lower bounds array of the x0 axis
        x0hi : array_like
            the higher bounds array of the xhi axis
        x1lo : array_like
            the lower bounds array of the x0 axis
        x1hi : array_like
            the higher bounds array of the xhi axis
        """
        self.filter = filter
        self.x0lo = _check(_check_indep(x0lo))
        self.x1lo = _check(_check_indep(x1lo))
        self.x0hi = _check(_check_indep(x0hi))
        self.x1hi = _check(_check_indep(x1hi))

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        IntegratedDataSpace2D
        """
        filter = bool_cast(filter)

        if not filter:
            return self

        data = self.grid

        data = tuple(self.filter.apply(axis) for axis in data)
        return IntegratedDataSpace2D(self.filter, *data)

    @property
    def grid(self):
        """
        Return the grid representation of this dataset.

        The x0 and x1 arrays in the grid are one-dimensional representations of the meshgrid obtained
        from the x and y axis arrays, as in `numpy.meshgrid(x, y)[0].ravel()`

        Returns
        -------
        tuple
            A tuple representing the x and y axes. The tuple will contain four arrays.
        """
        return self.x0lo, self.x1lo, self.x0hi, self.x1hi


class DataSpaceND():
    """
    Class for representing arbitray N-Dimensional data domains
    """
    def __init__(self, filter, indep):
        """
        Parameters
        ----------
        filter : Filter
            a filter object that initialized this data space
        indep : tuple of array_like
            the tuple of independent axes.
        """
        self.filter = filter
        self.indep = _check_indep(indep)

    def get(self, filter=False):
        """
        Get a filtered representation of this data set. If `filter` is `False` this object is returned.

        Parameters
        ----------
        filter : bool
            whether the data set should be filtered before being returned

        Returns
        -------
        DataSpaceND
        """
        filter = bool_cast(filter)

        if not filter:
            return self

        data = tuple(self.filter.apply(axis) for axis in self.indep)
        return DataSpaceND(self.filter, data)

    @property
    def grid(self):
        """
        Return the grid representation of this dataset.

        The independent arrays are returned unchanged, i.e. unlike the DataSpace2D class they are not meshed

        Returns
        -------
        tuple
            A tuple representing the independent axes.
        """
        return self.indep


class Filter():
    """
    A class for representing filters of N-Dimentional datasets.
    """
    def __init__(self):
        self._mask = True

    @property
    def mask(self):
        """
        Mask array for dependent variable

        Returns
        -------
        mask : bool or numpy.ndarray
        """
        return self._mask

    @mask.setter
    def mask(self, val):
        if (val is True) or (val is False):
            self._mask = val
        # if val is of type np.bool_ and True, it was failed the previous test because
        # "is True" compares with Python "True" singelton.
        # Yet, we do not want to allow arbitrary values that evaluate as True.
        elif (val is numpy.ma.nomask) or (numpy.isscalar(val) and val and isinstance(val, numpy.bool_)):
            self._mask = True
        elif (val is None) or numpy.isscalar(val):
            raise DataErr('ismask')
        else:
            self._mask = numpy.asarray(val, numpy.bool_)

    def apply(self, array):
        """
        Apply this filter to an array

        Parameters
        ----------
        array : array_like
            Array to be filtered

        Returns
        -------
        array_like : filtered array
        """
        if array is None:
            return

        if self.mask is False:
            raise DataErr('notmask')

        if self.mask is not True:  # mask is not False and not True, so it's something else we'll try to use as an array
            array = numpy.asarray(array)
            if array.shape != self.mask.shape:
                raise DataErr('mismatch', 'mask', 'data array')
            return array[self.mask]

        return array

    def notice(self, mins, maxes, axislist, ignore=False):
        ignore = bool_cast(ignore)
        if str in [type(min) for min in mins]:
            raise DataErr('typecheck', 'lower bound')
        elif str in [type(max) for max in maxes]:
            raise DataErr('typecheck', 'upper bound')
        elif str in [type(axis) for axis in axislist]:
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


class BaseData(metaclass=ABCMeta):
    """
    Base class for all data classes. Left for compatibility with older versions.
    """
    pass


# DATA-NOTE: ND Data cannot be plotted
class Data(NoNewAttributesAfterInit, BaseData):
    """
    Data class for generic, N-Dimensional data sets, where N depends on the number of independent axes passed during
    initialization.

    A data class is the collection of a data space and a number of data array for the dependent variable and
    associated errors.

    This class can be extended by classes definining data sets of specific dimensionality. Extending classes should
    override the `_init_data_space` method.

    This classe provides most of the infrastructure for extending classes for free.
    """
    _fields = ("name", "indep", "dep", "staterror", "syserror")

    def __init__(self, name, indep, y, staterror=None, syserror=None):
        """
        Parameters
        ----------
        name : basestring
            name of this dataset
        indep: tuple of array_like
            the tuple of independent arrays. If this is an numpy masked array,
            the mask will used.
        y : array_like
            the values of the dependent observable
        staterror : array_like
            the statistical error associated with the data
        syserror : array_like
            the systematic error associated with the data
        """
        self.name = name
        self._data_space = self._init_data_space(Filter(), *indep)
        self.y, self.mask = _check_dep(y)
        self.staterror = staterror
        self.syserror = syserror
        NoNewAttributesAfterInit.__init__(self)

    def _init_data_space(self, filter, *data):
        """
        Extending classes should implement this method to provide the proper data space construction.

        Parameters
        ----------
        filter : Filter
            a filter object passed by the initializer upon initialization of extending classes.
        data : tuple of array_like
            the tuple of independent arrays used to build the data space.

        Returns
        -------
        object
            an instance of the dataspace associated with this data set.

        """
        return DataSpaceND(filter, data)

    @property
    def dep(self):
        """
        Left for compatibility with older versions
        """
        return self.y

    @dep.setter
    def dep(self, val):
        self.y = val

    @property
    def mask(self):
        """
        Mask array for dependent variable

        Returns
        -------
        mask : bool or numpy.ndarray
        """
        return self._data_space.filter.mask

    @mask.setter
    def mask(self, val):
        self._data_space.filter.mask = val

    def get_dims(self):
        """
        Return the dimensions of this data space as a tuple of tuples.
        The first element in the tuple is a tuple with the dimensions of the data space, while the second element
        provides the size of the dependent array.
        Returns
        -------
        tuple
        """
        indep_size = tuple(indep.size for indep in self.indep)
        return indep_size, self.dep.size

    @property
    def indep(self):
        """
        Return the grid of the data space associated with this data set.
        Returns
        -------
        tuple of array_like
        """
        return self._data_space.get().grid

    @indep.setter
    def indep(self, val):
        self._data_space = self._init_data_space(self._data_space.filter, *val)

    def get_indep(self, filter=False):
        """Return the independent axes of a data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis: tuple of arrays
           The independent axis values for the data set. This gives
           the coordinates of each point in the data set.

        See Also
        --------
        get_dep : Return the dependent axis of a data set.

        """
        data_space = self._data_space.get(filter)
        return data_space.grid

    def set_indep(self, val):
        self.indep = val

    def get_dep(self, filter=False):
        """Return the dependent axis of a data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis: array
           The dependent axis values for the data set. This gives
           the value of each point in the data set.

        See Also
        --------
        get_indep : Return the independent axis of a data set.
        get_error : Return the errors on the dependent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.

        """
        dep = self.dep
        filter = bool_cast(filter)
        if filter:
            dep = self.apply_filter(dep)
        return dep

    def set_dep(self, val):
        """
        Set the dependent variable values"

        Parameters
        ----------
        val

        Returns
        -------

        """
        if numpy.iterable(val):
            dep = numpy.asarray(val, SherpaFloat)
        else:
            val = SherpaFloat(val)
            dep = numpy.array([val] * len(self.get_indep()[0]))
        self.y = dep

    def get_y(self, filter=False, yfunc=None, use_evaluation_space=False):
        """
        Return dependent axis in N-D view of dependent variable

        Parameters
        ----------
        filter
        yfunc
        use_evaluation_space

        Returns
        -------

        """
        y = self.get_dep(filter)

        if yfunc is not None:
            if filter:
                yfunc = self.eval_model_to_fit(yfunc)
            else:
                yfunc = self.eval_model(yfunc)
            y = (y, yfunc)

        return y

    def get_staterror(self, filter=False, staterrfunc=None):
        """Return the statistical error on the dependent axis of a data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        staterrfunc : function
           If no statistical error has been set, the errors will
           be calculated by applying this function to the
           dependent axis of the data set.

        Returns
        -------
        axis : array or `None`
           The statistical error for each data point. A value of
           `None` is returned if the data set has no statistical error
           array and `staterrfunc` is `None`.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.

        """
        staterror = getattr(self, 'staterror', None)
        filter = bool_cast(filter)
        if filter:
            staterror = self.apply_filter(staterror)

        if (staterror is None) and (staterrfunc is not None):
            dep = self.get_dep()
            if filter:
                dep = self.apply_filter(dep)
            staterror = staterrfunc(dep)
        return staterror

    def get_syserror(self, filter=False):
        """Return the statistical error on the dependent axis of a data set.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.

        Returns
        -------
        axis : array or `None`
           The systematic error for each data point. A value of
           `None` is returned if the data set has no systematic
           errors.

        See Also
        --------
        get_error : Return the errors on the dependent axis of a data set.
        get_indep : Return the independent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.

        """
        syserr = getattr(self, 'syserror', None)
        filter = bool_cast(filter)
        if filter:
            syserr = self.apply_filter(syserr)
        return syserr

    def get_error(self, filter=False, staterrfunc=None):
        """Return the total error on the dependent variable.

        Parameters
        ----------
        filter : bool, optional
           Should the filter attached to the data set be applied to
           the return value or not. The default is `False`.
        staterrfunc : function
           If no statistical error has been set, the errors will
           be calculated by applying this function to the
           dependent axis of the data set.

        Returns
        -------
        axis : array or `None`
           The error for each data point, formed by adding the
           statistical and systematic errors in quadrature.

        See Also
        --------
        get_dep : Return the independent axis of a data set.
        get_staterror : Return the statistical errors on the dependent axis of a data set.
        get_syserror : Return the systematic errors on the dependent axis of a data set.

        """
        return calc_total_error(self.get_staterror(filter, staterrfunc),
                                self.get_syserror(filter))

    def get_yerr(self, filter=False, staterrfunc=None):
        """
        Return errors in dependent axis in N-D view of dependent variable

        Parameters
        ----------
        filter
        staterrfunc

        Returns
        -------

        """
        return self.get_error(filter, staterrfunc)

    def get_ylabel(self, yfunc=None):
        """
        Return label for dependent axis in N-D view of dependent variable"

        Parameters
        ----------
        yfunc

        Returns
        -------

        """
        return 'y'

    def apply_filter(self, data):
        return self._data_space.filter.apply(data)

    def notice(self, mins, maxes, ignore=False):
        self._data_space.filter.notice(mins, maxes, self.get_indep(), ignore)

    def ignore(self, *args, **kwargs):
        kwargs['ignore'] = True
        self.notice(*args, **kwargs)

    def eval_model(self, modelfunc):
        return modelfunc(*self.get_indep())

    def eval_model_to_fit(self, modelfunc):
        return modelfunc(*self.get_indep(filter=True))

    def to_guess(self):
        arrays = [self.get_y(True)]
        arrays.extend(self.get_indep(True))
        return tuple(arrays)

    def to_fit(self, staterrfunc=None):
        return (self.get_dep(True),
                self.get_staterror(True, staterrfunc),
                self.get_syserror(True))

    def __str__(self):
        """
        Return a listing of the attributes listed in self._fields and,
        if present, self._extra_fields.
        """

        fields = self._fields + getattr(self, '_extra_fields', ())
        fdict = dict(zip(fields, [getattr(self, f) for f in fields]))
        return print_fields(fields, fdict)

    def __repr__(self):
        r = '<%s data set instance' % type(self).__name__
        if hasattr(self, 'name'):
            r += " '%s'" % self.name
        r += '>'
        return r


class DataSimulFit(NoNewAttributesAfterInit):
    """Store multiple data sets.

    This class lets multiple data sets be treated as a single
    dataset by concatenation. That is, if two data sets have lengths
    n1 and n2 then they can be considered as a single data set of
    length n1 + n2.

    Parameters
    ----------
    name : str
        The name for the collection of data.
    datasets : sequence of Data objects
        The datasets to be stored; there must be at least one. They are
        assumed to behave as sherpa.data.Data objects, but there is no
        check for this condition.

    Attributes
    ----------
    datasets : sequence of Data

    See Also
    --------
    sherpa.models.model.SimulFitModel

    Examples
    --------

    >>> d1 = Data1D('d1', [1, 2, 3], [10, 12, 15])
    >>> d2 = Data1D('d2', [1, 2, 5, 7], [4, 15, 9, 24])
    >>> dall = DataSimulFit('comb', (d1, d2))
    >>> yvals, _, _ = dall.to_fit()
    >>> print(yvals)
    [10 12 15  4 15  9 24]

    """

    def __init__(self, name, datasets, numcores=1):
        if len(datasets) == 0:
            raise DataErr('zerodatasimulfit', type(self).__name__)
        self.name = name
        self.datasets = tuple(datasets)
        self.numcores = numcores
        NoNewAttributesAfterInit.__init__(self)

    def eval_model_to_fit(self, modelfuncs):
        if self.numcores == 1:
            total_model = []
            for func, data in zip(modelfuncs, self.datasets):
                tmp_model = data.eval_model_to_fit(func)
                total_model.append(tmp_model)
            return numpy.concatenate(total_model)
        else:
        # best to make this a different derived class
            funcs = []
            datasets = []
            for func, data in zip(modelfuncs, self.datasets):
                funcs.append(func)
                datasets.append(data.get_indep(filter=False))
            total_model = parallel_map_funcs(funcs, datasets, self.numcores)
            all_model = []
            for model, data in zip(total_model, self.datasets):
                all_model.append(data.apply_filter(model))
            return numpy.concatenate(all_model)

    def to_fit(self, staterrfunc=None):
        total_dep = []
        total_staterror = []
        total_syserror = []

        no_staterror = True
        no_syserror = True

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
        elif numpy.any([numpy.equal(array, None).any()
                        for array in total_staterror]):
            raise DataErr('staterrsimulfit')
        else:
            total_staterror = numpy.concatenate(total_staterror)

        if no_syserror:
            total_syserror = None
        else:
            total_syserror = numpy.concatenate(total_syserror)

        return total_dep, total_staterror, total_syserror

    # DATA-NOTE: this implementation is weird. Is this even used?
    def to_plot(self, yfunc=None, staterrfunc=None):
        return self.datasets[0].to_plot(yfunc.parts[0], staterrfunc)


class Data1D(Data):
    _fields = ("name", "x", "y", "staterror", "syserror")

    def __init__(self, name, x, y, staterror=None, syserror=None):
        Data.__init__(self, name, (x, ), y, staterror, syserror)

    def _init_data_space(self, filter, *data):
        return DataSpace1D(filter, *data)

    def get_x(self, filter=False, model=None, use_evaluation_space=False):
        return self.get_evaluation_indep(filter, model, use_evaluation_space)[0]

    def get_xerr(self, filter=False, yfunc=None):
        """
        Return linear view of bin size in independent axis/axes"

        Parameters
        ----------
        filter
        yfunc

        Returns
        -------

        """
        return None

    def get_xlabel(self):
        """
        Return label for linear view of independent axis/axes
        Returns
        -------

        """
        return 'x'

    def get_dims(self, filter=False):
        return len(self.get_x(filter)),

    def get_y(self, filter=False, yfunc=None, use_evaluation_space=False):
        """
        Return dependent axis in N-D view of dependent variable"

        Parameters
        ----------
        filter
        yfunc
        use_evaluation_space

        Returns
        -------

        """
        y = self.get_dep(filter)

        if yfunc is not None:
            model_evaluation = yfunc(*self.get_evaluation_indep(filter, yfunc, use_evaluation_space))
            y = y, model_evaluation

        return y

    def get_bounding_mask(self):
        mask = self.mask
        size = None
        if numpy.iterable(self.mask):
            # create bounding box around noticed image regions
            mask = numpy.array(self.mask)
            size = (mask.size,)
        return mask, size

    def get_img(self, yfunc=None):
        """
        Return 1D dependent variable as a 1 x N image

        Parameters
        ----------
        yfunc

        Returns
        -------

        """
        y_img = self.get_y(False, yfunc)
        if yfunc is not None:
            y_img = (y_img[0].reshape(1, y_img[0].size),
                     y_img[1].reshape(1, y_img[1].size))
        else:
            y_img = y_img.reshape(1, y_img.size)
        return y_img

    def get_imgerr(self):
        err = self.get_error()
        if err is not None:
            err = err.reshape(1, err.size)
        return err

    def get_filter(self, format='%.4f', delim=':'):
        # for derived intergrated classes, this will return values in center of
        # bin.
        x = self.get_x(filter=True)
        mask = numpy.ones(len(x), dtype=bool)
        if numpy.iterable(self.mask):
            mask = self.mask
        return create_expr(x, mask, format, delim)

    def get_filter_expr(self):
        return self.get_filter(delim='-') + ' ' + self.get_xlabel()

    def to_plot(self, yfunc=None, staterrfunc=None):
        # As we introduced models defined on arbitrary grids, the x array can also depend on the
        # model function, at least in principle.
        return (self.get_x(True, yfunc),
                self.get_y(True, yfunc),
                self.get_yerr(True, staterrfunc),
                self.get_xerr(True, yfunc),
                self.get_xlabel(),
                self.get_ylabel())

    def to_component_plot(self, yfunc=None, staterrfunc=None):
        # As we introduced models defined on arbitrary grids, the x array can also depend on the
        # model function, at least in principle.
        return (self.get_x(True, yfunc, use_evaluation_space=True),
                self.get_y(True, yfunc, use_evaluation_space=True),
                self.get_yerr(True, staterrfunc),
                self.get_xerr(True, yfunc),
                self.get_xlabel(),
                self.get_ylabel())

    def get_evaluation_indep(self, filter=False, model=None, use_evaluation_space=False):
        data_space = self._data_space.get(filter)
        if use_evaluation_space:
            return data_space.for_model(model).grid
        else:
            return data_space.grid

    def notice(self, xlo=None, xhi=None, ignore=False):
        Data.notice(self, (xlo,), (xhi,), ignore)

    @property
    def x(self):
        """
        Used for compatibility, in particular for __str__ and __repr__
        """
        return self.get_x()


class Data1DAsymmetricErrs(Data1D):
    """1-D data set with asymmetric errors
    Note: elo and ehi shall be stored as delta values from y"""

    _fields = ("name", "x", "y", "staterror", "syserror", "elo", "ehi")

    def __init__(self, name, x, y, elo, ehi, staterror=None, syserror=None):
        self.elo = elo
        self.ehi = ehi
        Data1D.__init__(self, name, x, y, staterror=staterror, syserror=syserror)

    def get_yerr(self, filter=False, staterrfunc=None):
        return self.elo, self.ehi


class Data1DInt(Data1D):
    """
    1-D integrated data set
    """
    _fields = ("name", "xlo", "xhi", "y", "staterror", "syserror")

    def __init__(self, name, xlo, xhi, y, staterror=None, syserror=None):
        Data.__init__(self, name, (xlo, xhi), y, staterror, syserror)

    def _init_data_space(self, filter, *data):
        return IntegratedDataSpace1D(filter, *data)

    def get_x(self, filter=False, model=None, use_evaluation_space=False):
        indep = self.get_evaluation_indep(filter, model, use_evaluation_space)
        return (indep[0] + indep[1]) / 2.0

    def get_xerr(self, filter=False, model=None):
        xlo, xhi = self.get_evaluation_indep(filter, model)
        return xhi - xlo

    def notice(self, xlo=None, xhi=None, ignore=False):
        Data.notice(self, (None, xlo), (xhi, None), ignore)

    @property
    def xlo(self):
        """
        Property kept for compatibility
        """
        return self._data_space.x_axis.lo

    @property
    def xhi(self):
        """
        Property kept for compatibility
        """
        return self._data_space.x_axis.hi


class Data2D(Data):
    _fields = ("name", "x0", "x1", "y", "shape", "staterror", "syserror")

    def __init__(self, name, x0, x1, y, shape=None, staterror=None, syserror=None):
        self.shape = shape
        Data.__init__(self, name, (x0, x1), y, staterror, syserror)

    def _init_data_space(self, filter, *data):
        return DataSpace2D(filter, *data)

    def get_x0(self, filter=False):
        return self._data_space.get(filter).x0

    def get_x1(self, filter=False):
        return self._data_space.get(filter).x1

    def get_x0label(self):
        """
        Return label for first dimension in 2-D view of independent axis/axes

        Returns
        -------

        """
        return 'x0'

    def get_x1label(self):
        """
        Return label for second dimension in 2-D view of independent axis/axes
        """
        return 'x1'

    def get_axes(self):
        self._check_shape()
        # FIXME: how to filter an axis when self.mask is size of self.y?
        return (numpy.arange(self.shape[1]) + 1,
                numpy.arange(self.shape[0]) + 1)

    def get_dims(self, filter=False):
        # self._check_shape()
        if self.shape is not None:
            return self.shape[::-1]
        return len(self.get_x0(filter)), len(self.get_x1(filter))

    def get_filter_expr(self):
        return ''

    def get_filter(self):
        return ''

    def get_max_pos(self, dep=None):
        if dep is None:
            dep = self.get_dep(True)
        x0 = self.get_x0(True)
        x1 = self.get_x1(True)

        pos = numpy.asarray(numpy.where(dep == dep.max())).squeeze()
        if pos.ndim == 0:  # DATA-NOTE: Could this ever be False?!
            pos = int(pos)
            return x0[pos], x1[pos]

        return [(x0[index], x1[index]) for index in pos]

    # For images, only need y-array
    # Also, we do not filter, as imager needs M x N (or
    # L x M x N) array
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

    def to_contour(self, yfunc=None):
        return (self.get_x0(True),
                self.get_x1(True),
                self.get_y(True, yfunc),
                self.get_x0label(),
                self.get_x1label())

    def _check_shape(self):
        if self.shape is None:
            raise DataErr('shape', self.name)

    @property
    def x0(self):
        """
        kept for compatibility
        """
        return self.get_x0()

    @property
    def x1(self):
        """
        kept for compatibility
        """
        return self.get_x1()

    def notice(self, x0lo=None, x0hi=None, x1lo=None, x1hi=None, ignore=False):
        Data.notice(self, (x0lo, x1lo), (x0hi, x1hi), ignore)


class Data2DInt(Data2D):
    """
    2-D integrated data set
    """
    _fields = ("name", "x0lo", "x1lo", "x0hi", "x1hi", "y", "shape", "staterror", "syserror")

    def __init__(self, name, x0lo, x1lo, x0hi, x1hi, y, shape=None, staterror=None, syserror=None):
        self.shape = shape
        Data.__init__(self, name, (x0lo, x1lo, x0hi, x1hi), y, staterror, syserror)

    def _init_data_space(self, filter, *data):
        return IntegratedDataSpace2D(filter, *data)

    def get_x0(self, filter=False):
        indep = self._data_space.get(filter)
        return (indep.x0lo + indep.x0hi) / 2.0

    def get_x1(self, filter=False):
        indep = self._data_space.get(filter)
        return (indep.x1lo + indep.x1hi) / 2.0

    def notice(self, x0lo=None, x0hi=None, x1lo=None, x1hi=None, ignore=False):
        Data.notice(self, (None, None, x0lo, x1lo), (x0hi, x1hi, None, None), ignore)
