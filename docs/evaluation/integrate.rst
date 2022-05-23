*************************
Binned and Unbinned grids
*************************

Sherpa supports models for both
:ref:`unbinned <data_unbinned>` and
:ref:`binned <data_binned>` data
sets. The output of a model depends on how it is called
(is it sent just the grid points or the bin edges), how
the :py:attr:`~sherpa.models.model.Model.integrate` flag of
the model component is set, and whether the model supports
both or just one case.

The :py:class:`~sherpa.models.basic.Const1D` model represents
a constant value, which means that for an unbinned
dataset the model evaluates to a single value (the
:py:attr:`~sherpa.models.basic.Const1D.c0` parameter)::

    >>> from sherpa.models.basic import Const1D
    >>> mdl = Const1D()
    >>> mdl.c0 = 0.1
    >>> mdl([1, 2, 3])
    array([ 0.1,  0.1,  0.1])
    >>> mdl([-4000, 12000])
    array([ 0.1,  0.1])

The default value for its
:py:attr:`~sherpa.models.basic.Const1D.integrate` flag is
``True``::

    >>> mdl.integrate
    True

which means that this value is multiplied by the bin width when
given a binned grid (i.e. when sent in the low and high
edges of each bin):

    >>> mdl([10, 20, 30], [20, 30, 50])
    array([ 1.,  1.,  2.])

When the ``integrate`` flag is unset, the model no longer
multiplies by the bin width, and so acts similarly to the
unbinned case:

    >>> mdl.integrate = False
    >>> mdl([10, 20, 30], [20, 30, 50])
    array([ 0.1,  0.1,  0.1])

The behavior in all these three cases depends on the model - for
instance some models may raise an exception, ignore the high-edge
values in the binned case, or use the mid-point - and so the
model documentation should be reviewed.

The following example uses the
:py:class:`~sherpa.models.basic.Polynom1D` class to model the
linear relation
:math:`y = mx + c` with the origin at :math:`x = 1400`,
an offset of 2, and a gradient of 1::

    >>> from sherpa.models.basic import Polynom1D
    >>> poly = Polynom1D()
    >>> poly.offset = 1400
    >>> poly.c0 = 2
    >>> poly.c1 = 1
    >>> x = [1391, 1396, 1401, 1406, 1411]
    >>> poly(x)
    array([ -7.,  -2.,   3.,   8.,  13.])

As the integrate flag is set, the model is integrated across
each bin::
    
    >>> poly.integrate
    True
    >>> xlo, xhi = x[:-1], x[1:]
    >>> y = poly(xlo, xhi)
    >>> y
    array([-22.5,   2.5,  27.5,  52.5])

Thanks to the easy functional form chosen for this example,
it is easy to confirm that these are the values of the
integrated model::
  
    >>> (y[:-1] + y[1:]) * 5 / 2.0
    array([-50.,  75., 200.])

Turning off the ``integrate`` flag for this model shows that it
uses the low-edge of the bin when evaluating the model::
    
    >>> poly.integrate = False
    >>> poly(xlo, xhi)
    array([-7., -2.,  3.,  8.])
