****************************************
Evaluating the model on a different grid
****************************************

.. todo::

   write this. Move the api into its own page.

Sherpa now provides *experimental* support for evaluating a model
on a different grid to the independent axis. This can be used to
better model the underlying distribution (by use of a finer grid),
or to include features that lie outside the data range but, due
to the use of a convolution model, can affect the values within
the data range.

Existing Sherpa models take advantage of this support by inheriting
from the
:py:class:`~sherpa.models.model.RegriddableModel1D`
or
:py:class:`~sherpa.models.model.RegriddableModel2D` classes.
At present the only documentation is provided in the
:doc:`low-level API module <../model_classes/regrid>`.
