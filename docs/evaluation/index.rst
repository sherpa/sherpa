******************
Evaluating a model
******************

.. todo::

   I think the text on this page needs to be re-written, since
   using a data object to evaluate the model should be made
   more obvious. This has been started but needs thinking about.

.. toctree::
   :maxdepth: 2

   integrate
   combine
   convolution
   regrid
   simulate
   cache
   examples

.. _evaluation_direct:

Direct evaluation of the model
==============================

Normally Sherpa will handle model evaluation automatically, such as
during a fit or displaying the model results. However, the
models can be evalutated directly by passing in the
grid
(:ref:`the independent axis <independent-axis>`)
directly. If ``mdl`` is an instance of a Sherpa model - that
is it is derived from the
:py:class:`~sherpa.models.model.Model`
class - then there are two standard ways to perform this
evaluation:

#. Call the model with the grid directly - e.g. for a one-dimensional
   grid use one of::

       mdl(x)
       mdl(xlo, xhi)

#. Use the :py:meth:`~sherpa.models.model.Model.calc` method, which requires
   a sequence of parameter values and then the grid; for the
   one-dimensional case this would be::

       mdl.calc(pars, x)
       mdl.calc(pars, xlo, xhi)

   In this case the parameter values do *not* need to match the values
   stored in the model itself. This can be useful when a model is to
   be embedded within another one, as shown in the
   :ref:`two-dimensional user model <example-usermodel-2d>`
   example. 

.. _evaluation_data:

Evaluating a model with a data object
=====================================

It is also possible to pass a model to a :doc:`data object <../data/index>`
and evaluate the model on a grid appropriate for the data,
using the
:py:meth:`~sherpa.data.Data.eval_model` and
:py:meth:`~sherpa.data.Data.eval_model_to_fit` methods.
This can be useful when working in an environment where the mapping
between the "native" grids used to represent data and models is
not a simple one-to-one relation, such as when analyzing
astronomical X-ray spectral data with an associated response
(i.e. a :term:`RMF` file), or
:ref:`when using a PSF <convolution-psf2d-evaluate>`.
