.. _available-models:

****************
Available Models
****************

.. todo::
   
   I am not convinced separating the models out here from the
   API below is a good idea, as we end up with multiple pages
   describing the `sherpa.astro.xspec` module (for instance).
   It seems less of an issue than I was originally worried about,
   but this should still be reviewed (in particular when looking
   at the auto-generated index or module-list pages).

.. note::
   
   The models in :py:mod:`sherpa.astro.xspec` are only available if
   Sherpa was built with support for the
   `XSPEC model library
   <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixExternal.html>`_.

This section describes the classes that implement models used
to describe and fit data, while the
:ref:`Reference/API <models-reference-api>`
section below describes the classes used to create these models.

.. toctree::
   :maxdepth: 2

   usermodel
   basic
   astro_models
   astro_optical
   astro_xspec

.. _models-reference-api:

Reference/API
=============

This section describes the classes used to create models and
the :ref:`Available Models <available-models>` section above
contains the classes that implement various models.

.. toctree::
   :maxdepth: 2

   model
   parameters
   instrument
   template
   astro_instrument
   xspec_model
   
