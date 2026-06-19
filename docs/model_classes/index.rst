.. _available-models:

****************
Available Models
****************

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

Models included in Sherpa
=========================

.. automodapi:: sherpa.models.basic

.. automodapi:: sherpa.astro.models

.. automodapi:: sherpa.astro.optical

.. automodapi:: sherpa.astro.xspec


.. _models-reference-api:

Reference/API
=============

This section describes the classes used to create models and
the :ref:`Available Models <available-models>` section above
contains the classes that implement various models.

.. automodapi:: sherpa.models.model

.. automodapi:: sherpa.models.parameter

.. automodapi:: sherpa.models.op

.. automodapi:: sherpa.models.regrid

.. automodapi:: sherpa.instrument

.. automodapi:: sherpa.models.template

.. automodapi:: sherpa.astro.instrument

.. automodapi:: sherpa.astro.xspec

.. automodapi:: sherpa.astro.xspec.utils
