*****************************
The sherpa.astro.xspec module
*****************************

.. currentmodule:: sherpa.astro.xspec

.. automodule:: sherpa.astro.xspec

This document describes the base classes for XSPEC models, and
the utility routines - such as querying and retrieving
the abundance table information. The models provided by XSPEC
are described in :doc:`astro_xspec`.

   .. rubric:: Classes

   .. autosummary::
      :toctree: api

      XSModel
      XSAdditiveModel
      XSMultiplicativeModel
      XSConvolutionKernel
      XSConvolutionModel
      XSTableModel

   .. rubric:: Functions

   .. autosummary::
      :toctree: api

      get_xsabund
      get_xschatter
      get_xscosmo
      get_xspath_manager
      get_xspath_model
      get_xsstate
      get_xsversion
      get_xsxsect
      get_xsxset
      read_xstable_model
      set_xsabund
      set_xschatter
      set_xscosmo
      set_xspath_manager
      set_xsstate
      set_xsxsect
      set_xsxset

Class Inheritance Diagram
=========================

.. inheritance-diagram:: XSModel XSAdditiveModel XSMultiplicativeModel XSConvolutionKernel XSConvolutionModel XSTableModel
   :parts: 1
