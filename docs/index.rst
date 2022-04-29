
.. hide the title (based on astropy)
.. raw:: html

    <style media="screen" type="text/css">h1 { display: none; }</style>

Welcome to Sherpa's documentation
=================================

.. image:: _static/sherpa_logo.png
    :width: 350px
    :height: 132px
    :target: https://cxc.harvard.edu/contrib/sherpa/

Welcome to the Sherpa documentation.
`Sherpa <https://cxc.harvard.edu/contrib/sherpa/>`_
is a Python package for
modeling and fitting data. It was originally developed by the
`Chandra X-ray Center <https://cxc.harvard.edu/>`_ for use in
`analysing X-ray data (both spectral and imaging)
<https://cxc.harvard.edu/sherpa>`_
from the  Chandra X-ray telescope, but it is designed to be a
general-purpose package, which can be enhanced with domain-specific
tasks (such as X-ray Astronomy).
Sherpa contains an expressive and powerful
modeling language, coupled with a
range of statistics and robust optimisers.

.. seealso::

   If you are looking for the similarly named package
   "SHERPA" for hyperparameter tuning of machine learning models
   go here:
   https://parameter-sherpa.readthedocs.io/


Sherpa is released under the
`GNU General Public License v3.0
<https://github.com/sherpa/sherpa/blob/master/LICENSE>`_,
and is compatible with Python versions 3.7 to 3.10.
Information on recent releases and citation information for
Sherpa is available using the Digital Object Identifier (DOI)
`10.5281/zenodo.593753 <https://doi.org/10.5281/zenodo.593753>`_.

The last version of Sherpa compatible with Python 2.7 was the
`4.11.1 release <https://doi.org/10.5281/zenodo.3358134>`_.

.. toctree::
   :maxdepth: 2
   :caption: Introduction
   :name: intro

   install
   quick
   ciao

.. toctree::
   :maxdepth: 2
   :caption: User Documentation
   :name: user_docs

   data/index
   models/index
   evaluation/index
   model_classes/index
   statistics/index
   optimisers/index
   fit/index
   plots/index
   mcmc/index
   overview/utilities

.. toctree::
   :maxdepth: 2
   :caption: Worked Examples
   :name: worked_examples

   examples/simple_interpolation
   examples/fit_peaked_data
   examples/simple_user_model

.. toctree::
   :maxdepth: 2
   :caption: An interactive application
   :name: ui

   ui/index

.. toctree::
   :maxdepth: 2
   :caption: Notebooks
   :name: notebooks

   SherpaQuickStart
   NotebookSupport
   ExamplePlots

.. toctree::
   :maxdepth: 2
   :caption: Extra Functionality
   :name: extra

   extra/optimize_function

.. toctree::
   :maxdepth: 2
   :caption: Getting Help
   :name: help_docs

   bugs
   developer/index
   indices

At present there is no developer mailing list for Sherpa.
