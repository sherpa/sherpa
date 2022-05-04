![Build Status: Conda](https://github.com/sherpa/sherpa/workflows/Conda%20CI/badge.svg)
![Build Status: Pip](https://github.com/sherpa/sherpa/workflows/Pip%20CI/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/sherpa/badge/)](https://sherpa.readthedocs.io/)
[![DOI](https://zenodo.org/badge/683/sherpa/sherpa.svg)](https://zenodo.org/badge/latestdoi/683/sherpa/sherpa)
[![GPLv3+ License](https://img.shields.io/badge/license-GPLv3+-blue.svg)](https://www.gnu.org/copyleft/gpl.html)
![Python version](https://img.shields.io/badge/Python-3.7,3.8,3.9,3.10-green.svg?style=flat)

<!-- TOC *generated with [DocToc](https://github.com/thlorenz/doctoc)* -->
**Table of Contents**

- [Sherpa](#sherpa)
- [License](#license)
- [How To Install Sherpa](#how-to-install-sherpa)
  - [Using Anaconda](#using-anaconda)
  - [Using pip](#using-pip)
  - [Building from source](#building-from-source)
- [History](#history)
  - [Release History](#release-history)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


Sherpa
======

Sherpa is a modeling and fitting application for Python. It contains a
powerful language for combining simple models into complex expressions
that can be fit to the data using a variety of statistics and
optimization methods.  It is easily extensible to include user models,
statistics, and optimization methods.  It provides a high-level User
Interface for interactive data-analysis work, such as within a
Jupyter notebook, and it can also be used as a library component,
providing fitting and modeling capabilities to an application.

What can you do with Sherpa?

- fit 1D (multiple) data including: spectra, surface brightness profiles, light curves, general ASCII arrays
- fit 2D images/surfaces in Poisson/Gaussian regime
- build complex model expressions
- import and use your own models
- use appropriate statistics for modeling Poisson or Gaussian data
- import new statistics, with priors if required by analysis
- visualize the parameter space with simulations or using 1D/2D cuts of the parameter space
- calculate confidence levels on the best fit model parameters
- choose a robust optimization method for the fit: Levenberg-Marquardt, Nelder-Mead Simplex or Monte Carlo/Differential Evolution.

Documentation for Sherpa is available at
[Read The Docs](https://sherpa.readthedocs.io/)
and also for [Sherpa in CIAO](http://cxc.harvard.edu/sherpa/).

A [Quick Start Tutorial](http://nbviewer.ipython.org/github/sherpa/sherpa/tree/main/notebooks/SherpaQuickStart.ipynb)
is included in the `notebooks` folder and can be opened with an `ipython notebook`.

License
=======

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. A copy of the GNU General Public License can be found in the
`LICENSE` file provided with the source code, or from the
[Free Software Foundation](http://www.gnu.org/licenses/).

How To Install Sherpa
=====================

[Full installation instructions](https://sherpa.readthedocs.io/en/latest/install.html)
are part of the [Read The Docs](https://sherpa.readthedocs.io/)
documentation, and should be read if the following is not sufficient.

It is strongly recommended that some form of *virtual environment* is
used with Sherpa.

Sherpa is tested against Python versions 3.7, 3.8, 3.9, and (as a beta) 3.10.

The last version of Sherpa which supported Python 2.7 is
[Sherpa 4.11.1](https://doi.org/10.5281/zenodo.3358134).

Using Anaconda
--------------

Sherpa is provided for both Linux and macOS operating systems running
Python 3.7 to 3.10. It can be installed with the `conda`
package manager by saying

    $ conda install -c sherpa sherpa

Using pip
---------

Sherpa is also available
[on PyPI](https://pypi.python.org/pypi/sherpa) and so can be installed
with the following command (which requires that the NumPy package is
already installed).

    % pip install sherpa

Building from source
--------------------

Source installation is available for platforms incompatible with the
binary builds, or for when the default build options are not sufficient
(such as including support for the
[`XSPEC` model library](https://heasarc.gsfc.nasa.gov/xanadu/xspec/)).
The steps are described in the
[building from source](https://sherpa.readthedocs.io/en/latest/install.html#building-from-source)
documentation.

History
=======

Sherpa is developed by the [Chandra X-ray
Observatory](http://chandra.harvard.edu/) to provide fitting and modelling
capabilities to the [CIAO](http://cxc.harvard.edu/ciao/) analysis package. It
has been released onto [GitHub](https://github.com/sherpa/sherpa) for users to
extend (whether to other areas of Astronomy or in other domains).

Release History
---------------

4.14.0: 07 October 2021 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5554957.svg)](https://doi.org/10.5281/zenodo.5554957)

4.13.1: 18 May 2021 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4770623.svg)](https://doi.org/10.5281/zenodo.4770623)

4.13.0: 08 January 2021 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4428938.svg)](https://doi.org/10.5281/zenodo.4428938)

4.12.2: 27 October 2020 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4141888.svg)](https://doi.org/10.5281/zenodo.4141888)

4.12.1: 14 July 2020 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3944985.svg)](https://doi.org/10.5281/zenodo.3944985)

4.12.0: 30 January 2020 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3631574.svg)](https://doi.org/10.5281/zenodo.3631574)

4.11.1: 1 August 2019 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3358134.svg)](https://doi.org/10.5281/zenodo.3358134)

4.11.0: 20 February 2019 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2573885.svg)](https://doi.org/10.5281/zenodo.2573885)

4.10.2: 14 December 2018 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2275738.svg)](https://doi.org/10.5281/zenodo.2275738)

4.10.1: 16 October 2018 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1463962.svg)](https://doi.org/10.5281/zenodo.1463962)

4.10.0: 11 May 2018 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1245678.svg)](https://doi.org/10.5281/zenodo.1245678)

4.9.1: 01 August 2017 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.838686.svg)](https://doi.org/10.5281/zenodo.838686)

4.9.0: 27 January 2017 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.260416.svg)](https://doi.org/10.5281/zenodo.260416)

4.8.2: 23 September 2016 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.154744.svg)](https://doi.org/10.5281/zenodo.154744)

4.8.1: 15 April 2016 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.49832.svg)](https://doi.org/10.5281/zenodo.49832)

4.8.0: 27 January 2016 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45243.svg)](https://doi.org/10.5281/zenodo.45243)
