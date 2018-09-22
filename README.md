[![Build Status](https://travis-ci.org/sherpa/sherpa.svg?branch=master)](https://travis-ci.org/sherpa/sherpa)
[![DOI](https://zenodo.org/badge/683/sherpa/sherpa.svg)](https://zenodo.org/badge/latestdoi/683/sherpa/sherpa)
[![GPLv3+ License](https://img.shields.io/badge/license-GPLv3+-blue.svg)](https://www.gnu.org/copyleft/gpl.html)
![Python version](https://img.shields.io/badge/Python-2.7,3.5,3.6-green.svg?style=flat)

<!-- TOC *generated with [DocToc](https://github.com/thlorenz/doctoc)* -->
**Table of Contents**

- [Sherpa](#sherpa)
- [License](#license)
- [How To Install Sherpa](#how-to-install-sherpa)
  - [Binary installation using Anaconda](#binary-installation-using-anaconda)
    - [1a. Anaconda](#1a-anaconda)
    - [1b. Starting from scratch](#1b-starting-from-scratch)
    - [1c. Other packages](#1c-other-packages)
  - [Source Build](#source-build)
    - [2a. Extract the source tarball](#2a-extract-the-source-tarball)
    - [2b. Get the code from the GitHub repository](#2b-get-the-code-from-the-github-repository)
    - [2c. Build Sherpa](#2c-build-sherpa)
    - [2d. Testing the build](#2d-testing-the-build)
    - [2e. Development mode](#2e-development-mode)
  - [Testing Sherpa](#testing-sherpa)
    - [3a. Binary installation](#3a-binary-installation)
    - [3b. Built from source](#3b-built-from-source)
- [Custom source build](#custom-source-build)
  - [FFTW library](#fftw-library)
  - [XSPEC](#xspec)
  - [Other customization options](#other-customization-options)
- [Building the documentation](#building-the-documentation)
  - [Using the build_sphinx target](#using-the-build_sphinx-target)
  - [Using the sphinx-build tool](#using-the-sphinx-build-tool)
  - [Testing the documentation with Travis](#testing-the-documentation-with-travis)
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
[Read The Docs](https://sherpa.readthedocs.io/) (at present only for
the `master` branch)
and also for [Sherpa in CIAO](http://cxc.harvard.edu/sherpa/).

A [Quick Start Tutorial](http://nbviewer.ipython.org/github/sherpa/sherpa/tree/master/notebooks/SherpaQuickStart.ipynb)
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

Sherpa can be installed from a binary distribution or built from
sources. The 4.10.0 release is available for Python 2.7, 3.5, and 3.6.

The binary distribution is available for Linux and Mac OS X via conda installation 
described in sections [1a](#1a-anaconda) and [1b](#1b-starting-from-scratch). This is the fastest
way to start using Sherpa.

Source installation is available for platforms incompatible with the
binary builds. It also allows for customization.

1. Binary installation (Anaconda)

2. Source build (from a source tarball or the GitHub repository)

Source builds can be customized, for instance:

- to point to a local build of the [FFTW library](http://www.fftw.org/)

- to build the [`XSPEC`](https://heasarc.gsfc.nasa.gov/xanadu/xspec/)
  extension to provide many common Astronomical X-ray spectral models 

These and other customization options are described below.


Binary installation using Anaconda
----------------------------------

If you already have Anaconda installed on your system, you can just follow the
easy steps in section [1a](#1a-anaconda).

If you don't have Anaconda you can follow the instructions in section [1b](#1b-starting-from-scratch),
or you can install Anaconda from:
[https://store.continuum.io/cshop/anaconda/](https://store.continuum.io/cshop/anaconda/)
and then refer to section [1a](#1a-anaconda).

Notice that section [1b](#1b-starting-from-scratch). only provides instructions on how to install a minimal
Anaconda-powered environment, not the full Anaconda distribution.

The Sherpa 4.10.0 release - which is the latest binary release - is
compatible with Python 2.7, Python 3.5, and Python 3.6.


### 1a. Anaconda

If you have Anaconda already installed on your system you can use it
to seamlessly install Sherpa.

First you need to add the Sherpa channel to your configuration,
and then install Sherpa:

    $ conda config --add channels sherpa
    $ conda install sherpa

To update Sherpa:

    $ conda update sherpa


### 1b. Starting from scratch

Miniconda is a minimal distribution of Anaconda that allows users to create
isolated virtual environments in which they can mix and match specific
versions of software. In this case the instructions show you how to install
Sherpa and its dependencies. After that you will be able to add more packages
and make them work with Sherpa. You do not need root/administrator privileges
for installing Miniconda and Sherpa, and you can activate and deactivate the
Sherpa environment every time you want.

The instructions cover both TCSH and BASH. However, Anaconda only supports
BASH, so we recommend that you start a BASH session before installing and
using Sherpa through Anaconda.

Download the Miniconda (a minimal distribution of Anaconda) installer for your
platform:

- Linux 64 bit - [https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh](https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh)
- OS X 64 bit (10.7 and forward) - [https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh](https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh)

Decide where you are going to install Miniconda, e.g.:

    $ export MINICONDA=/home/miniconda # BASH
    $ setenv MINICONDA /home/miniconda # TCSH

Run the Miniconda installer. It is assumed that you have read and agree with
the [Miniconda End User License Agreement (EULA)](http://docs.continuum.io/anaconda/eula.html)

    $ bash <Miniconda file you downloaded> -b -p $MINICONDA # BASH AND TCSH

Add miniconda to your PATH:

    $ export PATH=$PATH:$MINICONDA/bin # BASH
    $ setenv PATH "${PATH}:${MINICONDA}/bin" # TCSH

You may add these lines to your shell's startup script, e.g. `$HOME/.bash_profile`
for BASH or `$HOME/.cshrc` for TCSH.

Add the Sherpa conda repositories to your configuration:

    $ conda config --add channels sherpa

Create a new environment and install Sherpa:

    $ conda create -n sherpa sherpa=4.10

The above command will download and install Sherpa and its dependencies in an
isolated environment, so that Sherpa will not interfere with your System's
Python and you will be able to install other packages in a safe environment.

You need to activate the new environment in order to start using Sherpa. This
is straightforward in BASH, but it takes a little bit more work in TCSH. We
recommend TCSH users to create aliases for activating and deactivating the
Sherpa environment.

    $ source activate sherpa # BASH
    $ setenv OLDPATH $PATH; setenv PATH ${MINICONDA}/envs/sherpa/bin:${PATH} #TCSH

BASH users will be reminded that they are running in the Sherpa environment
by pre-pending the string (sherpa) to their BASH prompt.

When you are done working with Sherpa you can either close the terminal
window you were working with, or you can deactivate the Sherpa environment and
restore your default environment:

    $ source deactivate # BASH
    $ setenv $PATH $OLDPATH # TCSH

### 1c. Other packages

You can start using Sherpa by starting a Python shell, or you can install
`ipython` and use it as a more convenient shell. We recommend that you also install
`ipython-notebook` and `matplotlib` so that you can use the nice `ipython` notebook
features and the seamless integration with `matplotlib` for plotting from
Sherpa. We also recommend that you install `astropy` for enabling FITS I/O.

    $ conda install ipython-notebook matplotlib astropy

The minimum required version of `astropy` is version 1.3, although only 
versions 2 and higher are used in testing. There are no known
incompatabilities with `matplotlib`, but there has only been limited
testing. Please [report any problems](https://github.com/sherpa/sherpa/issues/)
you find.

Source Build
------------

The prerequisites for building from source are:

 - Python: `setuptools`, `numpy`
 - System: `gcc`, `g++`, `gfortran`, `make`, `flex`, `bison`

The full test suite requires the `mock`, `pytest >=3.3.0`, and `pytest-xvfb` packages,
which should be installed automatically if needed.

The current Sherpa code base works with Python 2.7, 3.5, and 3.6; support for
versions 3.3 and 3.4 is possible but would require community support.

It is *highly* recommended that [`matplotlib`](http://matplotlib.org/)
be installed, as this is used to create graphical output (although the
code can be built and used without this package), and
[`ipython`](http://ipython.org/), which is for interactive analysis.
Data I/O requires a `FITS` I/O library; at present the supported
libraries are [`astropy`](http://www.astropy.org) and
Crates from CIAO.

The instructions on how to set up the prerequisites vary from system to system,
and even on the same system there may be multiple ways of setting up the requirements.

There are two ways to download the source code:

 - as a source tarball

 - directly from GitHub as a git repository

**NOTE:** it is possible to build Sherpa with `fortran` compilers other than
`gfortran`. While this is not supported, [PR #202](https://github.com/sherpa/sherpa/pull/202/files)
shows how this can been accomplished with `g95` on `OS X` in a specific setup.
Similar changes are probably required for other compilers or setups.
The `fortran` extensions are compiled by
[`f2py` via `numpy.distutils`](http://docs.scipy.org/doc/numpy-1.11.0/f2py/distutils.html).

### 2a. Extract the source tarball

If you downloaded the Sherpa source tarball, you can extract it by:

    $ tar xf sherpa-<version>.tar.gz
    $ cd sherpa-<version>

### 2b. Get the code from the GitHub repository

You can clone the Sherpa repository with:

    $ git clone https://github.com/sherpa/sherpa
    $ cd sherpa

The most stable code is available through the 4.10.0 tag. The main
development code, which is unstable, is available in the `master`
branch. New features and bug fixes or other, even less stable versions
of the code may be available in other branches.

The master branch supports Python 2.7, 3.5, and 3.6. Note the
4.8.1 tag and earlier are only compatible with Python 2.7.

### 2c. Build Sherpa

Once the code is available, it can be built and installed with:

    $ python setup.py install

### 2d. Testing the build

To test that your installation of Sherpa is working, type:

    $ sherpa_test

which will run a small test suite (the script may not be in your path,
depending on where the installation step chose to install Sherpa).

Note that the test may report several `SKIPPED` lines.  These messages
are expected - as some of the tests require optional packages to be
installed alongside Sherpa. These warnings may be ignored, as long as
the test ends with an `OK` message.

**NOTE:** the `sherpa_test` command requires `pytest >=3.3.0` to run. If `pytest`
is not installed `sherpa_test` will try to install it.

### 2e. Development mode

If you plan to edit the Sherpa code, it is more convenient
to work in development mode rather than using the `install` command.

When in developer mode changes are picked up by `python` without
having to run `install` after every change.

In developer mode it may also be more convenient to use the `test`
command rather than the `sherpa_test` script:

    $ python setup.py test

The `test` command is a wrapper that calls `pytest` under the hood,
and includes the `develop` command.

You can pass additional arguments to `pytest` with the `-a` or
`--pytest-args` arguments.  As an example, a single test can be run
using the syntax:

    $ python setup.py test -a sherpa/astro/datastack/tests/test_datastack.py::test_load::test_case3

**NOTE:** if you run both `install` and `develop` or `test` in the same
Python environment you end up with two competing installations of Sherpa
which result in unexpected behavior. If this happens, simply run
`pip uninstall sherpa` as many times as necessary, until you get an
error message that no more Sherpa installations are available. At this
point you can re-install Sherpa.

The same issue may occur if you install the Sherpa binary release and
then try to build Sherpa from source in the same environment.

When both the [DS9 image viewer](http://ds9.si.edu/site/Home.html) and
[XPA toolset](http://hea-www.harvard.edu/RD/xpa/) are installed, the
test suite will include tests that check that DS9 can be used from
Sherpa. This causes several copies of the DS9 viewer to be created,
which can be distracting, as it can cause loss of mouse focus (depending
on how X-windows is set up). This can be avoided by installing the 
[X virtual-frame buffer (Xvfb)](https://en.wikipedia.org/wiki/Xvfb).


Testing Sherpa
--------------

To test that your installation works, just type:

    $ sherpa_test

The number of tests run by `sherpa_test` depends on what Python
packages are installed (for example, `astropy` and `matplotlib`), on
external software (are DS9 and the XPA toolset installed), and whether
the external Sherpa test data set is installed.

The external Sherpa test data is large, and mostly useful when
developing Sherpa. The download method depends on how Sherpa
was installed.

### 3a. Binary installation

The external test data files can be
installed from GitHub channel by saying:

    $ pip install https://github.com/sherpa/sherpa-test-data/archive/4.10.0.tar.gz

At this point, `sherpa_test` will pick up the data and so run more
tests.

### 3b. Built from source

At the top level of the Sherpa distribution, used to [build
Sherpa](#2c-build-sherpa), use the following commands to add the test
data set into the `sherpa-test-data/` directory (this assumes that the
source code was installed with `git` and not unpacked from a tarball):

    $ git submodule init
    $ git submodule update

At this point, the data will be picked up automatically by either of
the following commands:

    $ sherpa_test
    $ python setup.py test

The data files are also available as a standard Python package with its
own [`git` repository](https://github.com/sherpa/sherpa-test-data).


Custom source build
===================

There are several options for customizing the Sherpa build.

FFTW library
------------

Sherpa ships with the [`fftw`](http://www.fftw.org/) library source
code and builds it as part of its own build process by
default. However, users might want to point Sherpa to their own
version of the library. This might be required because the performance
of this library can be significantly increased by compiling it with
optimization flags specific for some system or architecture.

In order to make Sherpa build its modules against a local `fftw` library, users
need to change the default Sherpa build configuration as follows.

First, make sure you download the Sherpa source tarball, or get the source
code from GitHub:

    $ git clone https://github.com/sherpa/sherpa.git
    $ cd sherpa

Then, you need to edit the setup.cfg configuration.
This file is documented, so it should be easy to follow the instructions
therein.

In particular, you need to make sure that you set the following configuration
options:

    fftw=local
    fftw-include-dirs=/usr/local/include
    fftw-lib-dirs=/usr/local/lib
    fftw-libraries=fftw3

You might need to change the `/usr/local` path in the above example to the
actual directories on your system that contain the header (.h) files and the
`libfftw3.so` shared object.

Make sure there are no leading spaces or the python configuration system will
not parse the configuration options correctly.

Then, build Sherpa in the standard way:

    $ python setup.py install

XSPEC
-----

Sherpa can be built with support for
[`XSPEC`](https://heasarc.gsfc.nasa.gov/xanadu/xspec/), although
support is not enabled by default. The current supported XSPEC versions
are 12.10.0, 12.9.1, and 12.9.0, and it is expected that it will build against
newer versions, but without support for new models or features.

To build the XSPEC support in Sherpa, the `xspec_config` section of the
`setup.cfg` file will need changing to point to the libraries, and to
turn on the extension. In all cases, set

    with-xspec=True

The remaining settings depend on how the XSPEC libraries have
been built. In the examples below, the `$HEADAS` environment
variables **must be replaced** by the actual path to the
HEADAS installation, and the versions of the libraries - such
as `CCfits` - may need to be changed:

 1. If the full XSPEC 12.10.0 system has been built, then use

        xspec_include_dirs = $HEADAS/include
        xspec_lib_dirs = $HEADAS/lib
        xspec_libraries = XSFunctions XSModel XSUtil XS hdsp_3.0
        cfitsio_libraries = cfitsio
        ccfits_libraries = CCfits_2.5
        wcs_libraries = wcs-5.16

 2. If the full XSPEC 12.9.x system has been built, then use

        xspec_include_dirs = $HEADAS/include
        xspec_lib_dirs = $HEADAS/lib
        xspec_libraries = XSFunctions XSModel XSUtil XS
        cfitsio_libraries = cfitsio
        ccfits_libraries = CCfits_2.5
        wcs_libraries = wcs-5.16

 3. Use the model-only build of XSPEC, which will also require
    building the
    [cfitsio](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html),
    [CCfits](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/ccfits/),
    and
    [WCSLIB](http://www.atnf.csiro.au/people/mcalabre/WCS/wcslib/)
    libraries. If all the libraries are installed
    into the same location ($HEADAS), then a similar set up to the
    full XSPEC build is used

        xspec_include_dirs = $HEADAS/include
        xspec_lib_dirs = $HEADAS/lib
        xspec_libraries = XSFunctions XSModel XSUtil XS
        cfitsio_libraries = cfitsio
        ccfits_libraries = CCfits
        wcs_libraries = wcs

    although check on whether version numbers are required for the
    `cfitiso`, `CCfits`, and `wcs` libraries. If placed in different
    directories then the `cfitsio_lib_dirs`, `ccfits_lib_dirs`,
    and (possibly) `gfortran_lib_dirs` values should be set
    appropriately.

   Note that XSPEC 12.10.0 has simplified the models-only build,
   which means that the same settings as the full 12.10.0 build
   can be used as a starting point. However, by default not all
   libraries needed by Sherpa are built (e.g. `libXSModel.so`).

In all cases, the same version of `gfortran` should be used to build
Sherpa and XSPEC, in order to avoid possible incompatibilities.

If there are problems building, or using, the module, then the other
options may need to be set - in particular the `gfortran_lib_dirs` and
`gfortran_libraries` settings.

In order for the module to work, the `HEADAS` environment variable has
to be set in the shell from which the Python session is started.

In order to check that the module is working, importing the
`sherpa.astro.ui` module will no-longer warn you that the
`sherpa.astro.xspec` module is not available and you can use routines
such as:

    % python -c 'from sherpa.astro import xspec; print(xspec.get_xsversion())'
    12.10.0c

Other customization options
---------------------------

Sherpa supports other build configuration options that are required to
support Sherpa build in specific environments, for instance when
building the [CIAO analysis system](http://cxc.harvard.edu/ciao/). 
These options include:

- building Sherpa against a local version of the `CIAO` `region` and `group`
  libraries
- specify additional `CFLAGS` options for the group library
- building Sherpa against a local build of the `wcssubs` routines
- change the default `./configure` command line

The `setup.cfg` file in the Sherpa source distribution contains more information
about these options.

Building the documentation
==========================

The Sphinx documentation is available in the `docs/` directory. It is designed
so that the documentation can be built without needing to build (or install)
Sherpa, but this requires Python 3.5 or higher. An example for setting
up the build environment is:

    % conda create -n=sherpa-sphinx python=3.5 'sphinx >= 1.3' sphinx_rtd_theme
    % source activate sherpa-sphinx

Using the build_sphinx target
-----------------------------

    % python setup.py build_sphinx

The location of this output depends on the version os Sphinx in use. With
version 1.4 it appears to be located at `build/sphinx/html/index.html`.

Using the sphinx-build tool
---------------------------

The following will create the documentation in `docs/build/html/index.html`.

    % cd docs
    % sphinx-build -b html . build/html

Testing the documentation with Travis
-------------------------------------

There is a documentation build included as part of the Travis-CI test suite,
but it is not set up to do much validation. That is, you need to do something
quite severe to break this build. Please see
[issue 491](https://github.com/sherpa/sherpa/issues/491)
for more information.

History
=======

Sherpa is developed by the [Chandra X-ray
Observatory](http://chandra.harvard.edu/) to provide fitting and modelling
capabilities to the [CIAO](http://cxc.harvard.edu/ciao/) analysis package. It
has been released onto [GitHub](https://github.com/sherpa/sherpa) for users to
extend (whether to other areas of Astronomy or in other domains).

Release History
---------------

4.9.1: 01 August 2017 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.838686.svg)](https://doi.org/10.5281/zenodo.838686)

4.9.0: 27 January 2017 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.260416.svg)](https://doi.org/10.5281/zenodo.260416)

4.8.2: 23 September 2016 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.154744.svg)](https://doi.org/10.5281/zenodo.154744)

4.8.1: 15 April 2016 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.49832.svg)](https://doi.org/10.5281/zenodo.49832)

4.8.0: 27 January 2016 [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45243.svg)](https://doi.org/10.5281/zenodo.45243)

