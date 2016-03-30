[![Build Status](https://travis-ci.org/sherpa/sherpa.svg?branch=master)](https://travis-ci.org/sherpa/sherpa)
[![DOI](https://zenodo.org/badge/683/sherpa/sherpa.svg)](https://zenodo.org/badge/latestdoi/683/sherpa/sherpa)

<!-- TOC *generated with [DocToc](https://github.com/thlorenz/doctoc)* -->
**Table of Contents**

- [Sherpa](#sherpa)
- [How To Install Sherpa](#how-to-install-sherpa)
  - [Binary installation](#binary-installation)
    - [1a. Anaconda](#1a-anaconda)
    - [1b. Starting from scratch](#1b-starting-from-scratch)
  - [Source Build](#source-build)
    - [2a. Extract the source tarball](#2a-extract-the-source-tarball)
    - [2b. Get the code from the GitHub repository](#2b-get-the-code-from-the-github-repository)
    - [2c. Build Sherpa](#2c-build-sherpa)
    - [2d. Testing the build](#2d-testing-the-build)
- [Custom source build](#custom-source-build)
  - [FFTW library](#fftw-library)
  - [XSPEC](#xspec)
  - [Other customization options](#other-customization-options)
- [History](#history)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


Sherpa
======

Sherpa is a modeling and fitting application for Python. It contains a powerful language for combining simple models
into complex expressions that can be fit to the data using a variety of statistics and optimization methods.
It is easily extensible to include user models, statistics, and optimization methods.
It provides a high-level User Interface for interactive data-analysis work,
such as within an IPython notebook, and it can also be used as a library
component, providing fitting and modeling capabilities to an application.

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
- For detailed documentation see: [http://cxc.harvard.edu/sherpa](http://cxc.harvard.edu/sherpa)

A [Quick Start Tutorial](http://nbviewer.ipython.org/github/sherpa/sherpa/tree/master/docs/SherpaQuickStart.ipynb)
is included in the `docs` folder and can be opened with an `ipython notebook`.

How To Install Sherpa
=====================

Sherpa can be installed from a binary distribution or built from sources.

The binary distribution is suited for people wanting to have Sherpa up and
running as soon as possible in its standard form.

Source installation is available for platforms incompatible with the binary
builds, or for users wanting to customize the way Sherpa is built and installed.

If you are in doubt about which installation to perform, you should try
with the Conda installation (sections [1a](#1a-anaconda) and [1b](#1b-starting-from-scratch)).


1. Binary installation (Anaconda)

2. Source build (from a source tarball or the GitHub repository)

Source builds can be customized, for instance:

- to point to a local build of the [FFTW library](http://www.fftw.org/)

- to build the [`XSPEC`](https://heasarc.gsfc.nasa.gov/xanadu/xspec/)
  extension to provide many common Astronomical X-ray spectral models 

These and other customization options are descibed below.


Binary installation
-------------------

If you already have Anaconda installed on your system, you can just follow the
easy steps in section [1a](#1a-anaconda).

If you don't have Anaconda you can follow the instructions in section [1b](#1b-starting-from-scratch),
or you can install Anaconda from:
[https://store.continuum.io/cshop/anaconda/](https://store.continuum.io/cshop/anaconda/)
and then refer to section [1a](#1a-anaconda).

Notice that section [1b](#1b-starting-from-scratch). only provides instructions on how to install a minimal
Anaconda-powered environment, not the full Anaconda distribution.


### 1a. Anaconda

If you have Anaconda already installed on your system you can use it to seamlessly
install Sherpa.

First you need to add the Sherpa channel to your configuration,
and then install Sherpa:

    $ conda config --add channels https://conda.anaconda.org/sherpa
    $ conda install sherpa

To test that your installation works, just type:

    $ sherpa_test

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

- Linux 32 bit - [http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86.sh](http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86.sh)
- Linux 64 bit - [http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh](http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh)
- OS X 64 bit (10.7 and forward) - [http://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh](http://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh)

Decide where you are going to install Miniconda, e.g.:

    $ export MINICONDA=/home/miniconda # BASH
    $ setenv MINICONDA /home/miniconda # TCSH

Run the Miniconda installer. It is assumed that you have read and agree with
the [Miniconda EULA](http://docs.continuum.io/anaconda/eula.html)

    $ bash <Miniconda file you downloaded> -b -p $MINICONDA # BASH AND TCSH

Add miniconda to your PATH:

    $ export PATH=$PATH:$MINICONDA/bin # BASH
    $ setenv PATH "${PATH}:${MINICONDA}/bin" # TCSH

You may add these lines to your shell's startup script, e.g. `$HOME/.bash_profile`
for BASH or `$HOME/.cshrc` for TCSH.

Add the Sherpa conda repositories to your configuration:

    $ conda config --add channels https://conda.anaconda.org/sherpa

Create a new environment and install Sherpa:

    $ conda create -n sherpa sherpa=4.8

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

You can start using Sherpa by starting a Python shell, or you can install
`ipython` and use it as a more convenient shell. We recommend that you also install
`ipython-notebook` and `matplotlib` so that you can use the nice `ipython` notebook
features and the seamless integration with `matplotlib` for plotting from
Sherpa. We also recommend that you install `astropy` for enabling FITS I/O
(Sherpa will look for `pyfits` if `astropy` is not present).

    $ conda install ipython-notebook matplotlib astropy

When you are done working with Sherpa you can either close the terminal
window you were working with, or you can deactivate the Sherpa environment and
restore your default environment:

    $ source deactivate # BASH
    $ setenv $PATH $OLDPATH # TCSH

Please remember that you need to activate the Sherpa environment every time
you want to work with it. After the installation, the following commands are
necessary to run Sherpa:

    $ source activate sherpa # BASH
    $ setenv OLDPATH $PATH; setenv PATH ${MINICONDA}/envs/sherpa/bin:${PATH} #TCSH
    ... Run your favourite Python shell ....

When you are done and you want to go back to your environment, close the
terminal or deactivate the Sherpa environment.

    $ source deactivate sherpa # BASH  
    $ setenv $PATH $OLDPATH


Source Build
------------

The prerequisites for building from source are:

 - Python: `setuptools`, `numpy`
 - System: `gcc`, `g++`, `gfortran`, `make`, `flex`, `bison`

The current Sherpa code base only works with Python 2.7.

It is *highly* recommended that [`matplotlib`](http://matplotlib.org/)
be installed, as this is used to create graphical output (although the
code can be built and used without this package), and
[`ipython`](http://ipython.org/), which is for interactive analysis.
Data I/O requires a `FITS` I/O library. Sherpa looks for
[`astropy`](http://www.astropy.org) by default,
and it falls back to [`pyfits`](http://www.stsci.edu/institute/software_hardware/pyfits) 
if `astropy` is not installed.

The instructions on how to set up the prerequisites vary from system to system,
and even on the same system there may be multiple ways of setting up the requirements.

There are two ways to download the source code:

 - as a source tarball

 - directly from GitHub as a git repository

### 2a. Extract the source tarball

If you donwloaded the Sherpa source tarball, you can extract it by:

    $ tar xf sherpa-<version>.tar.gz
    $ cd sherpa-<version>

### 2b. Get the code from the GitHub repository

You can clone the Sherpa repository with:

    $ git clone https://github.com/sherpa/sherpa
    $ cd sherpa

The most stable code is available through the 4.8.1 tag. The main development
code, which is unstable, is available in the `master` branch. New features
and bug fixes or other, even less
stable versions of the code may be available in other branches.

### 2c. Build Sherpa

Once the code is available, it can be built and installed with:

    $ python setup.py install

### 2d. Testing the build

To test that your installation of Sherpa is working, type:

    $ sherpa_test

which will run a small test suite (the script may not be in your path,
depending on where the installation step chose to install Sherpa).

Note that the test may report several `SKIPPED` lines.
These messages are expected - as some of the
tests require optional packages to be installed alongside
Sherpa. These warnings may be ignored, as long as the test ends with
an `OK` message.

**NOTE:** the `sherpa_test` command requires `pytest` to run. If `pytest`
is not installed `sherpa_test` will try to install it.

### 2e. Development mode
If you plan to edit the Sherpa code, it is more convenient
to work in development mode rather than using the `install` command.

When in developer mode changes are picked up by `python` without
having to run `install` after every change.

In developer mode it may also be more convenient to use the `test`
command rather than the `sherpa_test` script:

    $ python setup.py test

The `test` command is a wrapper that calls `pytest` under the hood.

You can pass additional arguments to `pytest` with `-a` or `--pytest-args`.
For instance,
you can run a single test, i.e. a single test method, with:

    $ python setup.py test -a sherpa/astro/datastack/tests/test_datastack.py::test_load::test_case3

**NOTE:** if you run both `install` and `develop` or `test` in the same
Python environment you end up with two competing installations of Sherpa
which result in unexpected behavior. If this happens, simply run
`pip uninstall sherpa` as many times as necessary, until you get an
error message as no more Sherpa installations are available
and then install Sherpa again.

Also note that the `test` command executes `develop`.

The same issue may occur if you installed both the Sherpa binaries
and build Sherpa from sources in the same environment.

### 2f. Download Test Data
The `sherpa_test` and `python setup.py test` commands only execute
a small number of tests to ensure that your installation of Sherpa
is functional. The full test suite requires data files that are
not included in the Sherpa distribution by default.

If you want you can download such data files and run the whole test suite.

Since this is mostly useful when developing for Sherpa, the instructions
below assume that you are using `git` and that you are located in the
top directory:

    $ git submodule init
    $ git submodule update

This will install the data files under `sherpa-test-data`.

The data will be picked up automatically by `python setup.py test`.

The data files are included in a standard Python package with its
own [`git` repository](https://github.com/sherpa/sherpa-test-data).

You can download and install the `sherpatest` package as usual
if you are not using `git`.

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

Sherpa does not support
[`XSPEC`](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) models by
default. However, it is possible to instruct Sherpa to build its
`XSPEC` extension module by changing the build configuration options.

The `xspec_config` section of the `setup.cfg` file will need
changing to point to the libraries, and to turn on the extension.
In all cases, set

    with-xspec=True

The remaining settings depend on how the XSPEC libraries have
been built (in the examples below, environment variables are
used, but the full path should be in your own copy of the file):

 1. If the full XSPEC system has been built, then use

        xspec_lib_dirs=$HEADAS/lib
        xspec_libraries=XSFunctions XSModel XSUtil XS wcs-4.20
        cfitsio_libraries=cfitsio_3.37
        ccfits_libraries=CCfits_2.4

    The environment variable `$HEADAS` should be expanded out, and the
    version numbers of the `wcs`, `cfitsio`, and `CCfits` libraries
    may need to be changed.

 2. Use the model-only build of XSPEC, which will also require
    building the
    [cfitsio](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html),
    [CCfits](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/ccfits/),
    and
    [WCSLIB](http://www.atnf.csiro.au/people/mcalabre/WCS/wcslib/)
    libraries (it is not clear if version 5 is supported, since
    XSPEC 12.8.2 uses version 4.20). If all the libraries are installed
    into the same location ($HEADAS/lib), then a similar set up to the
    full XSPEC build is used

        xspec_lib_dirs=$HEADAS/lib
        xspec_libraries=XSFunctions XSModel XSUtil XS wcs

    except that the library names (`cfitiso`, `CCfits`, and
    `wcs`) do not need version numbers. If placed in different
    directories then the `cfitsio_lib_dirs`, `ccfits_lib_dirs`,
    and (possibly) `gfortran_lib_dirs` values should be set
    appropriately.

 3. or point to the XSPEC libraries provided by
    [CIAO](http://cxc.harvard.edu/ciao/). In this case the
    `wcs` library does not need to be specified because of
    the way the XSPEC models-only version was built with
    CIAO 4.8.

        xspec_lib_dirs=$ASCDS_INSTALL/ots/lib
        xspec_libraries=XSFunctions XSModel XSUtil XS

In all cases, the same version of `gfortran` should be used to
build Sherpa and XSPEC, to avoid possible incompatabilities.

If there are problems building, or using, the module, then the other
options may need to be set - in particular the `gfortran_lib_dirs` and
`gfortran_libraries` settings.

The XSpec module is designed for use with XSpec version 12.8.2e. It
can be used with other versions - either patches to 12.8.2 or
different versions - but there may be build or user issues.

In order for the module to work, the `HEADAS` environment variable has
to be set in the shell from which the Python session is started.  For
the CIAO-XSPEC build, `HEADAS` should be set to
`$ASCDS_INSTALL/ots/spectral`, otherwise it is the parent directory of
the `xspec_lib_dirs` directory.

In order to check that the module is working, importing the
`sherpa.astro.ui` module will no-longer warn you that the
`sherpa.astro.xspec` module is not available and you can use routines
such as:

    >>> from sherpa.astro import xspec
    >>> xspec.get_xsversion()
    '12.8.2e'

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

History
=======

Sherpa is developed by the
[Chandra X-ray Observatory](http://chandra.harvard.edu/) to provide
fitting and modelling capabilities to the
[CIAO](http://cxc.harvard.edu/ciao/) analysis package. It has been
released onto
[GitHub](https://github.com/sherpa/sherpa) under the
GNU General Public Licence, version 3 for users to extend (whether
to other areas of Astronomy or in other domains).
