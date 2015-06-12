[![Build Status](https://travis-ci.org/sherpa/sherpa.svg?branch=master)](https://travis-ci.org/sherpa/sherpa)

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

    $ conda config --add channels https://conda.binstar.org/sherpa
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

    $ conda config --add channels https://conda.binstar.org/sherpa

Create a new environment and install Sherpa:

    $ conda create -n sherpa sherpa=4.7

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
Sherpa. We also recommend that you install `pyfits` for enabling FITS I/O.

    $ conda install ipython-notebook matplotlib pyfits

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
Data I/O requires the [`pyfits`](http://www.stsci.edu/institute/software_hardware/pyfits) 
package.

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

The most stable code is available in the `stable` branch. The main development
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
Note that the test may report several `WARNING` lines and failed
attempts to load modules. These messages are expected - as some of the
tests require optional packages to be installed alongside
Sherpa. These warnings may be ignored, as long as the test ends with
an `OK` message.

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

You may need to build `XSPEC` yourself, and in any case to  point Sherpa to existing
binary libraries for `XSPEC`, `cfitsio`, and `CCfits`.
Additionaly, you will need to point Sherpa to the `libgfortran` shared library.
These dependencies are required to build `XSPEC` itself, so they are assumed to
be present on the system where you want to build Sherpa with `XSPEC` support.

First, download the Sherpa source tarball or get the source code from GitHub:

    $ git clone https://github.com/sherpa/sherpa.git
    $ cd sherpa

Now, edit the `setup.cfg` file. Find the XSPEC configuration section in the
file, uncomment the relative options and make sure they point to the location
of the `XSPEC`, `cfitsio`, `CCfits`, and `gfortran` libraries. For instance:

    with-xspec=True
    xspec_lib_dirs=/opt/xspec/lib
    cfitsio_lib_dirs=/usr/local/lib
    ccfits_lib_dirs=/usr/local/lib
    gfortran_lib_dirs=/usr/local/lib

You may need to change the values in the above example to reflect the actual
directories where the libraries are to be found on your system.

Then, build Sherpa in the standard way:

    $ python setup.py install


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
