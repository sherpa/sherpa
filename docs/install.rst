************
Installation
************

Requirements
============

Sherpa has the following requirements:

* Python 2.7, 3.5, or 3.6
* NumPy (the exact lower limit has not been determined,
  but it is likely to be 1.7.0 or later)
* Linux or OS-X (patches to add Windows support are welcome)
  
Sherpa can take advantage of the following Python packages
if installed:

* :term:`astropy`: for reading and writing files in
  :term:`FITS` format.
* :term:`matplotlib`: for visualisation of
  one-dimensional data or models, one- or two- dimensional
  error analysis, and the results of Monte-Carlo Markov Chain
  runs.

The Sherpa build can be configured to create the
:py:mod:`sherpa.astro.xspec` module, which provides the models and utility
functions from the :term:`XSPEC`.
The supported versions of XSPEC are 12.10.0, 12.9.1, and 12.9.0.

Interactive display and manipulation of two-dimensional images
is available if the :term:`DS9` image viewer is installed.

Releases and version numbers
============================

The Sherpa release policy has a major release at the start of
the year, corresponding to the code that is released in the
previous December as part of the
`CIAO release <http://cxc.harvard.edu/ciao/>`_, followed by
several smaller releases throughout the year.

Information on the Sherpa releases is available from the
Zenodo page for Sherpa, using the Digital Object Identifier
(DOI) `10.5281/zenodo.593753 <https://doi.org/10.5281/zenodo.593753>`_.

What version of Sherpa is installed?
------------------------------------

The version number and git commit id of Sherpa can be retrieved from
the ``sherpa._version`` module using the following command::

    python -c 'import sherpa._version; print(sherpa._version.get_versions())'
    {'full': '8638ca0fe411693ea3b1eebe0df47512ec5bd542', 'version': '4.9.0'}

Citing Sherpa
-------------

Information on citing Sherpa can be found from the
`CITATION document <https://github.com/sherpa/sherpa/blob/master/CITATION>`_
in the Sherpa repository, or from the 
`Sherpa Zenodo page <https://doi.org/10.5281/zenodo.593753>`_.
    
Installing a pre-compiled version of Sherpa
===========================================

Additional useful Python packages include ``astropy``, ``matplotlib``,
and ``ipython-notebook``.

Using the Anaconda python distribution
--------------------------------------

The Chandra X-ray Center provides releases of Sherpa that can be
installed using
`Anaconda <https://www.continuum.io/anaconda-overview>`_
from the ``sherpa`` channel. First check
to see what the latest available version is by using::

    conda install -c sherpa sherpa --dry-run

and then, if there is a version available and there are no
significant upgrades to the dependencies, Sherpa can be installed
using::

    conda install -c sherpa sherpa
    
It is **strongly** suggested that Sherpa is installed into a named
`conda environment <http://conda.pydata.org/docs/using/envs.html>`_
(i.e. not the default environment).

Using pip
---------

Sherpa is also available from PyPI at
https://pypi.python.org/pypi/sherpa and can be installed with the
command::

    pip install sherpa

The NumPy package must already have been installed for this to work.    

Building from source
====================

Prerequisites
-------------

The prerequisites for building from source are:

* Python versions: 2.7, 3.5, 3.6
* Python packages: ``setuptools``, ``numpy``
* System: ``gcc``, ``g++``, ``gfortran``, ``make``, ``flex``,
  ``bison`` (the aim is to support recent versions of these
  tools; please report problems to the
  `Sherpa issue tracker <https://github.com/sherpa/sherpa/issues/>`_).

It is *highly* recommended that `matplotlib` and `astropy` be installed
before building Sherpa, to avoid skipping a number of tests in the
test suite.

The full Sherpa test suite requires the `mock` package (Python 2.7 only)
and `pytest`. These packages should be installed
automatically for you by the test suite if they do not already exist.

If :term:`DS9` is installed, along with the :term:`XPA` commands,
then the Sherpa test suite will cause multiple instances of DS9 to
appear and then disappear, almost immediately. To avoid this
disruption, install both the `pytest-xvfb` package and the
`X virtual frame buffer (Xvfb) <https://en.wikipedia.org/wiki/Xvfb>`_,
which will cause the DS9 tests to use the virtual framebuffer rather
than the existing X session for display.

Obtaining the source package
----------------------------

The source code can be obtained as a release package from
Zenodo - e.g.
`the CIAO 4.9 release <https://zenodo.org/record/207470>`_ -
or from
`the Sherpa repository on GitHub <https://github.com/sherpa/sherpa>`_,
either a release version,
such as the
`ciao4.9 <https://github.com/sherpa/sherpa/tree/ciao4.9>`_ tag,
or the ``master`` branch , which is not guaranteed to be stable.

For example::

    git clone git://github.com/sherpa/sherpa.git
    cd sherpa
    git checkout ciao4.9

will use the ``ciao4.9`` tag.

Configuring the build
---------------------

The Sherpa build is controlled by the ``setup.cfg`` file in the
root of the Sherpa source tree. These configuration options
include:

FFTW
^^^^

Sherpa ships with the `fftw library <http://www.fftw.org/>`_ source
code and builds it by default. To use a different version, change
the ``fftw`` options in the ``sherpa_config`` section of the
``setup.cfg`` file. The options to change are::

    fftw=local
    fftw-include_dirs=/usr/local/include
    fftw-lib-dirs=/use/local/lib
    fftw-libraries=fftw3

The ``fftw`` option must be set to ``local`` and then the remaining
options changed to match the location of the local installation.

XSPEC
^^^^^

Sherpa does not build support for
`XSPEC models <https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_
by default. This can be changed by options in the ``xspec_config``
section of the ``setup.cfg`` file::
  
    with-xspec=True
    xspec_lib_dirs=None
    xspec_include_dirs=None
    xspec_libraries=XSFunctions XSModel XSUtil XS
    cfitsio_lib_dirs=None
    cfitsio_libraries=cfitsio
    ccfits_lib_dirs=None
    ccfits_libraries=CCfits
    wcslib_lib_dirs=None
    wcslib_libraries=wcs
    gfortran_lib_dirs=None
    gfortran_libraries=gfortran

The ``with-xspec`` option must be set to ``True`` and then the
remaining options set based on whether just the
XSPEC model library or the full XSPEC system has been installed.

The remaining settings depend on how the XSPEC libraries have
been built. In the examples below, the ``$HEADAS`` environment
variables **must be replaced** by the actual path to the
HEADAS installation, and the versions of the libraries - such
as ``CCfits`` - may need to be changed:

1. If the full XSPEC 12.10.0 system has been built then use::

       xspec_include_dirs = $HEADAS/include
       xspec_lib_dirs = $HEADAS/lib
       xspec_libraries = XSFunctions XSModel XSUtil XS hdsp_3.0
       cfitsio_libraries = cfitsio
       ccfits_libraries = CCfits_2.5
       wcslib_libraries = wcs-5.16

2. If the full XSPEC 12.9.x system has been built then use::

       xspec_include_dirs = $HEADAS/include
       xspec_lib_dirs = $HEADAS/lib
       xspec_libraries = XSFunctions XSModel XSUtil XS
       cfitsio_libraries = cfitsio
       ccfits_libraries = CCfits_2.5
       wcslib_libraries = wcs-5.16

3. If the model-only build of XSPEC has been installed, then
   the configuration is similar, but the library names may
   not need version numbers and locations, depending on how the
   ``cfistio``, ``CCfits``, and ``wcs`` libraries were installed

   Note that XSPEC 12.10.0 has simplified the models-only build,
   which means that the same settings as the full 12.10.0 build
   can be used as a starting point. However, by default not all
   libraries needed by Sherpa are built (e.g. `libXSModel.so`).

A common problem is to set the `xspec_lib_dirs` option to the value
of `$HEADAS` instead of `$HEADAS/lib`. This will cause the build to
fail with errors about being unable to find the various XSPEC libraries,
such as ``XSFunctions`` and ``XSModel``.

The ``gfortran`` options should be adjusted if there are problems
building or using the XSPEC module.

In order for the XSPEC module to be used from Python, the
``HEADAS`` environment variable **must** be set before the
:py:mod:`sherpa.astro.xspec` module is imported.

The Sherpa test suite includes an extensive set of tests of this
module, but a quick check of an installed version can be done with
the following::

    % python -c 'from sherpa.astro import xspec; print(xspec.get_xsversion())'
    12.10.0c

Other options
^^^^^^^^^^^^^

The remaining options in the ``setup.cfg`` file allow Sherpa to be
built in specific environments, such as when it is built as part
of the `CIAO analysis system <http://cxc.harvard.edu/ciao/>`_. Please
see the comments in the ``setup.cfg`` file for more information on
these options.

Building and Installing
-----------------------

.. note::
   
   It is highly recommended that some form of virtual environment,
   such as a
   `conda environment <http://conda.pydata.org/docs/using/envs.html>`_
   or that provided by
   `Virtualenv <https://virtualenv.pypa.io/en/stable/>`_,
   be used when building and installing Sherpa.

A standard installation
^^^^^^^^^^^^^^^^^^^^^^^

From the root of the Sherpa source tree, Sherpa can be built by saying::

    python setup.py build

and installed with one of::

    python setup.py install
    python setup.py install --user

A development build
^^^^^^^^^^^^^^^^^^^

The ``develop`` option should be used when developing Sherpa (such as
adding new functionality or fixing a bug)::

    python setup.py develop

Tests can then be run with the ``test`` option::

    python setup.py test

The
`Sherpa test data suite <https://github.com/sherpa/sherpa-test-data>`_
can be installed to reduce the number of tests
that are skipped with the following (this is only for those builds
which used ``git`` to access the source code)::

    git submodule init
    git submodule update

.. note::

   Although the standard Python setuptools approach is used to build
   Sherpa, there may be issues when using some of the other build
   targets, such as ``build_ext``. Please report these to the
   `Sherpa issues page <https://github.com/sherpa/sherpa/issues/>`_.   
  
Building the documentation
--------------------------

.. warning::

   The documentation support is **highly experimental**. It is also
   **restricted** to Python 3 systems.

Building the documentation requires the Sherpa source code and several
additional packages:

* Python 3.5 or greater
* `Sphinx <http://sphinx.pocoo.org/>`_, version 1.3 or later
* The ``sphinx_rtd_theme``
* NumPy and six

With these installed, the documentation can be built with the
``build_sphinx`` target::

    python setup.py build_sphinx

This can be done **without** building Sherpa (either an installation
or development version), since Mock objects are used to represent
compiled and optional components.

The documentation should be placed in ``build/sphinx/html/index.html``,
although this may depend on what version of Sphinx is used.

It is also possible to build the documentation from within the ``docs/``
directory::

    cd docs
    make html

This places the documentation in ``_build/html/index.html``.
    
Testing the Sherpa installation
===============================

A very-brief "smoke" test can be run from the command-line with
the ``sherpa_smoke`` executable::

    sherpa_smoke
    WARNING: failed to import sherpa.astro.xspec; XSPEC models will not be available
    ----------------------------------------------------------------------
    Ran 7 tests in 0.456s

    OK (skipped=5)

or from the Python prompt::
  
    >>> import sherpa
    >>> sherpa.smoke()
    WARNING: failed to import sherpa.astro.xspec; XSPEC models will not be available
    ----------------------------------------------------------------------
    Ran 7 tests in 0.447s

    OK (skipped=5)
    
This provides basic validation that Sherpa has been installed
correctly, but does not run many functional tests. The screen output
will include additional warning messages if the ``astropy`` or
``matplotlib`` packages are not installed, or Sherpa was built
without support for the XSPEC model library.
    
The Sherpa installation also includes the ``sherpa_test`` command-line
tool which will run through the Sherpa test suite (the number of
tests depends on what optional packages are available and how
Sherpa was configured when built)::

    sherpa_test

The ``sherpa`` Anaconda channel contains the ``sherpatest`` package, which
provides a number of data files in ASCII and FITS formats. This is
only useful when developing Sherpa, since the package is large. It
will automatically be picked up by the ``sherpa_test`` script
once it is installed.
