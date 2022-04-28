************
Installation
************

Quick overview
==============

For those users who have already read this page, and need a quick
refresher (or prefer to act first, and read documentation later),
the following commands can be used to install Sherpa, depending on
your environment and set up.

::

    conda install -c sherpa sherpa

::

    pip install sherpa

::

    python setup.py install

Requirements
============

Sherpa has the following requirements:

* Python 3.7, 3.8, 3.9, or (as a beta) 3.10
* NumPy (the exact lower limit has not been determined,
  1.17.0 or later will work, earlier version may work)
* Linux or OS-X (patches to add Windows support are welcome)

Sherpa can take advantage of the following Python packages
if installed:

* :term:`Astropy`: for reading and writing files in
  :term:`FITS` format. The minimum required version of astropy
  is version 1.3, although only versions 2 and higher are used in testing
  (version 3.2 is known to cause problems, but version 3.2.1 is okay).
* :term:`matplotlib`: for visualisation of
  one-dimensional data or models, one- or two- dimensional
  error analysis, and the results of Monte-Carlo Markov Chain
  runs. There are no known incompatabilities with matplotlib, but there
  has only been limited testing. Please
  `report any problems <https://github.com/sherpa/sherpa/issues/>`_
  you find.

The Sherpa build can be configured to create the
:py:mod:`sherpa.astro.xspec` module, which provides the models and utility
functions from :term:`XSPEC`.

Interactive display and manipulation of two-dimensional images
is available if the :term:`DS9` image viewer and the :term:`XPA`
commands are installed. It is expected that any recent version of
DS9 can be used.

Releases and version numbers
============================

The Sherpa release policy has a major release at the start of
the year, corresponding to the code that is released in the
previous December as part of the
`CIAO release <https://cxc.harvard.edu/ciao/>`_, followed by
several smaller releases throughout the year.

Information on the Sherpa releases is available from the
Zenodo page for Sherpa, using the Digital Object Identifier
(DOI) `10.5281/zenodo.593753 <https://doi.org/10.5281/zenodo.593753>`_.

What version of Sherpa is installed?
------------------------------------

The version number and git commit id of Sherpa can be retrieved from
the ``sherpa._version`` module using the following command::

    % python -c 'import sherpa._version; print(sherpa._version.get_versions())'
    {'version': '4.10.0', 'full': 'c7732043124b08d5e949b9a95c2eb6833e009421'}

Citing Sherpa
-------------

Information on citing Sherpa can be found from the
`CITATION document <https://github.com/sherpa/sherpa/blob/main/CITATION>`_
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
`conda environment <https://conda.pydata.org/docs/using/envs.html>`_
(i.e. not the default environment).

Using pip
---------

Sherpa is also available from PyPI at
https://pypi.python.org/pypi/sherpa and can be installed with the
command::

    pip install sherpa

The NumPy package must already have been installed for this to work.

.. _build-from-source:

Building from source
====================

Prerequisites
-------------

The prerequisites for building from source are:

* Python versions: 3.7, 3.8, 3.9, 3.10
* Python packages: ``setuptools``, ``numpy``
* System: ``gcc`` and ``g++`` or ``clang`` and ``clang++``, ``make``, ``flex``,
  ``bison`` (the aim is to support recent versions of these
  tools; please report problems to the
  `Sherpa issue tracker <https://github.com/sherpa/sherpa/issues/>`_).

It is *highly* recommended that `matplotlib` and `astropy` be installed
before building Sherpa, to avoid skipping a number of tests in the
test suite.

The full Sherpa test suite requires `pytest` and `pytest-xvfb`. These
packages should be installed automatically for you by the test suite
if they do not already exist.

.. warning::

   Sherpa includes a number of compiled extensions that use the NumPy
   C API. Following the `advice from NumPy
   <https://numpy.org/devdocs/reference/distutils_status_migration.html#numpy-setuptools-interaction>`_,
   it is **strongly suggested** that `setuptools < 60` is used when
   building Sherpa.

.. note::

   As of the Sherpa 4.10.1 release, a Fortran compiler is no-longer
   required to build Sherpa.

Obtaining the source package
----------------------------

The source code can be obtained as a release package from
Zenodo - e.g.
`the Sherpa 4.10.0 release <https://zenodo.org/record/1245678>`_ -
or from
`the Sherpa repository on GitHub <https://github.com/sherpa/sherpa>`_,
either a release version,
such as the
`4.10.0 <https://github.com/sherpa/sherpa/tree/4.10.0>`_ tag,
or the ``main`` branch (which is not guaranteed to be stable).

For example::

    git clone git://github.com/sherpa/sherpa.git
    cd sherpa
    git checkout 4.10.0

will use the ``4.10.0`` tag (although we strongly suggest using a
newer release now!).

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

.. _build-xspec:

XSPEC
^^^^^

.. note::

   The version number of XSPEC **must** be specified using the
   ``xspec_version`` configuration option, as described below. This is
   a change from previous releases of Sherpa, but is required in order
   to support changes made in XSPEC 12.10.0.

Sherpa can be built to use the Astronomy models provided by
:term:`XSPEC`. To enable XSPEC support, several changes must be
made to the ``xspec_config`` section of the ``setup.cfg`` file. The
available options (with default values) are::

    with-xspec = False
    xspec_version = 12.9.0
    xspec_lib_dirs = None
    xspec_include_dirs = None
    xspec_libraries = XSFunctions XSModel XSUtil XS
    cfitsio_lib_dirs = None
    cfitsio_libraries = cfitsio
    ccfits_lib_dirs = None
    ccfits_libraries = CCfits
    wcslib_lib_dirs = None
    wcslib_libraries = wcs
    gfortran_lib_dirs = None
    gfortran_libraries = gfortran

To build the :py:mod:`sherpa.astro.xspec` module, the
``with-xspec`` option must be set to ``True`` **and** the
``xspec_version`` option set to the correct version string (the XSPEC
patch level must not be included), and then the
remaining options depend on the version of XSPEC and whether
the XSPEC model library or the full XSPEC system has been installed.

In the examples below, the ``$HEADAS`` value **must be replaced**
by the actual path to the HEADAS installation, and the versions of
the libraries - such as ``CCfits_2.6`` - may need to be changed to
match the contents of the XSPEC installation.

1. If the full XSPEC 12.12.1 system has been built then use::

       with-xspec = True
       xspec_version = 12.12.1
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.30
       ccfits_libraries = CCfits_2.6
       wcslib_libraries = wcs-7.7

   where the version numbers were taken from version 6.30.1 of HEASOFT and
   may need updating with a newer release.

2. If the full XSPEC 12.12.0 system has been built then use::

       with-xspec = True
       xspec_version = 12.12.0
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.29
       ccfits_libraries = CCfits_2.6
       wcslib_libraries = wcs-7.3.1

   where the version numbers were taken from version 6.29 of HEASOFT and
   may need updating with a newer release.

3. If the full XSPEC 12.11.1 system has been built then use::

       with-xspec = True
       xspec_version = 12.11.1
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.28
       ccfits_libraries = CCfits_2.5
       wcslib_libraries = wcs-5.19.1

   where the version numbers were taken from version 6.28 of HEASOFT.

4. If the full XSPEC 12.11.0 system has been built then use::

       with-xspec = True
       xspec_version = 12.11.0
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.27
       ccfits_libraries = CCfits_2.5
       wcslib_libraries = wcs-5.19.1

   where the version numbers were taken from version 6.27 of HEASOFT.

5. If the full XSPEC 12.10.1 system has been built then use::

       with-xspec = True
       xspec_version = 12.10.1
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.26
       ccfits_libraries = CCfits_2.5
       wcslib_libraries = wcs-5.19.1

   where the version numbers were taken from version 6.26.1 of HEASOFT.

6. If the full XSPEC 12.10.0 system has been built then use::

       with-xspec = True
       xspec_version = 12.10.0
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSModel XSUtil XS hdsp_3.0
       ccfits_libraries = CCfits_2.5
       wcslib_libraries = wcs-5.16

7. If the full XSPEC 12.9.x system has been built then use::

       with-xspec = True
       xspec_version = 12.9.1
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSModel XSUtil XS
       ccfits_libraries = CCfits_2.5
       wcslib_libraries = wcs-5.16

   changing ``12.9.1`` to ``12.9.0`` as appropriate.

8. If the model-only build of XSPEC has been installed, then
   the configuration is similar, but the library names may
   not need version numbers and locations, depending on how the
   ``cfitsio``, ``CCfits``, and ``wcs`` libraries were installed.

   Note that XSPEC 12.10.0 introduces a new ``--enable-xs-models-only``
   flag when building HEASOFT which simplifies the installation of
   these extra libraries, but can cause problems for the Sherpa build.

A common problem is to set one or both of the ``xspec_lib_dirs``
and ``xspec_lib_include`` options to the value of ``$HEADAS`` instead of
``$HEADAS/lib`` and ``$HEADAS/include`` (after expanding out the
environment variable). Doing so will cause the build to fail with
errors about being unable to find various XSPEC libraries such as
``XSFunctions`` and ``XSModel``.

The ``gfortran`` options should be adjusted if there are problems
using the XSPEC module.

In order for the XSPEC module to be used from Python, the
``HEADAS`` environment variable **must** be set before the
:py:mod:`sherpa.astro.xspec` module is imported.

The Sherpa test suite includes an extensive set of tests of this
module, but a quick check of an installed version can be made with
the following command::

    % python -c 'from sherpa.astro import xspec; print(xspec.get_xsversion())'
    12.12.0

Other options
^^^^^^^^^^^^^

The remaining options in the ``setup.cfg`` file allow Sherpa to be
built in specific environments, such as when it is built as part
of the `CIAO analysis system <https://cxc.harvard.edu/ciao/>`_. Please
see the comments in the ``setup.cfg`` file for more information on
these options.

Installing all dependencies with conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`source-install-with-conda` for details on how to set up all
dependencies for the Sherpa build with conda.

Building and Installing
-----------------------

It is highly recommended that some form of virtual environment,
such as a
`conda environment <https://conda.pydata.org/docs/using/envs.html>`_
or that provided by
`Virtualenv <https://virtualenv.pypa.io/en/stable/>`_,
be used when building and installing Sherpa.

The ``CC`` and ``CXX`` environment variables can be set to the C and
C++ compilers to use if not found by ``setup.py``.

.. warning::

   When building Sherpa on macOS within a conda environment, the following
   environment variable must be set otherwise importing Sherpa will
   crash Python::

     export PYTHON_LDFLAGS=' '

   That is, the variable is set to a space, not the empty string.

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

The ``test`` command is a wrapper that calls ``pytest`` under the hood,
and includes the ``develop`` command.

You can pass additional arguments to ``pytest`` with the ``-a`` or
``--pytest-args`` arguments.  As examples, the following two commands
run all the tests in ``test_data.py`` and then a single named
test in this file::

    python setup.py test -a sherpa/tests/test_data.py
    python setup.py test -a sherpa/tests/test_data.py::test_data_eval_model

The full set of options, including those added by the Sherpa test
suite - which are listed at the end of the ``custom options``
section - can be found with::

    python setup.py test -a "--pyargs sherpa --help"

and to pass an argument to the Sherpa test suite (there are currently
three options, namely ``--test-data``, ``--runslow``, and
``-runzenodo``)::

    python setup.py test -a "--pyargs sherpa --runslow"

.. note::

   If you run both ``install`` and ``develop`` or ``test`` in the same
   Python environment you end up with two competing installations of
   Sherpa which result in unexpected behavior. If this happens, simply
   run ``pip uninstall sherpa`` as many times as necessary, until you
   get an error message that no more Sherpa installations are
   available. At this point you can re-install Sherpa.

   The same issue may occur if you install a Sherpa binary release and
   then try to build Sherpa from source in the same environment.

The
`Sherpa test data suite <https://github.com/sherpa/sherpa-test-data>`_
can be installed to reduce the number of tests
that are skipped with the following (this is only for those builds
which used ``git`` to access the source code)::

    git submodule init
    git submodule update

When both the `DS9 image viewer <https://ds9.si.edu/>`_ and
`XPA toolset <https://hea-www.harvard.edu/RD/xpa/>`_ are installed, the
test suite will include tests that check that DS9 can be used from
Sherpa. This causes several copies of the DS9 viewer to be created,
which can be distracting, as it can cause loss of mouse focus (depending
on how X-windows is set up). This can be avoided by installing the
`X virtual-frame buffer (Xvfb) <https://en.wikipedia.org/wiki/Xvfb>`_.

.. note::

   Although the standard Python setuptools approach is used to build
   Sherpa, there may be issues when using some of the other build
   targets, such as ``build_ext``. Please report these to the
   `Sherpa issues page <https://github.com/sherpa/sherpa/issues/>`_.

Building the documentation
--------------------------

Building the documentation requires the Sherpa source code and several
additional packages:

* `Sphinx <https://sphinx.pocoo.org/>`_, version 1.8 or later
* The ``sphinx_rtd_theme``
* NumPy and `sphinx-astropy <https://github.com/astropy/sphinx-astropy/>`_
  (the latter can be installed with ``pip``)
* `nbsphinx <https://pypi.org/project/nbsphinx/>`_, ``ipykernel``, and ``pandoc``
  for including Jupyter notebooks
* `Graphviz <https://www.graphviz.org/>`_ (for the inheritance diagrams)

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
provides a number of data files in ASCII and :term:`FITS` formats. This is
only useful when developing Sherpa, since the package is large. It
will automatically be picked up by the ``sherpa_test`` script
once it is installed.
