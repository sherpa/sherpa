************
Installation
************

Quick overview
==============

For those users who have already read this page, and need a quick
refresher (or prefer to act first, and read documentation later),
the following commands can be used to install Sherpa, depending on
your environment and set up.

#. Using conda

   ::

     conda install -c https://cxc.cfa.harvard.edu/conda/sherpa -c conda-forge sherpa

#. Install Sherpa using pip

   ::

     pip install sherpa

#. Building from source

   ::

     pip install .

Requirements
============

Sherpa has the following requirements:

* Python 3.11 to 3.13
* NumPy
* Linux or OS-X (patches to add Windows support are welcome)

Sherpa can take advantage of the following Python packages
if installed:

* :term:`Astropy`: for reading and writing files in
  :term:`FITS` format.
* :term:`matplotlib`: for visualisation of
  one-dimensional data or models, one- or two- dimensional
  error analysis, and the results of Monte-Carlo Markov Chain
  runs. There are no known incompatibilities with matplotlib, but there
  has only been limited testing. Please
  `report any problems <https://github.com/sherpa/sherpa/issues/>`_
  you find.
* `ArviZ <https://python.arviz.org>`_ for visualisation and analysis of
  `sherpa.sim.MCMC` results.

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

Using the Conda python distribution
--------------------------------------

The Chandra X-ray Center provides releases of Sherpa that can be
installed using
`Miniforge <https://github.com/conda-forge/miniforge>`_.
First check to see what the latest available version is by using::

    conda install -c https://cxc.cfa.harvard.edu/conda/sherpa -c conda-forge sherpa --dry-run

and then, if there is a version available and there are no
significant upgrades to the dependencies, Sherpa can be installed
using::

    conda install -c https://cxc.cfa.harvard.edu/conda/sherpa -c conda-forge sherpa

It is **strongly** suggested that Sherpa is installed into a named
`conda environment <https://conda.pydata.org/docs/using/envs.html>`_
(i.e. not the default environment).

Using pip
---------

Sherpa is also available from PyPI at
https://pypi.python.org/pypi/sherpa and can be installed with the
command::

    pip install sherpa

.. _build-from-source:

Building from source
====================

Prerequisites
-------------

The prerequisites for building from source are:

* Python versions: 3.9 to 3.12
* Python packages: ``setuptools``, ``numpy`` (these should be
  automatically installed by ``pip``)
* System: ``gcc`` and ``g++`` or ``clang`` and ``clang++``, ``make``, ``flex``,
  ``bison``, ``ar`` (which may be provided by the ``binutils`` package), ``file``.

The aim is to support recent versions of these tools and libraries;
please report problems to the
`Sherpa issue tracker <https://github.com/sherpa/sherpa/issues/>`_.

It is *highly* recommended that `matplotlib` and `astropy` be installed
before building Sherpa, to avoid skipping a number of tests in the
test suite.

The full Sherpa test suite requires `pytest`, which is included when
using the ``.[test]`` option with ``pip``. The `pytest-xvfb` package
can be useful if :term:`DS9` is installed, as it hides the DS9 windows
created during the tests.

.. note::

   As of the Sherpa 4.10.1 release, a Fortran compiler is no-longer
   required to build Sherpa.

Obtaining the source package
----------------------------

The source code can be obtained as a release package from
Zenodo - e.g.
`the Sherpa 4.16.0 release <https://zenodo.org/record/825839>`_ -
or from
`the Sherpa repository on GitHub <https://github.com/sherpa/sherpa>`_,
either a release version,
such as the
`4.16.0 <https://github.com/sherpa/sherpa/tree/4.16.0>`_ tag,
or the ``main`` branch (which is not guaranteed to be stable).

For example::

    git clone https://github.com/sherpa/sherpa.git
    cd sherpa
    git checkout 4.16.0

will use the ``4.16.0`` tag (although we strongly suggest using a
newer release now!).

Configuring the build
---------------------

The Sherpa build is controlled by the ``setup.cfg`` file in the
root of the Sherpa source tree. These configuration options
include:

FFTW
^^^^

Sherpa ships with the `fftw library <https://www.fftw.org/>`_ source
code and builds it by default. To use a different version, change
the ``fftw`` options in the ``sherpa_config`` section of the
``setup.cfg`` file. The options to change are::

    fftw=local
    fftw_include_dirs=/usr/local/include
    fftw_lib_dirs=/use/local/lib
    fftw_libraries=fftw3

The ``fftw`` option must be set to ``local`` and then the remaining
options changed to match the location of the local installation.

.. _build-xspec:

XSPEC
^^^^^

Sherpa can be built to use the Astronomy models provided by
:term:`XSPEC`. To enable XSPEC support, several changes must be
made to the ``xspec_config`` section of the ``setup.cfg`` file. The
available options (with default values) are::

    with_xspec = False
    xspec_version = 12.12.0
    xspec_lib_dirs = None
    xspec_include_dirs = None
    xspec_libraries = XSFunctions XSUtil XS
    cfitsio_lib_dirs = None
    cfitsio_libraries =
    ccfits_lib_dirs = None
    ccfits_libraries =
    wcslib_lib_dirs = None
    wcslib_libraries =
    gfortran_lib_dirs = None
    gfortran_libraries =

To build the :py:mod:`sherpa.astro.xspec` module, the
``with_xspec`` option must be set to ``True`` **and** the
``xspec_version`` option set to the correct version string (the XSPEC
patch level must not be included), and then the
remaining options depend on the version of XSPEC and whether
the XSPEC model library or the full XSPEC system has been installed.

In the examples below, the ``$HEADAS`` value **must be replaced**
by the actual path to the HEADAS installation, and the versions of
the libraries - such as ``CCfits_2.6`` - may need to be changed to
match the contents of the XSPEC installation.

1. If the full XSPEC 12.15.0 system has been built then use::

       with_xspec = True
       xspec_version = 12.15.0
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.35
       ccfits_libraries = CCfits_2.7
       wcslib_libraries = wcs-8.3

   where the version numbers were taken from version 6.35 of HEASOFT and
   may need updating with a newer release.

2. If the full XSPEC 12.14.1 system has been built then use::

       with_xspec = True
       xspec_version = 12.14.1
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.34
       ccfits_libraries = CCfits_2.6
       wcslib_libraries = wcs-8.3

   where the version numbers were taken from version 6.34 of HEASOFT and
   may need updating with a newer release.

3. If the full XSPEC 12.14.0 system has been built then use::

       with_xspec = True
       xspec_version = 12.14.0
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.33
       ccfits_libraries = CCfits_2.6
       wcslib_libraries = wcs-8.2.1

   where the version numbers were taken from version 6.33.1 of HEASOFT and
   may need updating with a newer release.

4. If the full XSPEC 12.13.1 system has been built then use::

       with_xspec = True
       xspec_version = 12.13.1
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.32
       ccfits_libraries = CCfits_2.6
       wcslib_libraries = wcs-7.7

   where the version numbers were taken from version 6.32 of HEASOFT and
   may need updating with a newer release.

5. If the full XSPEC 12.13.0 system has been built then use::

       with_xspec = True
       xspec_version = 12.13.0
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.31
       ccfits_libraries = CCfits_2.6
       wcslib_libraries = wcs-7.7

6. If the full XSPEC 12.12.1 system has been built then use::

       with_xspec = True
       xspec_version = 12.12.1
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.30
       ccfits_libraries = CCfits_2.6
       wcslib_libraries = wcs-7.7

7. If the full XSPEC 12.12.0 system has been built then use::

       with_xspec = True
       xspec_version = 12.12.0
       xspec_lib_dirs = $HEADAS/lib
       xspec_include_dirs = $HEADAS/include
       xspec_libraries = XSFunctions XSUtil XS hdsp_6.29
       ccfits_libraries = CCfits_2.6
       wcslib_libraries = wcs-7.3.1

8. If the model-only build of XSPEC - created with the
   ``--enable-xs-models-only`` flag when building HEASOFT - has been
   installed, then the configuration is similar, but the library names
   may not need version numbers and locations, depending on how the
   ``cfitsio``, ``CCfits``, and ``wcs`` libraries were installed.

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
    12.15.0

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

.. _install-build:

A standard installation
^^^^^^^^^^^^^^^^^^^^^^^

From the root of the Sherpa source tree, Sherpa can be built with

::

  pip install .

Please report any problems to the
`Sherpa issues page <https://github.com/sherpa/sherpa/issues/>`_.

.. _developer-build:

A development build
^^^^^^^^^^^^^^^^^^^

The code can be built locally, which is useful when adding new
functionality or fixing a bug (the ``[test]`` term just ensures that
``pytest`` is also installed)::

  pip install -e .[test]

This will need to be re-run when any of the extension models - that is,
any compiled code - is changed.

The ``--verbose`` flag is useful when diagnosing problems when building Sherpa::

  pip install -e .[test] --verbose

Testing Sherpa
^^^^^^^^^^^^^^

Tests can be run directly for the development build with::

  pytest

You can pass additional arguments to ``pytest``. As examples, the
following two commands run all the tests in ``test_data.py`` and then
a single named test in the file::

  pytest sherpa/tests/test_data.py
  pytest sherpa/tests/test_data.py::test_data_eval_model

The full set of options, including those added by the Sherpa test
suite - which are listed at the end of the ``custom options``
section - can be found with::

  pytest --pyargs sherpa --help

and to pass an argument to the Sherpa test suite (there are currently
three options, namely ``--test-data``, ``--runslow``, and
``--runzenodo``)::

  pytest --pyargs sherpa --runslow

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
`X virtual-frame buffer (Xvfb) <https://en.wikipedia.org/wiki/Xvfb>`_
and ensuring that the ``pytest-xvfb`` Python package is installed.

Tests can be run in parallel with the `pytest-xdist
<https://pytest-xdist.readthedocs.io/>`_ package installed. The safest
way is to include the `--dist=loadgroup` option (although this is only
needed if the DS9 tests are run)::

    pip install pytest-xdist
    pytest --dist=loadgroup -n auto

Building the documentation
--------------------------

Building the documentation requires a Sherpa installation and several
additional packages:

* `Sphinx <https://sphinx.pocoo.org/>`_, version 1.8 or later
* The ``sphinx_rtd_theme``
* NumPy and `sphinx-astropy <https://github.com/astropy/sphinx-astropy/>`_
  (the latter can be installed with ``pip``)
* `nbsphinx <https://pypi.org/project/nbsphinx/>`_, ``ipykernel``, and ``pandoc``
  for including Jupyter notebooks
* `Graphviz <https://www.graphviz.org/>`_ (for the inheritance diagrams)

The easiest way to install the Python packages is to install the ``doc``
option with::

  pip install .[doc]

This also ensures that Sherpa has been built, as this is needed to
build the documentation.

If conda is being used then the other packages can be installed with::

  conda install -c conda-forge pandoc graphviz

With these installed, the documentation can be built::

  cd docs
  make html

Only very specific modules are mocked out because they are hard to
build and are not needed for the documentation build (currently ds9
and XSPEC).

The documentation should be placed in ``docs/_build/html/index.html``.

.. note::

   Prior to Sherpa 4.16.0 the documentation was built directly from the
   source - using mock objects to handle compiled code - rather than
   using a Sherpa installation. As of 4.16.0, mock objects are only
   handled for the XSPEC and DS9 modules.

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
tool which will run through the Sherpa test suite (the number of tests
depends on what optional packages are available and how Sherpa was
configured when built)::

    sherpa_test

The ``sherpa_test`` command supports the same optional arguments as
``pytest`` does (the ``--pyargs sherpa`` option is, however, not
needed).

The
`Sherpa test data suite <https://github.com/sherpa/sherpa-test-data>`_
contains the ``sherpatest`` package, which provides a number of
data files in ASCII and :term:`FITS` formats. This is
only useful when developing Sherpa, since the package is large.
A version of the test data is released for each `version of Sherpa <https://doi.org/10.5281/zenodo.593753>`_.

As an example, the 4.15.1 version of the test data can be installed with pip::

   pip install https://github.com/sherpa/sherpa-test-data/archive/4.15.1.zip

The test data will automatically be picked up by the ``sherpa_test``
script once it is installed.
