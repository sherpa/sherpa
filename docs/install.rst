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

* Python 3.11 to 3.14 (there is no support for free-threaded Python)
* NumPy
* `FFTW <https://www.fftw.org/>`_
* Linux or OS-X (patches to add Windows support are welcome)

Sherpa can take advantage of the following Python packages
if installed:

* :term:`Astropy`: for reading and writing files in
  :term:`FITS` format.
* :term:`matplotlib` or :term:`bokeh`: for visualisation of
  one-dimensional data or models, one- or two- dimensional
  error analysis, and the results of Monte-Carlo Markov Chain
  runs.
* `scipy <https://www.scipy.org>`_ for the `sherpa.optmethods.optscipy`
  module, which provides an interface to several optimizers from Scipy.
* `optimagic <https://optimagic.readthedocs.io>`_ for the
  `sherpa.optmethods.optoptimagic` module which provides an interface to `optimagic.minimize`.
  Optimagic provides a single interface to dozens of optimizers, many of which
  rely on other external packages that need to be `installed separately
  <https://optimagic.readthedocs.io/en/latest/installation.html>`_.
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

    % python -c 'import sherpa; print(sherpa._version.version)'
    4.19.0
    % python -c 'import sherpa; print(sherpa._version.git_revision)'
    ...

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

.. versionchanged:: 4.19.0

   The build backend was changed from ``setuptools`` to ``meson-python``
   in the Sherpa 4.19.0 release. There are therefore a number of
   changes to what is needed to install Sherpa, how to configure the
   build, and how to test Sherpa.

Prerequisites
-------------

The prerequisites for building from source are:

* Python versions: 3.11 to 3.14 (there is no support for free-threaded Python)
* The FFTW3 library and headers, accessible via `pkg-config`
* Python packages: ``setuptools``, ``numpy`` (these should be
  automatically installed by ``pip``)
* System: ``gcc`` and ``g++`` or ``clang`` and ``clang++``, ``flex``,
  ``bison``, ``ar`` (which may be provided by the ``binutils`` package).

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

.. _configure:

Configuring the build
---------------------

The configuration options are listed in the ``meson.options``
file. These options need to be passed to the Python build step using
the syntax (this needs to be done for each option)::

  -Csetup-args=-D<option name>=<option value>

.. versionchanged:: 4.19.0

   Prior to 4.19.0 the configuration options were set in the
   ``setup.cfg`` file. The names and options have been changed,
   and they are now specified when building Sherpa, and not read from
   a file.

.. _build_limited_api:

Python limited API
^^^^^^^^^^^^^^^^^^

The build can select the
`Python limited API <https://docs.python.org/3/c-api/stable.html>`_
by setting::

  -Csetup-args=-Dpython.allow_limited_api=true

This uses the minimum-supported Python version (currently Python 3.11)
and produces a wheel that can be installed with different Python
versions.

.. _build_group:

group
^^^^^

The ``build-group`` option determines whether the ``group`` module is
built, and defaults to ``true``. Set to ``false`` if the :term:`CIAO`
``group`` module is already installed, or if support for PHA grouping
is not required.

.. _build_stk:

stk
^^^

The ``build-stk`` option determines whether the ``stk`` module is
built, and defaults to ``true``. Set to ``false`` if the :term:`CIAO`
``stk`` module is already installed, or if support for stacks
(specifying multiple arguments via a file) is not required.

.. _build_region:

region
^^^^^^

Support for region files, as used in calls like ``notice2d``, is
provided via the region code. If ``build-region`` is ``true`` (its
default) then the region code will be built along with Sherpa. If the
setting is ``false`` but ``region-prefix`` is set, then the region
library from that location will be used (the ``region-use-cxc-parser``
option should be set to ``false`` if this is not a :term:`CIAO`
environment).

If neither option is set then Sherpa will not include support for
region files.

.. _build-wcssubs:

wcssubs
^^^^^^^

Support for :term:`WCS` information, used to allow the coordinate
setting to be changed when fitting image data, is provided by either
setting the ``build-wcssubs`` option to ``true``, its default, or by
setting the ``wcssubs-prefix`` option to point to an existing wcssubs
library.

If neither option is set then Sherpa will not include support for
coordinate conversions of image data.

.. _build-fftw:

FFTW
^^^^

Prior to Sherpa 4.19.0, Sherpa included a copy of the
`fftw library <https://www.fftw.org/>`_ source
code and would build it be default. This copy has now been removed
and access is expected to be provided by ``pkg-config``.

.. _build-xspec:

XSPEC
^^^^^

Sherpa can be built to use the Astronomy models provided by
:term:`XSPEC`. The ``xspec-prefix`` argument is set to the location of
the XSPEC include and library directories (expected to be
``$HEADAS``). The ``xspec-libraries`` argument can be used to specify
the link libraries for rare cases.

As an example::

    pip install . -Cxsetup-args=-Dxspec-prefix=$HEADAS

In order for the XSPEC module to be used from Python, the
``HEADAS`` environment variable **must** be set before the
:py:mod:`sherpa.astro.xspec` module is imported.

The Sherpa test suite includes an extensive set of tests of this
module, but a quick check of an installed version can be made with
the following command::

    % python -c 'from sherpa.astro import xspec; print(xspec.get_xsversion())'
    12.15.0

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
C++ compilers to use if the default values picked up by meson are not
correct.

In the following the `pip tool
<https://pip.pypa.io/en/stable/cli/pip_install/>`_ is used for
installing Sherpa, but any modern Python installation tool that
supports building Python extension modules should also work.

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

If any :ref:`options are set <configure>` then they are given with the::

  -Csetup_args=-D<name>=<value>

syntax, one per option. As an example, the following installation
call will build Sherpa with :ref:`support for XSPEC <build-xspec>`::

    pip install . -Csetup-args=-Dxspec-prefix=$HEADAS

.. _developer-build:

A development build
^^^^^^^^^^^^^^^^^^^

The code can be built locally, which is useful when adding new
functionality or fixing a bug. Using the ``--no-build-isolation``
flag means that the ``build-system`` requirements from the
``pyproject.toml`` need to be installed, which can be done with::

  pip install numpy ninja meson-python

and then the code can be built with the following (the ``[test]`` term
just ensures that ``pytest`` is also installed)::

  pip install -e .[test] --no-build-isolation

Any extension module will be automatically re-built when it is imported,
although it can be useful to re-run the ``pip install`` line if there is
a compilation error.

The ``--verbose`` flag is useful when diagnosing problems when building Sherpa::

  pip install -e .[test] --no-build-isolation --verbose

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

Building the documentation requires a Sherpa installation, *including*
the test data suite (either as a submodule or installed with the
``sherpa-test-data`` package), and several additional packages:

* `Sphinx <https://sphinx.pocoo.org/>`_, version 1.8 or later
* The ``sphinx_rtd_theme``
* NumPy and `sphinx-astropy <https://github.com/astropy/sphinx-astropy/>`_
  (the latter can be installed with ``pip``)
* `nbsphinx <https://pypi.org/project/nbsphinx/>`_, ``ipykernel``, and ``pandoc``
  for including Jupyter notebooks
* `Graphviz <https://www.graphviz.org/>`_ (for the inheritance diagrams)
* ``make``

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

.. _install-sherpa-test-data:

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
