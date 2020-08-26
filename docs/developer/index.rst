**********************************
Contributing to Sherpa development
**********************************

.. todo::

   Needs work.

Contributions to Sherpa - whether it be bug reports, documentation
updates, or new code - are highly encouraged.  Please `report any
problems or feature requests on github
<https://github.com/sherpa/sherpa/issues/>`_.

At present we do not have any explicit documentation on how
to contribute to Sherpa, but it is similar to other open-source
packages such as
`AstroPy <https://docs.astropy.org/en/stable/index.html#contributing>`_.

The developer documentation is also currently lacking.

To do code development, Sherpa needs to be installed from source so
that tests can run locally and the documentation can be build locally
to test out any additions to code or docs.  :ref:`build-from-source`
describes several ways to build Sherpa from source, but one
particularly comfortable way is described in detail in the next
section.

.. _source-install-with-conda:

Install from source in conda
============================

Conda can be used to install all the dependencies for Sherpa.

::

    conda create -n sherpaciao -c https://cxc.cfa.harvard.edu/conda/ciao ds9 astropy ciao
    conda install -n sherpaciao --only-deps -c https://cxc.cfa.harvard.edu/conda/ciao sherpa
    conda install -n sherpaciao -c anaconda -c astropy sphinx graphviz sphinx-astropy sphinx_rtd_theme

The first line installes the full `CIAO release
<https://cxc.harvard.edu/ciao/>`_ and astropy, required for building
and running tests locally. The last line adds all requirements for
building the documentation.  Sherpa can use either astropy or crates
as backend for reading and writing files. The default configuration in
Sherpa is to use astropy. However, if crates is installed (e.g. by
installing the `ciao` package) and selected as backend in `sherpa.rc`,
then astropy can be omitted from the install (but is still needed to
build the docs).

As described in :ref:`build-from-source`, the file ``setup.cfg`` in
the root directory of the sherpa source needs to be modified to
configure the build. This is particularly easy in this setup, where
all external dependencies are installed in conda and the enviroment
variable ``ASCDS_LIB`` is set to the include directory, when the conda
environment is activated. Thus, all that is needed is to disable the
build of external dependencies and to set directories. The following
lists the lines in ``setup.cfg`` that need to be modified (adjust
xspec version as needed)::

    # GROUP Python module
    disable-group=True

    # File Stack Python module
    disable-stk=True

    # FFTW Library
    fftw=local
    fftw-include_dirs=${ASCDS_LIB}/../include
    fftw-lib-dirs=${ASCDS_LIB}
    fftw-libraries=fftw3

    # Region Library
    region=local
    region-include_dirs=${ASCDS_LIB}/../include
    region-lib-dirs=${ASCDS_LIB}
    region-libraries=region ascdm

    # WCS Subroutines
    wcs=local
    wcs-include-dirs=${ASCDS_LIB}/../include
    wcs-lib-dirs=${ASCDS_LIB}
    wcs-libraries=wcs

    # XSPEC Models
    [xspec_config]
    with-xspec=True
    xspec_version = 12.10.1
    xspec_lib_dirs = ${ASCDS_LIB}
    xspec_include_dirs = ${ASCDS_LIB}/../include

To avoid accidentially commiting the modified ``setup.cfg`` into git,
the file can be marked as "assumed unchanged".

::

    git update-index --assume-unchanged setup.cfg

After these steps, the conda enviroment (here called ``sherpaciao``)
can be activated and Sherpa can be build from source.

::

    conda activate sherpaciao
    python setup.py develop


.. warning::

   Just like in the case of a normal source install, when building Sherpa
   on recent versions of macOS within a conda environment, the following
   environment variable must be set::

     export PYTHON_LDFLAGS=' '

   That is, the variable is set to a space, not the empty string.


How do I ...
============

Update the XSPEC bindings?
--------------------------

The :py:mod:`sherpa.astro.xspec` module currently supports
:term:`XSPEC` versions 12.11.0 down to 12.9.0. It may build against
newer versions, but if it does it will not provide access
to any new models in the release. The following steps are needed
to update to a newer version, and assume that you have the new version
of XSPEC, or its model library, available. The following
sections of the
`XSPEC manual <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XspecManual.html>`__
should be reviewed:
"Appendix F: Using the XSPEC Models Library in Other Programs",
and
"Appendix C: Adding Models to XSPEC"
(direct links are not provided as there are no obvious stable URIs for
them).

#. Add a new version define in ``helpers/xspec_config.py``.

   Current version: `helpers/xspec_config.py <https://github.com/sherpa/sherpa/blob/master/helpers/xspec_config.py>`_.

   When adding support for XSPEC 12.11.0, the code in the ``run``
   method was changed to include::

       if xspec_version >= LooseVersion("12.11.0"):
           macros += [('XSPEC_12_11_0', None)]

   and the version check to::

       # Since there are patches (e.g. 12.10.0c), look for the
       # "next highest version.
       if xspec_version >= LooseVersion("12.11.1"):
           self.warn("XSPEC Version is greater than 12.11.0, which is the latest supported version for Sherpa")

   The define should be named ``XSPEC_<a>_<b>_<c>`` for XSPEC release
   ``<a>.<b>.<c>`` (the XSPEC patch level is not included). This define
   is used when compiling the XSPEC model interface, to select which
   functions to include.

   .. note:: The Sherpa build system requires that the user indicate the
	     version of XSPEC being used, via the ``xspec_config.xspec_version``
	     setting in their ``setup.cfg`` file (as attempts to identify
	     this value automatically were not successful). This version is
	     the value used in the checks in ``helpers/xspec_config.py``.

#. Attempt to build the XSPEC interface with::

     python setup.py develop

   This requires that the ``xspec_config`` section of the ``setup.cfg``
   file has been set up correctly for the new XSPEC release. The exact
   settings depend on how XSPEC was built (e.g. model only or as a
   full application), and are described in the
   :ref:`building XSPEC <build-xspec>` documentation. The most-common
   changes are that the version numbers of the ``CCfits``, ``wcslib``,
   and ``hdsp`` libraries need updating, and these can be checked by
   looking in ``$HEADAS/lib``.

   If the build succeeds, you can check that it has worked by directly
   importing the XSPEC module, such as with the following, which should
   print out the correct version::

     python -c 'from sherpa.astro import xspec; print(xspec.get_xsversion())'

   It may however fail, due to changes in the XSPEC interface (unfortunately,
   such changes are often not included in the release notes).

#. Identify changes in the XSPEC models.

   A new XSPEC release can add models, change parameter settings in
   existing models, change how a model is called, or even delete a
   model (the last case is rare, and may require a discussion on
   how to proceed). The
   `XSPEC release notes <https://heasarc.gsfc.nasa.gov/xanadu/xspec/CHANGELOG.txt>`_
   page provides an overview, but the ``model.dat`` file - found
   in ``headas-<version>/Xspec/src/manager/model.dat`` (build) or
   ``$HEADAS/../spectral/manager/model.dat`` (install) - provides
   the details. It greatly simplifies things if you have a copy of
   this file from the previous XSPEC version, since then a command
   like::

     diff heasoft-6.26.1/spectral/manager/model.dat heasoft-6.27/spectral/manager/model.dat

   will tell you the differences. If you do not have the previous
   version then the release notes will tell you which models to
   look for in the ``model.dat`` file.

   The ``model.dat`` is an ASCII file which is described in
   Appendix C: Adding Models to XSPEC of the
   `XSPEC manual <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XspecManual.html>`_.
   The Sherpa interface to XSPEC only supports models labelled
   as ``add``, ``mul``, and ``con`` (additive, multiplicative,
   and convolution, respectively).

   Each model is represented by a set of consecutive lines in
   the file, and as of XSPEC 12.11.0, the file begins with::

     % head -5 heasoft-6.27/Xspec/src/manager/model.dat
     agauss         2   0.         1.e20          C_agauss  add  0
     LineE   A      10.0   0.      0.      1.e6      1.e6      0.01
     Sigma   A      1.0    0.      0.      1.e6      1.e6      0.01

     agnsed        15   0.03       1.e20          agnsed    add  0

   The important parts of the model definition are the first line,
   which give the XSPEC model name (first parameter), number of
   parameters (second parameter), two numbers which we ignore, the
   name of the function that evaluates the model, the type
   (e.g. ``add``), and then 1 or more values which we ignore. Then
   there are lines which define the model parameters (the number match
   the second argument of the first line), and then one or more blank
   lines. In the output above we see that the XSPEC ``agauss`` model
   has 2 parameters, is an additive model provided by the ``C_agauss``
   function, and that the parameters are ``LineE`` and ``Sigma``.
   The ``agnsed`` model is then defined (which uses the ``agnsed``
   routines), but its 15 parameters have been cut off from the output.

   The parameter lines will mostly look like this: parameter name,
   unit string (is often ``" "``), the default value, the hard and then
   soft minimum, then the soft ahd hard maximum, and then a value used
   by the XSPEC optimiser, but we only care about if it is negative
   (which indicates that the parameter should be frozen by default).
   The other common variant is the "flag" parameter - that is, a
   parameter that should never be thawed in a fit - which is indicated
   by starting the parameter name with a ``$`` symbol (although the
   documentation says these should only be followed by a single value,
   you'll see a variety of formats in the ``model.dat`` file). These
   parameters are marked by setting the ``alwaysfrozen`` argument of
   the :py:class:`~sherpa.models.parameter.Parameter` constructor
   to ``True``. Another option is the "scale" parameter, which is
   labelled with a ``*`` prefix, and these are treated as normal
   parameter values.

   a. ``sherpa/astro/xspec/src/_xspec.cc``

      Current version: `sherpa/astro/xspec/src/_xspec.cc <https://github.com/sherpa/sherpa/blob/master/sherpa/astro/xspec/src/_xspec.cc>`_.

      New functions are added to the ``XspecMethods`` array,
      using macros defined in ``sherpa/include/sherpa/astro/xspec_extension.hh``,
      and should be surrounded by a pre-processor check for the
      version symbol added to ``helpers/xspec_config.py``.

      As an example::

        #ifdef XSPEC_12_10_1
          XSPECMODELFCT_NORM( agnsed, 16 ),
        #endif

      adds support for the ``agnsed`` function, but only for XSPEC
      12.10.1 and later. Note that the symbol name used here is
      **not** the XSPEC model name (the first argument of the model
      definition from ``model.dat``), but the function name (the fifth
      argument of the model definition (although for the ``agnsed``
      example they are the same).

      Some models have changed the name of the function over time,
      so the pre-processor directive may need to be more complex, such
      as::

        #ifdef XSPEC_12_10_0
          XSPECMODELFCT_C_NORM( C_nsmaxg, 6 ),
        #else
          XSPECMODELFCT_NORM( nsmaxg, 6 ),
        #endif

      The remaining pieces are the choice of macro
      (e.g. ``XSPECMODELFCT_NORM`` or ``XSPECMODELFCT_C_NORM``) and
      the value for the second argument.  The macro depends on the
      model type and the name of the function (which defines the
      interface that XSPEC provides for the model, such as single- or
      double- precision, and Fortran- or C- style linking). Additive
      models use the suffix ``_NORM`` and convolution models use the
      suffix ``_CON``. Model functions which begin with ``C_`` use the
      ``_C`` variant, while those which begin with ``c_`` currently
      require treating them as if they have no prefix.

      The numeric argument to the template defines the number of
      parameters supported by the model once in Sherpa, and should
      equal the value given in the ``model.dat`` file for
      multiplicative and convolution style models, and one larger than
      this for additive models (i.e. those which use a macro that ends
      in ``_NORM``).

      As an example, the following three models from ``model.dat``::

        apec           3  0.         1.e20           C_apec    add  0
        phabs          1  0.03       1.e20           xsphab    mul  0
        gsmooth        2  0.         1.e20           C_gsmooth    con  0

      are encoded as (ignoring any pre-processor directives)::

        XSPECMODELFCT_C_NORM( C_apec, 4 ),
        XSPECMODELFCT( xsphab, 1 ),
        XSPECMODELFCT_CON(C_gsmooth, 2),

      Those models that do not use the ``_C`` version of the macro (or,
      for convolution-style models, have to use
      ``XSPECMODELFCT_CON_F77``), also have to declare the function
      within the ``extern "C" {}`` block. For FORTRAN models, the
      declaration should look like (replacing ``func`` with the
      function name, and note the trailing underscore)::

        void func_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

      and for model functions called ``c_func``, the prefixless
      version should be declared as::

        void func(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

      If you are unsure, do not add a declaration and then try to
      build Sherpa: the compiler should fail with an indication of
      what symbol names are missing.

      .. note:: Ideally we would have a sensible ordering for the declarations in this
		file, but at present it is ad-hoc.

   b. ``sherpa/astro/xspec/__init__.py``

      Current version: `sherpa/astro/xspec/__init__.py <https://github.com/sherpa/sherpa/blob/master/sherpa/astro/xspec/__init__.py>`_.

      This is where the Python classes are added for additive and multiplicative
      models. The code additions are defined by the model and parameter
      specifications from the ``model.dat`` file, and the existing classes
      should be used for inspiration. The model class should be called
      ``XS<name>``, where ``<name>`` is the XSPEC model name, and the
      ``name`` argument to its constructor be set to the XSPEC model name.

      The two main issues are:

      * Documentation: there is no machine-readable version of the text, and
	so the documentation for the XSPEC model is used for inspiration.

        The idea is to provide minimal documentation, such as the
	model name and parameter descriptions, and then to point users to
	the XSPEC model page for more information.

	One wrinkle is that the
	`XSPEC manual <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/manual/manual.html>`__
	does not provide a stable URI for a model (as it can change with XSPEC
	version). However, it appears that you can use the following pattern:

	  https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodel<Name>.html

	where ``<Name>`` is the capitalised version of the model name (e.g.
	``Agnsed``).

      * Models that are not in older versions of XSPEC should be marked with
	the ``version_at_least`` decorator (giving it the minimum supported
	XSPEC version as a string), and the function (added to ``_xspec.cc``)
	is specified as a string using the ``__function__`` attribute. The
	:py:class:`sherpa.astro.xspec.utils.ModelMeta` metaclass performs
	a runtime check to ensure that the model can be used.

   c. ``sherpa/astro/xspec/tests/test_xspec.py``

      Current version: `sherpa/astro/xspec/tests/test_xspec.py <https://github.com/sherpa/sherpa/blob/master/sherpa/astro/xspec/tests/test_xspec.py>`_.

      The ``XSPEC_MODELS_COUNT`` version should be increased by the number
      of models classes added to ``__init__.py``.

      Additive and multiplicative models will be run as part of the test
      suite - using a simple test which runs on a default grid and
      uses the default parameter values - whereas convolution models
      are not (since their pre-conditions are harder to set up
      automatically).

   d. ``docs/model_classes/astro_xspec.rst``

      Current version: `docs/model_classes/astro_xspec.rst <https://github.com/sherpa/sherpa/blob/master/docs/model_classes/astro_xspec.rst>`_.

      New models should be added to both the ``Classes`` rubric - sorted
      by addtive and then multiplicative models, using an alphabetical
      sorting - and to the appropriate ``inheritance-diagram`` rule.

Notes
=====

Notes on the design and changes to Sherpa.

.. _model_dimensions:

The dimensionality of models
----------------------------

Originally the Sherpa model class did not enforce any requirement on
the models, so it was possible to combine 1D and 2D models, even though
the results are unlikely to make sense. With the start of the regrid
support, added in `PR #469 <https://github.com/sherpa/sherpa/pull/469>`_,
the class hierarchy included 1D- and 2D- specific classes, but there
was still no check on model expressions. This section describes the
current way that models are checked:

- the :py:class:`sherpa.models.model.Model` class defines a
  :py:attr:`sherpa.models.model.Model.ndim` attribute, which is set
  to ``None`` by default.
- the :py:class:`sherpa.models.model.RegriddableModel1D` and
  :py:class:`sherpa.models.model.RegriddableModel2D` classes set
  this attribute to 1 or 2, respectively (most user-callable classes
  are derived from one of these two classes).
- the :py:class:`sherpa.models.model.CompositeModel` class checks
  the ``ndim`` attribute for the components it is given (the
  ``parts`` argument) and checks that they all have the same
  ``ndim`` value (ignoring those models whose dimensionality
  is set to ``None``). If there is a mis-match then a
  :py:class:`sherpa.utils.err.ModelErr` is raised.

An alternative approach would have been to introdude 1D and 2D
specific classes, from which all models derive, and then require the
parent classes to match. This was not attempted as it would require
significantly-larger changes to Sherpa (but this change could still be
made in the future).

.. _model_combination:

Combining model expressions
---------------------------

Models can be combined in several ways (for models derived from the
:py:class:`sherpa.models.model.ArithmeticModel` class):

- a unary operator, taking advantage of the ``__neg__`` and
  ``__abs__`` special methods of a class;
- a binary operator, using the ``__add__``, ``__sub__``, ``__mul__``,
  ``__div__``, ``__floordiv__``, ``__truediv__``, ``__mod__`` and ``__pow__``
  methods.

This allows models such as::

    sherpa.models.basic.Polynom1D('continuum') + sherpa.models.basic.Gauss1D('line')

to be created, and relies on the :py:class:`sherpa.models.model.UnaryOpModel`
and :py:class:`sherpa.models.model.BinaryOpModel` classes.

The :py:class:`~sherpa.models.model.BinaryOpModel` class has special-case handling
for values that are not a model expression (i.e. that do not derive
from the :py:class:`~sherpa.models.model.ArithmeticModel` class),
such as::

    32424.43 * sherpa.astro.xspec.XSpowerlaw('pl')

In this case the term ``32424.43`` is converted to an
:py:class:`~sherpa.models.model.ArithmeticConstantModel` instance and then
combined with the remaining model instance (``XSpowerlaw``).

For those models that require the
full set of elements, such as multiplication by a :term:RMF or a convolution
kernel, requires creating a model that can "wrap" another
model. The wrapping model will evaluate the wrapped model on
the requested grid, and then apply any modifications.
Examples include the
:py:class:`sherpa.instrument.PSFModel` class,
which creats :py:class:`sherpa.instrument.ConvolutionModel`
instances,
and the :py:class:`sherpa.astro.xspec.XSConvolutionKernel`
class, which creates :py:class:`sherpa.astro.xspec.XSConvolutionModel`
instances.
