==================================
Contributing to Sherpa development
==================================

.. todo::

   Needs work.
   
Contributions to Sherpa - whether it be bug reports, documentation
updates, or new code - are highly encouraged.  Please `report any
problems or feature requests on github
<https://github.com/sherpa/sherpa/issues/>`_.


At present we do not have any explicit documentation on how
to contribute to Sherpa, but it is similar to other open-source
packages such as
`AstroPy <http://docs.astropy.org/en/stable/index.html#contributing>`_.

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
<http://cxc.harvard.edu/ciao/>`_ and astropy, required for building
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
