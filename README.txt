Building and Installing Sherpa

Sherpa 4.4.0
=======================

    * Required dependencies:
        o Python-2.7.X (preferred) or Python-2.6.X
        o Numpy-1.5.1 or later

    * Optional dependencies:
        o DS9 5.6 or later (SAO imager)
        o XPA 2.1.9 or later (for DS9 communication)
        o PyFITS 1.3 or later (for FITS I/O)
        o matplotlib 0.98.5.2 or later (for plotting and imaging,
          be sure to set "interactive=True" in ~/.matplotlib/matplotlibrc)
        o IPython 0.10

Dependencies
============

The new Sherpa can be downloaded, built, installed, and used independently of
CIAO. If you use Sherpa this way, you won't have access to CRATES I/O, ChIPS
plotting, 2D region, or dynamic grouping but all other functionality will be
available. The current (CIAO 4.4) source tarball is sherpa-4.4.0.tar.gz. It
has the following prerequisites:

    * Required:
          o GCC compilers with g77/gfortran/g95 3.4.5 or later
          o or Sun Studio with f77/f95 WS 10 or later
          o or Xcode with g77/gfortran/g95 3.1.2 or later
          o Python 2.7.X or Python 2.6.X
          o NumPy 1.5.1 or later
          o FFTW3 3.3 or later (be sure to build with "--enable-float" also)
    * Optional:
          o DS9 5.6 or later (SAO imager)
          o XPA 2.1.9 or later (for DS9 communication)
          o XSPEC 12.7.0 (for the XSPEC model functions)
	  o CCfits 2.3 (for use in XSPEC models)
	  o CFITSIO 3.27 (for use in XSPEC models)
          o PyFITS 1.3 or later (for non-DM FITS I/O)
          o matplotlib 0.98.5.2 or later (for plotting and imaging,
	    be sure to set "interactive=True" in ~/.matplotlib/matplotlibrc)
          o IPython (for a nicer interactive interface; this is distributed 
            with CIAO)

Build and Install
=================

To build and install the package, do the following:

$ tar xzf sherpa-4.4.0.tar.gz
$ cd sherpa-4.4.0
$ python setup.py [config-vars] install --prefix=<dest-dir>

config-vars is an optional list of arguments in the format var=value that 
specifies where to find prerequisites required at build time. The following 
variables can be set:

Variable 	     Default
========             =======
fftw_library_dir     /usr/local/lib (only needed if in alternate location)
fftw_include_dir     /usr/local/include (only needed if in alternate location)
wcs_library_dir      None (if not given, World Coordinate System, WCS, module
                           is not built)
wcs_include_dir      None (if not given, World Coordinate System, WCS, module
                           is not built)
reg_library_dir      None (if not given, Region 2D filtering module is not
                           built)
reg_include_dir      None (if not given, Region 2D filtering module is not
                           built)
fortran_library_dir  None (may be needed to find libg2c.so or libgfortran.so or
                           libf95.a depending on how XSPEC was built)
fortran_lib          None ('g2c', 'gfortran', or 'f95' depending on fcompiler )

xspec_library_dir    None (if not given, XSPEC module is not built)
cfitsio_lib 	     cfitsio  ('cfitsio', or 'cfitsio_3.27' depending on
		     version)
cfitsio_library_dir  <xspec_library_dir>

For example, to use the FFTW in /soft/fftw and the XSPEC library in
/opt/local/headas/lib, you would do

$ python setup.py \
     fftw_library_dir=/soft/fftw/lib \
     fftw_include_dir=/soft/fftw/include \
     xspec_library_dir=/opt/local/headas/lib \
     cfitsio_library_dir=/opt/local/headas/lib \
     install ...

The setup.py script distributed with Sherpa uses the standard Python distutils
package. For more information on using it, see Installing Python Modules.

Configuration Files
===================

Sherpa comes with a configuration file "sherpa.rc". For stand-alone 
installation, this file should be copied from source in to the home directory 
as "~/.sherpa.rc". Be sure to indicate the IO and Plotting backends as "pyfits"
and "pylab" depending on configuration.

Matplotlib comes with a configuration file "matplotlibrc". For smooth behavior
with Sherpa, be sure to indicate "interactive=True" in
~/.matplotlib/matplotlibrc.

Environment Variables
=====================

If Sherpa is installed in a non-standard location, be sure to update the 
following variables accordingly.

    PYTHONPATH
