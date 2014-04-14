# DISCLAIMER
We are currently in the process of moving the Sherpa codebase to GitHub, but be advised that this repository is
__not stable__, __not supported__, and __might change name and location at any time__ until we finalize the relocation of the code and all our tests pass.

While you are free to take a look, and maybe *star* this project to get updates, please do not assume that this repository is stable.

Thanks.

# Sherpa
Sherpa is the CIAO modeling and fitting application made available by the Chandra X-ray Center (CXC). It can be used for analysis of images, spectra and time series from many telescopes, including optical telescopes such as Hubble. Sherpa is flexible, modular and extensible. It has an IPython user interface and it is also an importable Python module. Sherpa models, optimization and statistic functions are available via both C++ and Python for software developers wishing to link such functions directly to their own compiled code.

Sherpa supports fitting of 1-D X-ray spectra from Chandra and other X-ray missions, as well as 1-D non-X-ray data, including ASCII data arrays, radial profiles, and lightcurves. The options for grating data analysis include fitting the spectrum with multiple response files required for overlapping orders in LETG observations. Modeling of 2-D spatial data is fully supported, including the PSF and exposure maps. User specified models can be added to Sherpa with advanced "user model" functionality.


## Get "raw" code
One can clone the git repository to get the whole repo

	git clone username@newdevel17:/data/scialg/staff/olaurino/sas/sherpa.git

Make sure you are on the standalone branch (temporary):

	git checkout feature/standalone

## Make distribution
For a standard distribution, just type:

        python setup.py sdist

## If you want to include the xspec extension in the source distribution:

	python setup.py xspec_config --with-xspec sdist

The dist directory will contain the tar ball with the distributable source code
for Sherpa.

## Default build
In order to make a default build of Sherpa one can simply type:

	pip install sherpa-version.tar.gz

or

	tar xf sherpa-version.tar.gz
	cd sherpa-version
	python setup.py install


## Configure
In order to enable or disable features, one can pass commandline arguments to
the special sherpa_config command, e.g.:

	python setup.py sherpa_config --fftw=local --fftw-lib-dirs=/usr/lib/ install

One can also change the configuration by editing the setup.cfg file, following
the instructions therein. Help on the sherpa_config command can be obtained by
typing

	python setup.py sherpa_config --help

Users can change the location and name of the fftw, region, and wcs libraries.

## XSPEC
In order to build sherpa with the support for XSPEC, one needs to provide the
location and names of the XSPEC libraries and of the other dependencies the
XSPEC extension needs to be linked to.

## Example setup.cfg file
Here is an example of customized setup.cfg file for building Sherpa with the
shipped C/C++ libraries, and locally-built XSPEC libraries:

~~~~
[sherpa_config]

# GROUP Python module
#disable-group=False

# FFTW Library
# Uncomment to use a local installation
#fftw=local

# If fftw=local uncomment the following lines and
# change the default location of libraries and the name
# of the library to be linked (usually fftw3)
# (include multiple values by separating them with spaces)
#fftw-include_dirs=build/include
#fftw-lib-dirs=build/lib
#fftw-libraries=fftw3

# Region Library
# Uncomment to use a local installation
#region=local

# If region=local uncomment the following lines and
# change the default location of libraries and headers and the name
# of the library to be linked (usually region)
# (include multiple values by separating them with spaces)
#region-include_dirs=build/include
#region-lib-dirs=build/lib
#region-libraries=region

# WCS Subroutines
# Uncomment and change default location if needed
#wcs-include-dirs=build/include
#wcs-lib-dirs=build/lib
#wcs-libraries=wcs

# XSPEC Models
[xspec_config]
# Uncomment (set to True) to build XSPEC extension
with-xspec=True

# If with-xspec is True, make sure to point Sherpa to right
# XSPEC-related libraries.
# Library names may probably be left to their default values,
# unless you have a very different setup
xspec_lib_dirs = /proj/port_ots/Linux64/xspec-modelsonly-12.8.0.k.1.Linux64/lib
xspec_libraries = XSUtil XSFunctions XSModel XS
cfitsio_lib_dirs = /proj/port_ots/Linux64/cfitsio-3.360.Linux64/lib/
cfitsio_libraries = cfitsio
ccfits_lib_dirs = /proj/port_ots/Linux64/cfitsio-3.360.Linux64/lib/
ccfits_libraries = CCfits
gfortran_lib_dirs = /proj/port_ots/Linux64/cfitsio-3.360.Linux64/lib/
gfortran_libraries = gfortran
~~~~
