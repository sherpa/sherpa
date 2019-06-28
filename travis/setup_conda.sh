#!/usr/bin/env bash -e

# Environment
libgfortranver="3.0"
sherpa_channel=sherpa
xspec_channel=xspec/channel/dev
miniconda=$HOME/miniconda

if [[ ${TRAVIS_OS_NAME} == linux ]];
then
    miniconda_os=Linux
    compilers="gcc_linux-64 gxx_linux-64 gfortran_linux-64"
else  # osx
    miniconda_os=MacOSX
    compilers="clang_osx-64 clangxx_osx-64 gfortran_osx-64"

    # On macOS we also need the conda libx11 libraries used to build xspec
    # We also need to pin down ncurses, for now only on macos.
    xorg="xorg-libx11 ncurses=5"
fi

# Download and install conda
wget http://repo.continuum.io/miniconda/Miniconda3-latest-${miniconda_os}-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b -p $miniconda
export PATH=$miniconda/bin:$PATH

# update and add channels
conda update --yes conda

# Note the order of the channels matter. We built the xspec conda packages for macos using the conda-forge channel
# with the highest priority, so we add it with a higher priority than the default channels, but with less priority
# than our own channels, so that we don't accidentally get conda-forge packages for cfitsio/ccfits.
if [[ ${TRAVIS_OS_NAME} == osx ]];
then
    conda config --add channels conda-forge
fi

conda config --add channels ${sherpa_channel}
conda config --add channels ${xspec_channel}
conda config --add channels anaconda

# Figure out requested dependencies
if [ -n "${MATPLOTLIBVER}" ]; then MATPLOTLIB="matplotlib=${MATPLOTLIBVER}"; fi
if [ -n "${NUMPYVER}" ]; then NUMPY="numpy=${NUMPYVER}"; fi
if [ -n "${XSPECVER}" ];
 then export XSPEC="xspec-modelsonly=${XSPECVER} ${xorg}";
fi
if [ "${DOCS}" == true ];
 then export DOCSBUILD="sphinx sphinx_rtd_theme graphviz";
fi

echo "dependencies: ${MATPLOTLIB} ${NUMPY} ${FITS} ${XSPEC} ${DOCSBUILD}"

# Hack to force AstroPy < 3.2 (see https://github.com/sherpa/sherpa/issues/632)
# This is a short-term hack just to get travis tests to pass
# Should there be better a way of specifying versions for these dependencies?
#
if [ "${FITS}" == astropy ]; then FITSBUILD="'astropy<3.2'"; else FITSBUILD="${FITS}"; fi

# Create and activate conda build environment
# We create a new environment so we don't care about the python version in the root environment.
conda create --yes -n build python=${TRAVIS_PYTHON_VERSION} pip ${MATPLOTLIB} ${NUMPY} ${XSPEC} ${FITSBUILD} ${DOCSBUILD} ${compilers}\
  libgfortran=${libgfortranver}

source activate build

# It looks like on some systems (well, linux) the F90 variable needs to be set. Not sure why
export F90=${F77}

# This is required to make sure that the CIAO python extensions being built pick the correct flags
export PYTHON_LDFLAGS=" "
