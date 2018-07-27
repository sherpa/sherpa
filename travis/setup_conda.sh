#!/usr/bin/env bash -e

# Environment
libgfortranver="3.0"
sherpa_channel=sherpa
xspec_channel=cxc/channel/dev
miniconda=$HOME/miniconda

if [[ ${TRAVIS_OS_NAME} == linux ]];
then
    miniconda_os=Linux
    compilers="gcc_linux-64 gxx_linux-64 gfortran_linux-64"
    sed -i.orig "s|#extra-fortran-link-flags=|extra-fortran-link-flags=-shared|" setup.cfg
else  # osx
    miniconda_os=MacOSX
    compilers="clang_osx-64 clangxx_osx-64 gfortran_osx-64"
    # This is required on macOS when building with conda
    sed -i.orig "s|#extra-fortran-link-flags=|extra-fortran-link-flags=-undefined dynamic_lookup -bundle|" setup.cfg

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

# Figure out requested dependencies
if [ -n "${MATPLOTLIBVER}" ]; then MATPLOTLIB="matplotlib=${MATPLOTLIBVER}"; fi
if [ -n "${NUMPYVER}" ]; then NUMPY="numpy=${NUMPYVER}"; fi
if [ -n "${XSPECVER}" ];
 then export XSPEC="xspec-modelsonly=${XSPECVER} ${xorg}";
fi
echo "dependencies: ${MATPLOTLIB} ${NUMPY} ${FITS} ${XSPEC}"

# Create and activate conda build environment
# We create a new environment so we don't care about the python version in the root environment.
conda create --yes -n build python=${TRAVIS_PYTHON_VERSION} pip ${MATPLOTLIB} ${NUMPY} ${XSPEC} ${FITS} ${compilers}\
  libgfortran=${libgfortranver}

source activate build

# It looks like on some systems (well, linux) the F90 variable needs to be set. Not sure why
export F90=${F77}

# This is required to make sure that the CIAO python extensions being built pick the correct flags
export PYTHON_LDFLAGS=" "
