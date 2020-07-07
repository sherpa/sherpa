#!/usr/bin/env bash -e

# Environment
sherpa_channel=sherpa
xspec_channel=xspec/channel/dev
miniconda=$HOME/miniconda

if [[ ${TRAVIS_OS_NAME} == linux ]];
then
    miniconda_os=Linux
    compilers="gcc_linux-64 gxx_linux-64 gfortran_linux-64"
else  # osx
    miniconda_os=MacOSX

    #Conda compilers - Comment out to test non-conda compilers
    #Unset the Travis compiler variables
    unset CC CFLAGS CXXFLAGS
    compilers="clang_osx-64 clangxx_osx-64 gfortran_osx-64"

    #Download and set the location of the macOS 10.9 SDK for the Conda Compilers to work
    mkdir -p 10.9SDK
    wget https://github.com/phracker/MacOSX-SDKs/releases/download/10.13/MacOSX10.9.sdk.tar.xz -O MacOSX10.9.sdk.tar.xz
    if [[ $? -ne 0 ]]; then
      echo "macOS 10.9 SDK download failed"
    fi
    tar -C 10.9SDK -xf MacOSX10.9.sdk.tar.xz
    export CONDA_BUILD_SYSROOT=$(pwd)/10.9SDK/MacOSX10.9.sdk
    echo "CONDA_BUILD_SYROOT=${CONDA_BUILD_SYSROOT}"
    #End of Conda compilers section
fi

# Download and install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-${miniconda_os}-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b -p $miniconda

#Source the Conda profile
source $miniconda/etc/profile.d/conda.sh

# update and add channels
conda update --yes conda

# To avoid issues with non-XSPEC builds (e.g.
# https://github.com/sherpa/sherpa/pull/794#issuecomment-616570995 )
# the XSPEC-related channels are only added if needed
#
conda config --add channels ${sherpa_channel}
if [ -n "${XSPECVER}" ]; then conda config --add channels ${xspec_channel}; fi

# Figure out requested dependencies
if [ -n "${MATPLOTLIBVER}" ]; then MATPLOTLIB="matplotlib=${MATPLOTLIBVER}"; fi
if [ -n "${NUMPYVER}" ]; then NUMPY="numpy=${NUMPYVER}"; fi
# Xspec >=12.10.1n Conda package includes wcslib & CCfits and pulls in cfitsio & fftw
if [ -n "${XSPECVER}" ];
 then export XSPEC="xspec-modelsonly=${XSPECVER}";
fi
if [ "${DOCS}" == true ];
 then export DOCSBUILD="sphinx sphinx_rtd_theme graphviz";
fi

echo "dependencies: ${MATPLOTLIB} ${NUMPY} ${FITS} ${XSPEC} ${DOCSBUILD}"

# Tests fail with AstroPy 3.2 but not 3.2.1. Since 3.2.1 is now available
# on conda, we no longer force AstroPy < 3.2 - see
# https://github.com/sherpa/sherpa/issues/632
#
FITSBUILD="${FITS}"

# Create and activate conda build environment
# We create a new environment so we don't care about the python version in the root environment.
conda create --yes -n build python=${TRAVIS_PYTHON_VERSION} pip ${MATPLOTLIB} ${NUMPY} ${XSPEC} ${FITSBUILD} ${DOCSBUILD} ${compilers}

conda activate build

# It looks like on some systems (well, linux) the F90 variable needs to be set. Not sure why
export F90=${F77}

# This is required to make sure that the CIAO python extensions being built pick the correct flags
export PYTHON_LDFLAGS=" "

# Create a customized configuration file where the number of CPU cores
# can be fixed: on VMs like on Travis the CPU detection logic we use
# picks up the total number of CPUS on the host machine, not the
# virtual machine we are running, as described at
# https://github.com/travis-ci/travis-ci/issues/4696
#
# TRAVIS_NUMCORES is from https://github.com/travis-ci/travis-build/pull/1079
# but doesn't seem to exist. I am hard-coding to 2 based on
# https://docs.travis-ci.com/user/reference/overview/#virtualisation-environment-vs-operating-system
#
ncores=`python -c 'import multiprocessing; print(multiprocessing.cpu_count());'`
echo "# Python thinks we have ${ncores} cores"
echo "# Forcing the number of cores to be 2"
sed -e "s/^numcore : None/numcore : 2/" sherpa/sherpa-standalone.rc > .sherpa-standalone.rc
grep numcore .sherpa-standalone.rc  # DBG
