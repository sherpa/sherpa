#!/usr/bin/env bash -e

# Occasionally useful to know what these values are
echo "** uname -s: `uname -s`"
echo "** uname -m: `uname -m`"

if [ "`uname -s`" == "Darwin" ] ; then

    if [ "`uname -m`" == "x86_64" ]; then
	sys="64"
    else
	sys="arm64"
    fi
    compilers="clang_osx-${sys} clangxx_osx-${sys} gfortran_osx-${sys}"

    #Download the macOS 11.0 SDK to the CONDA_BUILD_SYSROOT location for the Conda Compilers to work
    mkdir -p ${GITHUB_WORKSPACE}/../11.0SDK
    wget https://github.com/phracker/MacOSX-SDKs/releases/download/11.3/MacOSX11.0.sdk.tar.xz -O MacOSX11.0.sdk.tar.xz
    if [[ $? -ne 0 ]]; then
      echo "macOS 11.0 SDK download failed"
    fi
    tar -C ${GITHUB_WORKSPACE}/../11.0SDK -xf MacOSX11.0.sdk.tar.xz

else
    compilers="gcc_linux-64=14.2 gxx_linux-64=14.2 gfortran_linux-64"

    if [ -n "${MATPLOTLIBVER}" ]; then
        #Installed for qt-main deps which is needed for the QtAgg backend to work for matplotlib
        #This is only an issue on stripped down systems. You can check for this issue by:
        #  ldd $CONDA_PREFIX/plugins/platforms/libqxcb.so | grep "not found"
        sudo apt-get install -q libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-render-util0
   fi
fi

# To avoid issues with non-XSPEC builds (e.g.
# https://github.com/sherpa/sherpa/pull/794#issuecomment-616570995 )
# the XSPEC-related channels are only added if needed
#
if [ -n "${XSPECVER}" ]; then
  conda config --add channels ${xspec_channel}
fi

# Figure out requested dependencies
if [ -n "${MATPLOTLIBVER}" ]; then MATPLOTLIB="matplotlib=${MATPLOTLIBVER}"; fi
if [ -n "${BOKEHVER}" ]; then BOKEH="bokeh=${BOKEHVER}"; fi
if [ -n "${NUMPYVER}" ]; then NUMPY="numpy=${NUMPYVER}"; fi
# Xspec >=12.10.1n Conda package includes wcslib & CCfits and pulls in cfitsio & fftw
if [ -n "${XSPECVER}" ];
 then export XSPEC="xspec-modelsonly=${XSPECVER}";
fi

echo "dependencies: ${MATPLOTLIB} ${BOKEH} ${NUMPY} ${XSPEC} ${FITSBUILD}"
echo "compilers:    ${compilers}"

# Create and activate conda build environment
# conda create --yes -n build python"=${PYTHONVER}.*=*cpython*" pip ${MATPLOTLIB} ${BOKEH} ${NUMPY} ${XSPEC} ${FITSBUILD} ${compilers}
conda create --yes -n build python"=${PYTHONVER}.*" pip ${MATPLOTLIB} ${BOKEH} ${NUMPY} ${XSPEC} ${FITSBUILD} ${compilers}

conda activate build
