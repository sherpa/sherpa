#!/usr/bin/env bash -e

if [ "`uname -s`" == "Darwin" ] ; then
    compilers="clang_osx-64 clangxx_osx-64 gfortran_osx-64"

    #Download the macOS 11.0 SDK to the CONDA_BUILD_SYSROOT location for the Conda Compilers to work
    mkdir -p ${GITHUB_WORKSPACE}/11.0SDK
    wget https://github.com/phracker/MacOSX-SDKs/releases/download/11.3/MacOSX11.0.sdk.tar.xz -O MacOSX11.0.sdk.tar.xz
    if [[ $? -ne 0 ]]; then
      echo "macOS 11.0 SDK download failed"
    fi
    tar -C ${GITHUB_WORKSPACE}/11.0SDK -xf MacOSX11.0.sdk.tar.xz
    #End of Conda compilers section
else
    compilers="gcc_linux-64 gxx_linux-64 gfortran_linux-64"

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

echo "dependencies: ${MATPLOTLIB} ${BOKEH} ${NUMPY} ${FITS} ${XSPEC}"

# Create and activate conda build environment
conda create --yes -n build python"=${PYTHONVER}.*=*cpython*" pip ${MATPLOTLIB} ${BOKEH} ${NUMPY} ${XSPEC} ${FITSBUILD} ${compilers}

conda activate build
