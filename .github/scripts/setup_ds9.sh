#!/usr/bin/env bash -e

# The download/ version is the "latest" release whereas archive/ includes
# older releases.
# ds9_base_url=https://ds9.si.edu/download
ds9_base_url=https://ds9.si.edu/archive

if [[ "x${CONDA_PREFIX}" == "x" ]];
then
    echo "Error: CONDA_PREFIX not set. This should be set for active Conda environments."
    exit 1
fi

if [ "`uname -s`" == "Darwin" ] ; then
    if [ "`uname -m`" != "arm64" ] ; then
	echo "* DS9 testing not supported on `uname -s` : `uname -m`"
	exit 0
    fi
    ds9_os=darwinsonoma
    ds9_platform=arm64
else
    echo "* installing dev environment"

    # install build dependencies
    sudo apt-get update
    sudo apt-get install -qq libx11-dev libsm-dev libxrender-dev

    # set os-specific variables
    ds9_os=ubuntu24
    ds9_platform=x86
fi

echo "* ds9_os=$ds9_os"
echo "* ds9_platform=$ds9_platform"

ds9_label=${ds9_os}${ds9_platform}

download () {
  echo "* downloading $1"
  wget --quiet $1
  if [[ $? -ne 0 ]]; then
    echo "\n*** Unable to download $1\n"
  fi
}

# Tarballs to fetch
ds9_tar=ds9.${ds9_label}.8.7.tar.gz
xpa_tar=xpa.${ds9_label}.2.1.20.tar.gz

# Fetch them
download $ds9_base_url/$ds9_label/$ds9_tar
download $ds9_base_url/$ds9_label/$xpa_tar

# untar them; assume $CONDA_PREFIX/bin is in the path
echo "* unpacking ds9/XPA"
start_dir=$(pwd)
cd ${CONDA_PREFIX}/bin
tar xf ${start_dir}/${ds9_tar}
tar xf ${start_dir}/${xpa_tar}
cd -
