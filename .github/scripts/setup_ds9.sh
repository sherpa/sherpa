#!/usr/bin/env bash -e

ds9_base_url=https://ds9.si.edu/download/

if [[ "x${CONDA_PREFIX}" == "x" ]];
then
    echo "Error: CONDA_PREFIX not set. This should be set for active Conda environments."
    exit 1
fi

if [ "`uname -s`" == "Darwin" ] ; then
    ds9_os=darwinbigsurx86
else
    echo "* installing dev environment"

    # install build dependencies
    sudo apt-get update
    sudo apt-get install -qq libx11-dev libsm-dev libxrender-dev

    # set os-specific variables
    ds9_os=ubuntu20
fi

echo "* ds9_os=$ds9_os"

download () {
  echo "* downloading $1"
  wget --quiet $1
  if [[ $? -ne 0 ]]; then
    echo "\n*** Unable to download $1\n"
  fi
}

# Tarballs to fetch
ds9_tar=ds9.${ds9_os}.8.5.tar.gz
xpa_tar=xpa.${ds9_os}.2.1.20.tar.gz

# Fetch them
download $ds9_base_url/$ds9_os/$ds9_tar
download $ds9_base_url/$ds9_os/$xpa_tar

# untar them; assume $CONDA_PREFIX/bin is in the path
echo "* unpacking ds9/XPA"
start_dir=$(pwd)
cd ${CONDA_PREFIX}/bin
tar xf ${start_dir}/${ds9_tar}
tar xf ${start_dir}/${xpa_tar}
cd -
