#!/usr/bin/env bash

# Tarballs to fetch
DS9_SITE=http://ds9.si.edu/download/centos6/;
XPA_SITE=http://ds9.si.edu/download/centos6/;
DS9_TAR=ds9.centos6.7.5.tar.gz;
XPA_TAR=xpa.centos6.2.1.17.tar.gz;

# Fetch them
wget --quiet $DS9_SITE$DS9_TAR;
wget --quiet $XPA_SITE$XPA_TAR;

# untar them
THIS_DIR=$(pwd);
cd $MINICONDA/bin; tar xf $THIS_DIR/$DS9_TAR; tar xf $THIS_DIR/$XPA_TAR; cd -;

# install build dependencies
sudo apt-get update
sudo apt-get install -qq libwcs4 wcslib-dev libx11-dev libsm-dev libxrender-dev;

# Some relevant environment variables
export HEADAS=$MINICONDA/envs/build/Xspec/spectral;
export XSPEC_LIBRARY_PATH=$MINICONDA/envs/build/Xspec/x86_64-unknown-linux-gnu-libc2.15-0/lib/;
export XSPEC_INCLUDE_PATH=$MINICONDA/envs/build/Xspec/x86_64-unknown-linux-gnu-libc2.15-0/include/;
export WCS_DIR_PATH=/usr/lib/x86_64-linux-gnu/

# Change build configuration
sed -i.orig s/#with-xspec=True/with-xspec=True/g setup.cfg;
sed -i.orig "s|#xspec_lib_dirs = None|xspec_lib_dirs=${XSPEC_LIBRARY_PATH}|g" setup.cfg;
sed -i.orig "s|#xspec_include_dirs = None|xspec_include_dirs=${XSPEC_INCLUDE_PATH}|g" setup.cfg;
sed -i.orig "s|#wcslib_lib_dirs = None|wcslib_lib_dirs=${WCS_DIR_PATH}|g" setup.cfg;
sed -i.orig "s|#gfortran_libraries = gfortran|gfortran_libraries= :libgfortran.so.3|g" setup.cfg;
