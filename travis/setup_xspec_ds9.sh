#!/usr/bin/env bash -e

ds9_base_url=http://ds9.si.edu/download/

# variable $CONDA_PREFIX should be defined by conda by using conda activate (in setup_conda.sh)
if [[ "x${CONDA_PREFIX}" == "x" ]];
then
    echo "Error: CONDA_PREFIX not set. This should be set for active Conda environments."
    exit 1
fi

if [[ -z ${HEADAS} ]];
then
    echo "Error: HEADAS not set. Check that Xspec is installed in the Conda environment."
    exit 1
fi
echo "HEADAS=${HEADAS}"
#Set the xspec_root as the top of the Conda environment
xspec_root=${CONDA_PREFIX}

if [[ ${TRAVIS_OS_NAME} == linux ]];
then
    # install build dependencies
    sudo apt-get update
    sudo apt-get install -qq libx11-dev libsm-dev libxrender-dev

    # set os-specific variables
    ds9_os=ubuntu14
else  # osx
    ds9_os=darwinsierra

    # It looks like xvfb doesn't "just work" on osx travis, so...
    sudo Xvfb :99 -ac -screen 0 1024x768x8 &
fi

download () {
  wget --quiet $1
  if [[ $? -ne 0 ]]; then
    echo "\n*** Unable to download $1\n"
  fi
}

### DS9 and XPA
# Tarballs to fetch
ds9_tar=ds9.${ds9_os}.8.1.tar.gz
xpa_tar=xpa.${ds9_os}.2.1.20.tar.gz

# Fetch them
download $ds9_base_url/$ds9_os/$ds9_tar
download  $ds9_base_url/$ds9_os/$xpa_tar

# untar them
start_dir=$(pwd)
cd ${CONDA_PREFIX}/bin
tar xf ${start_dir}/${ds9_tar}
tar xf ${start_dir}/${xpa_tar}
cd -


### XSPEC
xspec_library_path=${xspec_root}/lib/
xspec_include_path=${xspec_root}/include/

case "${XSPECVER}" in
  12.10.1*)
      xspec_version_string="12.10.1"
      ;;
  *)
      echo "Xspec version listed currently unsuported in Travis jobs."
      exit 1
      ;;
esac

# Change build configuration
sed -i.orig "s/#with-xspec=True/with-xspec=True/g" setup.cfg
sed -i.orig "s|#xspec_lib_dirs = None|xspec_lib_dirs=${xspec_library_path}|g" setup.cfg
sed -i.orig "s|#xspec_include_dirs = None|xspec_include_dirs=${xspec_include_path}|g" setup.cfg
sed -i.orig "s|#xspec_version = 12.9.0|xspec_version = ${xspec_version_string}|g" setup.cfg

#Xspec ~12.10.1n now requires fftw. This disables the fftw build from extern
sed -i.orig "s|#fftw=local|fftw=local|g" setup.cfg
sed -i.orig "s|#fftw-include_dirs=build/include|fftw-include_dirs=${CONDA_PREFIX}/include|g" setup.cfg
sed -i.orig "s|#fftw-lib-dirs=build/lib|fftw-lib-dirs=${CONDA_PREFIX}/lib|g" setup.cfg
sed -i.orig "s|#fftw-libraries=fftw3|fftw-libraries=fftw3|g" setup.cfg
