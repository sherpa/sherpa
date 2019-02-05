#!/usr/bin/env bash -e

# variable $miniconda is defined in a previous script (setup_conda.sh)

ds9_base_url=http://ds9.si.edu/download/

if [[ ${TRAVIS_OS_NAME} == linux ]];
then
    # install build dependencies
    sudo apt-get update
    sudo apt-get install -qq libwcs4 wcslib-dev libx11-dev libsm-dev libxrender-dev

    # set os-specific variables
    ds9_os=ubuntu14
    xspec_root=$miniconda/envs/build/Xspec/x86_64-unknown-linux-gnu-libc2.15-0

    # Newer versions (or conda builds) of numpy bring in a libgfortran-ng dependency, which our
    # xspec packages cannot link against, so we need to pass a specific library name to the linker.
    libgfortran_name=":libgfortran.so.3"

    wcs_library_path=/usr/lib/x86_64-linux-gnu/
    sed -i.orig "s|#wcslib_lib_dirs = None|wcslib_lib_dirs=${wcs_library_path}|g" setup.cfg

else  # osx
    ds9_os=darwinsierra
    xspec_root=$miniconda/envs/build/Xspec/x86_64-apple-darwin16.3.0
    export DYLD_LIBRARY_PATH=${xspec_root}/lib

    # It looks like on macs numpy does not bring in the dependency with libgfortran-ng,
    # so we can link as usual. Also, on macOS the :libgfortran.3.dylib syntax wouldn't work.
    libgfortran_name="gfortran"

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
ds9_tar=ds9.${ds9_os}.8.0.tar.gz
xpa_tar=xpa.${ds9_os}.2.1.18.tar.gz

# Fetch them
download $ds9_base_url/$ds9_os/$ds9_tar
download  $ds9_base_url/$ds9_os/$xpa_tar

# untar them
start_dir=$(pwd)
cd ${miniconda}/bin
tar xf ${start_dir}/${ds9_tar}
tar xf ${start_dir}/${xpa_tar}
cd -


### XSPEC
xspec_library_path=${xspec_root}/lib/
xspec_include_path=${xspec_root}/include/

case "${XSPECVER}" in
  12.10.1b)
      xspec_version_string="12.10.1"
      xspec_include_path="$miniconda/envs/build/include"
      xspec_library_path="$miniconda/envs/build/lib"
      ;;
  *)
      xspec_version_string="12.9.1"
      sed -i.orig "s/#cfitsio_libraries/cfitsio_libraries/g" setup.cfg
      sed -i.orig "s/#ccfits_libraries/ccfits_libraries/g" setup.cfg
      sed -i.orig "s/#wcslib_libraries/wcslib_libraries/g" setup.cfg
      ;;
esac

# Change build configuration
sed -i.orig "s/#with-xspec=True/with-xspec=True/g" setup.cfg
sed -i.orig "s|#xspec_lib_dirs = None|xspec_lib_dirs=${xspec_library_path}|g" setup.cfg
sed -i.orig "s|#xspec_include_dirs = None|xspec_include_dirs=${xspec_include_path}|g" setup.cfg
sed -i.orig "s|#gfortran_libraries = gfortran|gfortran_libraries= ${libgfortran_name}|g" setup.cfg
sed -i.orig "s|#xspec_version = 12.9.0|xspec_version = ${xspec_version_string}|g" setup.cfg


# Set HEADAS environment variables
export HEADAS=$miniconda/envs/build/Xspec/spectral
