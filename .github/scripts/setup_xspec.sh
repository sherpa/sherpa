#!/usr/bin/env bash -e

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
xspec_library_path=${xspec_root}/lib/
xspec_include_path=${xspec_root}/include/

case "${XSPECVER}" in
  12.1*.* )
      xspec_version_string="${XSPECVER:0:7}"
      ;;
  *)
      echo "Xspec version listed currently unsupported in GitHub Actions jobs."
      exit 1
      ;;
esac

echo "* configuring XSPEC"

# Change build configuration
sed -i.orig "s|#with_xspec=True|with_xspec=True|" setup.cfg
sed -i.orig "s|#xspec_lib_dirs = None|xspec_lib_dirs=${xspec_library_path}|" setup.cfg
sed -i.orig "s|#xspec_include_dirs = None|xspec_include_dirs=${xspec_include_path}|" setup.cfg
sed -i.orig "s|#xspec_version = .*|xspec_version = ${xspec_version_string}|" setup.cfg

# Xspec ~12.10.1n now requires fftw. This disables the fftw build from extern
sed -i.orig "s|#fftw=local|fftw=local|" setup.cfg
sed -i.orig "s|#fftw_include_dirs=build/include|fftw_include_dirs=${CONDA_PREFIX}/include|" setup.cfg
sed -i.orig "s|#fftw_lib_dirs=build/lib|fftw_lib_dirs=${CONDA_PREFIX}/lib|" setup.cfg
sed -i.orig "s|#fftw_libraries=fftw3|fftw_libraries=fftw3|" setup.cfg

echo "* START setup.cfg"
cat setup.cfg
echo "* END   setup.cfg"
