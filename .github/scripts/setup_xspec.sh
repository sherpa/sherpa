#!/usr/bin/env bash -e

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
xspec_library_path=${xspec_root}/lib/
xspec_include_path=${xspec_root}/include/

case "${XSPECVER}" in
  12.13.0*)
      xspec_version_string="12.13.0"
      ;;
  12.12.1*)
      xspec_version_string="12.12.1"
      ;;
  12.12.0*)
      xspec_version_string="12.12.0"
      ;;
  *)
      echo "Xspec version listed currently unsupported in GitHub Actions jobs."
      exit 1
      ;;
esac

echo "* configuring XSPEC"

# Change build configuration
sed -i.orig "s|#with-xspec=True|with-xspec=True|" setup.cfg
sed -i.orig "s|#xspec_lib_dirs = None|xspec_lib_dirs=${xspec_library_path}|" setup.cfg
sed -i.orig "s|#xspec_include_dirs = None|xspec_include_dirs=${xspec_include_path}|" setup.cfg
sed -i.orig "s|#xspec_version = .*|xspec_version = ${xspec_version_string}|" setup.cfg

# Xspec ~12.10.1n now requires fftw. This disables the fftw build from extern
sed -i.orig "s|#fftw=local|fftw=local|" setup.cfg
sed -i.orig "s|#fftw-include_dirs=build/include|fftw-include_dirs=${CONDA_PREFIX}/include|" setup.cfg
sed -i.orig "s|#fftw-lib-dirs=build/lib|fftw-lib-dirs=${CONDA_PREFIX}/lib|" setup.cfg
sed -i.orig "s|#fftw-libraries=fftw3|fftw-libraries=fftw3|" setup.cfg

echo "* START setup.cfg"
cat setup.cfg
echo "* END   setup.cfg"
