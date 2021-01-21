#!/usr/bin/env bash -e
#

# edit the setup.cfg file; break up into a lot of separate commands!
CFG=setup.cfg
sed -i -e 's|#disable-|disable-|' ${CFG}
sed -i -e 's|#fftw=local|fftw=local|' ${CFG}
sed -i -e 's|#fftw-include_dirs=build/include|fftw-include_dirs='${ASCDS_LIB}'/../include|' ${CFG}
sed -i -e 's|#fftw-lib-dirs=build/lib|fftw-lib-dirs='${ASCDS_LIB}'|' ${CFG}
sed -i -e 's|#fftw-libraries=fftw3|fftw-libraries=fftw3|' ${CFG}

sed -i -e 's|#region=local|region=local|' ${CFG}
sed -i -e 's|#region-include_dirs=build/include|region-include_dirs='${ASCDS_LIB}'/../include|' ${CFG}
sed -i -e 's|#region-lib-dirs=build/lib|region-lib-dirs='${ASCDS_LIB}'|' ${CFG}
sed -i -e 's|#region-libraries=region|region-libraries=region ascdm|' ${CFG}

sed -i -e 's|#wcs=local|wcs=local|' ${CFG}
sed -i -e 's|#wcs-include_dirs=build/include|wcs-include_dirs='${ASCDS_LIB}'/../include|' ${CFG}
sed -i -e 's|#wcs-lib-dirs=build/lib|wcs-lib-dirs='${ASCDS_LIB}'|' ${CFG}
sed -i -e 's|#wcs-libraries=wcs|wcs-libraries=wcs|' ${CFG}

sed -i -e 's|#with-xspec=True|with-xspec=True|' ${CFG}
sed -i -e 's|#xspec_version = 12.9.0|xspec_version = 12.10.1|' ${CFG}
sed -i -e 's|#xspec_lib_dirs = None|xspec_lib_dirs = '${ASCDS_LIB}'|' ${CFG}
sed -i -e 's|#xspec_include_dirs = None|xspec_include_dirs = '${ASCDS_LIB}'/../include|' ${CFG}

# Show the differences
git diff ${CFG}
