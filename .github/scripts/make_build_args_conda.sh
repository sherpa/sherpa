#!/usr/bin/env bash -e
#
# What are the build options for Sherpa? At present the only option is whether
# to turn on XSPEC support, which is triggered by
#
# - HEADAS environment variable is set
#
# and the XSPECVER variable is used to set the XSPEC version. The XSPEC
# library and include files are assumed to be in $CONDA_PREFIX.
#

if [[ "x${CONDA_PREFIX}" == "x" ]];
then
    echo "Error: CONDA_PREFIX not set. This should be set for active Conda environments."
    exit 1
fi

# If HEADAS is not set then there is nothing to do
if [[ -z ${HEADAS} ]];
then
    exit 0
fi

if [[ -z ${XSPECVER} ]];
then
    echo "Error: XSPECVER is not set."
    exit 1
fi

# The xspec-modelsonly package has the "non-XSPEC" libraries already
# included, so the build-xspec-libraries argument can be simplified.
#
echo -Csetup-args=-Dxspec-prefix=$CONDA_PREFIX -Csetup-args=-Dxspec-libraries=XSFunctions,XSUtil,XS

# End
