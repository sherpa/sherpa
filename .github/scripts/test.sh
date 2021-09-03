#!/usr/bin/env bash -e

echo "*** Test=${TEST} ***"

if [ "`uname -s`" == "Darwin" ] ; then
    export DISPLAY=":99"
    export PATH="${PATH}:/opt/X11/bin"
    # Run headless Xvfb
    sudo Xvfb :99 -ac -screen 0 1024x768x8 &
    if [ $? != 0 ] ; then
        exit 1
    fi
fi

# No test data, then remove submodule (git automatically clones recursively)
if [ ${TEST} == none ];
 then git submodule deinit -f .;
fi

# Install test data as a package, then remove the submodule
if [ ${TEST} == package ];
 then pip install ./sherpa-test-data;
 pip install pytest-xvfb;
 git submodule deinit -f .;
fi

# Build smoke test switches, to ensure requested dependencies are reachable
if [ -n "${XSPECVER}" ]; then XSPECTEST="-x -d"; fi
if [ -n "${FITS}" ] ; then FITSTEST="-f ${FITS}"; fi
smokevars="${XSPECTEST} ${FITSTEST} -v 3"

if [ ${TEST} == submodule ]; then
    echo "*** Runnng Sherpa tests ***"
    # We want to be able to say
    #   pytest --cov sherpa --cov-report term || exit 1;
    # but this fails for a source installation (pip install .)
    # due to problems with the compiled components (such as
    # sherpa.utils._utils). So stick with setuptools for now.
    #
    python setup.py -q test -a "--cov sherpa --cov-report term" || exit 1;
    codecov;
fi

# Run smoke test
echo "*** Running sherpa_smoke ***"
cd /home
sherpa_smoke ${smokevars} || exit 1

# Run regression tests using sherpa_test
if [ ${TEST} == package ] || [ ${TEST} == none ]; then
    echo "*** Running sherpa_test ***"
    cd $HOME;
    # This automatically picks up the sherpatest modile when TEST==package
    sherpa_test --cov sherpa --cov-report term || exit 1;
    codecov;
fi
