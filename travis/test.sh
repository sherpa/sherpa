#!/usr/bin/env bash -e

# For now do not run tests with the documentation build
if [ "${DOCS}" == true ]; then exit 0; fi

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

# Install coverage tooling and run tests using setuptools
if [ ${TEST} == submodule ]; then
    pip install pytest-cov; python setup.py -q test -a "--cov sherpa --cov-report term" || exit 1;
fi

# Run smoke test
cd /home
sherpa_smoke ${smokevars} || exit 1

# Run regression tests using sherpa_test
if [ ${TEST} == package ] || [ ${TEST} == none ];
    then cd $HOME
    sherpa_test || exit 1
fi
