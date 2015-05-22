#!/bin/bash -x

# CONDA
conda create --yes -n test -c astropy-ci-extras python=$PYTHON_VERSION pip
source activate test

# EGG_INFO
if [[ $SETUP_CMD == egg_info ]]
then
  exit  # no more dependencies needed
fi

# PEP8
if [[ $MAIN_CMD == pep8* ]]
then
  pip install pep8
  exit  # no more dependencies needed
fi

# CORE DEPENDENCIES
conda install --yes pytest Cython jinja2 psutil

# Install Sherpa
# TODO: really we only want to install all the dependencies. Is there a conda command for this?
conda config --add channels https://conda.binstar.org/cxc
conda install sherpa

# NUMPY
if [[ $NUMPY_VERSION == dev ]]
then
  pip install git+http://github.com/numpy/numpy.git
  export CONDA_INSTALL="conda install --yes python=$PYTHON_VERSION"
else
  conda install --yes numpy=$NUMPY_VERSION
  export CONDA_INSTALL="conda install --yes python=$PYTHON_VERSION numpy=$NUMPY_VERSION"
fi

# Now set up shortcut to conda install command to make sure the Python and Numpy
# versions are always explicitly specified.

# OPTIONAL DEPENDENCIES
if $OPTIONAL_DEPS
then
  $CONDA_INSTALL scipy h5py matplotlib pyyaml scikit-image pandas
  pip install beautifulsoup4
fi

# DOCUMENTATION DEPENDENCIES
# build_sphinx needs sphinx and matplotlib (for plot_directive). Note that
# this matplotlib will *not* work with py 3.x, but our sphinx build is
# currently 2.7, so that's fine
if [[ $SETUP_CMD == build_sphinx* ]]
then
  $CONDA_INSTALL Sphinx=1.2.2 Pygments matplotlib
fi

# COVERAGE DEPENDENCIES
if [[ $SETUP_CMD == 'test --coverage' ]]
then
  pip install coverage coveralls
fi
