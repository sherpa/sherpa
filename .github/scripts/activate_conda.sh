#!/usr/bin/env bash -e

# Activate conda and select the environment

if [ -z "${CONDA}" ]; then
    # set up for macOS
    CONDA="${HOME}/conda"
fi

source ${CONDA}/etc/profile.d/conda.sh
conda activate build
