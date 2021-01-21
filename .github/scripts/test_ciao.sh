#!/usr/bin/env bash -e

# We need to set up a Sherpa configuration file which uses the crates I/O
# backend.
#
cp sherpa/sherpa.rc $HOME/.sherpa-standalone.rc

echo "Sherpa is using the configuration file:"
python -c 'import sherpa; print(sherpa.get_config());'

pytest --cov sherpa --cov-report term || exit 1
codecov

# Run smoke test
cd $HOME
sherpa_smoke -x -d -f pycrates || exit 1
