#!/usr/bin/env python

#
#  Copyright (C) 2021, 2025
#  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""Usage:

  ./add_xspec_model.py version infile model

Aim:

Given an XSPEC model.dat file report and an XSPEC model name
(e.g. without the xs prefix), output the necessary code to create the
Python and C++ wrapper interface. You still need to add the code to
the necessary files and update them to add documentation and deal with
any possible version updates.

"""

import sys

from sherpa.astro.utils import xspec


Version = tuple[int, int, int]


def convert_version(smajor: str, sminor: str, smicro: str) -> Version:
    """Convert version tuple"""

    try:
        major = int(smajor)
        minor = int(sminor)
        micro = int(smicro)
    except ValueError as ve:
        raise ValueError(f"version={smajor}/{sminor}/{smicro} elements must be integers") from ve

    return major, minor, micro


def parse_version(ver: str) -> Version:
    """Convert version to its components."""

    match ver.split("."):
        case [smajor, sminor, smicro]:
            return convert_version(smajor, sminor, smicro)

        case _:
            raise ValueError(f"version={ver} should match 'X.Y.Z'")


def edit_python(version: Version,
                model: xspec.ModelDefinition,
                code: str
                ) -> str:
    """Update the Python code.

    The Python code is meant for user models, and not for the Sherpa
    XSPEC module, so we need to make a few simple changes to avoid the
    need for the caller to do this.

    """

    return code


def create_xspec_model(version: Version,
                       models: list[xspec.ModelDefinition],
                       modelname: str
                       ) -> None:
    """Write out the Python and C++ code needed to wrap the XSPEC model.

    Parameters
    ----------
    version
       The XSPEC version
    models
    modelname
        XSPEC model name (as given in the model.dat file)

    """

    for amodel in models:
        if amodel.name == modelname:
            model = amodel
            break
    else:
        raise ValueError(f"Unknown model name: {modelname}")

    code = xspec.create_xspec_code([model], internal=version)
    print("# C++ code for sherpa/astro/xspec/src/_xspec.cc\n")
    print(code.compiled)
    print("\n# Python code for sherpa/astro/xspec/__init__.py\n")
    print(edit_python(version, model, code.python))


if __name__ == "__main__":

    if len(sys.argv) != 4:
        sys.stderr.write(f"Usage: {sys.argv[0]} version infile model\n")
        sys.exit(1)

    version = parse_version(sys.argv[1])

    # This errors out in case of significant model-support issues
    # (e.g. a periodic model parameter) but can also just be an
    # issue with a poorly-specified interchange format (the model.dat
    # file) which may just need changes to this routine.
    #
    models = xspec.parse_xspec_model_description(sys.argv[2])
    create_xspec_model(version, models, sys.argv[3])
