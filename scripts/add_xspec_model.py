#!/usr/bin/env python

#
#  Copyright (C) 2021
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

  ./add_xspec_model.py infile model

Aim:

Given an XSPEC model.dat file report and an XSPEC model name
(e.g. without the xs prefix), output the necessary code to create the
Python and C++ wrapper interface. This does the minimum and you still
need to add the code to the necessary files and update them to add
documentation and deal with any possible version updates.

"""

import sys

from sherpa.astro.utils.xspec import parse_xspec_model_file, create_xspec_code


def create_xspec_model(models, modelname):
    """Write out the Python and C++ code needed to wrap the XSPEC model.

    Parameters
    ----------
    models : list of ModelDefinition
    modelname : str
        XSPEC model name (as given in the model.dat file)

    """

    model = None
    for amodel in models:
        if amodel.name == modelname:
            model = amodel
            break

    if model is None:
        raise ValueError(f"Unknown model name: {modelname}")

    code = create_xspec_code([model])
    print("# C++ code for sherpa/astro/xspec/src/_xspec.cc\n")
    print(code.compiled)
    print("\n# Python code for sherpa/astro/xspec/__init__.py\n")
    print(code.python)


if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} infile model\n")
        sys.exit(1)

    # This errors out in case of significant model-support issues
    # (e.g. a periodic model parameter) but can also just be an
    # issue with a poorly-specified interchange format (the model.dat
    # file) which may just need changes to this routine.
    #
    models = parse_xspec_model_file(sys.argv[1])
    create_xspec_model(models, sys.argv[2])
