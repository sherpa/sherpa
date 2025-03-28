#!/usr/bin/env python

#
#  Copyright (C) 2024, 2025
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

  ./update_xspec_functions.py version infile

Aim:

Given an XSPEC model.dat file produce the output needed to be
included in sherpa/astro/xspec/src/_xspec.cc to fill up the
XSpecMethods array. The output is to stdout.

This code should be placed in _xspec.cc between the

  // Start model definitions

  // End model definitions

marks. To try and support old versions the script takes the current
_xspec.cc file and uses that as a basis for this section, so that
existing version checks will remain, but this is only meant as a
starting point, so check the output.

The XSPEC version should be given in dotted form - e.g. 12.15.0 - and
without any patch level.

"""

from importlib import resources
from pathlib import Path
import sys
from typing import Literal

from sherpa.astro.utils import xspec


# Should really provide better typing
ModelName = str
ModelDefinition = str

# A version can be major, minor, micro or, to support '#else' lines,
# the string "older", but that's only for the "state" handling.
#
Version = tuple[int, int, int]
StateVersion = Version | Literal["older"]

State = tuple[ModelName, ModelDefinition, StateVersion | None] | None

# Store a set of #ifdef...#else lines by
# - a dictionary where the key is each "#" line and the
#   values are the lines guarded by that check
# - a checksum value which just gives a combination of
#   all the "#.." lines so we can easily look for
#   a compatible set of lines.
#
Checksum = str
Storage = tuple[dict[str, list[str]], Checksum]

NO_CHECKSUM: Checksum = "<none>"

START_MACRO = "// Start model definitions"
END_MACRO = "// End model definitions"


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


def grab_version(inline: str) -> Version:
    """Given an #ifdef/#else line, grab the version"""

    toks = inline.split()
    if len(toks) != 2:
        raise ValueError(f"Expected two tokens in '{inline}'")

    if not toks[1].startswith("XSPEC_"):
        raise ValueError(f"Expected XSPEC_X_Y_Z in '{inline}'")

    version = toks[1][6:]
    match version.split("_"):
        case [smajor, sminor, smicro]:
            return convert_version(smajor, sminor, smicro)

        case _:
            raise ValueError(f"version={toks[1]} should match 'XSPEC_X_Y_Z'")


def make_xspec_macro(version: Version) -> str:
    """Create the XSPEC macro for this version."""
    return "XSPEC_{0}_{1}_{2}".format(*version)


def find_file() -> Path:
    """Find the _xspec.cc file"""

    basepath = Path("sherpa") / "astro" / "xspec" / "src" / "_xspec.cc"

    # Is this run from the base directory?
    if basepath.exists():
        return basepath

    # Is it run from scripts?
    inpath = Path("..") / basepath
    if inpath.exists():
        return inpath

    # Can we find it from the installed location (which implies an
    # editable install)? This is not ideal, since it is somewhat
    # mis-using the importlib machinery. The files is labelled as
    # returning a Traversable object which is not quite a Path, so
    # explicitly convert it here (even though the actual return type
    # is a Path).
    #
    with resources.as_file(resources.files("sherpa.astro.xspec")) as fs:
        infile = str(fs)

    inpath = Path(infile) / "src" / "_xspec.cc"
    if inpath.exists():
        return inpath

    raise OSError("Unable to find sherpa/astro/xspec/src/_xspec.cc")


def get_current_version() -> list[str]:
    """Return the current model definitions.

    This looks for, in order
       sherpa/astro/xspec/src/_xspec.cc
       ../sherpa/astro/xspec/src/_xspec.cc
       <installed version>

    """

    inpath = find_file()
    store = []
    with open(inpath, "rt") as fh:
        state = 0
        for inline in fh.readlines():
            if state == 0:
                if inline.find(START_MACRO) > -1:
                    state = 1
                continue

            if inline.find(END_MACRO) > -1:
                break

            store.append(inline.strip())

    if len(store) == 0:
        raise OSError(f"No definitions found in {inpath}")

    return store


def process_definitions(invals: list[str]
                        ) -> tuple[dict[ModelName, ModelDefinition],
                                   dict[ModelName, dict[StateVersion, ModelDefinition]]
                                   ]:
    """Extract the 'meaning' of each line

    We have
       - blank lines
       - model definitins
       - start macro definition
       - change macro definition
       - end macro definition

    Drop the blank lines and identify the model name and any version
    constraints.

    """

    versioned = {}
    unversioned = {}
    version: StateVersion | None = None
    for inval in invals:
        if inval == "":
            continue

        if inval == "#endif":
            assert version is not None, "Unexpected #endif"
            version = None
            continue

        if inval.startswith("#ifdef"):
            assert version is None, "Unexpected #ifdef"
            version = grab_version(inval)
            continue

        if inval == "#else":
            assert version is not None
            version = "older"
            continue

        if inval == "#elif":
            assert version is not None
            nversion = grab_version(inval)
            assert version != "older" and nversion < version, f"Version mismatch: {version} {nversion}"
            version = nversion
            continue

        assert not inval.startswith("#"), f"Errr: {inval}"

        pos = inval.find("// ")
        if pos == -1:
            assert ValueError(f"Missing // <model name> in '{inval}'")

        # To allow extra info added to the model name - e.g. "// XSfoobar  look at Xsfoo"
        # grab just the first word
        modelname = inval[pos + 3:].split()[0]
        model = inval[:pos].rstrip()

        assert modelname not in unversioned, f"model={modelname} already seen"

        if version is None:
            assert modelname not in versioned, f"model={modelname} already seen"
            unversioned[modelname] = model
            continue

        if modelname in versioned:
            versioned[modelname][version] = model
            continue

        versioned[modelname] = {version: model}

    assert len(unversioned) > 0, "What did I do?"
    assert len(versioned) > 0, "What did I do?"
    assert version is None, "somehow version constraint checks failed"
    return unversioned, versioned


def get_definition(mname: ModelName, mdef: ModelDefinition) -> str:
    return f"  {mdef:48s} // {mname}"


def fake_storage(line: str) -> Storage:
    """There should be a better way of indicating this."""

    # Using "older" here is not ideal
    return {"older": [line]}, NO_CHECKSUM


def print_storage(stored: Storage) -> None:
    """Display the stored lines."""

    for keyline, lines in stored[0].items():
        print(keyline)
        for line in lines:
            print(line)

    print("#endif")


def get_multiple(versioned: dict[StateVersion, ModelDefinition],
                 name: ModelName
                 ) -> Storage:
    """Create the #ifdef/... chain."""

    out = {}
    checksum = ""
    for version, mdef in versioned.items():
        if len(out) == 0:
            assert version != "older"
            key = f"#ifdef {make_xspec_macro(version)}"
        elif version != "older":
            key = f"#elif {make_xspec_macro(version)}"
        else:
            key = "#else"

        # We need this to store a list, even though here the list is
        # only ever a singleton.
        #
        out[key] = [get_definition(name, mdef)]
        checksum += key

    return out, checksum


def get_unversioned(version: Version,
                    olddef: ModelDefinition,
                    name: ModelName,
                    newdef: ModelDefinition
                    ) -> Storage:
    """Do we need to add a version to the definition?"""

    if olddef == newdef:
        # Ugh: special casing the "no-version" constraint is not ideal
        line = get_definition(name, olddef)
        return fake_storage(line)

    return get_multiple({version: newdef, "older": olddef}, name)


def get_versioned(version: Version,
                  old: dict[StateVersion, ModelDefinition],
                  name: ModelName,
                  newdef: ModelDefinition
                  ) -> Storage:
    """Add the new definion to the existing ones."""

    assert len(old) > 0

    # If the first element is the same as new then we do not need to
    # change anything.
    #
    if list(old.values())[0] == newdef:
        return get_multiple(old, name)

    # Assume that version is not in old
    assert version not in old
    new = {version: newdef} | old
    return get_multiple(new, name)


def get_new(version: Version,
            name: ModelName,
            newdef: ModelDefinition
            ) -> Storage:
    """Add the new definition"""

    return get_multiple({version: newdef}, name)


def report_model_definitions(version: Version,
                             old_unversioned: dict[ModelName, ModelDefinition],
                             old_versioned: dict[ModelName, dict[StateVersion, ModelDefinition]],
                             models: list[xspec.ModelDefinition]
                             ) -> None:
    """Print to stdout the code needed to bind to the models.

    There is a limited attempt to combine multiple neighbouring
    conditions which have the same set of version constraints, but it
    is not intended to be perfect.

    """

    print(f"  {START_MACRO}")

    stored: Storage | None = None

    prev_mtype = None
    for model in models:

        # For the moment models that are calculated per-spectrum are unsupported.
        try:
            if model.flags[1] == 1:
                continue
        except IndexError:
            pass

        # Assume an unsupported model; should this happen?
        try:
            new_def = xspec.model_to_compiled(model)[0]

            # Remove the leading spaces, to make the handling of old
            # and new definitions easier.
            #
            if new_def.startswith("  "):
                new_def = new_def[2:]

        except ValueError:
            continue

        if prev_mtype is None or prev_mtype != model.modeltype:
            # Display any stored lines as this is the "end" of this
            # sequence.
            #
            if stored is not None:
                print_storage(stored)
                stored = None

            print("")

        prev_mtype = model.modeltype

        mname = f"XS{model.name}"
        if mname in old_unversioned:
            nstore = get_unversioned(version, old_unversioned[mname], mname, new_def)

        elif mname in old_versioned:
            nstore = get_versioned(version, old_versioned[mname], mname, new_def)

        else:
            nstore = get_new(version, mname, new_def)

        if nstore[1] == NO_CHECKSUM:
            # This is unversioned, so we can display any stored values
            # and then this line.
            if stored is not None:
                print_storage(stored)
                stored = None

            lines = nstore[0]["older"]
            assert len(lines) == 1
            print(lines[0])
            continue

        # Store them if there's no current storage.
        #
        if stored is None:
            stored = nstore
            continue

        # Are the checksums the same? If they are then we can add the
        # new lines into the current storage area. If they are not
        # then we display the storage and swap in the new lines.
        #
        if stored[1] == nstore[1]:
            # Add the various lines into each section
            for key, lines in nstore[0].items():
                assert len(lines) == 1, f"len={len(lines)} lines='{lines}'"
                stored[0][key].extend(lines)

            continue

        print_storage(stored)
        stored = nstore

    # Print out any remaining stored lines.
    #
    if stored is not None:
        print_storage(stored)
        stored = None  # not really needed

    print("")
    print(f"  {END_MACRO}")



if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} version infile\n")
        sys.exit(1)

    version = parse_version(sys.argv[1])
    current = process_definitions(get_current_version())

    # This errors out in case of significant model-support issues
    # (e.g. a periodic model parameter) but can also just be an
    # issue with a poorly-specified interchange format (the model.dat
    # file) which may just need changes to this routine.
    #
    mdefs = xspec.parse_xspec_model_description(sys.argv[2])
    report_model_definitions(version, current[0], current[1], mdefs)
