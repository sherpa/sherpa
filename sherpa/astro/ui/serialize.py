#
#  Copyright (C) 2015, 2016, 2019, 2021, 2023
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

"""Serialize the Sherpa session state.

This module is used by ``sherpa.astro.ui.utils`` and is not
intended for public use. The API and semantics of the
routines in this module are subject to change.
"""

import inspect
import logging
import sys

import numpy

import sherpa.utils
from sherpa.utils.err import ArgumentErr

from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.astro.data import DataIMG, DataPHA

logger = logging.getLogger(__name__)
warning = logger.warning

string_types = (str, )

# Note: a lot of the serialization logic should probably be moved into
#       the objects (or modules) being serialized.
#


def _output(msg, fh=None):
    """Display the message.

    Parameters
    ----------
    msg : str
       The message to output.
    fh : None or a file handle
       The file handle to write the message to. If fh is ``None``
       then the standard output is used.
    """

    if fh is None:
        fh = sys.stdout

    fh.write(msg + '\n')


def _output_nl(fh=None):
    """Add a new-line.

    Parameters
    ----------
    fh : None or a file handle
       The file handle to write the message to. If fh is ``None``
       then the standard output is used.

    """

    _output("", fh)


def _output_banner(msg, fh=None, indent=0):
    """Display the banner message.

    Parameters
    ----------
    msg : str
       The label to output.
    fh : None or a file handle
       The file handle to write the message to. If fh is ``None``
       then the standard output is used.
    indent : int, optional
       How many times to indent the comment (if set the leading and
       trailing newlines aren't added).

    """

    if indent == 0:
        _output_nl(fh)

    space = ' ' * (indent * 4)
    _output(f"{space}######### {msg}", fh)
    if indent == 0:
        _output_nl(fh)


def _id_to_str(id):
    """Convert a data set identifier to a string value.

    Parameters
    ----------
    id : int or str
       The data set identifier.

    Returns
    -------
    out : str
       A string representation of the identifier for use
       in the Python serialization.
    """

    if isinstance(id, string_types):
        return f'"{id}"'

    return str(id)


def _save_entries(store, tostatement, fh=None):
    """Iterate through entries in the store and serialize them.

    Write out the key, value pairs in store in lexical order,
    rather than rely on the ordering of the container (e.g.
    whatever hash is used). The idea is to come up with a
    repeatable ordering, primarily to make testing easier.

    Parameters
    ----------
    store
       A container with keys. The elements of the container are
       passed to tocommand to create the string that is then
       written to fh.
    tostatement : func
       A function which accepts two arguments, the key and value
       from store, and returns a string. The reason for the name
       is that it is expected that the returned string will be
       a Python statement to restore this setting.
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    keys = list(store)
    keys.sort()
    for key in keys:
        cmd = tostatement(key, store[key])
        _output(cmd, fh)


def _save_intro(fh=None):
    """The set-up for the serialized file (imports).

    Parameters
    ----------
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    # QUS: should numpy only be loaded if it is needed?
    _output("import numpy", fh)
    _output("from sherpa.astro.ui import *", fh)


def _save_response(label, respfile, id, rid, bid=None, fh=None):
    """Save the ARF or RMF

    Parameters
    ----------
    label : str
       Either ``arf`` or ``rmf``.
    respfile : str
       The name of the ARF or RMF.
    id : id or str
       The Sherpa data set identifier.
    rid
       The Sherpa response identifier for the data set.
    bid
       If not ``None`` then this indicates that this is the ARF for
       a background dataset, and which such data set to use.
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    id = _id_to_str(id)
    rid = _id_to_str(rid)

    cmd = f'load_{label}({id}, "{respfile}", resp_id={rid}'
    if bid is not None:
        cmd += f", bkg_id={_id_to_str(bid)}"

    cmd += ")"
    _output(cmd, fh)


def _save_arf_response(state, id, rid, bid=None, fh=None):
    """Save the ARF.

    Parameters
    ----------
    state
    id : id or str
       The Sherpa data set identifier.
    rid
       The Sherpa response identifier for the data set.
    bid
       If not ``None`` then this indicates that this is the ARF for
       a background dataset, and which such data set to use.
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    try:
        respfile = state.get_arf(id, resp_id=rid, bkg_id=bid).name
    except:
        return

    _save_response('arf', respfile, id, rid, bid=bid, fh=fh)


def _save_rmf_response(state, id, rid, bid=None, fh=None):
    """Save the RMF.

    Parameters
    ----------
    state
    id : id or str
       The Sherpa data set identifier.
    rid
       The Sherpa response identifier for the data set.
    bid
       If not ``None`` then this indicates that this is the RMF for
       a background dataset, and which such data set to use.
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    try:
        respfile = state.get_rmf(id, resp_id=rid, bkg_id=bid).name
    except:
        return

    _save_response('rmf', respfile, id, rid, bid=bid, fh=fh)


def _save_pha_array(state, label, id, bid=None, fh=None):
    """Save a grouping or quality array for a PHA data set.

    Parameters
    ----------
    state
    label : "grouping" or "quality"
    id : id or str
       The Sherpa data set identifier.
    bid
       If not ``None`` then this indicates that the background dataset
       is to be used.
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    # This is an internal routine, so just protect against accidents
    if label not in ['grouping', 'quality']:
        raise ValueError(f"Invalid label={label}")

    if bid is None:
        data = state.get_data(id)
        lbl = 'Data'
    else:
        data = state.get_bkg(id, bid)
        lbl = 'Background'

    vals = getattr(data, label)
    if vals is None:
        return

    _output_banner(f"{lbl} {label} flags", fh)

    # QUS: can we not use the vals variable here rather than
    #      reassign it? i.e. isn't the "quality" (or "grouping"
    #      field of get_data/get_bkg the same as state.get_grouping
    #      or state.get_quality?
    #
    func = getattr(state, f'get_{label}')
    vals = func(id, bkg_id=bid)

    # The OGIP standard is for quality and grouping to be 2-byte
    # integers -
    # http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node7.html
    # - then force that type here.
    vals = vals.astype(numpy.int16)

    cmd = f"set_{label}({_id_to_str(id)}, " + \
          "val=numpy.array(" + repr(vals.tolist()) + \
          f", numpy.{vals.dtype})"
    if bid is not None:
        cmd += f", bkg_id={_id_to_str(bid)}"

    cmd += ")"
    _output(cmd, fh)


def _save_pha_grouping(state, id, bid=None, fh=None):
    """Save the grouping column values for a PHA data set.

    Parameters
    ----------
    state
    id : id or str
       The Sherpa data set identifier.
    bid
       If not ``None`` then this indicates that the background dataset
       is to be used.
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    _save_pha_array(state, "grouping", id, bid=bid, fh=fh)


def _save_pha_quality(state, id, bid=None, fh=None):
    """Save the quality column values for a PHA data set.

    Parameters
    ----------
    state
    id : id or str
       The Sherpa data set identifier.
    bid
       If not ``None`` then this indicates that the background dataset
       is to be used.
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    _save_pha_array(state, "quality", id, bid=bid, fh=fh)


def _handle_filter(state, id, fh):
    """Set any filter expressions for source and background
    components for data set id.

    It is expected that there is a dataset with the given id
    in the Sherpa session object (state).
    """

    cmd_id = _id_to_str(id)
    _output_banner("Filter Data", fh)
    d = state.get_data(id)
    fvals = d.get_filter()
    ndims = len(d.get_dims())
    if ndims == 1:
        cmd = f'notice_id({cmd_id}, "{fvals}")'
        _output(cmd, fh)
    elif ndims == 2:
        # Only output the filter if it does anything
        if fvals != '':
            cmd = f'notice2d_id({cmd_id}, "{fvals}")'
            _output(cmd, fh)
    else:
        # just in case
        _output(f'print("Set notice range of id={cmd_id} to {fvals}")',
                fh)

    try:
        bids = state.list_bkg_ids(id)
    except ArgumentErr:
        # Not a PHA data set
        return

    # Only set the noticed range if the data set does not have
    # the background subtracted. It might be useful to keep any
    # noticed range the user may have previously set - if switching
    # between fitting and subtracting the background - but that is
    # probably beyond the use case of the serialization.
    #
    if d.subtracted:
        return

    # NOTE: have to clear out the source filter before applying the
    #       background one.
    for bid in bids:
        bkg_id = _id_to_str(bid)
        fvals = state.get_bkg(id, bkg_id=bid).get_filter()

        # We know this is a PHA dataset so we do not have to worry
        # about things like 2D data.
        #
        _output(f'notice_id({cmd_id}, bkg_id={bkg_id})', fh)
        cmd = f'notice_id({cmd_id}, "{fvals}", bkg_id={bkg_id})'
        _output(cmd, fh)


def _save_data(state, fh=None):
    """Save the data.

    This can just be references to files, or serialization of
    the data (or both).

    Parameters
    ----------
    state
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.

    Notes
    -----
    This does not - at least at present - save data to one or
    more files to be read in by the script. If any data needs
    to be serialized it is included in the script.
    """

    _output_banner("Load Data Sets", fh)

    cmd_id = ""
    cmd_bkg_id = ""

    for id in state.list_data_ids():
        # But if id is a string, then quote as a string
        # But what about the rest of any possible load_data() options;
        # how do we replicate the optional keywords that were possibly
        # used?  Store them with data object?
        cmd_id = _id_to_str(id)

        cmd = _save_dataset(state, id)
        _output(cmd, fh)

        # Set physical or WCS coordinates here if applicable
        # If can't be done, just pass to next
        try:
            # TODO: add a test of the following
            _output_banner("Set Image Coordinates", fh)
            cmd = f"set_coord({_id_to_str(id)}, '{state.get_coord(id)}')"
            _output(cmd, fh)
        except:
            pass

        # PHA attributes; group data if applicable
        try:
            # Only store group flags and quality flags if they were changed
            # from flags in the file
            if not state.get_data(id)._original_groups:
                _save_pha_grouping(state, id, fh=fh)
                _save_pha_quality(state, id, fh=fh)

            # End check for original groups and quality flags
            if state.get_data(id).grouped:
                cmd = f"if get_data({cmd_id}).grouping is not None " + \
                    f"and not get_data({cmd_id}).grouped:"
                _output(cmd, fh)
                _output_banner("Group Data", fh, indent=1)
                _output(f"    group({cmd_id})", fh)
        except:
            pass

        # Add responses and ARFs, if any
        try:
            _output_banner("Data Spectral Responses", fh)
            rids = state.list_response_ids(id)

            for rid in rids:
                _save_arf_response(state, id, rid, fh=fh)
                _save_rmf_response(state, id, rid, fh=fh)

        except:
            pass

        # Check if this data set has associated backgrounds
        try:
            _output_banner("Load Background Data Sets", fh)
            bids = state.list_bkg_ids(id)
            cmd_bkg_id = ""
            for bid in bids:
                cmd_bkg_id = _id_to_str(bid)

                bname = state.get_bkg(id, bid).name
                cmd = f'load_bkg({cmd_id}, "{bname}", bkg_id={cmd_bkg_id})'
                _output(cmd, fh)

                # Group data if applicable
                try:
                    # Only store group flags and quality flags if they were
                    # changed from flags in the file
                    if not state.get_bkg(id, bid)._original_groups:
                        if state.get_bkg(id, bid).grouping is not None:
                            _save_pha_grouping(state, id, bid, fh=fh)
                            _save_pha_quality(state, id, bid, fh=fh)

                    # End check for original groups and quality flags
                    if state.get_bkg(id, bid).grouped:
                        cmd = f"if get_bkg({cmd_id}, {cmd_bkg_id}).grouping is not None " + \
                            f"and not get_bkg({cmd_id}, {cmd_bkg_id}).grouped:"
                        _output(cmd, fh)
                        _output_banner("Group Background", fh, indent=1)
                        _output(f"    group({cmd_id}, {cmd_bkg_id})", fh)
                except:
                    pass

                # Load background response, ARFs if any
                _output_banner("Background Spectral Responses", fh)
                rids = state.list_response_ids(id, bid)
                for rid in rids:
                    _save_arf_response(state, id, rid, bid, fh=fh)
                    _save_rmf_response(state, id, rid, bid, fh=fh)

        except:
            pass

        # Set energy units if applicable
        # If can't be done, just pass to next
        _output_banner("Set Energy or Wave Units", fh)
        try:
            units = state.get_data(id).units
            if state.get_data(id).rate:
                rate = "rate"
            else:
                rate = "counts"
            factor = state.get_data(id).plot_fac
            cmd = f'set_analysis({cmd_id}, quantity="{units}", ' + \
                f'type="{rate}", factor={factor})'
            _output(cmd, fh)
        except:
            pass

        # Subtract background data if applicable
        try:
            if state.get_data(id).subtracted:
                cmd = f"if not get_data({cmd_id}).subtracted:"
                _output(cmd, fh)
                _output_banner("Subtract Background Data", fh, indent=1)
                _output(f"    subtract({cmd_id})", fh)
        except:
            pass

        _handle_filter(state, id, fh)


def _print_par(par):
    """Convert a Sherpa parameter to a string.

    Note that we have to be careful with XSParameter parameters,
    to see if the hard limits need updating.

    Parameters
    ----------
    par
       The Sherpa parameter object to serialize.

    Returns
    -------
    out_pars, out_link : (str, str)
       A multi-line string serializing the contents of the
       parameter and then any link setting.
    """

    linkstr = ""
    if par.link is not None:
        linkstr = f"\nlink({par.fullname}, {par.link.fullname})\n"

    unitstr = ""
    if isinstance(par.units, string_types):
        unitstr = f'"{par.units}"'

    # Do we have to worry about XSPEC parameters which have changed their
    # hard min/max ranges?
    #
    parstrs = []
    try:
        if par.hard_min_changed():
            parstrs.append(f'{par.fullname}.hard_min    = {repr(par.hard_min)}')
    except AttributeError:
        pass

    try:
        if par.hard_max_changed():
            parstrs.append(f'{par.fullname}.hard_max    = {repr(par.hard_max)}')
    except AttributeError:
        pass

    parstrs.extend([f'{par.fullname}.default_val = {repr(par.default_val)}',
                    f'{par.fullname}.default_min = {repr(par.default_min)}',
                    f'{par.fullname}.default_max = {repr(par.default_max)}',
                    f'{par.fullname}.val     = {repr(par.val)}',
                    f'{par.fullname}.min     = {repr(par.min)}',
                    f'{par.fullname}.max     = {repr(par.max)}',
                    f'{par.fullname}.units   = {unitstr}',
                    f'{par.fullname}.frozen  = {par.frozen}'])
    parstr = '\n'.join(parstrs) + '\n'
    return (parstr, linkstr)


def _save_statistic(state, fh=None):
    """Save the statistic settings.

    Parameters
    ----------
    state
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    _output_banner("Set Statistic", fh)
    _output(f'set_stat("{state.get_stat_name()}")', fh)
    _output_nl(fh)


def _save_fit_method(state, fh=None):
    """Save the fit method settings.

    Parameters
    ----------
    state
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    # Save fitting method

    _output_banner("Set Fitting Method", fh)
    _output(f'set_method("{state.get_method_name()}")', fh)
    _output_nl(fh)

    def tostatement(key, val):
        # TODO: Using .format() returns more decimal places, which
        # is probably what we want but is a change, so leave
        # for now.
        # return 'set_method_opt("{}", {})'.format(key, val)
        return 'set_method_opt("%s", %s)' % (key, val)

    _save_entries(state.get_method_opt(), tostatement, fh)
    _output_nl(fh)


def _save_iter_method(state, fh=None):
    """Save the iterated-fit method settings, if any.

    Parameters
    ----------
    state
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    if state.get_iter_method_name() == 'none':
        return

    _output_banner("Set Iterative Fitting Method", fh)
    meth = state.get_iter_method_name()
    cmd = f'set_iter_method("{meth}")'
    _output(cmd, fh)
    _output_nl(fh)

    def tostatement(key, val):
        # There was a discussion here about the use of
        # str(val) vs val - the number of decimal places - but it
        # turns out for the test we only output integer values
        # so it makes no difference.
        return f'set_iter_method_opt("{key}", {val})'

    _save_entries(state.get_iter_method_opt(), tostatement, fh)
    _output_nl(fh)


# Is there something in the standard libraries that does this?
def _reindent(code):
    """Try to remove leading spaces. Somewhat hacky."""

    # Assume the first line is 'def func()'
    nspaces = code.find('def')
    if nspaces < 1:
        return code

    # minimal safety checks (e.g. if there was an indented
    # comment line).
    out = []
    for line in code.split("\n"):
        if line[:nspaces].isspace():
            out.append(line[nspaces:])
        else:
            out.append(line)

    return "\n".join(out)


# for user models, try to access the function definition via
# the inspect module and then re-create it in the script.
# An alternative would be to use the marshal module, and
# store the bytecode for the function in the code (or use
# pickle), but this is less readable and not guaranteed to
# be compatible with different major versions of Python.
# The idea is not to support all use case, but to try and
# support the simple use case.
#
# The user warnings are displayed for each user model,
# which could get annoying if there are many such models,
# but probably better to tell the user about each one
# than only the first.
#
def _handle_usermodel(mod, modelname, fh=None):

    try:
        pycode = inspect.getsource(mod.calc)
    except IOError:
        pycode = None

    # in case getsource can return None, have check here
    if pycode is None:
        msg = "Unable to save Python code for user model " + \
              f"'{modelname}' function {mod.calc.name}"
        warning(msg)
        _output(f'print("{msg}")', fh)
        _output(f"def {mod.calc.name}(*args):", fh)
        _output("    raise NotImplementedError('User model was " +
                "not saved by save_all().'", fh)
        _output_nl(fh)
        return

    msg = f"Found user model '{modelname}'; " + \
          "please check it is saved correctly."
    warning(msg)

    # Ensure the message is also seen if the script is run.
    _output(f'print("{msg}")', fh)

    _output(_reindent(pycode), fh)
    cmd = f'load_user_model({mod.calc.__name__}, "{modelname}")'
    _output(cmd, fh)

    # Work out the add_user_pars call; this is explicit, i.e.
    # it does not include logic to work out what arguments
    # are not needed.
    #
    # Some of these values are over-written later on, but
    # needed to set up the number of parameters, and good
    # documentation (hopefully).
    #
    parnames = [p.name for p in mod.pars]
    parvals = [p.default_val for p in mod.pars]
    # parmins = [p.default_min for p in mod.pars]
    # parmaxs = [p.default_max for p in mod.pars]
    parmins = [p.min for p in mod.pars]
    parmaxs = [p.max for p in mod.pars]
    parunits = [p.units for p in mod.pars]
    parfrozen = [p.frozen for p in mod.pars]

    spaces = '              '
    _output(f'add_user_pars("{modelname}",', fh)
    _output(f"{spaces}parnames={parnames},", fh)
    _output(f"{spaces}parvals={parvals},", fh)
    _output(f"{spaces}parmins={parmins},", fh)
    _output(f"{spaces}parmaxs={parmaxs},", fh)
    _output(f"{spaces}parunits={parunits},", fh)
    _output(f"{spaces}parfrozen={parfrozen}", fh)
    _output(f"{spaces})\n", fh)


def _save_model_components(state, fh=None):
    """Save the model components.

    Parameters
    ----------
    state
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    # Have to call elements in list in reverse order (item at end of
    # list was model first created), so that we get links right, if any

    # To recreate attributes, print out dictionary as ordered pairs,
    # for each parameter

    _output_banner("Set Model Components and Parameters", fh)
    all_model_components = state.list_model_components()
    all_model_components.reverse()

    # If there are any links between parameters, store link commands here
    # Then, *after* processing all models in the for loop below, send
    # link commands to outfile -- *all* models need to be created before
    # *any* links between parameters can be established.
    linkstr = ""
    for mod in all_model_components:

        # get actual model instance from the name we are given
        # then get model type, and name of this instance.
        mod = eval(mod)
        typename = mod.type
        modelname = mod.name.split(".")[1]

        # Special cases:

        # account for PSF creation elsewhere (above, with load_data
        # commands);

        # for table models, use "load_table_model", to ensure we
        # add to lists of known models *and* separate list of
        # tabel models;

        if typename == "usermodel":
            _handle_usermodel(mod, modelname, fh)

        elif typename == "psfmodel":
            cmd = f'load_psf("{mod._name}", "{mod.kernel.name}")'
            _output(cmd, fh)

        elif typename == "tablemodel":
            # Create table model with load_table_model
            cmd = f'load_table_model("{modelname}", "{mod.filename}")'
            _output(cmd, fh)

        else:
            # Normal case:  create an instance of the model.
            cmd = f'create_model_component("{typename}", "{modelname}")'
            _output(cmd, fh)

        # QUS: should this be included in the above checks?
        #      @DougBurke doesn't think so, as the "normal
        #      case" above should probably be run , but there's
        #      no checks to verify this.
        #
        if typename == "convolutionkernel":
            # Create general convolution kernel with load_conv
            cmd = f'load_conv("{modelname}", "{mod.kernel.name}")'
            _output(cmd, fh)

        if hasattr(mod, "integrate"):
            cmd = f"{modelname}.integrate = {mod.integrate}"
            _output(cmd, fh)
            _output_nl(fh)

        # Write out the parameters in the order they are stored in
        # the model. The original version of the code iterated
        # through mod.__dict__.values() and picked out Parameter
        # values.
        #
        for par in mod.pars:
            par_attributes, par_linkstr = _print_par(par)
            _output(par_attributes, fh)
            linkstr = linkstr + par_linkstr

        # If the model is a PSFModel, could have special
        # attributes "size" and "center" -- if so, record them.
        if typename == "psfmodel":
            spacer = False
            if hasattr(mod, "size") and mod.size is not None:
                _output(f"{modelname}.size = {mod.size}", fh)
                spacer = True

            if hasattr(mod, "center") and mod.center is not None:
                _output(f"{modelname}.center = {mod.center}", fh)
                spacer = True

            if spacer:
                _output_nl(fh)

    # If there were any links made between parameters, send those
    # link commands to outfile now; else, linkstr is just an empty string
    _output(linkstr, fh)

    # Now associate any PSF models with their appropriate datasets.
    # This is done after creating all the models and datasets.
    #
    if len(state._psf) == 0:
        return

    _output_banner("Associate PSF models with the datasets", fh)
    for idval, psfmod in state._psf.items():
        cmd_id = _id_to_str(idval)
        _output(f"set_psf({cmd_id}, {psfmod._name})", fh)


def _save_models(state, fh=None):
    """Save the source, pileup, and background models.

    Parameters
    ----------
    state
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    # Save all source, pileup and background models

    _output_banner("Set Source, Pileup and Background Models", fh)
    for id in state.list_data_ids():
        cmd_id = _id_to_str(id)

        # If a data set has a source model associated with it,
        # set that here -- try to distinguish cases where
        # source model is different from whole model.
        # If not, just pass
        try:
            try:
                the_source = state.get_source(id)
            except:
                the_source = None

            try:
                the_full_model = state.get_model(id)
            except:
                the_full_model = None

            have_source = the_source is not None
            have_full_model = the_full_model is not None

            if have_source:
                if have_full_model:
                    # for now assume that set_full_model is only
                    # used by PHA data sets.
                    try:
                        is_pha = isinstance(state.get_data(id), DataPHA)
                    except:
                        is_pha = False

                    if is_pha and repr(the_source) == repr(the_full_model):
                        cmd = f"set_full_model({cmd_id}, {the_full_model.name})"
                    else:
                        cmd = f"set_source({cmd_id}, {the_source.name})"
                else:
                    cmd = f"set_source({cmd_id}, {the_source.name})"

            elif have_full_model:
                cmd = f"set_full_model({cmd_id}, {the_full_model.name})"

            else:
                cmd = ""

            _output(cmd, fh)
            _output_nl(fh)
        except:
            pass

        # If any pileup models, try to set them.  If not, just pass.
        try:
            pname = state.get_pileup_model(id).name
            cmd = f"set_pileup_model({cmd_id}, {pname})"
            _output(cmd, fh)
        except:
            pass

        # Set background models (if any) associated with backgrounds
        # tied to this data set -- if none, then pass.  Again, try
        # to distinguish cases where background "source" model is
        # different from whole background model.
        try:
            bids = state.list_bkg_ids(id)
            cmd_bkg_id = ""
            for bid in bids:
                cmd_bkg_id = _id_to_str(bid)

                try:
                    the_bkg_source = state.get_bkg_source(id, bkg_id=bid)
                except:
                    the_bkg_source = None

                try:
                    the_bkg_full_model = state.get_bkg_model(id, bkg_id=bid)
                except:
                    the_bkg_full_model = None

                have_source = the_bkg_source is not None
                have_full_model = the_bkg_full_model is not None

                if have_source:
                    # This does not check for the dataset being a DataPHA
                    # object, since (at present) it has to be, as it's the
                    # only one to support backgrounds
                    if have_full_model:
                        if repr(the_bkg_source) == repr(the_bkg_full_model):
                            cmd = f"set_bkg_full_model({cmd_id}, {the_bkg_full_model.name}, bkg_id={cmd_bkg_id})"
                        else:
                            cmd = f"set_bkg_source({cmd_id}, {the_bkg_source.name}, bkg_id={cmd_bkg_id})"
                    else:
                        cmd = f"set_bkg_source({cmd_id}, {the_bkg_source.name}, bkg_id={cmd_bkg_id})"

                elif have_full_model:
                    cmd = f"set_bkg_full_model({cmd_id}, {the_bkg_full_model.name}, bkg_id={cmd_bkg_id})"

                else:
                    cmd = ""

                _output(cmd, fh)
                _output_nl(fh)

        except:
            pass


def _save_xspec(fh=None):
    """Save the XSPEC settings, if the module is loaded.

    Parameters
    ----------
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    if not hasattr(sherpa.astro, "xspec"):
        return

    # TODO: should this make sure that the XSPEC module is in use?
    #       i.e. only write these out if an XSPEC model is being used?
    _output_banner("XSPEC Module Settings", fh)
    xspec_state = sherpa.astro.xspec.get_xsstate()

    chatter =  xspec_state["chatter"]
    abund = xspec_state["abund"]
    xsect = xspec_state["xsect"]
    cs = xspec_state["cosmo"]
    _output(f"set_xschatter({chatter})", fh)
    _output(f'set_xsabund("{abund}")', fh)
    _output(f"set_xscosmo({cs[0]:g}, {cs[1]:g}, {cs[2]:g})", fh)
    _output(f'set_xsxsect("{xsect}")', fh)

    def tostatement(key, val):
        return f'set_xsxset("{key}", "{val}")'

    _save_entries(xspec_state["modelstrings"], tostatement, fh)


def _save_dataset(state, id):
    """Given a dataset identifier, return the text needed to
    re-create it.

    The data set design does not make it easy to tell:

    - if the data was read in from a file, or by load_arrays
      (and the name field set by the user)

    - if the data has been modified - e.g. by a call to set_counts -
      after it was loaded.

    - if in the correct directory (so paths may be wrong)

    """

    idstr = _id_to_str(id)
    dset = state.get_data(id)

    # For now assume that a missing name indicates the data
    # was created by the user, otherwise it's a file name.
    # The checks and error messages do not cover all situations,
    # but hopefully the common ones.
    #
    # This logic should be moved into the DataXXX objects, since
    # this is more extensible, and also the data object can
    # retain knowledge of where the data came from.
    #
    if dset.name.strip() == '':
        # Do not attempt to recreate all data sets at this
        # time. They could be serialized, but leave for a
        # later time (it may also make sense to actually
        # save the data to an external file).
        #
        if isinstance(dset, DataPHA):
            msg = f"Unable to re-create PHA data set '{id}'"
            warning(msg)
            return f'print("{msg}")'

        if isinstance(dset, DataIMG):
            msg = "Unable to re-create image data set '{id}'"
            warning(msg)
            return f'print("{msg}")'

        # Fall back to load_arrays. As using isinstance,
        # need to order the checks, since Data1DInt is
        # a subclass of Data1D.
        #
        xs = dset.get_indep()
        ys = dset.get_dep()
        staterr = dset.get_staterror()
        syserr = dset.get_syserror()

        need_sys = syserr is not None
        if need_sys:
            syserr = f"{syserr.tolist()}"
        else:
            syserr = "None"

        need_stat = staterr is not None or need_sys
        if staterr is not None:
            staterr = f"{staterr.tolist()}"
        else:
            staterr = "None"

        ys = f"{ys.tolist()}"

        out = f'load_arrays({idstr},\n'
        if isinstance(dset, Data1DInt):
            out += f'            {xs[0].tolist()},\n'
            out += f'            {xs[1].tolist()},\n'
            out += f'            {ys},\n'
            if need_stat:
                out += f'            {staterr},\n'
            if need_sys:
                out += f'            {syserr},\n'
            out += '            Data1DInt)'

        elif isinstance(dset, Data1D):
            out += f'            {xs[0].tolist()},\n'
            out += f'            {ys},\n'
            if need_stat:
                out += f'            {staterr},\n'
            if need_sys:
                out += f'            {syserr},\n'
            out += '            Data1D)'

        elif isinstance(dset, Data2DInt):
            msg = f"Unable to re-create Data2DInt data set '{id}'"
            warning(msg)
            out = f'print("{msg}")'

        elif isinstance(dset, Data2D):
            out += f'            {xs[0].tolist()},\n'
            out += f'            {xs[1].tolist()},\n'
            out += f'            {ys},\n'
            out += f'            {dset.shape},\n'
            if need_stat:
                out += f'            {staterr},\n'
            if need_sys:
                out += f'            {syserr},\n'
            out += '            Data2D)'

        else:
            msg = f"Unable to re-create {dset.__class__} data set '{id}'"
            warning(msg)
            out = f'print("{msg}")'

        return out

    # TODO: this does not handle options like selecting the columns
    #       from a file, or the number of columns.
    #
    if isinstance(dset, DataPHA):
        dtype = 'pha'
    elif isinstance(dset, DataIMG):
        dtype = 'image'
    else:
        dtype = 'data'

    return f'load_{dtype}({idstr}, "{dset.name}")'


def save_all(state, fh=None):
    """Save the information about the current session to a file handle.

    This differs to the `save` command in that the output is human
    readable. Three consequences are:

     1. numeric values may not be recorded to their full precision

     2. data sets are not included in the file

     3. some settings and values may not be recorded.

    Parameters
    ----------
    fh : file_like, optional
       If not given the results are displayed to standard out,
       otherwise they are written to this file handle.

    See Also
    --------
    save : Save the current Sherpa session to a file.
    restore : Load in a Sherpa session from a file.

    Notes
    -----

    This command will create a series of commands that restores
    the current Sherpa set up. It does not save the set of commands
    used. Not all Sherpa settings are saved. Items not fully restored
    include:

    - data sets changed from the version on disk - e.g. by calls to
      `sherpa.astro.ui.set_counts`

    - any optional keywords to comands such as `load_data`
      or `load_pha`

    - user models may not be restored correctly

    - only a subset of Sherpa commands are saved.

    Examples
    --------

    Write the current Sherpa session to the standard output:

    >>> save_all()

    Save the session to a StringIO handle:

    >>> from io import StringIO
    >>> store = StringIO()
    >>> save_all(store)

    """

    _save_intro(fh)
    _save_data(state, fh)
    _output_nl(fh)
    _save_statistic(state, fh)
    _save_fit_method(state, fh)
    _save_iter_method(state, fh)
    _save_model_components(state, fh)
    _save_models(state, fh)
    _save_xspec(fh)
