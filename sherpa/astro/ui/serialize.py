#
#  Copyright (C) 2015  Smithsonian Astrophysical Observatory
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

import sys
import logging

import numpy

import sherpa.utils

logger = logging.getLogger(__name__)
warning = logger.warning

# Note: a lot of the serialization logic should probably be moved into
#       the objects (or modules) being serialized.
#


def _output(msg, fh=None):
    """Display the message.

    Parameters
    ----------
    msg : None or str
       The message to output. If ``None`` then the routine
       returns immediately (with no output).
    fh : None or a file handle
       The file handle to write the message to. If fh is ``None``
       then the standard output is used.
    """

    if msg is None:
        return

    if fh is None:
        fh = sys.stdout

    fh.write(msg + '\n')


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

    if isinstance(id, basestring):
        return '"{}"'.format(id)
    else:
        return str(id)


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

    cmd = 'load_{}({}, "{}", resp_id={}'.format(label, id, respfile, rid)
    if bid is not None:
        cmd += ", bkg_id={}".format(_id_to_str(bid))

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
        raise ValueError("Invalid label={}".format(label))

    if bid is None:
        data = state.get_data(id)
        lbl = 'Data'
    else:
        data = state.get_bkg(id, bid)
        lbl = 'Background'

    vals = getattr(data, label)
    if vals is None:
        return

    _output("\n######### {} {} flags\n".format(lbl, label), fh)

    # QUS: can we not use the vals variable here rather than
    #      reassign it? i.e. isn't the "quality" (or "grouping"
    #      field of get_data/get_bkg the same as state.get_grouping
    #      or state.get_quality?
    #
    func = getattr(state, 'get_{}'.format(label))
    vals = func(id, bkg_id=bid)

    # The OGIP standard is for quality and grouping to be 2-byte
    # integers -
    # http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node7.html
    # - then force that type here.
    vals = vals.astype(numpy.int16)

    cmd = "set_{}({}, ".format(label, _id_to_str(id)) + \
          "val=numpy.array(" + repr(vals.tolist()) + \
          ", numpy.{})".format(vals.dtype)
    if bid is not None:
        cmd += ", bkg_id={}".format(_id_to_str(bid))

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


def _save_data(state, funcs, fh=None):
    """Save the data.

    This can just be references to files, or serialization of
    the data (or both).

    Parameters
    ----------
    state
    funcs : dict
       A dictionary of function references with keys for ``load_data``
       and ``set_coord``, where the function accepts the data set
       identifier and returns a string representation of that command.
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.

    Notes
    -----
    This does not - at least at present - save data to one or
    more files to be read in by the script. If any data needs
    to be serialized it is included in the script.
    """

    _output("\n######### Load Data Sets\n", fh)

    cmd_id = ""
    cmd_bkg_id = ""

    for id in state.list_data_ids():
        # But if id is a string, then quote as a string
        # But what about the rest of any possible load_data() options;
        # how do we replicate the optional keywords that were possibly
        # used?  Store them with data object?
        cmd_id = _id_to_str(id)

        cmd = funcs['load_data'](id)
        _output(cmd, fh)

        # Set physical or WCS coordinates here if applicable
        # If can't be done, just pass to next
        try:
            _output("\n######### Set Image Coordinates\n", fh)
            cmd = funcs['set_coord'](id)
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
                cmd = "if get_data(%s).grouping is not None and not get_data(%s).grouped:" % (
                    cmd_id, cmd_id)
                _output(cmd, fh)
                _output("    ######### Group Data", fh)
                cmd = "    group(%s)" % cmd_id
                _output(cmd, fh)
        except:
            pass

        # Add responses and ARFs, if any
        try:
            _output(
                "\n######### Data Spectral Responses\n", fh)
            rids = state.list_response_ids(id)

            for rid in rids:
                _save_arf_response(state, id, rid, fh=fh)
                _save_rmf_response(state, id, rid, fh=fh)

        except:
            pass

        # Check if this data set has associated backgrounds
        try:
            _output(
                "\n######### Load Background Data Sets\n", fh)
            bids = state.list_bkg_ids(id)
            cmd_bkg_id = ""
            for bid in bids:
                cmd_bkg_id = _id_to_str(bid)

                cmd = 'load_bkg(%s, "%s", bkg_id=%s)' % (
                    cmd_id, state.get_bkg(id, bid).name, cmd_bkg_id)
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
                        cmd = "if get_bkg(%s, %s).grouping is not None and not get_bkg(%s, %s).grouped:" % (
                            cmd_id, cmd_bkg_id, cmd_id, cmd_bkg_id)
                        _output(cmd, fh)
                        _output(
                            "    ######### Group Background", fh)
                        cmd = "    group(%s, %s)" % (cmd_id, cmd_bkg_id)
                        _output(cmd, fh)
                except:
                    pass

                # Load background response, ARFs if any
                _output(
                    "\n######### Background Spectral Responses\n", fh)
                rids = state.list_response_ids(id, bid)
                for rid in rids:
                    _save_arf_response(state, id, rid, bid, fh=fh)
                    _save_rmf_response(state, id, rid, bid, fh=fh)

        except:
            pass

        # Set energy units if applicable
        # If can't be done, just pass to next
        try:
            _output(
                "\n######### Set Energy or Wave Units\n", fh)
            units = state.get_data(id).units
            rate = state.get_data(id).rate
            if rate:
                rate = '"rate"'
            else:
                rate = '"counts"'
            factor = state.get_data(id).plot_fac
            cmd = "set_analysis(%s, %s, %s, %s)" % (cmd_id,
                                                    repr(units),
                                                    rate,
                                                    repr(factor))
            _output(cmd, fh)
        except:
            pass

        # Subtract background data if applicable
        try:
            if state.get_data(id).subtracted:
                cmd = "if not get_data(%s).subtracted:" % cmd_id
                _output(cmd, fh)
                _output(
                    "    ######### Subtract Background Data", fh)
                cmd = "    subtract(%s)" % cmd_id
                _output(cmd, fh)
        except:
            pass

        # Set filter if applicable
        try:
            _output("\n######### Filter Data\n", fh)
            fvals = state.get_data(id).get_filter()
            ndims = len(state.get_data(id).get_dims())
            if ndims == 1:
                cmd = 'notice_id({}, "{}")'.format(cmd_id, fvals)
                _output(cmd, fh)
            elif ndims == 2:
                cmd = 'notice2d_id({}, "{}")'.format(cmd_id, fvals)
                _output(cmd, fh)
        except:
            pass


def _print_par(par):
    """Convert a Sherpa parameter to a string.

    Parameters
    ----------
    par
       The Sherpa parameter object to serialize.

    Returns
    -------
    out : str
       A multi-line string serializing the contents of the
       parameter.
    """

    linkstr = ""
    if par.link is not None:
        linkstr = "\nlink(%s, %s)\n" % (
            par.fullname, par.link.fullname)

    unitstr = ""
    if isinstance(par.units, basestring):
        unitstr = '"%s"' % par.units

    return ((('%s.default_val = %s\n' +
              '%s.default_min = %s\n' +
              '%s.default_max = %s\n' +
              '%s.val     = %s\n' +
              '%s.min     = %s\n' +
              '%s.max     = %s\n' +
              '%s.units   = %s\n' +
              '%s.frozen  = %s\n') %
             (par.fullname, repr(par.default_val),
              par.fullname, repr(par.default_min),
              par.fullname, repr(par.default_max),
              par.fullname, repr(par.val),
              par.fullname, repr(par.min),
              par.fullname, repr(par.max),
              par.fullname, unitstr,
              par.fullname, par.frozen)), linkstr)


def _save_statistic(state, fh=None):
    """Save the statistic settings.

    Parameters
    ----------
    state
    fh : None or file-like
       If ``None``, the information is printed to standard output,
       otherwise the information is added to the file handle.
    """

    _output("\n######### Set Statistic\n", fh)
    cmd = 'set_stat("%s")' % state.get_stat_name()
    _output(cmd, fh)
    _output("", fh)


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

    _output("\n######### Set Fitting Method\n", fh)
    cmd = 'set_method("%s")' % state.get_method_name()
    _output(cmd, fh)
    _output("", fh)

    mdict = state.get_method_opt()
    for key in mdict:
        val = mdict.get(key)
        cmd = 'set_method_opt("%s", %s)' % (key, val)
        _output(cmd, fh)

    _output("", fh)


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

    _output("\n######### Set Iterative Fitting Method\n", fh)
    cmd = 'set_iter_method("%s")' % state.get_iter_method_name()
    _output(cmd, fh)
    _output("", fh)

    mdict = state.get_iter_method_opt()
    for key in mdict:
        val = mdict.get(key)
        cmd = 'set_iter_method_opt("%s", %s)' % (key, val)
        _output(cmd, fh)

    _output("", fh)


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

    _output("\n######### Set Model Components and Parameters\n", fh)
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
        modelname = mod.name.partition(".")[2]

        # Special cases:

        # account for PSF creation elsewhere (above, with load_data
        # commands);

        # for table models, use "load_table_model", to ensure we
        # add to lists of known models *and* separate list of
        # tabel models;

        # skip user models entirely, as they require importation of
        # user modules, beyond scope of this script.

        if typename == "usermodel":
            # Skip user models -- don't create, don't set parameters
            #
            # At present this is reported for each component
            msg = "User model '{}'".format(mod.name) + \
                  " not saved, add any user model to save file manually"
            warning(msg)

            # This will cause a syntax error when the file is run,
            # but this is a "good thing" here.
            _output("WARNING: {}\n".format(msg), fh)

            # Go directly to next model in the model component list.
            continue

        elif typename == "psfmodel":
            cmd = 'load_psf("%s", "%s")' % (mod._name, mod.kernel.name)
            _output(cmd, fh)
            try:
                psfmod = state.get_psf(id)
                cmd_id = _id_to_str(id)
                cmd = "set_psf(%s, %s)" % (cmd_id, psfmod._name)
                _output(cmd, fh)
            except:
                pass

        elif typename == "tablemodel":
            # Create table model with load_table_model
            cmd = 'load_table_model("%s", "%s")' % (
                modelname, mod.filename)
            _output(cmd, fh)

        else:
            # Normal case:  create an instance of the model.
            cmd = 'create_model_component("{}", "{}")'.format(
                typename, modelname)
            _output(cmd, fh)

        # QUS: should this be included in the above checks?
        #      @DougBurke doesn't think so, as the "normal
        #      case" above should probably be run , but there's
        #      no checks to verify this.
        #
        if typename == "convolutionkernel":
            # Create general convolution kernel with load_conv
            cmd = 'load_conv("%s", "%s")' % (
                modelname, mod.kernel.name)
            _output(cmd, fh)

        if hasattr(mod, "integrate"):
            cmd = "%s.integrate = %s" % (modelname, mod.integrate)
            _output(cmd, fh)
            _output("", fh)

        from sherpa.models import Parameter
        for par in mod.__dict__.values():
            if isinstance(par, Parameter):
                par_attributes, par_linkstr = _print_par(par)
                _output(par_attributes, fh)
                linkstr = linkstr + par_linkstr

        # If the model is a PSFModel, could have special
        # attributes "size" and "center" -- if so, record them.
        if typename == "psfmodel":
            if hasattr(mod, "size"):
                cmd = "%s.size = %s" % (modelname, repr(mod.size))
                _output(cmd, fh)
                _output("", fh)
            if hasattr(mod, "center"):
                cmd = "%s.center = %s" % (modelname, repr(mod.center))
                _output(cmd, fh)
                _output("", fh)

    # If there were any links made between parameters, send those
    # link commands to outfile now; else, linkstr is just an empty string
    _output(linkstr, fh)


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

    _output("\n######### Set Source, Pileup and Background Models\n", fh)
    for id in state.list_data_ids():
        cmd_id = _id_to_str(id)

        # If a data set has a source model associated with it,
        # set that here -- try to distinguish cases where
        # source model is different from whole model.
        # If not, just pass
        try:
            the_source = None
            the_full_model = None
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
                    if repr(the_source) == repr(the_full_model):
                        cmd = "set_full_model(%s, %s)" % (
                            cmd_id, the_full_model.name)
                    else:
                        cmd = "set_source(%s, %s)" % (
                            cmd_id, the_source.name)
                else:
                    cmd = "set_source(%s, %s)" % (cmd_id, the_source.name)

            elif have_full_model:
                cmd = "set_full_model(%s, %s)" % (
                    cmd_id, the_full_model.name)

            else:
                cmd = ""

            _output(cmd, fh)
            _output("", fh)
        except:
            pass

        # If any pileup models, try to set them.  If not, just pass.
        try:
            cmd = "set_pileup_model(%s, %s)" % (
                cmd_id, state.get_pileup_model(id).name)
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

                the_bkg_source = None
                the_bkg_full_model = None
                try:
                    the_bkg_source = state.get_bkg_source(bid)
                except:
                    the_bkg_source = None

                try:
                    the_bkg_full_model = state.get_bkg_model(id)
                except:
                    the_bkg_full_model = None

                have_source = the_bkg_source is not None
                have_full_model = the_bkg_full_model is not None

                if have_source:
                    if have_full_model:
                        if repr(the_bkg_source) == repr(the_bkg_full_model):
                            cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (
                                cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                        else:
                            cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (
                                cmd_id, the_bkg_source.name, cmd_bkg_id)
                    else:
                        cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (
                            cmd_id, the_bkg_source.name, cmd_bkg_id)

                elif have_full_model:
                    cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (
                        cmd_id, the_bkg_full_model.name, cmd_bkg_id)

                else:
                    cmd = ""

                _output(cmd, fh)
                _output("", fh)

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

    # TODO: should this make sure that the XSPEC module is loaded?
    _output("\n######### XSPEC Module Settings\n", fh)
    xspec_state = sherpa.astro.xspec.get_xsstate()

    cmd = "set_xschatter(%d)" % xspec_state["chatter"]
    _output(cmd, fh)
    cmd = 'set_xsabund("%s")' % xspec_state["abund"]
    _output(cmd, fh)
    cmd = "set_xscosmo(%g, %g, %g)" % (xspec_state["cosmo"][0],
                                       xspec_state["cosmo"][1],
                                       xspec_state["cosmo"][2])
    _output(cmd, fh)
    cmd = 'set_xsxsect("%s")' % xspec_state["xsect"]
    _output(cmd, fh)
    for name in xspec_state["modelstrings"].keys():
        mstring = xspec_state["modelstrings"][name]
        cmd = 'set_xsxset("%s", "%s")' % (name, mstring)
        _output(cmd, fh)


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

    Items which are not saved include:

    - user models

    - any optional keywords to comands such as `load_data`
      or `load_pha`

    - only a subset of Sherpa commands are saved.

    Examples
    --------

    Write the current Sherpa session to the standard output:

    >>> save_all()

    Save the session to a StringIO handle:

    >>> import StringIO
    >>> store = StringIO.StringIO()
    >>> save_all(store)

    """

    # TODO:
    #
    #    1) Finish RMF, ARF settings for backgrounds    DONE
    #    2) Add PSF models, table models, kernels etc.  DONE
    #    2a) Account for multi-response model           DONE
    #    2b) And background model (set_bkg)             DONE
    #    2c) And pileup (set_pileup_model)              DONE
    #    3) Any way to deal with user models?           SKIP
    #    4) Energy, coord settings for every data set   DONE
    #    5) Filters for each data set                   DONE
    #    6) Group flags for each data set               DONE

    # 7) Save optional keyword arguments for load_data/load_bkg
    #    8) Set model integrate flags                   DONE
    #    9) Subtract flags                              DONE

    # Check output file can be written to

    funcs = {
        'load_data': lambda id:
        'load_data({}, "{}")'.format(_id_to_str(id), state.get_data(id).name),
        'set_coord': lambda id:
        'set_coord({}, {})'.format(_id_to_str(id), repr(state.get_coord(id))),
    }

    _save_intro(fh)
    _save_data(state, funcs, fh)
    _output("", fh)
    _save_statistic(state, fh)
    _save_fit_method(state, fh)
    _save_iter_method(state, fh)
    _save_model_components(state, fh)
    _save_models(state, fh)
    _save_xspec(fh)


# What is this routine used for, and how is it different to save_all?
# The only difference is that the load_data and set_coord values are
# retrieved from the "loggable" facility. Do we need to keep this?
#
def save_session(state, fh=None):

    def get_logged_call(call_name, id=None):
        try:
            if id is None:
                return state._calls_tracker[call_name]
            else:
                return state._calls_tracker[id][call_name]
        except KeyError:
            # QUS: should this error out or return None?
            return None

    funcs = {
        'load_data': lambda id: get_logged_call('load_data', id),
        'set_coord': lambda id: get_logged_call('set_coord', id),
    }

    _save_intro(fh)
    _save_data(state, funcs, fh)
    _output("", fh)
    _save_statistic(state, fh)
    _save_fit_method(state, fh)
    _save_iter_method(state, fh)
    _save_model_components(state, fh)
    _save_models(state, fh)
    _save_xspec(fh)
