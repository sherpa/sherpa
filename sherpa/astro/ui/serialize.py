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

import os

import sherpa.utils

from sherpa.utils.err import ArgumentTypeErr, IOErr

# Note: a lot of the serialization logic should probably be moved into
#       the objects (or modules) being serialized.
#


def _send_to_outfile(msg, filename=None):
    """Display the message.

    Parameters
    ----------
    msg : None or str
       The message to output. If ``None`` then the routine
       returns immediately (with no output).
    filename : None or str
       If ``None``, the message is printed to standard output,
       otherwise the file is opened (in append mode) and the
       message printed to it.
    """

    if msg is None:
        return

    try:
        if filename is None:
            print msg
        else:
            outfile = file(filename, 'a')
            print >> outfile, msg
    except:
        raise


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
    if type(par.units) == str:
        unitstr = "\"%s\"" % par.units

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

    if type(id) == str:
        return '"{}"'.format(id)
    else:
        return str(id)


def _save_statistic(state, outfile):
    """Save the statistic settings.

    Parameters
    ----------
    state
    outfile : None or str
       If ``None``, the message is printed to standard output,
       otherwise the file is opened (in append mode) and the
       statistic settings printed to it.
    """

    _send_to_outfile("\n######### Set Statistic\n", outfile)
    cmd = "set_stat(\"%s\")" % state.get_stat_name()
    _send_to_outfile(cmd, outfile)
    _send_to_outfile("", outfile)


def _save_fit_method(state, outfile):
    """Save the fit method settings.

    Parameters
    ----------
    state
    outfile : None or str
       If ``None``, the message is printed to standard output,
       otherwise the file is opened (in append mode) and the
       fitting-method settings printed to it.
    """

    # Save fitting method

    _send_to_outfile("\n######### Set Fitting Method\n", outfile)
    cmd = "set_method(\"%s\")" % state.get_method_name()
    _send_to_outfile(cmd, outfile)
    _send_to_outfile("", outfile)

    mdict = state.get_method_opt()
    for key in mdict:
        val = mdict.get(key)
        cmd = "set_method_opt(\"%s\", %s)" % (key, val)
        _send_to_outfile(cmd, outfile)

    _send_to_outfile("", outfile)


def _save_iter_method(state, outfile=None):
    """Save the iterated-fit method settings, if any.

    Parameters
    ----------
    state
    outfile : None or str
       If ``None``, the message is printed to standard output,
       otherwise the file is opened (in append mode) and the
       settings printed to it.
    """

    if state.get_iter_method_name() == 'none':
        return

    _send_to_outfile(
        "\n######### Set Iterative Fitting Method\n", outfile)
    cmd = "set_iter_method(\"%s\")" % state.get_iter_method_name()
    _send_to_outfile(cmd, outfile)
    _send_to_outfile("", outfile)

    mdict = state.get_iter_method_opt()
    for key in mdict:
        val = mdict.get(key)
        cmd = "set_iter_method_opt(\"%s\", %s)" % (key, val)
        _send_to_outfile(cmd, outfile)

    _send_to_outfile("", outfile)


def _save_model_components(state, outfile=None):
    """Save the model components.

    Parameters
    ----------
    state
    outfile : None or str
       If ``None``, the message is printed to standard output,
       otherwise the file is opened (in append mode) and the
       models printed to it.
    """

    # Have to call elements in list in reverse order (item at end of
    # list was model first created), so that we get links right, if any

    # To recreate attributes, print out dictionary as ordered pairs,
    # for each parameter

    _send_to_outfile(
        "\n######### Set Model Components and Parameters\n", outfile)
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

        if typename != "psfmodel" and typename != "tabelmodel" and \
           typename != "usermodel":
            # Normal case:  create an instance of the model.
            cmd = "eval(\"%s.%s\")" % (typename, modelname)
            _send_to_outfile(cmd, outfile)
        if typename == "psfmodel":
            cmd = "load_psf(\"%s\", \"%s\")" % (mod._name, mod.kernel.name)
            _send_to_outfile(cmd, outfile)
            try:
                psfmod = state.get_psf(id)
                cmd_id = _id_to_str(id)
                cmd = "set_psf(%s, %s)" % (cmd_id, psfmod._name)
                _send_to_outfile(cmd, outfile)
            except:
                pass
        if typename == "tablemodel":
            # Create table model with load_table_model
            cmd = "load_table_model(\"%s\", \"%s\")" % (
                modelname, mod.filename)
            _send_to_outfile(cmd, outfile)

        if typename == "convolutionkernel":
            # Create general convolution kernel with load_conv
            cmd = "load_conv(\"%s\", \"%s\")" % (
                modelname, mod.kernel.name)
            _send_to_outfile(cmd, outfile)

        if typename == "usermodel":
            # Skip user models -- don't create, don't set parameters
            # Go directly to next model in the model component list.
            _send_to_outfile(
                "WARNING: User model not saved, add any user model to save file manually\n", outfile)
            continue

        if hasattr(mod, "integrate"):
            cmd = "%s.integrate = %s" % (modelname, mod.integrate)
            _send_to_outfile(cmd, outfile)
            _send_to_outfile("", outfile)

        from sherpa.models import Parameter
        for par in mod.__dict__.values():
            if type(par) == Parameter or issubclass(Parameter, type(par)):
                par_attributes, par_linkstr = _print_par(par)
                _send_to_outfile(par_attributes, outfile)
                linkstr = linkstr + par_linkstr

        # If the model is a PSFModel, could have special
        # attributes "size" and "center" -- if so, record them.
        if typename == "psfmodel":
            if hasattr(mod, "size"):
                cmd = "%s.size = %s" % (modelname, repr(mod.size))
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("", outfile)
            if hasattr(mod, "center"):
                cmd = "%s.center = %s" % (modelname, repr(mod.center))
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("", outfile)

    # If there were any links made between parameters, send those
    # link commands to outfile now; else, linkstr is just an empty string
    _send_to_outfile(linkstr, outfile)


def _save_xspec(outfile=None):
    """Save the XSPEC settings, if the module is loaded.

    Parameters
    ----------
    outfile : None or str
       If ``None``, the message is printed to standard output,
       otherwise the file is opened (in append mode) and the
       XSPEC settings printed to it.
    """

    if not hasattr(sherpa.astro, "xspec"):
        return

    # TODO: should this make sure that the XSPEC module is loaded?
    _send_to_outfile("\n######### XSPEC Module Settings\n", outfile)
    xspec_state = sherpa.astro.xspec.get_xsstate()

    cmd = "set_xschatter(%d)" % xspec_state["chatter"]
    _send_to_outfile(cmd, outfile)
    cmd = "set_xsabund(\"%s\")" % xspec_state["abund"]
    _send_to_outfile(cmd, outfile)
    cmd = "set_xscosmo(%g, %g, %g)" % (xspec_state["cosmo"][0],
                                       xspec_state["cosmo"][1],
                                       xspec_state["cosmo"][2])
    _send_to_outfile(cmd, outfile)
    cmd = "set_xsxsect(\"%s\")" % xspec_state["xsect"]
    _send_to_outfile(cmd, outfile)
    for name in xspec_state["modelstrings"].keys():
        mstring = xspec_state["modelstrings"][name]
        cmd = "set_xsxset(\"%s\", \"%s\")" % (name, mstring)
        _send_to_outfile(cmd, outfile)


def save_all(state, outfile=None, clobber=False):
    """Save the information about the current session to a text file.

    This differs to the `save` command in that the output is human
    readable. Three consequences are:

     1. numeric values may not be recorded to their full precision

     2. data sets are not included in the file

     3. some settings and values may not be recorded.

    Parameters
    ----------
    outfile : str, optional
       If not given the results are displayed to the screen,
       otherwise it is taken to be the name of the file to
       write the results to.
    clobber : bool, optional
       If ``outfile`` is not ``None``, then this flag controls
       whether an existing file can be overwritten (``True``)
       or if it raises an exception (``False``, the default
       setting).

    Raises
    ------
    sherpa.utils.err.IOErr
       If ``outfile`` already exists and ``clobber`` is ``False``.

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

    Write the current Sherpa session to the screen:

    >>> save_all()

    Save the session to the file 'fit.sherpa', overwriting
    it if it already exists:

    >>> save_all('fit.sherpa', clobber=True)

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

    clobber = sherpa.utils.bool_cast(clobber)
    if type(outfile) == str:
        if os.path.isfile(outfile):
            if clobber:
                os.remove(outfile)
            else:
                raise IOErr('filefound', outfile)
    elif outfile is not None:
        raise ArgumentTypeErr('badarg', 'string or None')

    # Import numpy
    _send_to_outfile("import numpy", outfile)

    # Save data files

    _send_to_outfile("\n######### Load Data Sets\n", outfile)
    dids = state.list_data_ids()

    cmd_id = ""
    cmd_resp_id = ""
    cmd_bkg_id = ""

    for id in dids:
        # But if id is a string, then quote as a string
        # But what about the rest of any possible load_data() options;
        # how do we replicate the optional keywords that were possibly
        # used?  Store them with data object?
        cmd_id = _id_to_str(id)
        cmd = "load_data(%s,\"%s\")" % (cmd_id, state.get_data(id).name)
        _send_to_outfile(cmd, outfile)

        # Set physical or WCS coordinates here if applicable
        # If can't be done, just pass to next
        try:
            _send_to_outfile(
                "\n######### Set Image Coordinates \n", outfile)
            cmd = "set_coord(%s, %s)" % (cmd_id, repr(state.get_coord(id)))
            _send_to_outfile(cmd, outfile)
        except:
            pass

        # PHA attributes; group data if applicable
        try:
            # Only store group flags and quality flags if they were changed
            # from flags in the file
            if not state.get_data(id)._original_groups:
                if state.get_data(id).grouping is not None:
                    _send_to_outfile(
                        "\n######### Data Group Flags\n", outfile)
                    cmd = "set_grouping(%s, " % cmd_id
                    cmd = cmd + "val=numpy.array(" + repr(state.get_grouping(
                        id).tolist()) + ", numpy." + str(state.get_grouping(id).dtype) + "))"
                    _send_to_outfile(cmd, outfile)
                if state.get_data(id).quality is not None:
                    _send_to_outfile(
                        "\n######### Data Quality Flags\n", outfile)
                    cmd = "set_quality(%s, " % cmd_id
                    cmd = cmd + "val=numpy.array(" + repr(state.get_quality(
                        id).tolist()) + ", numpy." + str(state.get_quality(id).dtype) + "))"
                    _send_to_outfile(cmd, outfile)
            # End check for original groups and quality flags
            if state.get_data(id).grouped:
                cmd = "if get_data(%s).grouping is not None and not get_data(%s).grouped:" % (
                    cmd_id, cmd_id)
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("\t######### Group Data", outfile)
                cmd = "\tgroup(%s)" % cmd_id
                _send_to_outfile(cmd, outfile)
        except:
            pass

        # Add responses and ARFs, if any
        try:
            _send_to_outfile(
                "\n######### Data Spectral Responses\n", outfile)
            rids = state.list_response_ids(id)
            cmd_resp_id = ""

            for rid in rids:
                if type(rid) == str:
                    cmd_resp_id = "\"%s\"" % rid
                else:
                    cmd_resp_id = "%s" % rid

                try:
                    arf = state.get_arf(id, rid)
                    cmd = "load_arf(%s,\"%s\",%s)" % (
                        cmd_id, arf.name, cmd_resp_id)
                    _send_to_outfile(cmd, outfile)
                except:
                    pass

                try:
                    rmf = state.get_rmf(id, rid)
                    cmd = "load_rmf(%s,\"%s\",%s)" % (
                        cmd_id, rmf.name, cmd_resp_id)
                    _send_to_outfile(cmd, outfile)
                except:
                    pass
        except:
            pass

        # Check if this data set has associated backgrounds
        try:
            _send_to_outfile(
                "\n######### Load Background Data Sets\n", outfile)
            bids = state.list_bkg_ids(id)
            cmd_bkg_id = ""
            for bid in bids:
                if type(bid) == str:
                    cmd_bkg_id = "\"%s\"" % bid
                else:
                    cmd_bkg_id = "%s" % bid
                cmd = "load_bkg(%s,\"%s\", bkg_id=%s)" % (
                    cmd_id, state.get_bkg(id, bid).name, cmd_bkg_id)
                _send_to_outfile(cmd, outfile)

                # Group data if applicable
                try:
                    # Only store group flags and quality flags if they were changed
                    # from flags in the file
                    if not state.get_bkg(id, bid)._original_groups:
                        if state.get_bkg(id, bid).grouping is not None:
                            _send_to_outfile(
                                "\n######### Background Group Flags\n", outfile)
                            cmd = "set_grouping(%s, " % cmd_id
                            cmd = cmd + "val=numpy.array(" + repr(state.get_grouping(id).tolist()) + ", numpy." + str(
                                state.get_grouping(id, bid).dtype) + "), bkg_id=" + cmd_bkg_id + ")"
                            _send_to_outfile(cmd, outfile)
                        if state.get_bkg(id, bid).quality is not None:
                            _send_to_outfile(
                                "\n######### Background Quality Flags\n", outfile)
                            cmd = "set_quality(%s, " % cmd_id
                            cmd = cmd + "val=numpy.array(" + repr(state.get_quality(id).tolist()) + ", numpy." + str(
                                state.get_quality(id, bid).dtype) + "), bkg_id=" + cmd_bkg_id + ")"
                            _send_to_outfile(cmd, outfile)
                    # End check for original groups and quality flags
                    if state.get_bkg(id, bid).grouped:
                        cmd = "if get_bkg(%s,%s).grouping is not None and not get_bkg(%s,%s).grouped:" % (
                            cmd_id, cmd_bkg_id, cmd_id, cmd_bkg_id)
                        _send_to_outfile(cmd, outfile)
                        _send_to_outfile(
                            "\t######### Group Background", outfile)
                        cmd = "\tgroup(%s,%s)" % (cmd_id, cmd_bkg_id)
                        _send_to_outfile(cmd, outfile)
                except:
                    pass

                # Load background response, ARFs if any
                _send_to_outfile(
                    "\n######### Background Spectral Responses\n", outfile)
                rids = state.list_response_ids(id, bid)
                cmd_resp_id = ""
                for rid in rids:
                    if type(rid) == str:
                        cmd_resp_id = "\"%s\"" % rid
                    else:
                        cmd_resp_id = "%s" % rid

                    try:
                        arf = state.get_arf(id, rid, bid)
                        cmd = "load_arf(%s,\"%s\",%s,%s)" % (
                            cmd_id, arf.name, cmd_resp_id, cmd_bkg_id)
                        _send_to_outfile(cmd, outfile)
                    except:
                        pass

                    try:
                        rmf = state.get_rmf(id, rid, bid)
                        cmd = "load_rmf(%s,\"%s\",%s,%s)" % (
                            cmd_id, rmf.name, cmd_resp_id, cmd_bkg_id)
                        _send_to_outfile(cmd, outfile)
                    except:
                        pass

        except:
            pass

        # Set energy units if applicable
        # If can't be done, just pass to next
        try:
            _send_to_outfile(
                "\n######### Set Energy or Wave Units\n", outfile)
            units = state.get_data(id).units
            rate = state.get_data(id).rate
            if rate:
                rate = "\"rate\""
            else:
                rate = "\"counts\""
            factor = state.get_data(id).plot_fac
            cmd = "set_analysis(%s, %s, %s, %s)" % (cmd_id,
                                                    repr(units),
                                                    rate,
                                                    repr(factor))
            _send_to_outfile(cmd, outfile)
        except:
            pass

        # Subtract background data if applicable
        try:
            if state.get_data(id).subtracted:
                cmd = "if not get_data(%s).subtracted:" % cmd_id
                _send_to_outfile(cmd, outfile)
                _send_to_outfile(
                    "\t######### Subtract Background Data", outfile)
                cmd = "\tsubtract(%s)" % cmd_id
                _send_to_outfile(cmd, outfile)
        except:
            pass

        # Set filter if applicable
        try:
            _send_to_outfile("\n######### Filter Data\n", outfile)
            if len(state.get_data(id).get_filter()) > 0:
                filter = state.get_data(id).get_filter()
                if len(state.get_data(id).get_dims()) == 1:
                    cmd = "notice_id(%s,\"%s\")" % (cmd_id, filter)
                    _send_to_outfile(cmd, outfile)
                if len(state.get_data(id).get_dims()) == 2:
                    cmd = "notice2d_id(%s,\"%s\")" % (cmd_id, filter)
                    _send_to_outfile(cmd, outfile)
        except:
            pass

    _send_to_outfile("", outfile)

    _save_statistic(state, outfile)
    _save_fit_method(state, outfile)
    _save_iter_method(state, outfile)

    _save_model_components(state, outfile)

    # Save all source, pileup and background models

    _send_to_outfile(
        "\n######### Set Source, Pileup and Background Models\n", outfile)
    for id in dids:
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
                pass

            try:
                the_full_model = state.get_model(id)
            except:
                the_full_model = None
                pass

            if the_source is None and the_full_model is None:
                cmd = ""
                pass
            elif the_source is None and the_full_model is not None:
                cmd = "set_full_model(%s, %s)" % (
                    cmd_id, the_full_model.name)
            elif the_source is not None and the_full_model is None:
                cmd = "set_source(%s, %s)" % (cmd_id, state.the_source.name)
            elif the_source is not None and the_full_model is not None:
                if repr(the_source) == repr(the_full_model):
                    cmd = "set_full_model(%s, %s)" % (
                        cmd_id, the_full_model.name)
                else:
                    cmd = "set_source(%s, %s)" % (
                        cmd_id, state.the_source.name)
            else:
                # You can't actually get here
                cmd = ""
                pass
            _send_to_outfile(cmd, outfile)
            _send_to_outfile("", outfile)
        except:
            pass

        # If any pileup models, try to set them.  If not, just pass.
        try:
            cmd = "set_pileup_model(%s, %s)" % (
                cmd_id, state.get_pileup_model(id).name)
            _send_to_outfile(cmd, outfile)
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
                if type(bid) == str:
                    cmd_bkg_id = "\"%s\"" % bid
                else:
                    cmd_bkg_id = "%s" % bid

                the_bkg_source = None
                the_bkg_full_model = None
                try:
                    the_bkg_source = state.get_bkg_source(bid)
                except:
                    the_bkg_source = None
                    pass

                try:
                    the_bkg_full_model = state.get_bkg_model(id)
                except:
                    the_bkg_full_model = None
                    pass

                if the_bkg_source is None and the_bkg_full_model is None:
                    cmd = ""
                    pass
                elif the_bkg_source is None and the_bkg_full_model is not None:
                    cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (
                        cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                elif the_bkg_source is not None and the_bkg_full_model is None:
                    cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (
                        cmd_id, the_bkg_source.name, cmd_bkg_id)
                elif the_bkg_source is not None and \
                      the_bkg_full_model is not None:
                    if repr(the_bkg_source) == repr(the_bkg_full_model):
                        cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (
                            cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                    else:
                        cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (
                            cmd_id, the_bkg_source.name, cmd_bkg_id)
                else:
                    # You can't actually get here
                    cmd = ""
                    pass
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("", outfile)

        except:
            pass

    _save_xspec(outfile)


# What is this routine used for, and how is it different to save_all?
def save_session(state, outfile=None, clobber=False):

    # Check output file can be written to

    clobber = sherpa.utils.bool_cast(clobber)
    if type(outfile) == str:
        if os.path.isfile(outfile):
            if clobber:
                os.remove(outfile)
            else:
                raise IOErr('filefound', outfile)
    elif outfile is not None:
        raise ArgumentTypeErr('badarg', 'string or None')

    # Import numpy
    _send_to_outfile("import numpy", outfile)

    # Save data files

    _send_to_outfile("\n######### Load Data Sets\n", outfile)
    dids = state.list_data_ids()

    def get_logged_call(call_name, id=None):
        if id is not None:
            if state._calls_tracker.has_key(id) and state._calls_tracker[id].has_key(call_name):
                return state._calls_tracker[id][call_name]
        else:
            if state._calls_tracker.has_key(call_name):
                return state._calls_tracker[call_name]
    cmd_id = ""
    cmd_resp_id = ""
    cmd_bkg_id = ""

    for id in dids:
        # But if id is a string, then quote as a string
        # But what about the rest of any possible load_data() options;
        # how do we replicate the optional keywords that were possibly
        # used?  Store them with data object?
        cmd_id = _id_to_str(id)

        cmd = get_logged_call('load_data', id)
        _send_to_outfile(cmd, outfile)

        # Set physical or WCS coordinates here if applicable
        # If can't be done, just pass to next
        try:
            _send_to_outfile(
                "\n######### Set Image Coordinates \n", outfile)
            cmd = get_logged_call('set_coord', id)
            _send_to_outfile(cmd, outfile)
        except:
            pass

        # PHA attributes; group data if applicable
        try:
            # Only store group flags and quality flags if they were changed
            # from flags in the file
            if not state.get_data(id)._original_groups:
                if state.get_data(id).grouping is not None:
                    _send_to_outfile(
                        "\n######### Data Group Flags\n", outfile)
                    cmd = get_logged_call('set_grouping')
                    cmd = "set_grouping(%s, " % cmd_id
                    cmd = cmd + "val=numpy.array(" + repr(state.get_grouping(
                        id).tolist()) + ", numpy." + str(state.get_grouping(id).dtype) + "))"
                    _send_to_outfile(cmd, outfile)
                if state.get_data(id).quality is not None:
                    _send_to_outfile(
                        "\n######### Data Quality Flags\n", outfile)
                    cmd = "set_quality(%s, " % cmd_id
                    cmd = cmd + "val=numpy.array(" + repr(state.get_quality(
                        id).tolist()) + ", numpy." + str(state.get_quality(id).dtype) + "))"
                    _send_to_outfile(cmd, outfile)
            # End check for original groups and quality flags
            if state.get_data(id).grouped:
                cmd = "if get_data(%s).grouping is not None and not get_data(%s).grouped:" % (
                    cmd_id, cmd_id)
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("\t######### Group Data", outfile)
                cmd = "\tgroup(%s)" % cmd_id
                _send_to_outfile(cmd, outfile)
        except:
            pass

        # Add responses and ARFs, if any
        try:
            _send_to_outfile(
                "\n######### Data Spectral Responses\n", outfile)
            rids = state.list_response_ids(id)
            cmd_resp_id = ""

            for rid in rids:
                if type(rid) == str:
                    cmd_resp_id = "\"%s\"" % rid
                else:
                    cmd_resp_id = "%s" % rid

                try:
                    arf = state.get_arf(id, rid)
                    cmd = "load_arf(%s,\"%s\",%s)" % (
                        cmd_id, arf.name, cmd_resp_id)
                    _send_to_outfile(cmd, outfile)
                except:
                    pass

                try:
                    rmf = state.get_rmf(id, rid)
                    cmd = "load_rmf(%s,\"%s\",%s)" % (
                        cmd_id, rmf.name, cmd_resp_id)
                    _send_to_outfile(cmd, outfile)
                except:
                    pass
        except:
            pass

        # Check if this data set has associated backgrounds
        try:
            _send_to_outfile(
                "\n######### Load Background Data Sets\n", outfile)
            bids = state.list_bkg_ids(id)
            cmd_bkg_id = ""
            for bid in bids:
                if type(bid) == str:
                    cmd_bkg_id = "\"%s\"" % bid
                else:
                    cmd_bkg_id = "%s" % bid
                cmd = "load_bkg(%s,\"%s\", bkg_id=%s)" % (
                    cmd_id, state.get_bkg(id, bid).name, cmd_bkg_id)
                _send_to_outfile(cmd, outfile)

                # Group data if applicable
                try:
                    # Only store group flags and quality flags if they were changed
                    # from flags in the file
                    if not state.get_bkg(id, bid)._original_groups:
                        if state.get_bkg(id, bid).grouping is not None:
                            _send_to_outfile(
                                "\n######### Background Group Flags\n", outfile)
                            cmd = "set_grouping(%s, " % cmd_id
                            cmd = cmd + "val=numpy.array(" + repr(state.get_grouping(id).tolist()) + ", numpy." + str(
                                state.get_grouping(id, bid).dtype) + "), bkg_id=" + cmd_bkg_id + ")"
                            _send_to_outfile(cmd, outfile)
                        if state.get_bkg(id, bid).quality is not None:
                            _send_to_outfile(
                                "\n######### Background Quality Flags\n", outfile)
                            cmd = "set_quality(%s, " % cmd_id
                            cmd = cmd + "val=numpy.array(" + repr(state.get_quality(id).tolist()) + ", numpy." + str(
                                state.get_quality(id, bid).dtype) + "), bkg_id=" + cmd_bkg_id + ")"
                            _send_to_outfile(cmd, outfile)
                    # End check for original groups and quality flags
                    if state.get_bkg(id, bid).grouped:
                        cmd = "if get_bkg(%s,%s).grouping is not None and not get_bkg(%s,%s).grouped:" % (
                            cmd_id, cmd_bkg_id, cmd_id, cmd_bkg_id)
                        _send_to_outfile(cmd, outfile)
                        _send_to_outfile(
                            "\t######### Group Background", outfile)
                        cmd = "\tgroup(%s,%s)" % (cmd_id, cmd_bkg_id)
                        _send_to_outfile(cmd, outfile)
                except:
                    pass

                # Load background response, ARFs if any
                _send_to_outfile(
                    "\n######### Background Spectral Responses\n", outfile)
                rids = state.list_response_ids(id, bid)
                cmd_resp_id = ""
                for rid in rids:
                    if type(rid) == str:
                        cmd_resp_id = "\"%s\"" % rid
                    else:
                        cmd_resp_id = "%s" % rid

                    try:
                        arf = state.get_arf(id, rid, bid)
                        cmd = "load_arf(%s,\"%s\",%s,%s)" % (
                            cmd_id, arf.name, cmd_resp_id, cmd_bkg_id)
                        _send_to_outfile(cmd, outfile)
                    except:
                        pass

                    try:
                        rmf = state.get_rmf(id, rid, bid)
                        cmd = "load_rmf(%s,\"%s\",%s,%s)" % (
                            cmd_id, rmf.name, cmd_resp_id, cmd_bkg_id)
                        _send_to_outfile(cmd, outfile)
                    except:
                        pass

        except:
            pass

        # Set energy units if applicable
        # If can't be done, just pass to next
        try:
            _send_to_outfile(
                "\n######### Set Energy or Wave Units\n", outfile)
            units = state.get_data(id).units
            rate = state.get_data(id).rate
            if rate:
                rate = "\"rate\""
            else:
                rate = "\"counts\""
            factor = state.get_data(id).plot_fac
            cmd = "set_analysis(%s, %s, %s, %s)" % (cmd_id,
                                                    repr(units),
                                                    rate,
                                                    repr(factor))
            _send_to_outfile(cmd, outfile)
        except:
            pass

        # Subtract background data if applicable
        try:
            if state.get_data(id).subtracted:
                cmd = "if not get_data(%s).subtracted:" % cmd_id
                _send_to_outfile(cmd, outfile)
                _send_to_outfile(
                    "\t######### Subtract Background Data", outfile)
                cmd = "\tsubtract(%s)" % cmd_id
                _send_to_outfile(cmd, outfile)
        except:
            pass

        # Set filter if applicable
        try:
            _send_to_outfile("\n######### Filter Data\n", outfile)
            if len(state.get_data(id).get_filter()) > 0:
                filter = state.get_data(id).get_filter()
                if len(state.get_data(id).get_dims()) == 1:
                    cmd = "notice_id(%s,\"%s\")" % (cmd_id, filter)
                    _send_to_outfile(cmd, outfile)
                if len(state.get_data(id).get_dims()) == 2:
                    cmd = "notice2d_id(%s,\"%s\")" % (cmd_id, filter)
                    _send_to_outfile(cmd, outfile)
        except:
            pass

    _send_to_outfile("", outfile)

    _save_statistic(state, outfile)
    _save_fit_method(state, outfile)
    _save_iter_method(state, outfile)

    _save_model_components(state, outfile)

    # Save all source, pileup and background models

    _send_to_outfile(
        "\n######### Set Source, Pileup and Background Models\n", outfile)
    for id in dids:
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
                pass

            try:
                the_full_model = state.get_model(id)
            except:
                the_full_model = None
                pass

            if the_source is None and the_full_model is None:
                cmd = ""
                pass
            elif the_source is None and the_full_model is not None:
                cmd = "set_full_model(%s, %s)" % (
                    cmd_id, the_full_model.name)
            elif the_source is not None and the_full_model is None:
                cmd = "set_source(%s, %s)" % (cmd_id, state.the_source.name)
            elif the_source is not None and the_full_model is not None:
                if repr(the_source) == repr(the_full_model):
                    cmd = "set_full_model(%s, %s)" % (
                        cmd_id, the_full_model.name)
                else:
                    cmd = "set_source(%s, %s)" % (
                        cmd_id, state.the_source.name)
            else:
                # You can't actually get here
                cmd = ""
                pass
            _send_to_outfile(cmd, outfile)
            _send_to_outfile("", outfile)
        except:
            pass

        # If any pileup models, try to set them.  If not, just pass.
        try:
            cmd = "set_pileup_model(%s, %s)" % (
                cmd_id, state.get_pileup_model(id).name)
            _send_to_outfile(cmd, outfile)
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
                if type(bid) == str:
                    cmd_bkg_id = "\"%s\"" % bid
                else:
                    cmd_bkg_id = "%s" % bid

                the_bkg_source = None
                the_bkg_full_model = None
                try:
                    the_bkg_source = state.get_bkg_source(bid)
                except:
                    the_bkg_source = None
                    pass

                try:
                    the_bkg_full_model = state.get_bkg_model(id)
                except:
                    the_bkg_full_model = None
                    pass

                if the_bkg_source is None and the_bkg_full_model is None:
                    cmd = ""
                    pass
                elif the_bkg_source is None and the_bkg_full_model is not None:
                    cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (
                        cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                elif the_bkg_source is not None and \
                     the_bkg_full_model is None:
                    cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (
                        cmd_id, the_bkg_source.name, cmd_bkg_id)
                elif the_bkg_source is not None and \
                     the_bkg_full_model is not None:
                    if repr(the_bkg_source) == repr(the_bkg_full_model):
                        cmd = "set_bkg_full_model(%s, %s, bkg_id=%s)" % (
                            cmd_id, the_bkg_full_model.name, cmd_bkg_id)
                    else:
                        cmd = "set_bkg_source(%s, %s, bkg_id=%s)" % (
                            cmd_id, the_bkg_source.name, cmd_bkg_id)
                else:
                    # You can't actually get here
                    cmd = ""
                    pass
                _send_to_outfile(cmd, outfile)
                _send_to_outfile("", outfile)

        except:
            pass

    _save_xspec(outfile)
