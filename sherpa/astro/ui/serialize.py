#
#  Copyright (C) 2015-2016, 2019, 2021, 2023-2026
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

The aim is that each load_xxx call from the
sherpa.astro.ui.utils.Session object will fill a "FileStore"
dictionary which will indicate what the call was (so it can be
re-created). The FileStore dataclass is a simple way to store this
information.

There are three issues that complicate this:

1) For non-PHA datasets we can just index on the "idval" argument,
   but PHA files also add in responses and backrounds (which also
   use the IdType as an index). This means that there are some
   cases where we have multiple levels of access.

2) PHA responses can be automatically loaded (using the ANCRFILE,
   BACKFILE, and RESPFILE keywords)

   These are indicated by storing a None rather than a FileStore
   object.

3) PHA2 files can load in multiple datasets with a single call
   (and the id argument may not match the values associated
   with the dataset). So we need a way to associate those datasets
   read in from a single call, as well as knowing when they have
   been deleted.

"""

from collections.abc import Callable, Sequence
from dataclasses import KW_ONLY, dataclass
import inspect
import logging
import os
import sys
import textwrap
from types import ModuleType
from typing import TYPE_CHECKING, Any, Mapping, TextIO, TypedDict, \
    overload

import numpy

from sherpa.astro.data import DataIMG, DataPHA, DataARF, DataRMF
from sherpa.astro import io
from sherpa.astro.io.wcs import WCS

from sherpa.data import Data, Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.models.basic import UserModel
from sherpa.utils import get_keyword_defaults
from sherpa.utils.types import IdType

if TYPE_CHECKING:
    # Avoid an import cycle
    from sherpa.astro.ui.utils import Session

xspec: ModuleType | None
try:
    from sherpa.astro import xspec
except ImportError:
    xspec = None


logger = logging.getLogger(__name__)
warning = logger.warning

string_types = (str, )

# Note: a lot of the serialization logic should probably be moved into
#       the objects (or modules) being serialized.
#

OutType = TypedDict("OutType", {"imports": set[str], "main": list[str]})
MaybeIdType = IdType | None
DataType = Data1D | Data2D

# Typing the state argument is awkward since it causes an import
# error, so for runtime we default to Any and only when run under a
# type checker do we use the actual type.
#
if TYPE_CHECKING:
    SessionType = Session
else:
    SessionType = Any

# The par parameter is hard to type as it's Parameter |
# XSBaseParameter but the latter is only defined if XSPEC support is
# enabled. One way around this is to create a Protocol for defining
# the parameter interface rather than have it be based on a class.
#
# For now we use Any to skip this, at the expense of not type checking
# the module property
#
# ParameterType = Parameter
ParameterType = Any


def _id_to_str(id: IdType) -> str:
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


def showval(v: Any) -> str:
    """How to display a "value" for a function call."""

    # Special case certain values.
    # Is this check too generic?
    #
    try:
        return v.__name__
    except AttributeError:
        pass

    if isinstance(v, str):
        return f'"{v}"'

    # Should there be some attempt to replace 'np.int64(0)' by '0'?
    return str(v)


def remove_default_args(func: Callable,
                        kwargs: Mapping[str, Any]
                        ) -> dict[str, Any]:
    """Remove elements from kwargs that match the defaults for func."""

    out = {}
    defaults = get_keyword_defaults(func)
    for k, v in kwargs.items():
        try:
            if v == defaults[k]:
                continue
        except KeyError:
            # The function may have **kwargs in its signature, so we
            # want to just pass those on.
            pass

        out[k] = v

    return out


@dataclass
class FileStore:
    """Record how a dataset was read in."""

    loadfunc: Callable
    """The function to load the data."""

    idval: IdType
    """The dataset identifier."""

    filename: str
    """The location of the file."""

    _: KW_ONLY

    args: Sequence[Any]
    """Positional arguments used to read in the file."""

    kwargs: Mapping[str, Any]
    """Named arguments used to read in the file."""

    autoloaded: bool = False
    """Was the file auto-loaded?"""

    def show(self) -> str:
        """How to load the data, ignoring the autoloaded setting."""

        idstr = _id_to_str(self.idval)
        out = f'{self.loadfunc.__name__}({idstr}, '
        out += f'"{self.filename}"'
        for arg in self.args:
            out += f", {showval(arg)}"

        kwargs = remove_default_args(self.loadfunc, self.kwargs)
        for k, v in kwargs.items():
            out += f", {k}={showval(v)}"

        return f"{out})"


def _output(out: OutType,
            msg: str,
            indent: int = 0
            ) -> None:
    """Output the line (if it exists)."""

    space = ' ' * (indent * 4)
    out["main"].append(textwrap.indent(msg, space))


def _output_nl(out: OutType) -> None:
    """Add a new-line."""
    _output(out, "")


def _output_banner(out: OutType, msg: str) -> None:
    """Display the banner message.

    Parameters
    ----------
    out : dict
       The output state
    msg : str
       The label to output.

    """

    _output_nl(out)
    _output(out, f"######### {msg}")
    _output_nl(out)


def _get_out_pos(out: OutType) -> int:
    """Return the current position.

    This is to make it easier to remove a banner call. The code
    should be re-written so we only add text once we know we have
    anything, but that is a relatively large change.
    """
    return len(out["main"])


def _remove_banner(out: OutType, orig_pos: int) -> None:
    """Remove the banner message if no text has been added.

    Parameters
    ----------
    out : dict
       The output state
    orig_pos : int
       The position after the banner was added.

    """

    if _get_out_pos(out) != orig_pos:
        return

    out["main"].pop()
    out["main"].pop()
    out["main"].pop()


def _save_entries(out: OutType,
                  store: Mapping[str, Any],
                  tostatement: Callable[[str, Any], str]) -> None:
    """Iterate through entries in the store and serialize them.

    Write out the key, value pairs in store in lexical order,
    rather than rely on the ordering of the container (e.g.
    whatever hash is used). The idea is to come up with a
    repeatable ordering, primarily to make testing easier.

    Parameters
    ----------
    out : dict
       The output state
    store
       A container with keys. The elements of the container are
       passed to tostatement to create the string that is then
       saved in out.
    tostatement : func
       A function which accepts two arguments, the key and value
       from store, and returns a string. The reason for the name
       is that it is expected that the returned string will be
       a Python statement to restore this setting.
    """

    keys = list(store)
    keys.sort()
    for key in keys:
        cmd = tostatement(key, store[key])
        _output(out, cmd)


def _show_response(state: SessionType,
                   label: str,
                   respfile: str,
                   idval: IdType,
                   rid: IdType,
                   bid: MaybeIdType = None
                   ) -> str:
    """How to load the ARF or RMF

    Parameters
    ----------
    state
    label : str
       Either ``arf`` or ``rmf``.
    respfile : str
       The name of the ARF or RMF.
    idval : id or str
       The Sherpa data set identifier.
    rid
       The Sherpa response identifier for the data set.
    bid
       If not ``None`` then this indicates that this is the ARF for
       a background dataset, and which such data set to use.
    """

    kwargs = {"resp_id": rid}
    if bid is not None:
        kwargs["bkg_id"] = bid

    func = getattr(state, f'load_{label}')
    store = FileStore(func, idval, respfile, args=(),
                      kwargs=kwargs)
    return store.show()


def _show_arf_response(state: SessionType,
                       arf: DataARF,
                       id: IdType,
                       rid: IdType,
                       bid: MaybeIdType = None
                       ) -> str:
    """How to load the ARF.

    Parameters
    ----------
    state
    arf : DataARF
       The ARF object
    id : id or str
       The Sherpa data set identifier.
    rid
       The Sherpa response identifier for the data set.
    bid
       If not ``None`` then this indicates that this is the ARF for
       a background dataset, and which such data set to use.
    """

    return _show_response(state, 'arf', arf.name, id, rid, bid=bid)


def _show_rmf_response(state: SessionType,
                       rmf: DataRMF,
                       id: IdType,
                       rid: IdType,
                       bid: MaybeIdType = None
                       ) -> str:
    """How to load the RMF.

    Parameters
    ----------
    state
    rmf : DataRMF
       The RMF object
    id : id or str
       The Sherpa data set identifier.
    rid
       The Sherpa response identifier for the data set.
    bid
       If not ``None`` then this indicates that this is the RMF for
       a background dataset, and which such data set to use.
    """

    return _show_response(state, 'rmf', rmf.name, id, rid, bid=bid)


def _save_pha_array(out: OutType,
                    state: SessionType,
                    pha: DataPHA,
                    label: str,
                    id: IdType,
                    bid: MaybeIdType = None) -> None:
    """Save a grouping or quality array for a PHA data set.

    Parameters
    ----------
    out : dict
       The output state
    state
    pha : DataPHA
       The PHA object
    label : "grouping" or "quality"
    id : id or str
       The Sherpa data set identifier.
    bid
       If not ``None`` then this indicates that the background dataset
       is to be used.
    """

    vals = getattr(pha, label)
    if vals is None:
        return

    # The OGIP standard is for quality and grouping to be 2-byte
    # integers -
    # https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node7.html
    # - then force that type here.
    vals = vals.astype(numpy.int16)

    lbl = 'Data' if bid is None else 'Background'
    _output_banner(out, f"{lbl} {label} flags")

    cmd = f"set_{label}({_id_to_str(id)}, " + \
          "val=numpy.array(" + repr(vals.tolist()) + \
          f", numpy.{vals.dtype})"
    if bid is not None:
        cmd += f", bkg_id={_id_to_str(bid)}"

    cmd += ")"
    _output(out, cmd)


def _save_pha_grouping(out: OutType,
                       state: SessionType,
                       pha: DataPHA,
                       id: IdType,
                       bid: MaybeIdType = None) -> None:
    """Save the grouping column values for a PHA data set.

    Parameters
    ----------
    out : dict
       The output state
    state
    pha : DataPHA
       The PHA object
    id : id or str
       The Sherpa data set identifier.
    bid
       If not ``None`` then this indicates that the background dataset
       is to be used.
    """

    _save_pha_array(out, state, pha, "grouping", id, bid=bid)


def _save_pha_quality(out: OutType,
                      state: SessionType,
                      pha: DataPHA,
                      id: IdType,
                      bid: MaybeIdType = None) -> None:
    """Save the quality column values for a PHA data set.

    Parameters
    ----------
    out : dict
       The output state
    state
    pha : DataPHA
       The PHA object
    id : id or str
       The Sherpa data set identifier.
    bid
       If not ``None`` then this indicates that the background dataset
       is to be used.
    """

    _save_pha_array(out, state, pha, "quality", id, bid=bid)


def _add_filter(out: OutType,
                data: DataType,
                cmd_id: str) -> None:
    """Add any filter applied to the source."""

    # We have mask being
    #  - True or an array of Trues
    #    => no filter is needed
    #
    #  - False or an array of Falses
    #    => everything is filtered out
    #
    #  - array of True and False
    #    => want to filter the data
    #
    # We only report a filter if it actually excludes some data.  The
    # following relies on numpy.all/any behaving the same when given
    # boolval as [boolval], so we do not have to explicitly check if
    # the mask is a boolean.
    #
    if numpy.all(data.mask):
        return

    extra = "2d" if len(data.get_dims()) == 2 else ""
    if numpy.any(data.mask):
        fvals = data.get_filter()
        _output(out, f'notice{extra}_id({cmd_id}, "{fvals}")')
        return

    # All data has been ignored.
    #
    _output(out, f'ignore{extra}_id({cmd_id})')


def _add_bkg_filter(out: OutType,
                    bpha: DataPHA,
                    cmd_id: str,
                    bkg_id: int) -> None:
    """Add any filter applied to the background."""

    if numpy.all(bpha.mask):
        return

    if numpy.any(bpha.mask):
        # We need to clear any existing background filter set by
        # the source.
        #
        _output(out, f'notice_id({cmd_id}, bkg_id={bkg_id})')
        fvals = bpha.get_filter()
        _output(out, f'notice_id({cmd_id}, "{fvals}", bkg_id={bkg_id})')
        return

    _output(out, f'ignore_id({cmd_id}, bkg_id={bkg_id})')


def _handle_filter(out: OutType,
                   state: SessionType,  # currently unused
                   data: DataType,
                   id: IdType) -> None:
    """Set any filter expressions for source and background
    components for data set id.

    It is expected that there is a dataset with the given id
    in the Sherpa session object (state).
    """

    _output_banner(out, "Filter Data")
    orig_pos = _get_out_pos(out)

    cmd_id = _id_to_str(id)

    _add_filter(out, data, cmd_id)

    # Since the state doesn't have type annotations, help the system out.
    # Technically this is only true after the try call, not here.
    #
    if TYPE_CHECKING:
        assert isinstance(data, DataPHA)

    try:
        bids = data.background_ids
    except AttributeError:
        # Not a PHA data set
        _remove_banner(out, orig_pos)
        return

    # Only set the noticed range if the data set does not have
    # the background subtracted. It might be useful to keep any
    # noticed range the user may have previously set - if switching
    # between fitting and subtracting the background - but that is
    # probably beyond the use case of the serialization.
    #
    if data.subtracted:
        _remove_banner(out, orig_pos)
        return

    # NOTE: have to clear out the source filter before applying the
    #       background one.
    for bid in bids:
        bpha = data.get_background(bid)

        if TYPE_CHECKING:
            assert isinstance(bpha, DataPHA)

        _add_bkg_filter(out, bpha, cmd_id, bid)

    _remove_banner(out, orig_pos)


def _load_from_store(out: OutType,
                     idval: IdType,
                     store: dict[IdType, FileStore],
                     auto_load: bool
                     ) -> bool:
    """If the file is in the store, display the load call if wanted."""

    try:
        val = store[idval]
    except KeyError:
        return False

    if not auto_load or not val.autoloaded:
        _output(out, val.show())

    return True


def _load_bkg(out: OutType,
              idval: IdType,
              bid: IdType,
              bpha: DataPHA,
              store: dict[IdType, FileStore],
              auto_load: bool):
    """How to load the background?"""

    if _load_from_store(out, bid, store, auto_load):
        return

    # Manually create the background.
    raise RuntimeError("what to do")


def _load_response(state: SessionType,
                   out: OutType,
                   resp: DataARF | DataRMF | None,
                   idval: IdType,
                   rid: IdType,
                   store: dict[IdType, FileStore],
                   show: Callable[..., str],
                   auto_load: bool
                   ) -> None:
    """Do we need to load in a response?"""

    if resp is None:
        return

    if _load_from_store(out, rid, store, auto_load):
        return

    # TODO: This should really create the response from first principles
    # rather than reading it in from a file.
    #
    cmd = show(state, resp, idval, rid)
    _output(out, cmd)


def _load_bkg_response(state: SessionType,
                       out: OutType,
                       resp: DataARF | DataRMF | None,
                       idval: IdType,
                       bid: IdType,
                       rid: IdType,
                       store: dict[IdType, FileStore],
                       show: Callable[..., str],
                       auto_load: bool
                       ) -> None:
    """Do we need to load a background response."""

    if resp is None:
        return

    if _load_from_store(out, rid, store, auto_load):
        return

    # TODO: This should really create the response from first principles
    # rather than reading it in from a file.
    #
    cmd = show(state, resp, idval, rid, bid=bid)
    _output(out, cmd)


def _save_dataset_settings_pha(out: OutType,
                               state: SessionType,
                               pha: DataPHA,
                               idval: IdType,
                               auto_load: bool
                               ) -> None:
    """What settings need to be set for DataPHA"""

    cmd_id = _id_to_str(idval)

    # Only store group flags and quality flags if they were changed
    # from flags in the file.
    #
    if not pha._original_groups:
        _save_pha_grouping(out, state, pha, idval)
        _save_pha_quality(out, state, pha, idval)

    if pha.grouped:
        _output(out, f"group({cmd_id})")

    # Do we need to load any responses (ARF, RMF, background)?
    #
    arf_store = state._load_arf_store.get(idval, {})
    rmf_store = state._load_rmf_store.get(idval, {})

    rids = pha.response_ids
    if len(rids) > 0:
        _output_banner(out, "Data Spectral Responses")

        for rid in rids:
            arf, rmf = pha.get_response(rid)
            _load_response(state, out, arf, idval, rid, store=arf_store,
                           show=_show_arf_response, auto_load=auto_load)
            _load_response(state, out, rmf, idval, rid, store=rmf_store,
                           show=_show_rmf_response, auto_load=auto_load)

    bkg_store = state._load_bkg_store.get(idval, {})
    bkg_arf_store = state._load_bkg_arf_store.get(idval, {})
    bkg_rmf_store = state._load_bkg_rmf_store.get(idval, {})

    bids = pha.background_ids
    if len(bids) > 0:

        # Although there is support for PHA2 background files, do not
        # make use of this knowledge here as we do not have any such
        # files to test on.
        #
        _output_banner(out, "Load Background Data Sets")
        for bid in bids:
            cmd_bkg_id = _id_to_str(bid)

            bpha = pha.get_background(bid)
            if TYPE_CHECKING:
                assert isinstance(bpha, DataPHA)

            _load_bkg(out, idval, bid, bpha, store=bkg_store, auto_load=auto_load)

            # Only store group flags and quality flags if they were
            # changed from flags in the file
            #
            if not bpha._original_groups:
                _save_pha_grouping(out, state, bpha, idval, bid=bid)
                _save_pha_quality(out, state, bpha, idval, bid=bid)

            if bpha.grouped:
                _output(out, f"group({cmd_id}, bkg_id={cmd_bkg_id})")

            # Load background response, ARFs if any.
            #
            rids = bpha.response_ids
            if len(rids) > 0:
                # Can we only print this banner when we need to load
                # the response?
                #
                _output_banner(out, "Background Spectral Responses")

                bkg_arf_store_bid = bkg_arf_store.get(bid, {})
                bkg_rmf_store_bid = bkg_rmf_store.get(bid, {})
                for rid in rids:
                    bkg_arf, bkg_rmf = bpha.get_response(rid)
                    _load_bkg_response(state, out, bkg_arf, idval, bid, rid,
                                       store=bkg_arf_store_bid, show=_show_arf_response,
                                       auto_load=auto_load)
                    _load_bkg_response(state, out, bkg_rmf, idval, bid, rid,
                                       store=bkg_rmf_store_bid, show=_show_rmf_response,
                                       auto_load=auto_load)

    # Set energy units if applicable
    #
    _output_banner(out, "Set Energy or Wave Units")
    rate = "rate" if pha.rate else "counts"
    cmd = f'set_analysis({cmd_id}, quantity="{pha.units}", ' + \
        f'type="{rate}", factor={pha.plot_fac})'
    _output(out, cmd)

    if pha.subtracted:
        _output(out, f"subtract({cmd_id})")


def _save_dataset_settings_2d(out: OutType,
                              state: SessionType,  # currently unused
                              data: DataIMG,
                              id: IdType) -> None:
    """What settings need to be set for Data2D/IMG?"""

    _output_banner(out, "Set Image Coordinates")
    _output(out, f"set_coord({_id_to_str(id)}, '{data.coord}')")


def _save_data(out: OutType,
               state: SessionType,
               auto_load: bool
               ) -> None:
    """Save the data.

    This can just be references to files, or serialization of
    the data (or both).

    Parameters
    ----------
    out : dict
       The output state
    state
    auto_load : bool
       If ``False`` then the output will contain `load_arf`,
       `load_rmf`, and `load_bkg` calls for ancillary PHA files.

    Notes
    -----
    This does not - at least at present - save data to one or
    more files to be read in by the script. If any data needs
    to be serialized it is included in the script.
    """

    ids = state.list_data_ids()
    if len(ids) == 0:
        return

    _output_banner(out, "Load Data Sets")

    # There are two sets of data ids:
    #   a) state.list_data_ids()
    #   b) the keys of state._load_data_store
    #
    # Entries in b can be converted to a simple 'load_xxx' call [*]
    # whereas those only in set a require manual recreation.
    #
    # [*] This **assumes** that the data values have not been modified
    #     after being read in.
    #
    # PHA2 datasets also complicate this, as a single load call will
    # create multiple datasets (with a potentially different id
    # value), and then some of these may get deleted.  We process
    # these first, using state._multi_data_store.
    #
    for store, mids in state._multi_data_store.values():

        # Write out a banner line indicating all the ids
        # (we only do this once, which also means we only get the
        # load_pha/data call once).
        #
        msg = f"# Load PHA2 into: {mids}"
        if msg in out["main"]:
            continue

        _output(out, msg)
        _save_dataset_store(out, store)

        # Do we need to delete any of the PHA2 datasets?
        #
        for idval in mids:
            if idval not in ids:
                idstr = _id_to_str(idval)
                _output(out, f"delete_data({idstr})")

    for idval in ids:
        data = state.get_data(idval)
        if TYPE_CHECKING:
            # Assert an actual type rather than the base type of Data
            assert isinstance(data, (Data1D, Data2D))

        _save_dataset(out, state, data, idval)

        if isinstance(data, DataPHA):
            _save_dataset_settings_pha(out, state, data, idval,
                                       auto_load=auto_load)

        elif isinstance(data, DataIMG):
            _save_dataset_settings_2d(out, state, data, idval)

        _handle_filter(out, state, data, idval)


def _print_par(par: ParameterType) -> tuple[str, str]:
    """Convert a Sherpa parameter to a string.

    Note that we have to be careful with XSParameter parameters,
    to see if the hard limits need updating.

    Parameters
    ----------
    par : Parameter
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
            parstrs.append(f'{par.fullname}.hard_min    = {par.hard_min}')
    except AttributeError:
        pass

    try:
        if par.hard_max_changed():
            parstrs.append(f'{par.fullname}.hard_max    = {par.hard_max}')
    except AttributeError:
        pass

    parstrs.extend([f'{par.fullname}.default_val = {par.default_val}',
                    f'{par.fullname}.default_min = {par.default_min}',
                    f'{par.fullname}.default_max = {par.default_max}',
                    f'{par.fullname}.val     = {par.val}',
                    f'{par.fullname}.min     = {par.min}',
                    f'{par.fullname}.max     = {par.max}',
                    f'{par.fullname}.units   = {unitstr}',
                    f'{par.fullname}.frozen  = {par.frozen}'])
    parstr = '\n'.join(parstrs) + '\n'
    return (parstr, linkstr)


def _save_statistic(out: OutType, state: SessionType) -> None:
    """Save the statistic settings.

    Parameters
    ----------
    out : dict
       The output state
    state
    """

    _output_banner(out, "Set Statistic")
    _output(out, f'set_stat("{state.get_stat_name()}")')
    _output_nl(out)


def _save_fit_method(out: OutType, state: SessionType) -> None:
    """Save the fit method settings.

    Parameters
    ----------
    out : dict
       The output state
    state
    """

    # Save fitting method

    _output_banner(out, "Set Fitting Method")
    _output(out, f'set_method("{state.get_method_name()}")')
    _output_nl(out)

    def tostatement(key, val):
        # TODO: Using .format() returns more decimal places, which
        # is probably what we want but is a change, so leave
        # for now.
        # return 'set_method_opt("{}", {})'.format(key, val)
        return 'set_method_opt("%s", %s)' % (key, val)

    _save_entries(out, state.get_method_opt(), tostatement)
    _output_nl(out)


def _save_iter_method(out: OutType, state: SessionType) -> None:
    """Save the iterated-fit method settings, if any.

    Parameters
    ----------
    out : dict
       The output state
    state
    """

    iname = state.get_iter_method_name()
    if iname == 'none':
        return

    _output_banner(out, "Set Iterative Fitting Method")
    cmd = f'set_iter_method("{iname}")'
    _output(out, cmd)
    _output_nl(out)

    def tostatement(key, val):
        # There was a discussion here about the use of
        # str(val) vs val - the number of decimal places - but it
        # turns out for the test we only output integer values
        # so it makes no difference.
        return f'set_iter_method_opt("{key}", {val})'

    _save_entries(out, state.get_iter_method_opt(), tostatement)
    _output_nl(out)


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
def _handle_usermodel(out: OutType,
                      mod: UserModel,
                      modelname: str) -> None:

    try:
        pycode = inspect.getsource(mod.calc)
    except IOError:
        pycode = None

    # in case getsource can return None, have check here
    if pycode is None:
        msg = "Unable to save Python code for user model " + \
              f"'{modelname}' function {mod.calc.__name__}"
        warning(msg)
        _output(out, f'print("{msg}")')
        _output(out, f"def {mod.calc.__name__}(*args):")
        _output(out, "raise NotImplementedError('User model was "
                "not saved by save_all().'", indent=1)
        _output_nl(out)
        return

    msg = f"Found user model '{modelname}'; " + \
          "please check it is saved correctly."
    warning(msg)

    # Ensure the message is also seen if the script is run.
    _output(out, f'print("{msg}")')

    _output(out, textwrap.dedent(pycode))
    cmd = f'load_user_model({mod.calc.__name__}, "{modelname}")'
    _output(out, cmd)

    # Work out the add_user_pars call; this is explicit, i.e.  it does
    # not include logic to work out what arguments are not needed.
    #
    # Some of these values are over-written later on, but needed to
    # set up the number of parameters, and good documentation
    # (hopefully).
    #
    # Explicitly cast to float to avoid "np.float64(...)" output.
    #
    parnames = [p.name for p in mod.pars]
    parvals = [float(p.default_val) for p in mod.pars]
    parmins = [float(p.min) for p in mod.pars]
    parmaxs = [float(p.max) for p in mod.pars]
    parunits = [p.units for p in mod.pars]
    parfrozen = [p.frozen for p in mod.pars]

    spaces = '              '
    _output(out, f'add_user_pars("{modelname}",')
    _output(out, f"{spaces}parnames={parnames},")
    _output(out, f"{spaces}parvals={parvals},")
    _output(out, f"{spaces}parmins={parmins},")
    _output(out, f"{spaces}parmaxs={parmaxs},")
    _output(out, f"{spaces}parunits={parunits},")
    _output(out, f"{spaces}parfrozen={parfrozen}")
    _output(out, f"{spaces})\n")


def _save_model_components(out: OutType, state: SessionType) -> bool:
    """Save the model components.

    Parameters
    ----------
    out : dict
       The output state
    state

    Returns
    -------
    xspec_state : bool
       True if any XSPEC models are found.

    """

    found_xspec = False

    # Try to be careful about the ordering of the components here.
    #
    all_model_components = state.list_model_components()
    if len(all_model_components) == 0:
        return found_xspec

    all_model_components.reverse()
    _output_banner(out, "Set Model Components and Parameters")

    # If there are any links between parameters, store link commands here
    # Then, *after* processing all models in the for loop below, send
    # link commands to outfile -- *all* models need to be created before
    # *any* links between parameters can be established.
    linkstr = ""
    for modval in all_model_components:

        # get actual model instance from the name we are given
        # then get model type, and name of this instance.
        mod = eval(modval)
        typename = mod.type
        modelname = mod.name.split(".")[1]

        if typename == "usermodel":
            _handle_usermodel(out, mod, modelname)

        elif typename == "psfmodel":
            cmd = f'load_psf("{mod._name}", "{mod.kernel.name}")'
            _output(out, cmd)

        elif typename == "tablemodel":
            cmd = f'load_table_model("{modelname}", "{mod.filename}")'
            _output(out, cmd)

        elif typename == "xstablemodel":
            cmd = f'load_xstable_model("{modelname}", "{mod.filename}"'
            if mod.etable:
                cmd += ', etable=True'

            cmd += ')'
            _output(out, cmd)
            found_xspec = True

        else:
            # Normal case:  create an instance of the model.
            cmd = f'create_model_component("{typename}", "{modelname}")'
            _output(out, cmd)

            # Is this an XSPEC model? We only care if it's the first
            # one we find.
            #
            if not found_xspec and xspec is not None:
                found_xspec = isinstance(mod, xspec.XSModel)

        # QUS: should this be included in the above checks?
        #      @DougBurke doesn't think so, as the "normal
        #      case" above should probably be run , but there's
        #      no checks to verify this.
        #
        if typename == "convolutionkernel":
            # Create general convolution kernel with load_conv
            cmd = f'load_conv("{modelname}", "{mod.kernel.name}")'
            _output(out, cmd)

        if hasattr(mod, "integrate"):
            cmd = f"{modelname}.integrate = {mod.integrate}"
            _output(out, cmd)
            _output_nl(out)

        # Write out the parameters in the order they are stored in
        # the model. The original version of the code iterated
        # through mod.__dict__.values() and picked out Parameter
        # values.
        #
        for par in mod.pars:
            par_attributes, par_linkstr = _print_par(par)
            _output(out, par_attributes)
            linkstr = linkstr + par_linkstr

        # If the model is a PSFModel, could have special
        # attributes "size" and "center" -- if so, record them.
        if typename == "psfmodel":
            spacer = False
            if hasattr(mod, "size") and mod.size is not None:
                _output(out, f"{modelname}.size = {mod.size}")
                spacer = True

            if hasattr(mod, "center") and mod.center is not None:
                _output(out, f"{modelname}.center = {mod.center}")
                spacer = True

            if spacer:
                _output_nl(out)

    # If there were any links made between parameters, send those
    # link commands to outfile now; else, linkstr is just an empty string
    _output(out, linkstr)
    return found_xspec


def _save_psf_components(out: OutType, state: SessionType) -> None:
    """Save the PSF components.

    This must be done after setting up the data and model components.

    Parameters
    ----------
    out : dict
       The output state
    state

    """

    # This is done after creating all the models and datasets.
    #
    if len(state._psf) == 0:
        return

    # Associate any PSF models with their appropriate datasets.
    #
    _output_banner(out, "Associate PSF models with the datasets")
    for idval, psfmod in state._psf.items():
        cmd_id = _id_to_str(idval)
        _output(out, f"set_psf({cmd_id}, {psfmod._name})")


def _save_models(out: OutType, state: SessionType) -> None:
    """Save the source, pileup, and background models.

    Parameters
    ----------
    out : dict
       The output state
    state
    """

    ids = state.list_data_ids()
    if len(ids) == 0:
        return

    _output_banner(out, "Set Source, Pileup and Background Models")
    orig_pos = _get_out_pos(out)

    for id in ids:
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
                continue

            _output(out, cmd)
            _output_nl(out)
        except:
            pass

        # Set background models (if any) associated with backgrounds
        # tied to this data set -- if none, then pass.  Again, try
        # to distinguish cases where background "source" model is
        # different from whole background model.
        try:
            bids = state.list_bkg_ids(id)
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
                    continue

                _output(out, cmd)
                _output_nl(out)

        except:
            pass

    # separate out the pileup models from the source models
    #
    for id in ids:
        try:
            pname = state.get_pileup_model(id).name
        except:
            continue

        cmd_id = _id_to_str(id)
        cmd = f"set_pileup_model({cmd_id}, {pname})"
        _output(out, cmd)

    # Do want to remove the initial banner?
    #
    _remove_banner(out, orig_pos)


def _save_xspec(out: OutType) -> None:
    """Save the XSPEC settings.

    Parameters
    ----------
    out : dict
       The output state
    """

    # This case should not happen, but just in case.
    #
    if xspec is None:
        return

    _output_banner(out, "XSPEC Module Settings")
    xspec_state = xspec.get_xsstate()

    chatter =  xspec_state["chatter"]
    abund = xspec_state["abund"]
    xsect = xspec_state["xsect"]
    cs = xspec_state["cosmo"]
    _output(out, f"set_xschatter({chatter})")
    _output(out, f'set_xsabund("{abund}")')
    _output(out, f"set_xscosmo({cs[0]:g}, {cs[1]:g}, {cs[2]:g})")
    _output(out, f'set_xsxsect("{xsect}")')

    def tostatement(key, val):
        return f'set_xsxset("{key}", "{val}")'

    _save_entries(out, xspec_state["modelstrings"], tostatement)


def _save_dataset_store(out: OutType,
                        store: FileStore
                        ) -> None:
    """The data can be read in from a file."""

    outstr = store.show()
    if outstr in out["main"]:
        # This indicates this is from a PHA2 file which has already
        # been loaded.
        return

    _output(out, outstr)


def _output_wcs_import(out: OutType) -> None:
    """Import the WCS symbol if not done already."""

    out["imports"].add("from sherpa.astro.io.wcs import WCS")


def _output_add_wcs(out: OutType,
                    idstr: str,
                    label: str,
                    wcs: WCS) -> None:
    """Copy over the WCS"""

    _output(out, f"get_data({idstr}).{label} = WCS('{wcs.name}', '{wcs.type}',")
    _output(out, f"crval={wcs.crval.tolist()},", indent=1)
    _output(out, f"crpix={wcs.crpix.tolist()},", indent=1)
    if wcs.type == "LINEAR":
        _output(out, f"cdelt={wcs.cdelt.tolist()})", indent=1)
    else:
        _output(out, f"cdelt={wcs.cdelt.tolist()},", indent=1)
        _output(out, f"crota={wcs.crota}, epoch={wcs.epoch}, equinox={wcs.equinox})", indent=1)


def _save_dataset_pha_manual(out: OutType, idstr: str, pha: DataPHA) -> None:
    """Try to recreate the PHA"""

    spacer = "            "
    _output(out, f'load_arrays({idstr},')
    _output(out, f"{spacer}{pha.channel.tolist()},")
    _output(out, f"{spacer}{pha.counts.tolist()},")
    _output(out, f"{spacer}DataPHA)")

    def setval(key):
        val = getattr(pha, key)
        if val is None:
            return

        _output(out, f"set_{key}({idstr}, {val})")

    def setarray(key):
        val = getattr(pha, key)
        if val is None:
            return

        _output(out, f"set_{key}({idstr}, {val.tolist()})")

    setval("exposure")
    setval("backscal")
    setval("areascal")

    setarray("staterror")
    setarray("syserror")

    # These are set later so do not need to be set here.
    #
    # setarray("quality")
    # setarray("grouping")


def _save_dataset(out: OutType,
                  state: SessionType,
                  data: Data,
                  id: IdType) -> None:
    """Given a dataset identifier, return the text needed to
    re-create it.

    The state._load_data_store dictionary is intended to record how
    a file was loaded, but there is no tracking of the data values
    (e.g. the independent and dependent axes) to know if they have
    been changed after loading. If there is no information in this
    dictionary then the dataset is assumed to be manually created.

    """

    store = state._load_data_store.get(id, None)
    if store is not None:
        _save_dataset_store(out, store)
        return

    idstr = _id_to_str(id)

    # We could use dataspace1d/2d but easiest to just use load_arrays.
    # The isinstance checks have to pick the more-specific classes
    # first (e.g. DataPHA and Data1DInt before Data1D, DataIMG before
    # Data2D).
    #
    # Note that for DataIMG classes we use the "hidden"
    # _orig_indep_axis attribute rather than get_indep(), since the
    # latter will change when the coord system is not logical, and
    # using this would cause downstream problems. See PR #1414 for
    # more details. This should also be done for DataIMGInt, but this
    # class is poorly tested, documented, or used.
    #
    if isinstance(data, DataPHA):
        _save_dataset_pha_manual(out, idstr, data)
        return

    if isinstance(data, DataIMG):
        # This assumes that orig[0] == "logical"
        orig = data._orig_indep_axis
        xs = (orig[1], orig[2])
    else:
        xs = data.get_indep()

    ys = f"{data.get_dep().tolist()}"

    spacer = "            "

    if isinstance(data, Data1DInt):
        _output(out, f'load_arrays({idstr},')
        _output(out, f"{spacer}{xs[0].tolist()},")
        _output(out, f"{spacer}{xs[1].tolist()},")
        _output(out, f"{spacer}{ys},")
        _output(out, f"{spacer}Data1DInt)")

    elif isinstance(data, Data1D):
        _output(out, f'load_arrays({idstr},')
        _output(out, f"{spacer}{xs[0].tolist()},")
        _output(out, f"{spacer}{ys},")
        _output(out, f"{spacer}Data1D)")

    elif isinstance(data, DataIMG):

        # Note that we set transforms outside the load_arrays call
        if data.sky is not None or data.eqpos is not None:
            _output_wcs_import(out)

        # This does not save the header information in the file
        _output(out, f'load_arrays({idstr},')
        _output(out, f"{spacer}{xs[0].tolist()},")
        _output(out, f"{spacer}{xs[1].tolist()},")
        _output(out, f"{spacer}{ys},")
        if data.shape is not None:
            _output(out, f"{spacer}({data.shape[0]}, {data.shape[1]}),")

        _output(out, f"{spacer}DataIMG)")

        # Note that we set transforms outside the load_arrays call
        if data.sky is not None:
            _output_add_wcs(out, idstr, "sky", data.sky)

        if data.eqpos is not None:
            _output_add_wcs(out, idstr, "eqpos", data.eqpos)

    elif isinstance(data, Data2DInt):
        msg = f"Unable to re-create Data2DInt data set '{idstr}'"
        warning(msg)
        _output(out, f'print("{msg}")')

    elif isinstance(data, Data2D):
        _output(out, f'load_arrays({idstr},')
        _output(out, f"{spacer}{xs[0].tolist()},")
        _output(out, f"{spacer}{xs[1].tolist()},")
        _output(out, f"{spacer}{ys},")
        if data.shape is not None:
            _output(out, f"{spacer}{tuple(data.shape)},")

        _output(out, f"{spacer}Data2D)")

    else:
        msg = f"Unable to re-create {data.__class__} data set '{idstr}'"
        warning(msg)
        _output(out, f'print("{msg}")')
        return

    staterr = data.get_staterror()
    syserr = data.get_syserror()

    if staterr is not None:
        _output(out, f"set_staterror({idstr}, {staterr.tolist()})")

    if syserr is not None:
        _output(out, f"set_syserror({idstr}, {syserr.tolist()})")


def _save_id(out: OutType, state: SessionType) -> None:
    "Does the default id need changing?"

    defid = state.get_default_id()
    if defid == 1:
        return

    _output_banner(out, "Default identifier")
    _output(out, f'set_default_id({repr(defid)})')
    _output_nl(out)


def save_all(state: SessionType,
             fh: TextIO | None = None,
             *,
             auto_load: bool = True
             ) -> None:
    """Save the information about the current session to a file handle.

    This differs to the `save` command in that the output is human
    readable. Three consequences are:

     1. numeric values may not be recorded to their full precision

     2. data sets are not included in the file

     3. some settings and values may not be recorded.

     .. versionchanged:: 4.18.0
        Handling of PHA data has been improved, and the output now
        defaults to not including automatically-loaded ancillary files
        (such as background and responses). Set the ``auto_load``
        flag to False to add these commands back to the save file (and
        match previous versions of Sherpa).

    Parameters
    ----------
    state : sherpa.astro.ui.utils.Session
    fh : file_like, optional
       If not given the results are displayed to standard out,
       otherwise they are written to this file handle.
    auto_load : bool, optional
       If ``False`` then the output will contain `load_arf`,
       `load_rmf`, and `load_bkg` calls for ancillary PHA files.

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

    - any optional keywords to commands such as `load_data`

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

    # Record the various elements of the output.
    #
    out: OutType
    out = { "imports": set(),
            "main": []
           }

    _save_data(out, state, auto_load=auto_load)
    _output_nl(out)
    _save_statistic(out, state)
    _save_fit_method(out, state)
    _save_iter_method(out, state)
    req_xspec = _save_model_components(out, state)
    _save_psf_components(out, state)
    _save_models(out, state)
    if req_xspec:
        _save_xspec(out)

    _save_id(out, state)

    outfh = sys.stdout if fh is None else fh

    def write(msg: str) -> None:
        outfh.write(f"{msg}\n")

    # Hard code the required imports.
    #
    write("import numpy")
    write("from sherpa.astro.ui import *")

    # Add additional imports.
    #
    for iline in out["imports"]:
        write(iline)

    # Write out the contents.
    #
    for txt in out["main"]:
        write(txt)
