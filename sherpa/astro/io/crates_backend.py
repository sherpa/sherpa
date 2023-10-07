#
#  Copyright (C) 2011, 2015, 2016, 2019 - 2024
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

from collections import defaultdict
import logging
import os
from typing import Any, Optional, Sequence, Union

import numpy as np

from pycrates import CrateDataset, IMAGECrate, TABLECrate  # type: ignore
import pycrates  # type: ignore
import pytransform  # type: ignore

from sherpa.utils import is_binary_file
from sherpa.utils.err import ArgumentTypeErr, IOErr
from sherpa.utils.numeric_types import SherpaFloat

from .types import ColumnsType, DataType, DataTypeArg, \
    HdrType, HdrTypeArg, KeyType, NamesType, BlockType, \
    HeaderItem, Header, Column, ImageBlock, TableBlock, BlockList

warning = logging.getLogger(__name__).warning
error = logging.getLogger(__name__).error
info = logging.getLogger(__name__).info

try:
    from .wcs import WCS
    HAS_TRANSFORM = True
except ImportError:
    HAS_TRANSFORM = False
    warning('failed to import WCS module; WCS routines will not be '
            'available')

__all__ = ('get_table_data', 'get_header_data', 'get_image_data',
           'get_column_data', 'get_ascii_data',
           'get_arf_data', 'get_rmf_data', 'get_pha_data',
           'pack_table_data', 'pack_image_data', 'pack_pha_data',
           'pack_arf_data', 'pack_rmf_data', 'pack_hdus',
           'set_table_data', 'set_image_data', 'set_pha_data',
           'set_arf_data', 'set_rmf_data', 'set_hdus')


name: str = "crates"
"""The name of the I/O backend."""


CrateType = Union[TABLECrate, IMAGECrate]
DatasetType = Union[str, CrateDataset]


def open_crate(filename: str,
               mode: str = "r") -> CrateType:
    """Get the "most interesting" block of the file.

    Parameters
    ----------
    filename : str
    mode : str

    Returns
    -------
    crate : IMAGECrate or TABLECrate

    Notes
    -----
    This routine could just be inlined where necessary, but it is
    useful because it allows the code to provide unified behaviour
    when crates can not read a file because of an unknown datatype.

    """
    try:
        dataset = CrateDataset(filename, mode=mode)
    except OSError as oe:
        raise IOErr('openfailed', f"unable to open {filename}: {oe}") from oe

    current = dataset.get_current_crate()
    try:
        return dataset.get_crate(current)
    except Exception as exc:
        # crates can fail with unsupported data types - e.g.  see
        # #1898 but also long-long columns for CIAO 4.15.
        #
        emsg = f"crates is unable to read '{filename}': {str(exc)}"
        raise IOErr("openfailed", emsg) from exc


def get_filename_from_dmsyntax(filename: str) -> str:
    """What is the filename to use?

    When given a non-binary file the routine will attempt to determine
    what the column names are, to decide whether the CIAO DataModel
    syntax "[opt colnames=none]" needs to be added to the file name.
    Other parts of the CIAO virtual file specifier in filename will
    be dropped.

    .. versionchanged:: 4.17.0
       The code used to check whether the `colnames=none` option is
       needed has been revamped and may behave differently with a file
       using '%' as the comment character (support for such a file has
       not been documented so it is not clear what the behavior should
       be).

    """

    # Is this routine still useful? For example, the column-name check
    # doesn't work for TEXT/DTF or TEXT/TSV file formats.
    #
    out = filename
    dmsyn = ''

    if '[' in filename and ']' in filename:
        parts = filename.split('[')
        filename = parts.pop(0)
        if parts:
            dmsyn = parts.pop(0).lower()

    if not is_binary_file(filename):
        colnames = True

        # Crates only guarantees to support ASCII-like variants but it
        # should be okay to use a less-restrictive encoding here.
        #
        with open(filename, mode='r', encoding="UTF8") as fh:
            last = None
            line = fh.readline().strip()
            while len(line) > 0 and line.startswith("#"):
                last = line
                line = fh.readline().strip()

        # Check whether we have a header line and an actual data line
        # and, if so, whether they match.
        #
        if last is not None:
            cols_header = last.strip('#').strip().split(' ')
            cols_data = line.strip().split(' ')
            if len(cols_header) != len(cols_data):
                colnames = False

        if not colnames and 'cols' not in dmsyn:
            out += "[opt colnames=none]"

    return out


def _try_key(crate: CrateType,
             name: str,
             dtype: type = str) -> Optional[KeyType]:
    """Access the key from the crate returning None if it does not exist.

    There's no way to differentiate between a key that does not exist
    and a key that is set but has been set to a value of None (which
    is unlikely as the ASCII and FITS serializations do not support
    this).

    Parameters
    ----------
    crate : TABLECrate or IMAGECrate
    name : str
       The name is case insensitive.
    dtype

    Returns
    -------
    value
       The key value or None.

    """

    # The get_key method will raise an IndexError if sent an integer
    # outside the valid range, but here name is a string, where it
    # returns None. This is inconsistent behaviour in crates circa
    # CIAO 4.15.
    #
    key = crate.get_key(name)
    if key is None:
        return None

    value = key.value

    # Note that this check seems a bit too permissive, since a value
    # of "This key is not set to none" will be matched.
    #
    if str(value).find('none') != -1:
        return None

    # byte strings are to be decoded to strings
    #
    if dtype == str and isinstance(value, bytes):
        return dtype(value, "utf-8")

    return dtype(value)


# This is being replaced by _get_meta_data_header
#
def _get_meta_data(crate: CrateType) -> HdrType:
    """Retrieve the header information from the crate.

    This loses the "specialized" keywords like HISTORY and COMMENT.

    """

    meta = {}
    for name in crate.get_keynames():
        val = crate.get_key(name).value
        meta[name] = val

    return meta


def _get_meta_data_header(crate: CrateType) -> Header:
    """Retrieve the header information from the crate.

    This loses the "specialized" keywords like HISTORY and COMMENT.

    """

    out = []
    for name in crate.get_keynames():
        key = crate.get_key(name)
        item = HeaderItem(name=name, value=key.value)
        item.unit = key.unit if key.unit != '' else None
        item.desc = key.desc if key.desc != '' else None
        out.append(item)

    return Header(out)


def _set_key(crate: CrateType,
             name: str,
             val: KeyType,
             fix_type: bool = False,
             dtype: type = str) -> None:
    """Add a key + value to the header."""

    key = pycrates.CrateKey()
    key.name = name.upper()
    key.value = dtype(val) if fix_type else val
    crate.add_key(key)


# The typing for this is not ideal.
#
def _set_column(crate: TABLECrate,
                name: str,
                val: Any) -> None:
    """Add a column to the table crate."""

    col = pycrates.CrateData()
    col.name = name.upper()
    col.values = np.array(val)
    crate.add_column(col)


def _require_key(crate: CrateType,
                 name: str,
                 dtype: type = str) -> KeyType:
    """Access the key from the crate, erroring out if it does not exist.

    There's no way to differentiate between a key that does not exist
    and a key that is set but has been set to a value of None (which
    is unlikely as the ASCII and FITS serializations do not support
    this).

    Parameters
    ----------
    crate : TABLECrate or IMAGECrate
    name : str
       The name is case insensitive.
    dtype

    Returns
    -------
    value
       The key value

    Raises
    ------
    IOErr
       The key does not exist.

    """
    key = _try_key(crate, name, dtype)
    if key is None:
        raise IOErr('nokeyword', crate.get_filename(), name)
    return key


def _try_col(crate: TABLECrate,
             colname: str,
             make_copy: bool = False,
             fix_type: bool = False,
             dtype: type = SherpaFloat) -> Optional[np.ndarray]:
    """Return the column data or None (if it does not exist)."""

    try:
        return _require_col(crate, colname, make_copy=make_copy,
                            fix_type=fix_type, dtype=dtype)
    except IOErr:
        return None


# This is surprisingly hard to type (because of column_stack) so skip
# it.
#
def _require_tbl_col(crate, colname, cnames, make_copy=False,
                     fix_type=False):
    """
    checked for new crates
    """

    try:
        col = crate.get_column(colname)
    except ValueError as ve:
        raise IOErr('reqcol', colname, cnames) from ve

    if make_copy:
        # Make a copy if a filename passed in
        data = np.array(col.values)
    else:
        # Use a reference if a crate passed in
        data = np.asarray(col.values)

    if fix_type:
        data = data.astype(SherpaFloat)

    return np.column_stack(data)


def _try_key_list(crate: CrateType,
                  keyname: str,
                  num: int,
                  dtype: type = SherpaFloat,
                  fix_type: bool = False) -> Optional[np.ndarray]:
    """
    checked for new crates
    """
    if not crate.key_exists(keyname):
        return None

    key = crate.get_key(keyname)
    keys = np.full(num, key.value)
    if fix_type:
        keys = keys.astype(dtype)

    return keys


def _try_col_list(crate: TABLECrate,
                  colname: str,
                  num: int,
                  make_copy: bool = False,
                  fix_type: bool = False) -> np.ndarray:
    """
    checked for the new crates
    """

    try:
        return _require_col_list(crate, colname,
                                 make_copy=make_copy,
                                 fix_type=fix_type)
    except IOErr:
        return np.full(num, None)


def _require_col_list(crate: TABLECrate,
                      colname: str,
                      make_copy: bool = False,
                      fix_type: bool = False) -> np.ndarray:
    """
    check for new crates
    """

    try:
        col = crate.get_column(colname)
    except ValueError:
        raise IOErr("reqcol", colname, crate.get_filename()) from None

    if col.is_varlen():
        values = col.get_fixed_length_array()
    else:
        values = col.values

    if make_copy:
        # Make a copy if a filename passed in
        col = np.array(values)
    else:
        # Use a reference if a crate passed in
        col = np.asarray(values)

    if fix_type:
        col = col.astype(SherpaFloat)

    return col


def _require_col(crate: TABLECrate,
                 colname: str,
                 make_copy: bool = False,
                 fix_type: bool = False,
                 label: Optional[str] = None,
                 dtype: type = SherpaFloat) -> np.ndarray:
    """Return the column data or error out."""

    try:
        col = crate.get_column(colname)
    except ValueError:
        msg = colname if label is None else label
        raise IOErr("reqcol", msg, crate.get_filename()) from None

    if col.is_varlen():
        values = col.get_fixed_length_array()
    else:
        values = col.values

    kwargs = {}
    if fix_type:
        kwargs['dtype'] = dtype

    return np.array(values, copy=make_copy, **kwargs).ravel()



def _require_image(crate: IMAGECrate,
                   filename: str,
                   make_copy: bool = False,
                   fix_type: bool = False) -> np.ndarray:
    """Return the image pixel data.

    Note that "extra" dimensions - e.g. radio cube data where all but
    the first two dimensions are 1 - will be removed by this call.

    """

    try:
        img = crate.get_image()
    except LookupError:
        raise IOErr("badimg", filename) from None

    if make_copy:
        dat = img.values.copy()
    else:
        dat = img.values

    # This is presumably to convert [n, m, 1, 1] to [n, m] but it's
    # not-at-all obvious why we need it (or, rather, we should
    # probably handle extra dimensions better).
    #
    out = dat.squeeze()
    if fix_type:
        out = out.astype(SherpaFloat)

    return out


def _get_crate_by_blockname(dataset: CrateDataset,
                            blkname: str) -> Optional[CrateType]:
    """Select the given block name.

    Parameters
    ----------
    dataset : CrateDataset
    blkname : str

    Returns
    -------
    mblock : TABLECrate, IMAGECrate, or None
    """

    # Access used to only be by number but it can now
    # also be by name, but we need to hide the error
    # for unknown blocks.
    #
    try:
        return dataset.get_crate(blkname)
    except IndexError:
        return None
    except Exception as exc:
        # crates can fail with unsupported data types - e.g.  see
        # #1898 but also long-long columns for CIAO 4.15.
        #
        filename = dataset.get_filename()
        emsg = f"crates is unable to read '{filename}': {str(exc)}"
        raise IOErr("openfailed", emsg) from exc


# Read Functions #

def read_table_blocks(arg: Union[str, CrateDataset, TABLECrate],
                      make_copy: bool = False
                      ) -> tuple[str,
                                 dict[int, ColumnsType],
                                 dict[int, HdrType]
                                 ]:
    """Read in tabular data with no restrictions on the columns."""

    dataset = None
    if isinstance(arg, TABLECrate):
        filename = arg.get_filename()
        dataset = arg.get_dataset()
    elif isinstance(arg, CrateDataset):
        filename = arg.get_filename()
        dataset = arg
    elif isinstance(arg, str):
        filename = arg
        try:
            dataset = CrateDataset(arg)
        except OSError as oe:
            raise IOErr('openfailed', f"unable to open {arg}: {oe}") from oe

    else:
        raise IOErr('badfile', arg, "CrateDataset obj")

    cols: dict[int, ColumnsType] = {}
    hdr: dict[int, HdrType] = {}
    for idx in range(1, dataset.get_ncrates() + 1):
        crate = dataset.get_crate(idx)
        hdr[idx] = {}
        names = crate.get_keynames()
        for name in names:
            hdr[idx][name] = _try_key(crate, name)

        cols[idx] = {}
        # skip over primary
        if crate.name == 'PRIMARY':
            continue

        names = crate.get_colnames()
        for name in names:
            cols[idx][name] = crate.get_column(name).values

    return filename, cols, hdr


def get_header_data(arg: Union[str, TABLECrate],
                    blockname: Optional[str] = None,
                    hdrkeys: Optional[NamesType] = None
                    ) -> HdrType:
    """Read the metadata."""

    if isinstance(arg, str):
        arg = get_filename_from_dmsyntax(arg)
        tbl = open_crate(arg)
    elif isinstance(arg, TABLECrate):
        tbl = arg
    else:
        raise IOErr('badfile', arg, 'TABLECrate obj')

    # Crates "caches" open files by their filename in memory.  If you try
    # to open a file multiple times (with DM syntax) it corrupts the Crate
    # in memory.  This is a work-around to open the CrateDataset without
    # DM syntax and iterate through the crates looking for the block
    # name that matches.
    if blockname is not None:
        crate = _get_crate_by_blockname(tbl.get_dataset(), blockname)
        tbl = crate or tbl

    hdr = {}
    if hdrkeys is None:
        hdrkeys = tbl.get_keynames()

    for key in hdrkeys:
        hdr[key] = _require_key(tbl, key)

    return hdr


def get_column_data(*args) -> list[np.ndarray]:
    """Extract the column data."""

    # args is passed as type list
    if len(args) == 0:
        raise IOErr('noarrays')

    cols = []
    for arg in args:
        if isinstance(arg, pycrates.CrateData):
            vals = arg.values
        elif arg is None or isinstance(arg, (np.ndarray, list, tuple)):
            vals = arg
        else:
            raise IOErr('badarray', arg)

        if arg is not None:
            vals = np.asanyarray(vals)
            for col in np.atleast_2d(vals.T):
                cols.append(col)
        else:
            cols.append(vals)

    return cols


def get_ascii_data(filename: str,
                   ncols: int = 2,
                   colkeys: Optional[NamesType] = None,
                   **kwargs
                   ) -> tuple[list[str], list[np.ndarray], str]:
    """Read columns from an ASCII file"""
    return get_table_data(filename, ncols, colkeys)[:3]


def get_table_data(arg: Union[str, TABLECrate],
                   ncols: int = 1,
                   colkeys: Optional[NamesType] = None,
                   make_copy: bool = True,
                   fix_type: bool = True,
                   blockname: Optional[str] = None,
                   hdrkeys: Optional[NamesType] = None
                   ) -> tuple[list[str], list[np.ndarray], str, HdrType]:
    """Read columns from a file or crate."""

    if isinstance(arg, str):
        arg = get_filename_from_dmsyntax(arg)
        tbl = open_crate(arg)
        if not isinstance(tbl, TABLECrate):
            raise IOErr('badfile', arg, 'TABLECrate obj')

        filename = tbl.get_filename()

    elif isinstance(arg, TABLECrate):
        tbl = arg
        filename = arg.get_filename()
        make_copy = False
    else:
        raise IOErr('badfile', arg, 'TABLECrate obj')

    # Crates "caches" open files by their filename in memory.  If you try
    # to open a file multiple times (with DM syntax) it corrupts the Crate
    # in memory.  This is a work-around to open the CrateDataset without
    # DM syntax and iterate through the crates looking for the block
    # name that matches.
    #
    # Is this still a valid issue with crates?
    if blockname is not None:
        crate = _get_crate_by_blockname(tbl.get_dataset(), blockname)
        tbl = crate or tbl

    cnames = list(pycrates.get_col_names(tbl, vectors=False, rawonly=True))

    if colkeys is not None:
        colkeys = [str(name).strip() for name in list(colkeys)]

    elif (isinstance(arg, str) and (not os.path.isfile(arg))
          and '[' in arg and ']' in arg):
        colkeys = cnames

    # Try Channel, Counts or X,Y before defaulting to first two table cols
    elif 'CHANNEL' in cnames and 'COUNTS' in cnames:
        colkeys = ['CHANNEL', 'COUNTS']

    elif 'X' in cnames and 'Y' in cnames:
        colkeys = ['X', 'Y']

    else:
        colkeys = cnames[:ncols]

    cols = []
    for name in colkeys:
        for col in _require_tbl_col(tbl, name, cnames,
                                    make_copy=make_copy,
                                    fix_type=fix_type):
            cols.append(col)

    hdr = {}
    if hdrkeys is not None:
        for key in hdrkeys:
            hdr[key] = _require_key(tbl, key)

    return colkeys, cols, filename, hdr


def get_image_data(arg: Union[str, IMAGECrate],
                   make_copy: bool = True,
                   fix_type: bool = True
                   ) -> tuple[DataType, str]:
    """Read image data from a file or crate"""

    if isinstance(arg, str):
        img = open_crate(arg)
        if not isinstance(img, IMAGECrate):
            raise IOErr('badfile', arg, "IMAGECrate obj")

        filename = arg

    elif isinstance(arg, IMAGECrate):
        img = arg
        filename = arg.get_filename()
        make_copy = False

    else:
        raise IOErr('badfile', arg, "IMAGECrate obj")

    data: DataType = {}

    data['y'] = _require_image(img, filename, make_copy=make_copy,
                               fix_type=fix_type)

    if HAS_TRANSFORM:
        sky = None
        skynames = ['SKY', 'sky', 'pos', 'POS']
        names = img.get_axisnames()

        # find the SKY name using the set intersection
        inter = list(set(names) & set(skynames))
        if inter:
            sky = img.get_transform(inter[0])

        wcs = None
        if 'EQPOS' in names:
            wcs = img.get_transform('EQPOS')

        if sky is not None:
            linear = pytransform.WCSTANTransform()
            linear.set_name("LINEAR")
            linear.set_transform_matrix(sky.get_transform_matrix())
            cdelt = np.array(linear.get_parameter_value('CDELT'))
            crpix = np.array(linear.get_parameter_value('CRPIX'))
            crval = np.array(linear.get_parameter_value('CRVAL'))
            data['sky'] = WCS('physical', 'LINEAR', crval, crpix, cdelt)

        if wcs is not None:
            cdelt = np.array(wcs.get_parameter_value('CDELT'))
            crpix = np.array(wcs.get_parameter_value('CRPIX'))
            crval = np.array(wcs.get_parameter_value('CRVAL'))
            crota = SherpaFloat(wcs.get_parameter_value('CROTA'))
            equin = SherpaFloat(wcs.get_parameter_value('EQUINOX'))
            epoch = SherpaFloat(wcs.get_parameter_value('EPOCH'))
            data['eqpos'] = WCS('world', 'WCS', crval, crpix, cdelt,
                                crota, epoch, equin)

    data['header'] = _get_meta_data(img)
    for key in ['CTYPE1P', 'CTYPE2P', 'WCSNAMEP', 'CDELT1P',
                'CDELT2P', 'CRPIX1P', 'CRPIX2P', 'CRVAL1P', 'CRVAL2P',
                'EQUINOX']:
        data['header'].pop(key, None)

    return data, filename


def get_arf_data(arg: Union[str, TABLECrate],
                 make_copy: bool = True
                 ) -> tuple[TableBlock, str]:
    """Read an ARF from a file or crate"""

    if isinstance(arg, str):
        arf = open_crate(arg)
        if not isinstance(arf, TABLECrate):
            raise IOErr('badfile', arg, "ARFCrate obj")
        filename = arg
    elif isinstance(arg, TABLECrate):
        arf = arg
        filename = arg.get_filename()
        make_copy = False

    else:
        raise IOErr('badfile', arg, "ARFCrate obj")

    headers = _get_meta_data_header(arf)
    cols = [copycol(arf, filename, name, dtype=SherpaFloat)
            for name in ["ENERG_LO", "ENERG_HI", "SPECRESP"]]

    # Optional columns: either both given or none
    #
    try:
        bin_lo = copycol(arf, filename, "BIN_LO", dtype=SherpaFloat)
        bin_hi = copycol(arf, filename, "BIN_HI", dtype=SherpaFloat)
        cols.extend([bin_lo, bin_hi])
    except IOErr:
        pass

    return TableBlock(arf.name, header=headers, columns=cols), filename


# Commonly-used block names for the MATRIX block. Only the first two
# are given in the OGIP standard.
#
RMF_BLOCK_NAMES = ["MATRIX", "SPECRESP MATRIX", "AXAF_RMF", "RSP_MATRIX"]


def _find_matrix_blocks(filename: str,
                        ds: CrateDataset) -> list[str]:
    """Report the block names that contain MATRIX data in a RMF.

    Arguments
    ---------
    filename : str
    ds : CrateDataset
       This is expected to be a RMFCrateDataset.

    Returns
    -------
    blnames : list of str
       The block names that contain the MATRIX block. It will not be
       empty.

    Raises
    ------
    IOErr
       No matrix block was found

    """

    # The naming of the matrix block can be complicated, and perhaps
    # we should be looking for
    #
    #    HDUCLASS = OGIP
    #    HDUCLAS1 = RESPONSE
    #    HDUCLAS2 = RSP_MATRIX
    #
    # and check the HDUVERS keyword, but it's easier just to go on the
    # name (as there's no guarantee that these keywords will be any
    # cleaner to use). As of December 2020 there is now the
    # possibility of RMF files with multiple MATRIX blocks (where the
    # EXTVER starts at 1 and then increases).
    #
    blnames = []
    for crnum in range(1, ds.get_ncrates() + 1):
        cr = ds.get_crate(crnum)
        if cr.name in RMF_BLOCK_NAMES or \
           (cr.get_key_value("HDUCLAS1") == "RESPONSE" and
            cr.get_key_value("HDUCLAS2") == "RSP_MATRIX"):
            blnames.append(cr.name)

    if not blnames:
        raise IOErr('notrsp', filename, 'an RMF')

    return blnames


def copycol(cr: TABLECrate,
            filename: str,
            name: str,
            dtype: Optional[type] = None) -> Column:
    """Copy the column data (which must exist).

    Parameters
    ----------
    cr : TableCrate
       The crate to search.
    filename : str
       The filename from which the crate was created. This is used for
       error messages only.
    name : str
       The column name (case insensitive).
    dtype : None or type
       The data type to convert the column values to.

    Returns
    -------
    col : Column
       The column data.

    Raises
    ------
    IOErr
       The column does not exist.

    Notes
    -----

    Historically Sherpa has converted RMF columns to double-precision
    values, even though the standard has them as single-precision.
    Using the on-disk type (e.g. float rather than double) has some
    annoynig consequences on the test suite, so for now we retain this
    behaviour.

    """

    try:
        col = cr.get_column(name)
    except ValueError:
        raise IOErr("reqcol", name, filename) from None

    vals = col.values
    if dtype is not None:
        vals = vals.astype(dtype)

    # do not expand out variable-length arrays for now.
    #
    out = Column(col.name, values=vals)
    out.desc = col.desc if col.desc != '' else None
    out.unit = col.unit if col.unit != '' else None

    # Only set the min/max if they are informative (ie not just
    # the min/max for that datatype). Assume we do not have to
    # worry about bit columns (for CIAO 4.16 it seems hard to find
    # this information out, hence the logic below).
    #
    if col.is_varlen():
        # Assume this works even if the first row has no elements
        coltype = col.values[0].dtype
    else:
        coltype = col.values.dtype

    tlmin = col.get_tlmin()
    tlmax = col.get_tlmax()

    # This is a bit-more verbose than normal to appease mypy.
    #
    minval: Any
    maxval: Any
    try:
        ftypeinfo = np.finfo(coltype)
        minval = ftypeinfo.min
        maxval = ftypeinfo.max
    except ValueError:
        dtypeinfo = np.iinfo(coltype)
        minval = dtypeinfo.min
        maxval = dtypeinfo.max

    out.minval = tlmin if tlmin != minval else None
    out.maxval = tlmax if tlmax != maxval else None
    return out


def _read_rmf_matrix(rmfdataset: pycrates.RMFCrateDataset,
                     filename: str,
                     blname: str) -> TableBlock:
    """Read in the given MATRIX block"""

    # By design we know this block exists but we are not guaranteed we
    # can read it in.
    #
    rmf = _get_crate_by_blockname(rmfdataset, blname)
    if rmf is None:
        raise IOErr(f"Unable to read {blname} from {filename}")

    # Ensure we have a DETCHANS keyword
    if rmf.get_key("DETCHANS") is None:
        raise IOErr("nokeyword", filename, "DETCHANS")

    mheaders = _get_meta_data_header(rmf)
    mcols = [copycol(rmf, filename, name, dtype=dtype) for name, dtype in
             [("ENERG_LO", SherpaFloat),
              ("ENERG_HI", SherpaFloat),
              ("N_GRP", None),
              ("F_CHAN", None),
              ("N_CHAN", None),
              ("MATRIX", None)]]
    return TableBlock(blname, header=mheaders, columns=mcols)


def _read_rmf_ebounds(rmfdataset: pycrates.RMFCrateDataset,
                      filename: str,
                      blname: str) -> TableBlock:
    """Read in the given EBOUNDS block"""

    ebounds = _get_crate_by_blockname(rmfdataset, "EBOUNDS")
    if ebounds is None:
        raise IOErr(f"Unable to read {blname} from {filename}")

    eheaders = _get_meta_data_header(ebounds)
    ecols = [copycol(ebounds, filename, name, dtype=dtype) for name, dtype in
             [("CHANNEL", None),
              ("E_MIN", SherpaFloat),
              ("E_MAX", SherpaFloat)]]
    return TableBlock("EBOUNDS", header=eheaders, columns=ecols)


def get_rmf_data(arg: Union[str, pycrates.RMFCrateDataset],
                 make_copy: bool = True
                 ) -> tuple[list[TableBlock], TableBlock, str]:
    """Read a RMF from a file or crate.

    This drops any information in the PRIMARY block.
    """

    if isinstance(arg, str):
        try:
            rmfdataset = pycrates.RMFCrateDataset(arg, mode="r")
        except OSError as oe:
            raise IOErr('openfailed', f"unable to open {arg}: {oe}") from oe

        # Is this check needed?
        if pycrates.is_rmf(rmfdataset) != 1:
            raise IOErr('badfile', arg, "RMFCrateDataset obj")

        filename = arg

    elif pycrates.is_rmf(arg) == 1:
        rmfdataset = arg
        filename = arg.get_filename()
        make_copy = False

    else:
        raise IOErr('badfile', arg, "RMFCrateDataset obj")

    # Find all the potential matrix blocks.
    #
    blnames = _find_matrix_blocks(filename, rmfdataset)
    matrixes = [_read_rmf_matrix(rmfdataset, filename, blname)
                for blname in blnames]
    ebounds = _read_rmf_ebounds(rmfdataset, filename, "EBOUNDS")
    return matrixes, ebounds, filename


def get_pha_data(arg: Union[str, pycrates.PHACrateDataset],
                 make_copy: bool = True,
                 use_background: bool = False
                 ) -> tuple[list[DataType], str]:
    """Read PHA data from a file or crate"""

    if isinstance(arg, str):
        try:
            phadataset = pycrates.PHACrateDataset(arg, mode="r")
        except OSError as oe:
            raise IOErr('openfailed', f"unable to open {arg}: {oe}") from oe

        if pycrates.is_pha(phadataset) != 1:
            raise IOErr('badfile', arg, "PHACrateDataset obj")

        filename = arg

    elif pycrates.is_pha(arg) == 1:
        phadataset = arg
        filename = arg.get_filename()
        make_copy = False

    else:
        raise IOErr('badfile', arg, "PHACrateDataset obj")

    pha = _get_crate_by_blockname(phadataset, "SPECTRUM")

    if pha is None:
        pha = phadataset.get_crate(phadataset.get_current_crate())
        if (_try_key(pha, 'HDUCLAS1') == 'SPECTRUM' or
                _try_key(pha, 'HDUCLAS2') == 'SPECTRUM'):
            pass
        else:
            pha = phadataset.get_crate(1)
            if (_try_key(pha, 'HDUCLAS1') == 'SPECTRUM' or
                    _try_key(pha, 'HDUCLAS2') == 'SPECTRUM'):
                pass
            else:
                # If background maybe better to go on to next block?
                pha = None

    if use_background:

        # Used to read BKGs found in an additional block of
        # Chandra Level 3 PHA files
        for idx in range(phadataset.get_ncrates()):
            block = phadataset.get_crate(idx + 1)
            if _try_key(block, 'HDUCLAS2') == 'BKG':
                pha = block

    if pha is None or pha.get_colnames() is None:
        raise IOErr('filenotfound', arg)

    keys = ['BACKFILE', 'ANCRFILE', 'RESPFILE',
            'BACKSCAL', 'AREASCAL', 'EXPOSURE']

    keys_or_cols = ['BACKSCAL', 'BACKSCUP', 'BACKSCDN', 'AREASCAL']

    datasets = []

    # Calling phadataset.is_pha_type1() is unreliable when
    # both TYPE:I and TYPE:II keywords are in the header.
    # Here, I instead test for a column, SPEC_NUM, that can
    # *only* be present in Type II. SMD 05/15/13
    if _try_col(pha, 'SPEC_NUM') is None:
        data: DataType = {}

        # Keywords
        data['exposure'] = _try_key(pha, 'EXPOSURE', SherpaFloat)
        # data['poisserr'] = _try_key(pha, 'POISSERR', bool)
        data['backfile'] = _try_key(pha, 'BACKFILE')
        data['arffile'] = _try_key(pha, 'ANCRFILE')
        data['rmffile'] = _try_key(pha, 'RESPFILE')

        # Keywords or columns
        for name in keys_or_cols:
            key = name.lower()
            data[key] = _try_key(pha, name, SherpaFloat)
            if data[key] is None:
                data[key] = _try_col(pha, name, make_copy)

        data['header'] = _get_meta_data(pha)
        for key in keys:
            data['header'].pop(key, None)

        # Columns

        data['channel'] = _require_col(pha, 'CHANNEL',
                                       make_copy=make_copy,
                                       fix_type=True)
        # Make sure channel numbers, not indices
        if int(data['channel'][0]) == 0 or pha.get_column('CHANNEL').get_tlmin() == 0:
            data['channel'] = data['channel'] + 1

        data['counts'] = None
        staterror = _try_col(pha, 'STAT_ERR', make_copy)
        if pha.column_exists('COUNTS'):
            data['counts'] = _require_col(pha, 'COUNTS',
                                          make_copy=make_copy,
                                          fix_type=True)
        else:
            data['counts'] = _require_col(pha, 'RATE', label="COUNTS or RATE",
                                          make_copy=make_copy,
                                          fix_type=True)
            data['counts'] *= data['exposure']
            if staterror is not None:
                staterror *= data['exposure']

        data['staterror'] = staterror
        data['syserror'] = _try_col(pha, 'SYS_ERR', make_copy)
        data['background_up'] = _try_col(
            pha, 'BACKGROUND_UP', make_copy, fix_type=True)
        data['background_down'] = _try_col(
            pha, 'BACKGROUND_DOWN', make_copy, fix_type=True)
        data['bin_lo'] = _try_col(pha, 'BIN_LO', make_copy, fix_type=True)
        data['bin_hi'] = _try_col(pha, 'BIN_HI', make_copy, fix_type=True)
        data['grouping'] = _try_col(pha, 'GROUPING', make_copy)
        data['quality'] = _try_col(pha, 'QUALITY', make_copy)

        datasets.append(data)

    else:
        # Type 2 PHA file support
        num = pha.get_nrows()

        # Keywords
        exposure = _try_key(pha, 'EXPOSURE', SherpaFloat)
        # poisserr = _try_key(pha, 'POISSERR', bool)
        backfile = _try_key(pha, 'BACKFILE')
        arffile = _try_key(pha, 'ANCRFILE')
        rmffile = _try_key(pha, 'RESPFILE')

        # Keywords or columns
        backscal = _try_key_list(pha, 'BACKSCAL', num)
        if backscal is None:
            backscal = _try_col_list(pha, 'BACKSCAL', num, make_copy)

        backscup = _try_key_list(pha, 'BACKSCUP', num)
        if backscup is None:
            backscup = _try_col_list(pha, 'BACKSCUP', num, make_copy)

        backscdn = _try_key_list(pha, 'BACKSCDN', num)
        if backscdn is None:
            backscdn = _try_col_list(pha, 'BACKSCDN', num, make_copy)

        areascal = _try_key_list(pha, 'AREASCAL', num)
        if areascal is None:
            areascal = _try_col_list(pha, 'AREASCAL', num, make_copy)

        # Columns

        channel = _require_col_list(pha, 'CHANNEL',
                                    make_copy=make_copy,
                                    fix_type=True)
        # Make sure channel numbers, not indices
        for idx in range(num):
            if int(channel[idx][0]) == 0:
                channel[idx] += 1

        staterror = _try_col_list(pha, 'STAT_ERR', num=num,
                                  make_copy=make_copy)

        counts = None
        if pha.column_exists('COUNTS'):
            counts = _require_col_list(pha, 'COUNTS',
                                       make_copy=make_copy, fix_type=True)

        else:
            if not pha.column_exists('RATE'):
                raise IOErr('reqcol', 'COUNTS or RATE', filename)
            counts = _require_col_list(pha, 'RATE',
                                       make_copy=make_copy, fix_type=True)
            counts *= exposure
            if staterror is not None:
                staterror *= exposure

        syserror = _try_col_list(pha, 'SYS_ERR', num, make_copy)
        background_up = _try_col_list(
            pha, 'BACKGROUND_UP', num, make_copy, fix_type=True)
        background_down = _try_col_list(
            pha, 'BACKGROUND_DOWN', num, make_copy, fix_type=True)
        bin_lo = _try_col_list(pha, 'BIN_LO', num, make_copy, fix_type=True)
        bin_hi = _try_col_list(pha, 'BIN_HI', num, make_copy, fix_type=True)
        grouping = _try_col_list(pha, 'GROUPING', num, make_copy)
        quality = _try_col_list(pha, 'QUALITY', num, make_copy)

        orders = _try_key_list(pha, 'TG_M', num)
        if orders is None:
            orders = _try_col_list(pha, 'TG_M', num, make_copy)

        parts = _try_key_list(pha, 'TG_PART', num)
        if parts is None:
            parts = _try_col_list(pha, 'TG_PART', num, make_copy)

        specnums = _try_col_list(pha, 'SPEC_NUM', num, make_copy)
        srcids = _try_col_list(pha, 'TG_SRCID', num, make_copy)

        # Iterate over all rows of channels, counts, errors, etc
        # Populate a list of dictionaries containing individual dataset info
        for (bscal, bscup, bscdn, arsc, chan, cnt, staterr, syserr,
             backup, backdown, binlo, binhi, grp, qual, ordr, prt,
             specnum, srcid
             ) in zip(backscal, backscup, backscdn, areascal, channel,
                      counts, staterror, syserror, background_up,
                      background_down, bin_lo, bin_hi, grouping, quality,
                      orders, parts, specnums, srcids):

            idata: DataType = {}

            idata['exposure'] = exposure
            # data['poisserr'] = poisserr
            idata['backfile'] = backfile
            idata['arffile'] = arffile
            idata['rmffile'] = rmffile

            idata['backscal'] = bscal
            idata['backscup'] = bscup
            idata['backscdn'] = bscdn
            idata['areascal'] = arsc

            idata['channel'] = chan
            idata['counts'] = cnt
            idata['staterror'] = staterr
            idata['syserror'] = syserr
            idata['background_up'] = backup
            idata['background_down'] = backdown
            idata['bin_lo'] = binlo
            idata['bin_hi'] = binhi
            idata['grouping'] = grp
            idata['quality'] = qual
            idata['header'] = _get_meta_data(pha)
            idata['header']['TG_M'] = ordr
            idata['header']['TG_PART'] = prt
            idata['header']['SPEC_NUM'] = specnum
            idata['header']['TG_SRCID'] = srcid

            for key in keys:
                idata['header'].pop(key, None)

            datasets.append(idata)

    return datasets, filename


#
# Write/Pack Functions #
#

def write_dataset(dataset: Union[TABLECrate, IMAGECrate, CrateDataset],
                  filename: str,
                  *,
                  ascii: bool,
                  clobber: bool) -> None:
    """Write out the data."""

    if not clobber and os.path.isfile(filename):
        raise IOErr("filefound", filename)

    if ascii and '[' not in filename and ']' not in filename:
        filename += "[opt kernel=text/simple]"

        # For CIAO 4.15-era crates, if this is a CrateDataset and the
        # first block is a "empty" image then the code will fail when
        # writing out as an ASCII file. So this will remove such a
        # block.
        #
        try:
            cr = dataset.get_crate(1)
        except AttributeError:
            cr = None

        if isinstance(cr, IMAGECrate) and cr.get_image().values.size == 0:
            dataset.delete_crate(1)

    dataset.write(filename, clobber=True)


def pack_image_data(data: DataTypeArg,
                    header: HdrTypeArg) -> IMAGECrate:
    """Pack up the image data."""

    img = IMAGECrate()

    # Write Image Header Keys
    _update_header(img, header, skip_if_known=False)

    # Write Image WCS Header Keys
    if data['eqpos'] is not None:
        cdeltw = data['eqpos'].cdelt
        crvalw = data['eqpos'].crval
        crpixw = data['eqpos'].crpix
        equin = data['eqpos'].equinox

    if data['sky'] is not None:
        cdeltp = data['sky'].cdelt
        crvalp = data['sky'].crval
        crpixp = data['sky'].crpix

        _set_key(img, 'MTYPE1', 'sky     ')
        _set_key(img, 'MFORM1', 'x,y     ')
        _set_key(img, 'CTYPE1P', 'x       ')
        _set_key(img, 'CTYPE2P', 'y       ')
        _set_key(img, 'WCSNAMEP', 'PHYSICAL')
        _set_key(img, 'CDELT1P', cdeltp[0])
        _set_key(img, 'CDELT2P', cdeltp[1])
        _set_key(img, 'CRPIX1P', crpixp[0])
        _set_key(img, 'CRPIX2P', crpixp[1])
        _set_key(img, 'CRVAL1P', crvalp[0])
        _set_key(img, 'CRVAL2P', crvalp[1])

        if data['eqpos'] is not None:
            # Simply the inverse of read transformations in get_image_data
            cdeltw = cdeltw * cdeltp
            crpixw = (crpixw - crvalp) / cdeltp + crpixp

    if data['eqpos'] is not None:
        _set_key(img, 'MTYPE2', 'EQPOS   ')
        _set_key(img, 'MFORM2', 'RA,DEC  ')
        _set_key(img, 'CTYPE1', 'RA---TAN')
        _set_key(img, 'CTYPE2', 'DEC--TAN')
        _set_key(img, 'CDELT1', cdeltw[0])
        _set_key(img, 'CDELT2', cdeltw[1])
        _set_key(img, 'CRPIX1', crpixw[0])
        _set_key(img, 'CRPIX2', crpixw[1])
        _set_key(img, 'CRVAL1', crvalw[0])
        _set_key(img, 'CRVAL2', crvalw[1])
        _set_key(img, 'EQUINOX', equin)

    # Write Image pixel values
    pix_col = pycrates.CrateData()
    pix_col.values = data['pixels']
    img.add_image(pix_col)
    return img


def set_image_data(filename: str,
                   data: DataTypeArg,
                   header: HdrTypeArg,
                   ascii: bool = False,
                   clobber: bool = False) -> None:
    """Write out the image data."""

    if ascii and '[' not in filename and ']' not in filename:
        raise IOErr('writenoimg')

    img = pack_image_data(data, header)
    write_dataset(img, filename, ascii=ascii, clobber=clobber)


def pack_table_data(data: ColumnsType,
                    col_names: NamesType,
                    header: Optional[HdrTypeArg] = None) -> TABLECrate:
    """Pack up the table data."""

    tbl = TABLECrate()
    hdr = {} if header is None else header
    for name in col_names:
        _set_column(tbl, name, data[name])

    _update_header(tbl, hdr, skip_if_known=False)
    return tbl


def set_table_data(filename: str,
                   data: ColumnsType,
                   col_names: NamesType,
                   header: Optional[HdrTypeArg] = None,
                   ascii: bool = False,
                   clobber: bool = False) -> None:
    """Write out the table data."""

    tbl = pack_table_data(data, col_names, header=header)
    write_dataset(tbl, filename, ascii=ascii, clobber=clobber)


def pack_arf_data(blocks: BlockList) -> CrateDataset:
    """Pack the ARF"""

    return pack_hdus(blocks)


def set_arf_data(filename: str,
                 blocks: BlockList,
                 ascii: bool = False,
                 clobber: bool = False) -> None:
    """Write out the ARF"""

    arf = pack_arf_data(blocks)
    write_dataset(arf, filename, ascii=ascii, clobber=clobber)


def pack_pha_data(data: ColumnsType,
                  col_names: NamesType,
                  header: Optional[HdrTypeArg] = None
                  ) -> pycrates.PHACrateDataset:
    """Pack the PHA data."""

    if header is None:
        raise ArgumentTypeErr("badarg", "header", "set")

    phadataset = pycrates.PHACrateDataset()

    # FIXME: Placeholder for pycrates2 bug
    phadataset.set_rw_mode('rw')

    pha = TABLECrate()
    pha.name = "SPECTRUM"

    _update_header(pha, header, skip_if_known=False)

    # Write column values using CrateData objects
    for name in col_names:
        if data[name] is None:
            continue
        _set_column(pha, name, data[name])

    phadataset.add_crate(pha)
    return phadataset


def set_pha_data(filename: str,
                 data: ColumnsType,
                 col_names: NamesType,
                 header: Optional[HdrTypeArg] = None,
                 ascii: bool = False,
                 clobber: bool = False) -> None:
    """Create a PHA dataset/file"""

    pha = pack_pha_data(data, col_names, header)
    write_dataset(pha, filename, ascii=ascii, clobber=clobber)


def _update_header(cr: CrateType,
                   header: HdrTypeArg,
                   skip_if_known: bool = True) -> None:
    """Update the header of the crate."""

    for key, value in header.items():
        if skip_if_known and cr.key_exists(key):
            continue

        # value should not be None, but just in case
        if value is None:
            continue

        # Do we need to worry about keys like TTYPE1 which crates
        # hides from us?
        _set_key(cr, key, value)


def pack_rmf_data(blocks: BlockList) -> pycrates.RMFCrateDataset:
    """Pack up the RMF data."""

    _validate_block_names(blocks.blocks)

    rmf = pycrates.RMFCrateDataset()
    if blocks.header is not None:
        rmf.add_crate(_create_primary_crate(blocks.header))

    # Should only be tables in the list but allow images to make the
    # code simpler.
    #
    for block in blocks.blocks:
        if isinstance(block, TableBlock):
            cr = _create_table_crate(block)
        elif isinstance(block, ImageBlock):
            cr = _create_image_crate(block)
        else:
            raise RuntimeError(f"Unsupported block: {block}")

        rmf.add_crate(cr)

    return rmf


def set_rmf_data(filename: str,
                 blocks: BlockList,
                 clobber: bool = False) -> None:
    """Save the RMF data to disk.

    Unlike the other save_*_data calls this does not support the ascii
    argument.

    """

    rmf = pack_rmf_data(blocks)
    write_dataset(rmf, filename, ascii=False, clobber=clobber)


def set_arrays(filename: str,
               args: Sequence[np.ndarray],
               fields: Optional[NamesType] = None,
               ascii: bool = True,
               clobber: bool = False) -> None:
    """Write out the columns."""

    # Historically the clobber command has been checked before
    # processing the data, so do so here.
    #
    if not clobber and os.path.isfile(filename):
        raise IOErr("filefound", filename)

    # Check args is a sequence of sequences (although not a complete
    # check).
    #
    try:
        size = len(args[0])
    except (TypeError, IndexError) as exc:
        raise IOErr('noarrayswrite') from exc

    for arg in args[1:]:
        try:
            argsize = len(arg)
        except (TypeError, IndexError) as exc:
            raise IOErr('noarrayswrite') from exc

        if argsize != size:
            raise IOErr('arraysnoteq')

    nargs = len(args)
    if fields is None:
        fieldnames = [f'col{idx + 1}' for idx in range(nargs)]
    elif nargs == len(fields):
        fieldnames = list(fields)
    else:
        raise IOErr('wrongnumcols', nargs, len(fields))

    tbl = TABLECrate()
    for val, name in zip(args, fieldnames):
        _set_column(tbl, name, val)

    write_dataset(tbl, filename, ascii=ascii, clobber=clobber)


def _add_header(cr: CrateType,
                header: Header) -> None:
    """Add the header keywords to the crate."""

    for item in header.values:
        key = pycrates.CrateKey()
        key.name = item.name
        key.value = item.value
        key.unit = item.unit
        key.desc = item.desc
        cr.add_key(key)


def _create_primary_crate(hdr: Header) -> IMAGECrate:
    """Create the primary block."""

    out = IMAGECrate()
    out.name = "PRIMARY"

    # For some reason we need to add an empty image
    # (CIAO 4.15).
    #
    null = pycrates.CrateData()
    null.values = np.asarray([], dtype=np.uint8)
    out.add_image(null)

    _add_header(out, hdr)
    return out


def _create_image_crate(hdu: ImageBlock) -> IMAGECrate:
    """Create an image block."""

    out = IMAGECrate()
    out.name = hdu.name

    img = pycrates.CrateData()
    img.values = hdu.image
    out.add_image(img)

    _add_header(out, hdu.header)
    return out


def _create_table_crate(hdu: TableBlock) -> TABLECrate:
    """Create a table block."""

    out = TABLECrate()
    out.name = hdu.name

    _add_header(out, hdu.header)

    for idx, col in enumerate(hdu.columns, 1):
        cdata = pycrates.CrateData()
        cdata.name = col.name
        cdata.desc = col.desc
        cdata.unit = col.unit
        cdata.values = col.values
        out.add_column(cdata)

        # For CIAO 4.16 crates does not have a "nice" way to set
        # TLMIN/MAX, but it appears we can just set the TLMIN/MAX
        # values directly in the header and it works.
        #
        if col.minval is not None:
            key = pycrates.CrateKey()
            key.name = f"TLMIN{idx}"
            key.value = col.minval
            out.add_key(key)

        if col.maxval is not None:
            key = pycrates.CrateKey()
            key.name = f"TLMAX{idx}"
            key.value = col.maxval
            out.add_key(key)

    return out


def _validate_block_names(hdulist: Sequence[BlockType]) -> None:
    """Ensure the block names are "independent".

    Crates will correctly handle the case when blocks are numbered
    MATRIX1, MATRIX2 (setting EXTNAME to MATRIX and EXTVER to 1 or 2),
    but let's automatically add this if the user hasn't already done
    so (since AstroPy works differently).

    This may change the input list contents.
    """

    blnames: dict[str, int] = defaultdict(int)
    for hdu in hdulist:
        blnames[hdu.name.upper()] += 1

    multi = [n for n,v in blnames.items() if v > 1]

    # If the names are unique do nothing
    #
    if len(multi) == 0:
        return

    extvers: dict[str, int] = defaultdict(int)
    for hdu in hdulist:
        if hdu.name not in multi:
            continue

        # Either create or update the EXTVER value
        extvers[hdu.name] += 1
        extver = extvers[hdu.name]

        # Change the name field
        hdu.name = f"{hdu.name}{extver}"


def pack_hdus(blocks: BlockList) -> CrateDataset:
    """Create a dataset."""

    _validate_block_names(blocks.blocks)

    dset = CrateDataset()
    if blocks.header is not None:
        dset.add_crate(_create_primary_crate(blocks.header))

    for block in blocks.blocks:
        if isinstance(block, TableBlock):
            cr = _create_table_crate(block)
        elif isinstance(block, ImageBlock):
            cr = _create_image_crate(block)
        else:
            raise RuntimeError(f"Unsupported block: {block}")

        dset.add_crate(cr)

    return dset


def set_hdus(filename: str,
             blocks: BlockList,
             clobber: bool = False) -> None:
    """Write out a dataset."""

    dset = pack_hdus(blocks)
    write_dataset(dset, filename, ascii=False, clobber=clobber)
