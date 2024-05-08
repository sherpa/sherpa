#
#  Copyright (C) 2011, 2015 - 2024
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

"""
Read and write FITS [1]_ files using the ``astropy.io.fits`` module [2]_.

Notes
-----

Support for PyFITS [3]_ was dropped in Sherpa 4.10.1.

References
----------

.. [1] https://en.wikipedia.org/wiki/FITS

.. [2] http://astropy.readthedocs.org/en/latest/io/fits/

.. [3] http://www.stsci.edu/institute/software_hardware/pyfits

"""

import logging
import os
from typing import Any, Mapping, Optional, Sequence, Union
import warnings

import numpy as np

from astropy.io import fits  # type: ignore
from astropy.io.fits.column import _VLF  # type: ignore
from astropy.table import Table  # type: ignore


import sherpa.utils
from sherpa.utils.err import ArgumentTypeErr, IOErr
from sherpa.utils.numeric_types import SherpaInt, SherpaUInt, \
    SherpaFloat
from sherpa.io import get_ascii_data, write_arrays

from .types import ColumnsType, DataType, DataTypeArg, \
    HdrType, HdrTypeArg, KeyType, NamesType
from .xstable import HeaderItem, TableHDU

warning = logging.getLogger(__name__).warning
error = logging.getLogger(__name__).error

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
           'pack_table_data', 'pack_arf_data', 'pack_rmf_data',
           'pack_image_data', 'pack_hdus',
           'set_table_data', 'set_image_data', 'set_pha_data',
           'set_arf_data', 'set_rmf_data', 'set_hdus')


DatasetType = Union[str, fits.HDUList]
HDUType = Union[fits.PrimaryHDU, fits.BinTableHDU]


def _try_key(hdu: HDUType,
             name: str,
             *,
             fix_type: bool = False,
             dtype: type = SherpaFloat) -> Optional[KeyType]:

    value = hdu.header.get(name, None)
    if value is None:
        return None

    # TODO: As noted with the crates backend, this test is probably
    # overly broad, and we should look at simplifying it.
    #
    if str(value).find('none') != -1:
        return None

    if not fix_type:
        return value

    return dtype(value)


def _require_key(hdu: HDUType,
                 name: str,
                 *,
                 fix_type: bool = False,
                 dtype: type = SherpaFloat) -> KeyType:
    value = _try_key(hdu, name, fix_type=fix_type, dtype=dtype)
    if value is None:
        raise IOErr('nokeyword',
                    hdu._file.name if hdu._file is not None else "unknown",
                    name)
    return value


def _get_meta_data(hdu: HDUType) -> HdrType:
    # If the header keywords are not specified correctly then
    # astropy will error out when we try to access it. Since
    # this is not an uncommon problem, there is a verify method
    # that can be used to fix up the data to avoid this: the
    # "silentfix" option is used so as not to complain too much.
    #
    hdu.verify('silentfix')

    meta = {}
    for key, val in hdu.header.items():
        meta[key] = val

    return meta


def _add_keyword(hdrlist: fits.Header,
                 name: str,
                 val: Any) -> None:
    """Add the name,val pair to hdulist.

    Ensures that the name is upper-case.
    """

    name = name.upper()
    hdrlist.append(fits.Card(name, val))


def _try_col(hdu, name, *, dtype=SherpaFloat, fix_type=False):
    try:
        col = hdu.data.field(name)
    except KeyError:
        return None

    if isinstance(col, _VLF):
        col = np.concatenate([np.asarray(row) for row in col])
    else:
        col = np.asarray(col).ravel()

    if fix_type:
        col = col.astype(dtype)

    return col


def _try_tbl_col(hdu, name, *, dtype=SherpaFloat, fix_type=False):
    try:
        col = hdu.data.field(name)
    except KeyError:
        return (None, )

    if isinstance(col, _VLF):
        col = np.concatenate([np.asarray(row) for row in col])
    else:
        col = np.asarray(col)

    if fix_type:
        col = col.astype(dtype)

    return np.column_stack(col)


def _try_vec(hdu, name, *, size=2, dtype=SherpaFloat, fix_type=False):
    try:
        col = hdu.data.field(name)
    except KeyError:
        return np.full(size, None)

    if isinstance(col, _VLF):
        col = np.concatenate([np.asarray(row) for row in col])
    else:
        col = np.asarray(col)

    if fix_type:
        return col.astype(dtype)

    return col


def _require_col(hdu, name, *, dtype=SherpaFloat, fix_type=False):
    col = _try_col(hdu, name, dtype=dtype, fix_type=fix_type)
    if col is None:
        raise IOErr('reqcol', name, hdu._file.name)
    return col


def _require_tbl_col(hdu, name, *, dtype=SherpaFloat, fix_type=False):
    col = _try_tbl_col(hdu, name, dtype=dtype, fix_type=fix_type)
    if len(col) > 0 and col[0] is None:
        raise IOErr('reqcol', name, hdu._file.name)
    return col


def _require_vec(hdu, name, *, size=2, dtype=SherpaFloat, fix_type=False):
    col = _try_vec(hdu, name, size=size, dtype=dtype, fix_type=fix_type)
    if np.equal(col, None).any():
        raise IOErr('reqcol', name, hdu._file.name)
    return col


def _try_col_or_key(hdu, name, *, dtype=SherpaFloat, fix_type=False):
    col = _try_col(hdu, name, dtype=dtype, fix_type=fix_type)
    if col is not None:
        return col
    return _try_key(hdu, name, fix_type=fix_type, dtype=dtype)


def _try_vec_or_key(hdu, name, size, *, dtype=SherpaFloat, fix_type=False):
    col = _try_col(hdu, name, dtype=dtype, fix_type=fix_type)
    if col is not None:
        return col

    got = _try_key(hdu, name, fix_type=fix_type, dtype=dtype)
    return np.full(size, got)


# Read Functions

# Crates supports reading gzipped files both when the '.gz' suffix
# is given and even when it is not (i.e. if you ask to load 'foo.fits'
# and 'foo.fits.gz' exists, then it will use this). This requires some
# effort to emulate with the astropy backend.

def _infer_and_check_filename(filename: str) -> str:
    if os.path.exists(filename):
        return filename

    gzname = f"{filename}.gz"
    if os.path.exists(gzname):
        return gzname

    # error message reports the original filename requested
    raise IOErr('filenotfound', filename)


def is_binary_file(filename: str) -> bool:
    """Do we think the file contains binary data?

    Notes
    -----
    This is a wrapper around the version in sherpa.utils which
    attempts to transparently support reading from a `.gz`
    version if the input file name does not exist.
    """
    fname = _infer_and_check_filename(filename)
    return sherpa.utils.is_binary_file(fname)


def open_fits(filename: str) -> fits.HDUList:
    """Try and open filename as a FITS file.

    Parameters
    ----------
    filename : str
       The name of the FITS file to open. If the file
       can not be opened then '.gz' is appended to the
       name and the attempt is tried again.

    Returns
    -------
    out
       The return value from the `astropy.io.fits.open`
       function.
    """
    fname = _infer_and_check_filename(filename)

    # With AstroPy v4.2.1 we start to see warnings when this call
    # fails (i.e. when the argument is not a FITS file).  This leads
    # to spurious messages to the user, so we just remove all
    # warnings. This may need tweaking at some point: ideally we would
    # only want to hide the warnings when it is not a FITS file but
    # show them when it is a FITS file.
    #
    # Note that this is not thread safe.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', module='astropy.io.fits')
        out = fits.open(fname)

    return out


def read_table_blocks(arg: DatasetType,
                      make_copy: bool = False
                      ) -> tuple[str,
                                 dict[int, ColumnsType],
                                 dict[int, HdrType]]:
    """Read in tabular data with no restrictions on the columns."""

    # This sets nobinary=True to match the original version of
    # the code, which did not use _get_file_contents.
    #
    hdus, filename, close = _get_file_contents(arg,
                                               exptype="BinTableHDU",
                                               nobinary=True)

    cols: dict[int, ColumnsType] = {}
    hdr: dict[int, HdrType] = {}
    try:
        for blockidx, hdu in enumerate(hdus, 1):
            hdr[blockidx] = {}
            header = hdu.header
            if header is not None:
                for key in header.keys():
                    hdr[blockidx][key] = header[key]

            # skip over primary, hdu.data is None

            cols[blockidx] = {}
            recarray = hdu.data
            if recarray is not None:
                for colname in recarray.names:
                    cols[blockidx][colname] = recarray[colname]

    finally:
        if close:
            hdus.close()

    return filename, cols, hdr


def _get_file_contents(arg: DatasetType,
                       exptype: str = "PrimaryHDU",
                       nobinary: bool = False
                       ) -> tuple[fits.HDUList, str, bool]:
    """Read in the contents if needed.

    Set nobinary to True to avoid checking that the input
    file is a binary file (via the is_binary_file routine).
    Is this needed?
    """

    if isinstance(arg, str) and (not nobinary or is_binary_file(arg)):
        tbl = open_fits(arg)
        filename = arg
        close = True
    elif isinstance(arg, fits.HDUList) and len(arg) > 0 and \
            isinstance(arg[0], fits.PrimaryHDU):
        tbl = arg
        filename = arg._file.name if arg._file is not None else "unknown"
        close = False
    else:
        msg = f"a binary FITS table or a {exptype} list"
        raise IOErr('badfile', arg, msg)

    return (tbl, filename, close)


def _find_binary_table(tbl: fits.HDUList,
                       filename: str,
                       blockname: Optional[str] = None) -> fits.BinTableHDU:
    """Return the first binary table extension we find. If blockname
    is not None then the name of the block has to match (case-insensitive
    match), and any spaces are removed from blockname before checking.

    Throws an exception if there isn't a matching HDU.
    """

    # This behaviour appears odd:
    # - if blockname is not set then find the first table block
    #   (this seems okay)
    # - if blockname is set then loop through the tables and
    #   find the first block which is either a table or matches
    #   blockname
    #   (why the "is table" check; why not only match the blockname?)
    #
    if blockname is None:
        for hdu in tbl:
            if isinstance(hdu, fits.BinTableHDU):
                return hdu

    else:
        blockname = blockname.strip().lower()
        for hdu in tbl:
            if hdu.name.lower() == blockname or \
                    isinstance(hdu, fits.BinTableHDU):
                return hdu

    raise IOErr('badext', filename)


def get_header_data(arg: DatasetType,
                    blockname: Optional[str] = None,
                    hdrkeys: Optional[NamesType] = None
                    ) -> HdrType:
    """Read in the header data."""

    tbl, filename, close = _get_file_contents(arg, exptype="BinTableHDU")

    hdr = {}
    try:
        hdu = _find_binary_table(tbl, filename, blockname)

        if hdrkeys is None:
            hdrkeys = hdu.header.keys()

        for key in hdrkeys:
            # TODO: should this set require_type=true or remove dtype,
            # as currently it does not change the value to str.
            hdr[key] = _require_key(hdu, key, dtype=str)

    finally:
        if close:
            tbl.close()

    return hdr


def get_column_data(*args) -> list[np.ndarray]:
    """Extract the column data."""

    # args is passed as type list
    if len(args) == 0:
        raise IOErr('noarrays')

    cols = []
    for arg in args:
        if arg is not None and not isinstance(arg, (np.ndarray, list, tuple)):
            raise IOErr('badarray', arg)
        if arg is not None:
            vals = np.asanyarray(arg)
            for col in np.atleast_2d(vals.T):
                cols.append(col)
        else:
            cols.append(arg)

    return cols


def get_table_data(arg: DatasetType,
                   ncols: int = 1,
                   colkeys: Optional[NamesType] = None,
                   make_copy: bool = False,
                   fix_type: bool = False,
                   blockname: Optional[str] = None,
                   hdrkeys: Optional[NamesType] = None
                   ) -> tuple[list[str], list[np.ndarray], str, HdrType]:
    """Read columns."""

    tbl, filename, close = _get_file_contents(arg, exptype="BinTableHDU")

    try:
        hdu = _find_binary_table(tbl, filename, blockname)
        cnames = list(hdu.columns.names)

        # Try Channel, Counts or X,Y before defaulting to the first
        # ncols columns in cnames (when colkeys is not given).
        #
        if colkeys is not None:
            colkeys = [name.strip().upper() for name in list(colkeys)]
        elif 'CHANNEL' in cnames and 'COUNTS' in cnames:
            colkeys = ['CHANNEL', 'COUNTS']
        elif 'X' in cnames and 'Y' in cnames:
            colkeys = ['X', 'Y']
        else:
            colkeys = cnames[:ncols]

        cols = []
        for name in colkeys:
            for col in _require_tbl_col(hdu, name, fix_type=fix_type):
                cols.append(col)

        hdr = {}
        if hdrkeys is not None:
            for key in hdrkeys:
                hdr[key] = _require_key(hdu, key)

    finally:
        if close:
            tbl.close()

    return colkeys, cols, filename, hdr


def get_image_data(arg: DatasetType,
                   make_copy: bool = False,
                   fix_type: bool = True  # Unused
                   ) -> tuple[DataType, str]:
    """Read image data."""

    hdus, filename, close = _get_file_contents(arg)

    #   FITS uses logical-to-world where we use physical-to-world.
    #   For all transforms, update their physical-to-world
    #   values from their logical-to-world values.
    #   Find the matching physical transform
    #      (same axis no, but sub = 'P' )
    #   and use it for the update.
    #   Physical tfms themselves do not get updated.
    #
    #  Fill the physical-to-world transform given the
    #  logical-to-world and the associated logical-to-physical.
    #      W = wv + wd * ( P - wp )
    #      P = pv + pd * ( L - pp )
    #      W = lv + ld * ( L - lp )
    # Then
    #      L = pp + ( P - pv ) / pd
    # so   W = lv + ld * ( pp + (P-pv)/pd - lp )
    #        = lv + ( ld / pd ) * ( P - [ pv +  (lp-pp)*pd ] )
    # Hence
    #      wv = lv
    #      wd = ld / pd
    #      wp = pv + ( lp - pp ) * pd

    #  EG suppose phys-to-world is
    #         W =  1000 + 2.0 * ( P - 4.0 )
    #  and we bin and scale to generate a logical-to-phys of
    #         P =  20 + 4.0 * ( L - 10 )
    #  Then
    #         W = 1000 + 2.0 * ( (20-4) - 4 * 10 ) + 2 * 4 $
    #

    try:
        data: DataType = {}

        # Look for data in the primary or first block.
        #
        img = hdus[0]
        if img.data is None:
            img = hdus[1]
            if img.data is None:
                raise IOErr('badimg', filename)

        if not isinstance(img, (fits.PrimaryHDU, fits.ImageHDU)):
            raise IOErr("badimg", filename)

        data['y'] = np.asarray(img.data)

        if HAS_TRANSFORM:

            def _get_wcs_key(key, suffix=""):
                val1 = _try_key(img, f"{key}1{suffix}")
                if val1 is None:
                    return None

                val2 = _try_key(img, f"{key}2{suffix}")
                if val2 is None:
                    return None

                return np.array([val1, val2], dtype=SherpaFloat)

            cdeltp = _get_wcs_key('CDELT', 'P')
            crpixp = _get_wcs_key('CRPIX', 'P')
            crvalp = _get_wcs_key('CRVAL', 'P')
            cdeltw = _get_wcs_key('CDELT')
            crpixw = _get_wcs_key('CRPIX')
            crvalw = _get_wcs_key('CRVAL')

            # proper calculation of cdelt wrt PHYSICAL coords
            if cdeltw is not None and cdeltp is not None:
                cdeltw = cdeltw / cdeltp

            # proper calculation of crpix wrt PHYSICAL coords
            if (crpixw is not None and crvalp is not None and
                cdeltp is not None and crpixp is not None):
                crpixw = crvalp + (crpixw - crpixp) * cdeltp

            if (cdeltp is not None and crpixp is not None and
                crvalp is not None):
                data['sky'] = WCS('physical', 'LINEAR', crvalp,
                                  crpixp, cdeltp)

            if (cdeltw is not None and crpixw is not None and
                crvalw is not None):
                data['eqpos'] = WCS('world', 'WCS', crvalw, crpixw,
                                    cdeltw)

        data['header'] = _get_meta_data(img)
        for key in ['CTYPE1P', 'CTYPE2P', 'WCSNAMEP', 'CDELT1P',
                    'CDELT2P', 'CRPIX1P', 'CRPIX2P', 'CRVAL1P', 'CRVAL2P',
                    'EQUINOX']:
            data['header'].pop(key, None)

    finally:
        if close:
            hdus.close()

    return data, filename


def _is_ogip_block(hdu: HDUType,
                   bltype1: str,
                   bltype2: Optional[str] = None) -> bool:
    """Does the block contain the expected HDUCLAS1 or HDUCLAS2 values?

    If given, we need both HDUCLAS1 and 2 to be set correctly.
    """

    if _try_key(hdu, 'HDUCLAS1') != bltype1:
        return False

    if bltype2 is None:
        return True

    return  _try_key(hdu, 'HDUCLAS2') == bltype2


def _has_ogip_type(hdus: fits.HDUList,
                   bltype: str,
                   bltype2: Optional[str] = None) -> Optional[HDUType]:
    """Return True if hdus[1] exists and has
    the given type (as determined by the HDUCLAS1 or HDUCLAS2
    keywords). If bltype2 is None then bltype is used for
    both checks, otherwise bltype2 is used for HDUCLAS2 and
    bltype is for HDUCLAS1.
    """

    try:
        hdu = hdus[1]
    except IndexError:
        return None

    if _is_ogip_block(hdu, bltype, bltype2):
        return hdu

    return None


def get_arf_data(arg: DatasetType,
                 make_copy: bool = False
                 ) -> tuple[DataType, str]:
    """Read in the ARF."""

    arf, filename, close = _get_file_contents(arg,
                                              exptype="BinTableHDU",
                                              nobinary=True)

    try:
        # Should we do the _has_ogip_type check first and only if that
        # fails fall back to the "common" block names? Given the lack
        # of consistency in following the OGIP standards by missions
        # let's keep this ordering.
        #
        for blname in ["SPECRESP", "AXAF_ARF"]:
            if blname in arf:
                hdu = arf[blname]
                break
        else:
            hdu = _has_ogip_type(arf, 'RESPONSE', 'SPECRESP')
            if hdu is None:
                raise IOErr('notrsp', filename, 'an ARF')

        data: DataType = {}

        data['exposure'] = _try_key(hdu, 'EXPOSURE', fix_type=True)

        data['energ_lo'] = _require_col(hdu, 'ENERG_LO', fix_type=True)
        data['energ_hi'] = _require_col(hdu, 'ENERG_HI', fix_type=True)
        data['specresp'] = _require_col(hdu, 'SPECRESP', fix_type=True)
        data['bin_lo'] = _try_col(hdu, 'BIN_LO', fix_type=True)
        data['bin_hi'] = _try_col(hdu, 'BIN_HI', fix_type=True)
        data['header'] = _get_meta_data(hdu)
        data['header'].pop('EXPOSURE', None)

    finally:
        if close:
            arf.close()

    return data, filename


def _read_col(hdu: fits.BinTableHDU, name: str) -> np.ndarray:
    """A specialized form of _require_col

    There is no attempt to convert from a variable-length field
    to any other form.

    """

    try:
        return hdu.data[name]
    except KeyError:
        raise IOErr("reqcol", name, hdu._file.name) from None


# Commonly-used block names for the MATRIX block. Only the first two
# are given in the OGIP standard.
#
RMF_BLOCK_NAMES = ["MATRIX", "SPECRESP MATRIX", "AXAF_RMF", "RSP_MATRIX"]


def _find_matrix_blocks(filename: str,
                        hdus: fits.HDUList) -> list[fits.BinTableHDU]:
    """Report the block names that contain MATRIX data in a RMF.

    Parameters
    ----------
    filename : str
    hdus : fits.HDUList

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
    blocks = []
    for hdu in hdus:
        if hdu.name in RMF_BLOCK_NAMES or \
           _is_ogip_block(hdu, "RESPONSE", "RSP_MATRIX"):
            blocks.append(hdu)

    if not blocks:
        raise IOErr('notrsp', filename, 'an RMF')

    return blocks


def _read_rmf_data(arg: DatasetType
                   ) -> tuple[DataType, str]:
    """Read in the data from the RMF."""

    rmf, filename, close = _get_file_contents(arg,
                                              exptype="BinTableHDU",
                                              nobinary=True)

    try:

        # Find all the potential matrix blocks.
        #
        mblocks = _find_matrix_blocks(filename, rmf)
        nmat = len(mblocks)
        if nmat > 1:
            # Warn the user that the multi-matrix RMF is not supported.
            #
            error("RMF in %s contains %d MATRIX blocks; "
                  "Sherpa only uses the first block!",
                  filename, nmat)

        hdu = mblocks[0]

        # The comments note the type we want the column to be, but
        # this conversion needs to be done after cleaning up the
        # data.
        #
        data: DataType = {}
        data['detchans'] = SherpaUInt(_require_key(hdu, 'DETCHANS'))
        data['energ_lo'] = _read_col(hdu, 'ENERG_LO')  # SherpaFloat
        data['energ_hi'] = _read_col(hdu, 'ENERG_HI')  # SherpaFloat
        data['n_grp'] = _read_col(hdu, 'N_GRP')        # SherpaUInt
        data['f_chan'] = _read_col(hdu, 'F_CHAN')      # SherpaUInt
        data['n_chan'] = _read_col(hdu, 'N_CHAN')      # SherpaUInt
        data['matrix'] = _read_col(hdu, "MATRIX")

        data['header'] = _get_meta_data(hdu)
        data['header'].pop('DETCHANS', None)

        # Beginning of non-Chandra RMF support
        fchan_col = list(hdu.columns.names).index('F_CHAN') + 1
        tlmin = _try_key(hdu, f"TLMIN{fchan_col}", fix_type=True,
                         dtype=SherpaUInt)

        if tlmin is not None:
            data['offset'] = tlmin
        else:
            error("Failed to locate TLMIN keyword for F_CHAN "
                  "column in RMF file '%s'; "
                  'Update the offset value in the RMF data set to '
                  'appropriate TLMIN value prior to fitting', filename)

        try:
            hdu = rmf['EBOUNDS']
        except KeyError:
            hdu = None

        if hdu is not None:
            data['e_min'] = _try_col(hdu, 'E_MIN', fix_type=True)
            data['e_max'] = _try_col(hdu, 'E_MAX', fix_type=True)

            # Beginning of non-Chandra RMF support
            chan_col = list(hdu.columns.names).index('CHANNEL') + 1
            tlmin = _try_key(hdu, f"TLMIN{chan_col}", fix_type=True,
                             dtype=SherpaUInt)
            if tlmin is not None:
                data['offset'] = tlmin

        else:
            data['e_min'] = None
            data['e_max'] = None
    finally:
        if close:
            rmf.close()

    return data, filename


def get_rmf_data(arg: DatasetType,
                 make_copy: bool = False
                 ) -> tuple[DataType, str]:
    """Read in the RMF.

    Notes
    -----
    For more information on the RMF format see the `OGIP Calibration
    Memo CAL/GEN/92-002
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html>`_.

    """

    # Read in the columns from the MATRIX and EBOUNDS extensions.
    #
    data, filename = _read_rmf_data(arg)

    # This could be taken from the NUMELT keyword, but this is not
    # a required value and it is not obvious how trustworthy it is.
    #
    # Since N_CHAN can be a VLF we can not just use sum().
    #
    numelt = 0
    for row in data["n_chan"]:
        try:
            numelt += sum(row)
        except TypeError:
            # assumed to be a scalar
            numelt += row

    # Remove unwanted data
    #  - rows where N_GRP=0
    #  - for each row, any excess data (beyond the sum(N_CHAN) for that row)
    #
    # The columns may be 1D or 2D vectors, or variable-field arrays,
    # where each row is an object containing a number of elements.
    #
    # Note that crates uses the sherpa.astro.utils.resp_init routine,
    # but it's not clear why, so it is not used here for now.
    #
    good = data['n_grp'] > 0

    matrix = data['matrix'][good]
    n_grp = data['n_grp'][good]
    n_chan = data['n_chan'][good]
    f_chan = data['f_chan'][good]

    # Flatten the array. There are four cases here:
    #
    # a) a variable-length field with no "padding"
    # b) a variable-length field with padding
    # c) a rectangular matrix is given, with no "padding"
    # d) a rectangular matrix is given, but a row can contain
    #    unused data
    #
    if numelt == matrix.size:
        if isinstance(matrix, _VLF):
            # case a
            matrix = np.concatenate([np.asarray(row)
                                     for row in matrix])
        else:
            # case c
            matrix = matrix.flatten()

    else:
        # cases b or d; assume we can use the same logic
        rowdata = []
        for mrow, ncs in zip(matrix, n_chan):
            # Need a RMF which ng>1 to test this with.
            if np.isscalar(ncs):
                ncs = [ncs]

            start = 0
            for nc in ncs:
                # n_chan can be an unsigned integer. Adding a Python
                # integer to a NumPy unsigned integer appears to return
                # a float.
                end = start + int(nc)

                # "perfect" RMFs may have mrow as a scalar
                try:
                    rdata = mrow[start:end]
                except IndexError as ie:
                    if start != 0 or end != 1:
                        raise IOErr('bad', 'format', 'MATRIX column formatting') from ie

                    rdata = [mrow]

                rowdata.append(rdata)
                start = end

        matrix = np.concatenate(rowdata)

    data['matrix'] = matrix.astype(SherpaFloat)

    # Flatten f_chan and n_chan vectors into 1D arrays as crates does
    # according to group. This is not always needed, but annoying to
    # check for so we always do it.
    #
    xf_chan = []
    xn_chan = []
    for grp, fch, nch, in zip(n_grp, f_chan, n_chan):
        if np.isscalar(fch):
            # The assumption is that grp==1 here
            xf_chan.append(fch)
            xn_chan.append(nch)
        else:
            # The arrays may contain extra elements (which
            # should be 0).
            #
            xf_chan.extend(fch[:grp])
            xn_chan.extend(nch[:grp])

    data['f_chan'] = np.asarray(xf_chan, SherpaUInt)
    data['n_chan'] = np.asarray(xn_chan, SherpaUInt)

    # Not all fields are "flattened" by this routine. In particular we
    # need to keep knowledge of these rows with 0 groups. This is why
    # we set the n_grp entry to data['n_grp'] rather than the n_grp
    # variable (which has the 0 elements removed).
    #
    data['n_grp'] = np.asarray(data['n_grp'], SherpaUInt)

    data['energ_lo'] = np.asarray(data['energ_lo'], SherpaFloat)
    data['energ_hi'] = np.asarray(data['energ_hi'], SherpaFloat)

    return data, filename


def get_pha_data(arg: DatasetType,
                 make_copy: bool = False,
                 use_background: bool = False
                 ) -> tuple[list[DataType], str]:
    """Read in the PHA."""

    pha, filename, close = _get_file_contents(arg,
                                              exptype="BinTableHDU")

    try:
        try:
            hdu = pha["SPECTRUM"]
        except KeyError:
            hdu = _has_ogip_type(pha, 'SPECTRUM')
            if hdu is None:
                raise IOErr('notrsp', filename, "a PHA spectrum") from None

        if use_background:
            for block in pha:
                if _try_key(block, 'HDUCLAS2') == 'BKG':
                    hdu = block

        keys = ['BACKFILE', 'ANCRFILE', 'RESPFILE',
                'BACKSCAL', 'AREASCAL', 'EXPOSURE']
        datasets = []

        if _try_col(hdu, 'SPEC_NUM') is None:
            data: DataType = {}

            # Create local versions of the "try" routines.
            #
            def try_any_sfloat(key):
                "Get col/key and return a SherpaFloat"
                return _try_col_or_key(hdu, key, fix_type=True)

            def try_sfloat(key):
                "Get col and return a SherpaFloat"
                return _try_col(hdu, key, fix_type=True)

            def try_sint(key):
                """Get col and return a SherpaInt

                Or, it looks like it should do this, but at the moment
                it does not force the data type.

                """
                # return _try_col(hdu, key, fix_type=True, dtype=SherpaInt)
                return _try_col(hdu, key, dtype=SherpaInt)

            def req_sfloat(key):
                "Get col and return a SherpaFloat"
                return _require_col(hdu, key, fix_type=True)

            # Keywords
            data['exposure'] = _try_key(hdu, 'EXPOSURE', fix_type=True)
            # data['poisserr'] = _try_key(hdu, 'POISSERR', True, bool)
            data['backfile'] = _try_key(hdu, 'BACKFILE')
            data['arffile'] = _try_key(hdu, 'ANCRFILE')
            data['rmffile'] = _try_key(hdu, 'RESPFILE')

            # Keywords or columns
            data['backscal'] = try_any_sfloat('BACKSCAL')
            data['backscup'] = try_any_sfloat('BACKSCUP')
            data['backscdn'] = try_any_sfloat('BACKSCDN')
            data['areascal'] = try_any_sfloat('AREASCAL')

            # Columns
            data['channel'] = req_sfloat("CHANNEL")

            # Make sure channel numbers not indices
            chan = list(hdu.columns.names).index('CHANNEL') + 1
            tlmin = _try_key(hdu, f"TLMIN{chan}", fix_type=True,
                             dtype=SherpaUInt)
            if int(data['channel'][0]) == 0 or tlmin == 0:
                data['channel'] = data['channel'] + 1

            data['counts'] = try_sfloat("COUNTS")
            data['staterror'] = _try_col(hdu, 'STAT_ERR')

            # The following assumes that EXPOSURE is set
            if data['counts'] is None:
                data['counts'] = req_sfloat("RATE") * data['exposure']
                if data['staterror'] is not None:
                    data['staterror'] = data['staterror'] * data['exposure']

            data['syserror'] = _try_col(hdu, 'SYS_ERR')
            data['background_up'] = try_sfloat('BACKGROUND_UP')
            data['background_down'] = try_sfloat('BACKGROUND_DOWN')
            data['bin_lo'] = try_sfloat('BIN_LO')
            data['bin_hi'] = try_sfloat('BIN_HI')
            data['grouping'] = try_sint('GROUPING')
            data['quality'] = try_sint('QUALITY')
            data['header'] = _get_meta_data(hdu)
            for key in keys:
                data['header'].pop(key, None)

            if data['syserror'] is not None:
                # SYS_ERR is the fractional systematic error
                data['syserror'] = data['syserror'] * data['counts']

            datasets.append(data)

        else:
            # Type 2 PHA file support
            specnum = _try_col_or_key(hdu, 'SPEC_NUM')
            num = len(specnum)

            # Create local versions of the "try" routines set for this
            # size=num value.
            #
            def try_any_sfloat(key):
                "Get col/key and return a SherpaFloat"
                return _try_vec_or_key(hdu, key, size=num, fix_type=True)

            def try_sfloat(key):
                "Get col and return a SherpaFloat"
                return _try_vec(hdu, key, size=num, fix_type=True)

            def try_sint(key):
                """Get col and return a SherpaInt

                Or, it looks like it should do this, but at the moment
                it does not force the data type.

                """
                # return _try_vec(hdu, key, size=num, fix_type=True, dtype=SherpaInt)
                return _try_vec(hdu, key, size=num, dtype=SherpaInt)

            def req_sfloat(key):
                "Get col and return a SherpaFloat"
                return _require_vec(hdu, key, size=num, fix_type=True)

            # Keywords
            exposure = _try_key(hdu, 'EXPOSURE', fix_type=True)
            # poisserr = _try_key(hdu, 'POISSERR', True, bool)
            backfile = _try_key(hdu, 'BACKFILE')
            arffile = _try_key(hdu, 'ANCRFILE')
            rmffile = _try_key(hdu, 'RESPFILE')

            # Keywords or columns
            backscal = try_any_sfloat('BACKSCAL')
            backscup = try_any_sfloat('BACKSCUP')
            backscdn = try_any_sfloat('BACKSCDN')
            areascal = try_any_sfloat('AREASCAL')

            # Columns
            # Why does this convert to SherpaFloat?
            channel = req_sfloat('CHANNEL')

            # Make sure channel numbers not indices
            chan = list(hdu.columns.names).index('CHANNEL') + 1
            tlmin = _try_key(hdu, f"TLMIN{chan}", fix_type=True,
                             dtype=SherpaUInt)

            for idx in range(num):
                if int(channel[idx][0]) == 0:
                    channel[idx] += 1

            # Why does this convert to SherpaFloat?
            counts = try_sfloat("COUNTS")
            staterror = _try_vec(hdu, 'STAT_ERR', size=num)
            if np.equal(staterror, None).any():
                counts = req_sfloat("RATE") * exposure
                if not np.equal(staterror, None).any():
                    staterror *= exposure

            syserror = _try_vec(hdu, 'SYS_ERR', size=num)
            background_up = try_sfloat('BACKGROUND_UP')
            background_down = try_sfloat('BACKGROUND_DOWN')
            bin_lo = try_sfloat('BIN_LO')
            bin_hi = try_sfloat('BIN_HI')
            grouping = try_sint('GROUPING')
            quality = try_sint('QUALITY')

            orders = try_sint('TG_M')
            parts = try_sint('TG_PART')
            specnums = try_sint('SPEC_NUM')
            srcids = try_sint('TG_SRCID')

            # Iterate over all rows of channels, counts, errors, etc
            # Populate a list of dictionaries containing
            # individual dataset info
            for (bscal, bscup, bscdn, arsc, chan, cnt, staterr, syserr,
                 backup, backdown, binlo, binhi, group, qual, ordr, prt,
                 specnum, srcid
                 ) in zip(backscal, backscup, backscdn, areascal, channel,
                          counts, staterror, syserror, background_up,
                          background_down, bin_lo, bin_hi, grouping, quality,
                          orders, parts, specnums, srcids):

                idata: DataType = {}

                idata['exposure'] = exposure
                # idata['poisserr'] = poisserr
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
                idata['grouping'] = group
                idata['quality'] = qual
                idata['header'] = _get_meta_data(hdu)
                idata['header']['TG_M'] = ordr
                idata['header']['TG_PART'] = prt
                idata['header']['SPEC_NUM'] = specnum
                idata['header']['TG_SRCID'] = srcid

                for key in keys:
                    idata['header'].pop(key, None)

                if syserr is not None:
                    # SYS_ERR is the fractional systematic error
                    idata['syserror'] = syserr * cnt

                datasets.append(idata)

    finally:
        if close:
            pha.close()

    return datasets, filename


# Write Functions

def _create_table(names: NamesType,
                  data: Mapping[str, Optional[np.ndarray]]) -> Table:
    """Create a Table.

    The idea is that by going via a Table we let the AstroPy code deal
    with all the conversions (e.g. to get the correct FITS data types
    on columns).

    Parameters
    ----------
    names : list of str
        The order of the column names (must exist in data).
    data : dict[str, ndarray or None]
        Any None values are dropped.

    Returns
    -------
    tbl : Table

    """

    store = []
    colnames = []
    for name in names:
        coldata = data[name]
        if coldata is None:
            continue

        store.append(coldata)
        colnames.append(name)

    return Table(names=colnames, data=store)


def check_clobber(filename: str, clobber: bool) -> None:
    """Error out if the file exists and clobber is not set."""

    if clobber or not os.path.isfile(filename):
        return

    raise IOErr("filefound", filename)


def pack_table_data(data: ColumnsType,
                    col_names: NamesType,
                    header: Optional[HdrTypeArg] = None) -> fits.BinTableHDU:
    """Pack up the table data."""

    tbl = _create_table(col_names, data)
    hdu = fits.table_to_hdu(tbl)

    # Add in the header. We should special case the HISTORY/COMMENT
    # keywords but at the moment we do nothing.
    #
    # There is the issue of conflicts between the existing and sent-in
    # header - fortunately we can just drop those keys from the sent-in
    # header, and also just whether the sent-in header keys (e.g. if
    # TUNIT1 and TLMIN1 are set are they valid)? For the latter we rely
    # on validation done by the code calling set_table_data.
    #
    if header is not None:
        _update_header(hdu, header)

    return hdu


def set_table_data(filename: str,
                   data: ColumnsType,
                   col_names: NamesType,
                   header: Optional[HdrTypeArg] = None,
                   ascii: bool = False,
                   clobber: bool = False) -> None:
    """Write out the table data."""

    check_clobber(filename, clobber)

    if ascii:
        tbl = _create_table(col_names, data)
        tbl.write(filename, format='ascii.commented_header',
                  overwrite=clobber)
        return

    hdu = pack_table_data(data, col_names, header)
    hdu.writeto(filename, overwrite=True)


def _create_header(header: HdrTypeArg) -> fits.Header:
    """Create a FITS header with the contents of header,
    the Sherpa representation of the key,value store.
    """

    hdrlist = fits.Header()
    for key, value in header.items():
        _add_keyword(hdrlist, key, value)

    return hdrlist


def _update_header(hdu: HDUType,
                   header: Mapping[str, Optional[KeyType]]) -> None:
    """Update the header of the HDU.

    Unlike the dict update method, this is left biased, in that
    it prefers the keys in the existing header (this is so that
    structural keywords are not over-written by invalid data from
    a previous FITS file), and to drop elements which are set
    to None.

    Parameters
    ----------
    hdu : HDU
    header : dict

    """

    for key, value in header.items():
        if key in hdu.header:
            continue

        # The value is not expected to be None but check anyway.
        if value is None:
            continue

        hdu.header[key] = value


def pack_arf_data(data: ColumnsType,
                  col_names: NamesType,
                  header: Optional[HdrTypeArg] = None) -> fits.BinTableHDU:
    """Pack the ARF"""

    if header is None:
        raise ArgumentTypeErr("badarg", "header", "set")

    return pack_table_data(data, col_names, header)


def set_arf_data(filename: str,
                 data: ColumnsType,
                 col_names: NamesType,
                 header: Optional[HdrTypeArg] = None,
                 ascii: bool = False,
                 clobber: bool = False) -> None:
    """Write out the ARF"""

    if header is None:
        raise ArgumentTypeErr("badarg", "header", "set")

    # This does not use pack_arf_data as we need to deal with ASCII
    # support.
    #
    set_table_data(filename, data, col_names, header=header,
                   ascii=ascii, clobber=clobber)


def pack_pha_data(data: ColumnsType,
                  col_names: NamesType,
                  header: Optional[HdrTypeArg] = None) -> fits.BinTableHDU:
    """Pack the PHA data."""

    if header is None:
        raise ArgumentTypeErr("badarg", "header", "set")

    return pack_table_data(data, col_names, header)


def set_pha_data(filename: str,
                 data: ColumnsType,
                 col_names: NamesType,
                 header: Optional[HdrTypeArg] = None,
                 ascii: bool = False,
                 clobber: bool = False) -> None:
    """ Write out the PHA."""

    if header is None:
        raise ArgumentTypeErr("badarg", "header", "set")

    # This does not use pack_pha_data as we need to deal with ASCII
    # support.
    #
    set_table_data(filename, data, col_names, header=header,
                   ascii=ascii, clobber=clobber)


def pack_rmf_data(blocks) -> fits.HDUList:
    """Pack up the RMF data."""

    # For now assume only two blocks:
    #    MATRIX
    #    EBOUNDS
    #
    matrix_data, matrix_header = blocks[0]
    ebounds_data, ebounds_header = blocks[1]

    # Extract the data:
    #   MATRIX:
    #     ENERG_LO
    #     ENERG_HI
    #     N_GRP
    #     F_CHAN
    #     N_CHAN
    #     MATRIX
    #
    #   EBOUNDS:
    #     CHANNEL
    #     E_MIN
    #     E_MAX
    #
    # We may need to convert F_CHAN/N_CHAN/MATRIX to Variable-Length
    # Fields. This is only needed if the ndarray type is object.
    #
    # It looks like we can not use the go-via-a-table approach
    # here, as that does not support Variable Length Fields.
    #
    def get_format(val):
        """We only bother with formats we expect"""
        if np.issubdtype(val, np.integer):
            if isinstance(val, np.int16):
                return "I"

            # default to Int32. AstroPy does support Int64/K but
            # this should not be needed here.
            return "J"

        if not np.issubdtype(val, np.floating):
            raise ValueError(f"Unexpected value '{val}' with type {val.dtype}")

        if isinstance(val,np.float32):
            return "E"

        # This should not be reached
        return "D"

    def get_full_format(name, vals):
        if vals.dtype == object:
            # We need to get the format from the first non-empty row.
            bformat = None
            ny = 0
            for v in vals:
                nv = len(v)
                ny = max(ny, nv)
                if nv == 0:
                    continue

                if bformat is not None:
                    continue

                bformat = get_format(v[0])
                break

            if bformat is None:
                # This should not happen
                raise ValueError(f"Unable to find data for column '{name}'")

            return f"P{bformat}({ny})"

        bformat = get_format(vals[0][0])
        ny = vals.shape[1]
        return f"{ny}{bformat}"

    def arraycol(name):
        vals = matrix_data[name]
        formatval = get_full_format(name, vals)
        return fits.Column(name=name, format=formatval, array=vals)

    n_grp = matrix_data["N_GRP"]
    col1 = fits.Column(name="ENERG_LO", format="E",
                       array=matrix_data["ENERG_LO"], unit="keV")
    col2 = fits.Column(name="ENERG_HI", format="E",
                       array=matrix_data["ENERG_HI"], unit="keV")
    col3 = fits.Column(name="N_GRP", format=get_format(n_grp[0]),
                       array=n_grp)
    col4 = arraycol("F_CHAN")
    col5 = arraycol("N_CHAN")
    col6 = arraycol("MATRIX")
    cols = [col1, col2, col3, col4, col5, col6]
    matrix_hdu = fits.BinTableHDU.from_columns(cols)
    _update_header(matrix_hdu, matrix_header)

    # Ensure the TLMIN4 value (for F_CHAN) is set.
    matrix_hdu.header["TLMIN4"] = matrix_data["OFFSET"]

    channel = ebounds_data["CHANNEL"]
    col1 = fits.Column(name="CHANNEL", format=get_format(channel[0]),
                       array=channel)
    col2 = fits.Column(name="E_MIN", format="E",
                       array=ebounds_data["E_MIN"], unit="keV")
    col3 = fits.Column(name="E_MAX", format="E",
                       array=ebounds_data["E_MAX"], unit="keV")
    ebounds_hdu = fits.BinTableHDU.from_columns([col1, col2, col3])
    _update_header(ebounds_hdu, ebounds_header)

    primary_hdu = fits.PrimaryHDU()
    return fits.HDUList([primary_hdu, matrix_hdu, ebounds_hdu])


def set_rmf_data(filename: str,
                 blocks,
                 clobber: bool = False) -> None:
    """Save the RMF data to disk.

    Unlike the other save_*_data calls this does not support the ascii
    argument. It also relies on the caller to have set up the headers
    and columns correctly apart for variable-length fields, which are
    limited to F_CHAN, N_CHAN, and MATRIX.

    """

    check_clobber(filename, clobber)
    hdus = pack_rmf_data(blocks)
    hdus.writeto(filename, overwrite=True)


def pack_image_data(data: DataTypeArg,
                    header: HdrTypeArg) -> fits.PrimaryHDU:
    """Pack up the image data."""

    hdrlist = _create_header(header)

    # Write Image WCS Header Keys
    if data['eqpos'] is not None:
        cdeltw = data['eqpos'].cdelt
        crpixw = data['eqpos'].crpix
        crvalw = data['eqpos'].crval
        equin = data['eqpos'].equinox

    if data['sky'] is not None:
        cdeltp = data['sky'].cdelt
        crpixp = data['sky'].crpix
        crvalp = data['sky'].crval

        _add_keyword(hdrlist, 'MTYPE1', 'sky     ')
        _add_keyword(hdrlist, 'MFORM1', 'x,y     ')
        _add_keyword(hdrlist, 'CTYPE1P', 'x      ')
        _add_keyword(hdrlist, 'CTYPE2P', 'y      ')
        _add_keyword(hdrlist, 'WCSNAMEP', 'PHYSICAL')
        _add_keyword(hdrlist, 'CDELT1P', cdeltp[0])
        _add_keyword(hdrlist, 'CDELT2P', cdeltp[1])
        _add_keyword(hdrlist, 'CRPIX1P', crpixp[0])
        _add_keyword(hdrlist, 'CRPIX2P', crpixp[1])
        _add_keyword(hdrlist, 'CRVAL1P', crvalp[0])
        _add_keyword(hdrlist, 'CRVAL2P', crvalp[1])

    if data['eqpos'] is not None:
        if data['sky'] is not None:
            # Simply the inverse of read transformations in get_image_data
            cdeltw = cdeltw * cdeltp
            crpixw = (crpixw - crvalp) / cdeltp + crpixp

        _add_keyword(hdrlist, 'MTYPE2', 'EQPOS   ')
        _add_keyword(hdrlist, 'MFORM2', 'RA,DEC  ')
        _add_keyword(hdrlist, 'CTYPE1', 'RA---TAN')
        _add_keyword(hdrlist, 'CTYPE2', 'DEC--TAN')
        _add_keyword(hdrlist, 'CDELT1', cdeltw[0])
        _add_keyword(hdrlist, 'CDELT2', cdeltw[1])
        _add_keyword(hdrlist, 'CRPIX1', crpixw[0])
        _add_keyword(hdrlist, 'CRPIX2', crpixw[1])
        _add_keyword(hdrlist, 'CRVAL1', crvalw[0])
        _add_keyword(hdrlist, 'CRVAL2', crvalw[1])
        _add_keyword(hdrlist, 'CUNIT1', 'deg     ')
        _add_keyword(hdrlist, 'CUNIT2', 'deg     ')
        _add_keyword(hdrlist, 'EQUINOX', equin)

    return fits.PrimaryHDU(data['pixels'], header=fits.Header(hdrlist))


def set_image_data(filename: str,
                   data: DataTypeArg,
                   header: HdrTypeArg,
                   ascii: bool = False,
                   clobber: bool = False) -> None:
    """Write out the image data."""

    check_clobber(filename, clobber)

    if ascii:
        set_arrays(filename, [data['pixels'].ravel()],
                   ascii=True, clobber=clobber)
        return

    img = pack_image_data(data, header)
    img.writeto(filename, overwrite=True)


def set_arrays(filename: str,
               args: Sequence[np.ndarray],
               fields: Optional[NamesType] = None,
               ascii: bool = True,
               clobber: bool = False) -> None:

    # Historically the clobber command has been checked before
    # processing the data, so do so here.
    #
    check_clobber(filename, clobber)

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
        fieldnames = [f'COL{idx + 1}' for idx in range(nargs)]
    elif nargs == len(fields):
        fieldnames = list(fields)
    else:
        raise IOErr("wrongnumcols", nargs, len(fields))

    if ascii:
        # Historically, the serialization doesn't quite match the
        # AstroPy Table version, so stay with write_arrays. We do
        # change the form for the comment line (used when fields is
        # set) to be "# <col1> .." rather than "#<col1> ..." to match
        # the AstroPy table output used for other "ASCII table" forms
        # in this module.
        #
        # The fields setting can be None here, which means that
        # write_arrays will not write out a header line.
        #
        write_arrays(filename, args, fields,
                     comment="# ", clobber=clobber)
        return

    data = dict(zip(fieldnames, args))
    tbl = _create_table(fieldnames, data)
    hdu = fits.table_to_hdu(tbl)
    hdu.name = 'TABLE'
    hdu.writeto(filename, overwrite=True)


def _add_header(hdu: HDUType,
                header: Sequence[HeaderItem]) -> None:
    """Add the header items to the HDU."""

    for hdr in header:
        card = [hdr.name, hdr.value]
        if hdr.desc is not None or hdr.unit is not None:
            comment = "" if hdr.unit is None else f"[{hdr.unit}] "
            if hdr.desc is not None:
                comment += hdr.desc
            card.append(comment)

        hdu.header.append(tuple(card))


def _create_primary_hdu(hdu: TableHDU) -> fits.PrimaryHDU:
    """Create a PRIMARY HDU."""

    out = fits.PrimaryHDU()
    _add_header(out, hdu.header)
    return out


def _create_table_hdu(hdu: TableHDU) -> fits.BinTableHDU:
    """Create a Table HDU."""

    if hdu.data is None:
        raise ValueError("No column data to write out")

    # First create a Table which handles the FITS column settings
    # correctly.
    #
    store = []
    colnames = []
    for col in hdu.data:
        colnames.append(col.name)
        store.append(col.values)

    out = fits.table_to_hdu(Table(names=colnames, data=store))
    out.name = hdu.name
    _add_header(out, hdu.header)

    # Add any column metadata
    #
    for idx, col in enumerate(hdu.data, 1):
        if col.unit is not None:
            out.header.append((f"TUNIT{idx}", col.unit))

        if col.desc is None:
            continue

        key = f"TTYPE{idx}"
        out.header[key] = (col.name, col.desc)

    return out


def pack_hdus(blocks: Sequence[TableHDU]) -> fits.HDUList:
    """Create a dataset.

    At present we are restricted to tables only.
    """

    out = fits.HDUList()
    out.append(_create_primary_hdu(blocks[0]))
    for hdu in blocks[1:]:
        out.append(_create_table_hdu(hdu))

    return out


def set_hdus(filename: str,
             blocks: Sequence[TableHDU],
             clobber: bool = False) -> None:
    """Write out multiple HDUS to a single file.

    At present we are restricted to tables only.
    """

    check_clobber(filename, clobber)
    hdus = pack_hdus(blocks)
    hdus.writeto(filename, overwrite=True)
