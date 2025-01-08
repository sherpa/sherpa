#
#  Copyright (C) 2007, 2015 - 2019, 2021 - 2024
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
"""Provide Astronomy-specific I/O routines for Sherpa.

This module contains read and write routines for handling `FITS
<https://en.wikipedia.org/wiki/FITS>`_ and ASCII format data. The
actual support is provided by the selected I/O backend package
(currently Crates, provided by CIAO, or the FITS support in the
AstroPy package, if installed).

Which backend is used?
----------------------

When this module is first imported, Sherpa tries to import the
backends installed with Sherpa in the order listed in the
``.sherpa.rc`` or ``.sherpa-standalone.rc`` file. The first module that imports
successfully is set as the active backend. The following command prints the
name of the backend:

   >>> from sherpa.astro import io
   >>> print(io.backend.name)

Change the backend
------------------

After the initial import, the backend can be changed by loading one of the
I/O backends shipped with sherpa (or any other module that provides the same
interface):

  >>> import sherpa.astro.io.pyfits_backend
  >>> io.backend = sherpa.astro.io.pyfits_backend

"""

from configparser import ConfigParser
from contextlib import suppress
import importlib
import importlib.metadata
import logging
import os
from pathlib import Path
import re
from typing import TYPE_CHECKING, Any, Callable, Mapping, Optional, \
    Sequence, Type, TypeVar, Union

import numpy as np

from sherpa import get_config
from sherpa.astro.data import DataIMG, DataIMGInt, DataARF, DataRMF, \
    DataPHA, DataRosatRMF
from sherpa.astro.utils import reshape_2d_arrays
from sherpa.data import Data, Data1D, Data2D, Data2DInt
from sherpa.io import _check_args
from sherpa.utils import is_subclass
from sherpa.utils.err import ArgumentErr, DataErr, IOErr
from sherpa.utils.numeric_types import SherpaFloat, SherpaUInt

from .types import KeyType, NamesType, HdrTypeArg, HdrType, DataType, \
    Header, HeaderItem, Column, Block, TableBlock, ImageBlock, \
    BlockList, BlockType, SpecrespBlock, MatrixBlock, EboundsBlock

# Responses can often be send as the Data object or the instrument
# version, so support it would be good to support this with types
# like
#
#  RMFType = Union[DataRMF, RMF1D]
#
# but that is annoying thanks to circular dependencies.
#
if TYPE_CHECKING:
    from sherpa.astro.instrument import ARF1D, RMF1D
    ARFType = Union[DataARF, ARF1D]
    RMFType = Union[DataRMF, RMF1D]

else:
    ARFType = DataARF
    RMFType = DataRMF


config = ConfigParser()
config.read(get_config())

# What should we use for the default package?
io_opt = [o.strip().lower() + '_backend' for o in
          config.get('options', 'io_pkg', fallback='dummy').split()]

ogip_emin_str = config.get('ogip', 'minimum_energy', fallback='1.0e-10')

ogip_emin : Optional[float]
if ogip_emin_str.upper() == 'NONE':
    ogip_emin = None
else:
    emsg = "Invalid value for [ogip] minimum_energy config value; " + \
           "it must be None or a float > 0"
    try:
        ogip_emin = float(ogip_emin_str)
    except ValueError as ve:
        raise ValueError(emsg) from ve

    if ogip_emin <= 0.0:
        raise ValueError(emsg)

for iotry in io_opt:
    with suppress(ImportError):
        backend = importlib.import_module('.' + iotry,
                                          package='sherpa.astro.io')
        break

else:
    # None of the options in the rc file work, e.g. because it's an
    # old file that does not have dummy listed
    import sherpa.astro.io.dummy_backend as backend


error = logging.getLogger(__name__).error
warning = logging.getLogger(__name__).warning
info = logging.getLogger(__name__).info

T = TypeVar('T')


__all__ = ('backend',
           'read_table', 'read_image', 'read_arf', 'read_rmf', 'read_arrays',
           'read_pha', 'write_image', 'write_pha', 'write_table',
           'write_arf', 'write_rmf',
           'pack_table', 'pack_image', 'pack_pha', 'read_table_blocks')


# Note: write_arrays is not included in __all__, so don't add to the
#       See Also section.
#
def read_arrays(*args) -> Data:
    """Create a dataset from multiple arrays.

    The return value defaults to a `sherpa.data.Data1D` instance,
    but this can be changed by supplying the required class
    as the last argument (anything that is derived from
    `sherpa.data.Data`).

    Parameters
    ----------
    *args
        There must be at least one argument. The number of
        arguments depends on the data type being read in.
        The supported argument types depends on the I/O backend
        in use (as supported by the ``get_column_data`` routine
        provided by the backend).

    Returns
    -------
    data : a sherpa.data.Data derived object

    Examples
    --------
    The following examples do not contain argument types specific to
    a particular I/O backend.

    Create a `sherpa.data.Data1D` instance from the data in the
    arrays ``x`` and ``y`` (taken to be the independent and
    dependent axes respectively):

    >>> d = read_arrays(x, y)

    As in the previous example, but explicitly declaring the data type:

    >>> d = read_arrays(x, y, sherpa.data.Data1D)

    Create a `sherpa.data.Data2D` instance with the independent
    axes ``x0`` and ``x1``, and dependent axis ``y``:

    >>> d = read_arrays(x0, x1, y, sherpa.data.Data2D)

    """
    largs = list(args)
    if len(largs) == 0:
        raise IOErr('noarrays')

    if is_subclass(largs[-1], Data):
        dstype = largs.pop()
    else:
        dstype = Data1D

    dargs = backend.get_column_data(*largs)

    # Determine max number of args for dataset constructor
    _check_args(len(dargs), dstype)

    return dstype('', *dargs)


def read_table(arg,
               ncols: int = 2,
               colkeys: Optional[NamesType] = None,
               dstype: Type[Data1D] = Data1D) -> Data1D:
    """Create a dataset from a tabular file.

    The supported file types (e.g. ASCII or FITS) depends on the
    selected I/O backend.

    Parameters
    ----------
    arg
        The name of the file or a representation of the file
        (the type depends on the I/O backend).
    ncols : int, optional
        The first ncols columns are read from the file.
    colkeys : None or list of strings, optional
        If given, select these columns from the file.
    dstype : optional
        The data type to create (it is expected to follow the
        `sherpa.data.Data` interface).

    Returns
    -------
    data : a sherpa.data.Data derived object

    See Also
    --------
    write_table

    Examples
    --------
    The following examples do not contain argument types specific to
    a particular I/O backend.

    Create a `sherpa.data.Data1D` object from the first two
    columns in the file ``src.fits``:

    >>> d = read_table('src.fits')

    Create a `sherpa.data.Data1DInt` object from the first three
    columns in the file ``src.fits``:

    >>> d = read_table('src.fits', ncols=3, dstype=Data1DInt)

    Create a `sherpa.data.Data1D` data set from the specified
    columns in ``tbl.fits``, where ``WLEN`` is used for the
    independent axis, ``FLUX`` the dependent axis, and
    ``FLUXERR`` for the statistical error on the dependent axis:

    >>> d = read_table('tbl.fits', colkeys=['WLEN', 'FLUX', 'FLUXERR'])

    """

    hdu, name = backend.get_table_data(arg, ncols, colkeys)

    cols = [col.values for col in hdu.columns]

    # Determine max number of args for dataset constructor
    _check_args(len(cols), dstype)
    return dstype(name, *cols)


# TODO: should this be exported?
#
def read_ascii(filename: str,
               ncols: int = 2,
               colkeys: Optional[NamesType] = None,
               dstype: Type[Data1D] = Data1D,
               **kwargs) -> Data1D:
    """Create a dataset from an ASCII tabular file.

    Parameters
    ----------
    filename : str
        The name of the file or a representation of the file
        (the type depends on the I/O backend).
    ncols : int, optional
        The first ncols columns are read from the file.
    colkeys : None or list of strings, optional
        If given, select these columns from the file.
    dstype : optional
        The data type to create (it is expected to follow the
        `sherpa.data.Data` interface).
    **kwargs
        The remaining arguments are passed through to the
        ``get_ascii_data`` routine of the I/O backend. It is
        expected that these include ``sep`` and ``comment``,
        which describe the column separators and indicator of
        a comment line, respectively.

    Returns
    -------
    data : a sherpa.data.Data derived object

    Examples
    --------
    The following examples do not contain argument types specific to
    a particular I/O backend.

    Create a `sherpa.data.Data1D` object from the first two
    columns in the file ``src.dat``:

    >>> d = read_ascii('src.dat')

    Create A `sherpa.data.Data1DInt` object from the first three
    columns in the file ``src.dat``:

    >>> d = read_ascii('src.dat', ncols=3, dstype=Data1DInt)

    Create a `sherpa.data.Data1D` data set from the specified
    columns in ``tbl.fits``, where ``WLEN`` is used for the
    independent axis, ``FLUX`` the dependent axis, and
    ``FLUXERR`` for the statistical error on the dependent axis:

    >>> d = read_ascii('tbl.fits', colkeys=['WLEN', 'FLUX', 'FLUXERR'])

    """

    hdu, name = backend.get_ascii_data(filename, ncols=ncols, colkeys=colkeys,
                                       dstype=dstype, **kwargs)

    cols = [col.values for col in hdu.columns]

    # Determine max number of args for dataset constructor
    _check_args(len(cols), dstype)
    return dstype(name, *cols)


# TODO: Can this read in Data2D or only DataIMG?
def read_image(arg,
               coord: str = 'logical',
               dstype: Type[Data2D] = DataIMG) -> Data2D:
    """Create an image dataset from a file.

    .. versionchanged:: 4.16.0
       Setting coord to a value other than 'logical' will now
       correctly change the coordinate setting for `DataIMG` datasets.

    Parameters
    ----------
    arg
        The name of the file or a representation of the file
        (the type depends on the I/O backend).
    coord : {'logical', 'physical', 'world'}, optional
        The type of axis coordinates to use. An error is raised if the
        file does not contain the necessary metadata for the selected
        coordinate system or the selected dstype does not support it.
    dstype : optional
        The data type to create (it is expected to follow the
        `sherpa.data.Data` interface).

    Returns
    -------
    data : a sherpa.data.Data derived object

    See Also
    --------
    write_image

    Examples
    --------
    The following examples do not contain argument types specific to
    a particular I/O backend.

    Create a `sherpa.astro.data.DataIMG` object from the FITS file
    ``img.fits``:

    >>> d = read_image('img.fits')

    Select the physical coordinate system from the file:

    >>> d = read_image('img.fits', coord='physical')

    """

    block, filename = backend.get_image_data(arg)

    img = block.image
    ny, nx = img.shape
    x0 = np.arange(nx, dtype=SherpaFloat) + 1.
    x1 = np.arange(ny, dtype=SherpaFloat) + 1.
    x0, x1 = reshape_2d_arrays(x0, x1)

    data = {"y": img.ravel(), "shape": img.shape}

    # We have to be careful with issubclass checks because
    #
    #   base class |   sub class
    #    Data2D        Data2DInt
    #    Data2D        DataIMG
    #    DataIMG       DataIMGInt
    #
    if issubclass(dstype, DataIMG):
        header = {}
        for item in _remove_structural_items(block.header):
            header[item.name] = item.value

        data["header"] = header
        data["sky"] = block.sky
        data["eqpos"] = block.eqpos

    # What are the independent axes?
    #
    if issubclass(dstype, (Data2DInt, DataIMGInt)):
        data['x0lo'] = x0 - 0.5
        data['x1lo'] = x1 - 0.5
        data['x0hi'] = x0 + 0.5
        data['x1hi'] = x1 + 0.5
    elif issubclass(dstype, Data2D):
        data['x0'] = x0
        data['x1'] = x1
    else:
        # backend.get_image_data should have caught this case
        #
        raise ArgumentErr("bad", "dstype argument",
                          "dstype is not derived from Data2D")

    # Note that we only set the coordinates after creating the
    # dataset, which assumes that the independent axes are always in
    # logical units. They may not be, but in this case we need to drop
    # the logical/physical/world conversion system (or improve it to
    # be more flexible). This is an attempt to resolve issue #1762,
    # where sherpa.astro.ui.load_image(infile, coord="physical") did
    # not behave sensibly (likely due to #1414 which was to address
    # issue #1380).
    #
    dataset = dstype(filename, **data)
    if isinstance(dataset, DataIMG):
        dataset.set_coord(coord)

    return dataset


def read_arf(arg) -> DataARF:
    """Create a DataARF object.

    Parameters
    ----------
    arg
        The name of the file or a representation of the file
        (the type depends on the I/O backend) containing the
        ARF data.

    Returns
    -------
    data : sherpa.astro.data.DataARF

    """

    block, filename = backend.get_arf_data(arg)

    header = {}
    exposure = None
    for item in _remove_structural_items(block.header):
        if item.name == "EXPOSURE":
            exposure = item.value
        else:
            header[item.name] = item.value

    data = {"energ_lo": block.rget("ENERG_LO").values,
            "energ_hi": block.rget("ENERG_HI").values,
            "specresp": block.rget("SPECRESP").values,
            "exposure": exposure,
            "header": header,
            "ethresh": ogip_emin}

    blo = block.get("BIN_LO")
    bhi = block.get("BIN_HI")
    if blo is not None and bhi is not None:
        data["bin_lo"] = blo.values
        data["bin_hi"] = bhi.values

    return DataARF(filename, **data)


def _extract_rmf(matrix: MatrixBlock,
                 ebounds: EboundsBlock,
                 filename: str) -> dict[str, Any]:
    """Extract the matrix data from the table.

    Parameters
    ----------
    matrix, ebounds : TableBlock
        The MATRIX and EBOUNDS blocks.
    filename : str
        The file name (used for error/warning messages).

    Returns
    -------
    data : dict
        The arguments needed for creating a RMF.

    """

    header = {}
    detchans = 0
    for item in _remove_structural_items(matrix.header):
        if item.name == "DETCHANS":
            detchans = int(item.value)
        else:
            header[item.name] = item.value

    # Validate the data:
    #  - check that OFFSET is set in either F_CHAN column of MATRIX
    #    or CHANNEL column of EBOUNDS
    #  - flatten out the F_CHAN, N_CHAN, and MATRIX columns
    #
    f_chan = matrix.rget("F_CHAN")
    channel = ebounds.rget("CHANNEL")

    if f_chan.minval is not None:
        offset = f_chan.minval
    elif channel.minval is not None:
        offset = channel.minval
    else:
        offset = 1
        error("Failed to locate TLMIN keyword for F_CHAN column "
              "in RMF file '%s'; Update the offset value in the "
              "RMF data set to the appropriate TLMIN value prior "
              "to fitting", filename)

    # The n_grp column remains as is but the other columns
    # - remove any rows where n_grp for that row is 0
    # - flatten out any 2D structure or variable-length
    #   data
    #
    # This is not designed to be the fastest code.
    #
    n_grp_raw = matrix.rget("N_GRP").values
    f_chan_raw = f_chan.values
    n_chan_raw = matrix.rget("N_CHAN").values
    mat_raw = matrix.rget("MATRIX").values

    good = n_grp_raw > 0
    ng_raw = n_grp_raw[good]
    fc_raw = f_chan_raw[good]
    nc_raw = n_chan_raw[good]
    mx_raw = mat_raw[good]

    # Since these values are sent to the fold_rmf routine in
    # sherpa.astro.utils._utils, we want the types to match up to
    # avoid conversion time in that routine. The n_grp, f_chan, and
    # n_chan arrays should be SherpaUIntArray and the matrix should be
    # given by sherpa.utils.numeric_types. This conversion should
    # probably be done by DataRMF itself though.
    #
    f_chan_l = []
    n_chan_l = []
    mat_l = []
    for ng, fc, nc, mxs in zip(ng_raw, fc_raw, nc_raw, mx_raw):
        # fc and nc may be scalars rather than an array (this is known
        # before the loop but for now we check for it each row).
        #
        try:
            f_chan_l.append(fc[:ng])
            n_chan_l.append(nc[:ng])
        except IndexError:
            f_chan_l.append([fc])
            n_chan_l.append([nc])

        # Loop through all the matrix elements, restricting to nchan
        # (in case we do not have a VLF matrix array).  The matrix may
        # also be a scalar.
        #
        start = 0
        for nchan in n_chan_l[-1]:
            # nchan may be an unsigned int, which could lead to type casting
            # errors, so explicitly convert to an int.
            end = start + int(nchan)
            try:
                mat_l.append(mxs[start:end])
            except IndexError:
                # Assume the scalar case
                mat_l.append([mxs])

            start = end

    def flatten(xs, dtype):
        return np.concatenate(xs).astype(dtype)

    return {"detchans": detchans,
            "energ_lo": matrix.rget("ENERG_LO").values,
            "energ_hi": matrix.rget("ENERG_HI").values,
            "n_grp": np.asarray(n_grp_raw, dtype=SherpaUInt),
            "f_chan": flatten(f_chan_l, dtype=SherpaUInt),
            "n_chan": flatten(n_chan_l, dtype=SherpaUInt),
            "matrix": flatten(mat_l, dtype=SherpaFloat),
            "offset": offset,
            "e_min": ebounds.rget("E_MIN").values,
            "e_max": ebounds.rget("E_MAX").values,
            "header": header,
            "ethresh": ogip_emin  # should we take this from the header?
            }


def read_rmf(arg) -> DataRMF:
    """Create a DataRMF object.

    Parameters
    ----------
    arg
        The name of the file or a representation of the file
        (the type depends on the I/O backend) containing the
        RMF data.

    Returns
    -------
    data : sherpa.astro.data.DataRMF or subclass

    """

    matrixes, ebounds, filename = backend.get_rmf_data(arg)

    # matrixes is guaranteed not to be empty
    nmat = len(matrixes)
    if nmat > 1:
        # Warn the user that the multi-matrix RMF is not supported.
        #
        error("RMF in %s contains %d MATRIX blocks; "
              "Sherpa only uses the first block!",
              filename, nmat)

    # Use the metadata from the first matrix
    #
    matrix = matrixes[0]
    data = _extract_rmf(matrix, ebounds, filename)
    return _rmf_factory(filename, data)


def _rmf_factory(filename: str,
                 data: Mapping[str, Any]) -> DataRMF:
    response_map = {
        'ROSAT': DataRosatRMF,
        'DEFAULT': DataRMF,
    }

    rmf_telescop = data['header'].get('TELESCOP', 'DEFAULT')
    rmf_class = response_map.get(rmf_telescop, DataRMF)

    return rmf_class(filename, **data)


def _read_ancillary(header: Header,
                    key: str,
                    label: str,
                    dpath: Path,
                    read_func: Callable[[str], T],
                    output_once: bool) -> Optional[T]:
    """Read in a file if the keyword is set.

    Parameters
    ----------
    header : Header
        The header information. This may be updated.
    key : str
        The key to look for in data. If not set, or its value compares
        case insensitively to "none" then nothing is read in.
    label : str
        This is used to identify the file being read in. It is only
        used if output_once is set.
    dname : Path
        The location of the file containing the metadata. This is prepended
        to  the value read from the data argument if it is not an
        absolute path.
    read_func
        The function used to read in the file: it expects one argument,
        the file to read in, and returns the data object.
    output_once : bool, optional
        If set then the routine uses the Sherpa logger to display
        informational or warning information.

    Returns
    -------
    data : None or a sherpa.data.Data derived object

    """

    val = header.get(key)
    if val is None:
        return None

    if not isinstance(val.value, str):
        return None

    if val.value.lower() == "none":
        return None

    # Rely on pathlib to handle the path combinations.
    #
    filename = str(dpath / Path(val.value))
    out = None
    try:
        out = read_func(filename)
        if output_once:
            info("read %s file %s", label, filename)

    except Exception as exc:
        if output_once:
            warning("unable to read %s: %s", label, str(exc))

    return out


def _process_pha_block(filename: str,
                       table: TableBlock,
                       *,
                       output_once: bool,
                       use_errors: bool,
                       use_background: bool) -> DataPHA:
    """Convert a PHA block into a DataPHA object.

    This is assumed to be a "type I" block.
    """

    # The assumption is that the table requires a CHANNEL column,
    # e.g. as guaranteed by read_pha, but leave this check in.
    #
    chancol = table.get("CHANNEL")
    if chancol is None:
        raise IOErr("reqcol", "CHANNEL", filename)

    channel = chancol.values
    nrows = channel.size

    expval = table.header.get("EXPOSURE")
    exposure = None if expval is None else expval.value

    def get(name: str,
            expand: bool = False) -> Optional[Union[np.ndarray, int, float, np.integer, np.floating]]:
        """Return the column values if they exist.

        This checks for columns and then the header. For the header
        case we can keep the scalar case (expand=False) or convert
        into a ndarray.

        """

        col = table.get(name)
        if col is not None:
            return col.values

        colval = table.header.get(name)
        if colval is None:
            return None

        if isinstance(colval.value, (str, np.floating)):
            val = float(colval.value)
        elif isinstance(colval.value, (bool, np.bool_, np.integer)):
            val = int(colval.value)
        else:
            val = colval.value

        if expand:
            return np.full(nrows, val)

        return val

    def getcol(name: str) -> Optional[np.ndarray]:
        """Return the column values if they exist.

        This does not check for a matching keyword to support the
        original code (even though the OGIP standard allows for these
        names).

        """

        col = table.get(name)
        if col is None:
            return None

        return col.values

    # Why do we remove the *FILE keywords? See #1885
    #
    unwanted = ["ANCRFILE", "BACKFILE", "RESPFILE",
                "BACKSCAL", "AREASCAL", "EXPOSURE"]
    header = _remove_structural_items(table.header)
    basic_header = dict((k.name, k.value) for k in header
                        if k.name not in unwanted)

    # The statistical and systematic errors can be
    # - not set
    # - a column
    # - a keyword
    #
    # However, we can have values of 0, which is not useful for
    # Sherpa, so we remove them (this is likely pointing at the
    # handling of the POISERR keyword being suboptimal).
    #
    staterr = get("STAT_ERR", expand=True)
    syserr = get("SYS_ERR", expand=True)
    if staterr is not None and np.all(staterr == 0.0):
        staterr = None

    if syserr is not None and np.all(syserr == 0.0):
        syserr = None

    # Do we warn the user that the errors read in from the file are
    # being ignored?
    #
    if not use_errors:
        if output_once and staterr is not None or syserr is not None:
            if staterr is None:
                msg = 'systematic'
            elif syserr is None:
                msg = 'statistical'

                # Extra warning
                warning("systematic errors were not found in "
                        "file '%s'", filename)

            else:
                msg = 'statistical and systematic'

            # TODO: remove extra space at end of first line
            info("%s errors were found in file '%s'\nbut not used; "
                 "to use them, re-read with use_errors=True", msg, filename)

        staterr = None
        syserr = None

    # For relative path we use the location of the PHA file as the
    # base.
    #
    dpath = Path(filename).parent
    albl = 'ARF'
    rlbl = 'RMF'
    if use_background:
        albl = albl + ' (background)'
        rlbl = rlbl + ' (background)'

    arf = _read_ancillary(table.header, "ANCRFILE", albl, dpath,
                          read_arf, output_once)
    rmf = _read_ancillary(table.header, "RESPFILE", rlbl, dpath,
                          read_rmf, output_once)

    # Do we create the backgrounds directly?
    #
    kwargs = {"channel": channel,
              "bin_lo": get("BIN_LO"),
              "bin_hi": get("BIN_HI"),
              # "grouping": get("GROUPING", expand=True),
              # "quality": get("QUALITY", expand=True),
              "grouping": getcol("GROUPING"),
              "quality": getcol("QUALITY"),
              "exposure": exposure,
              "header": basic_header}

    pha = DataPHA(filename, counts=get("COUNTS"),
                  staterror=staterr, syserror=syserr,
                  backscal=get("BACKSCAL"), areascal=get("AREASCAL"),
                  **kwargs)
    pha.set_response(arf, rmf)

    # Are there any backgrounds to process?
    #
    backgrounds = []

    # Do we read in the background data?
    # TODO: The use_background logic doesn't seem to make sense here.
    #
    backfile = table.header.get("BACKFILE")
    if backfile is not None and isinstance(backfile.value, str) \
       and backfile.value.lower() != "none" and \
           not use_background:

        bfilename = str(dpath / Path(backfile.value))
        try:
            backgrounds = _process_pha_backfile(bfilename, arf, rmf,
                                                output_once=output_once,
                                                use_errors=use_errors)
        except Exception as exc:
            if output_once:
                warning("unable to read background: %s", str(exc))

    for colname, scalname in [('BACKGROUND_UP', 'BACKSCUP'),
                              ('BACKGROUND_DOWN', 'BACKSCDN')]:
        col = table.get(colname)
        if col is None:
            continue

        bkg = DataPHA(filename, counts=col.values,
                      backscal=get(scalname), **kwargs)
        bkg.set_response(arf, rmf)
        if output_once:
            info("read %s into a dataset from file %s",
                 colname.lower(), filename)

        backgrounds.append(bkg)

    for idx, bkg in enumerate(backgrounds, 1):
        if bkg.grouping is None:
            bkg.grouping = pha.grouping
            bkg.grouped = bkg.grouping is not None

        if bkg.quality is None:
            bkg.quality = pha.quality

        pha.set_background(bkg, idx)

    # set units *after* bkgs have been set
    pha._set_initial_quantity()

    return pha


def _process_pha_backfile(filename: str,
                          arf: Optional[DataARF],
                          rmf: Optional[DataRMF],
                          *,
                          output_once: bool,
                          use_errors: bool) -> list[DataPHA]:

    bkgs = read_pha(filename, use_errors=use_errors, use_background=True)
    if output_once:
        info("read background file %s", filename)

    # read_pha can return a single item or a list.
    if isinstance(bkgs, DataPHA):
        blist = [bkgs]
    else:
        blist = list(bkgs)

    for bkg in blist:
        if rmf is not None and bkg.get_response() == (None, None):
            bkg.set_response(arf, rmf)

    return blist


def read_pha(arg,
             use_errors: bool = False,
             use_background: bool = False
             ) -> Union[DataPHA, list[DataPHA]]:
    """Create a DataPHA object.

    .. versionchanged:: 4.17.0
       Channel numbers that start at 0 are now left as is rather than
       be renumbered to start at 1.

    Parameters
    ----------
    arg
        The name of the file or a representation of the file
        (the type depends on the I/O backend) containing the
        PHA data.
    use_errors : bool, optional
        If the PHA file contains statistical error values for the
        count (or count rate) column, should it be read in. This
        defaults to ``False``.
    use_background : bool, optional
        Should the background PHA data (and optional responses) also
        be read in and associated with the data set?

    Returns
    -------
    data : sherpa.astro.data.DataPHA

    """

    pha, filename = backend.get_pha_data(arg,
                                         use_background=use_background)

    # Check all the columns have the same shape (e.g. all type I or
    # II).  The CHANNEL column is guaranteed to exist at this point
    # and be 1 or 2D.
    #
    chancol = pha.rget("CHANNEL")
    channel = chancol.values
    if len(channel.shape) == 1:
        ptype = 1
    elif len(channel.shape) == 2:
        ptype = 2
    else:
        # This should not be possible
        raise IOErr(f"PHA {filename}: unable to handle CHANNEL shape")

    # Prior to 4.17.0 the code tried to ensure that the first channel
    # was numbered 1. The approach is to now use the values read in
    # from file / memory.
    #
    # Convert from RATE to COUNT. Ideally this would be passed to
    # DataPHA so it knew (and could then enforce a non-Poisson
    # statistic).
    #
    # This is done before any II to I conversion. The error checks
    # should not be needed, because SpectrumBlock should enforce these
    # invariants, but leave in for now.
    #
    if pha.get("COUNTS") is None:
        rate = pha.get("RATE")
        if rate is None:
            raise IOErr("reqcol", "COUNTS or RATE", filename)

        expval = pha.header.get("EXPOSURE")
        if expval is None:
            raise IOErr("nokeyword", filename, "EXPOSURE")

        exposure = expval.value

        # Change the RATE column into a "COUNTS" column.
        #
        rate.values *= exposure
        rate.name = "COUNTS"

        # What other columns need to be converted?
        #
        # How about SYS_ERR, BACKGROUND_UP/DOWN?
        #
        # TODO: this likely fails if STAT_ERR is a keyword
        # and not a column.
        #
        stat_err = pha.get("STAT_ERR")
        if stat_err is not None:
            stat_err.values *= exposure

    if ptype == 1:
        phas = [pha]
    else:
        special = ["TG_M", "TG_PART", "TG_SRCID", "SPEC_NUM"]
        tg_m = pha.get("TG_M")
        tg_part = pha.get("TG_PART")
        tg_srcid = pha.get("TG_SRCID")
        spec_num = pha.get("SPEC_NUM")

        phas = []
        for idx in range(channel.shape[0]):
            header = Header([item for item in pha.header.values
                             if item.name not in special])

            if tg_m is not None:
                header.add(HeaderItem("TG_M", tg_m.values[idx]))

            if tg_part is not None:
                header.add(HeaderItem("TG_PART", tg_part.values[idx]))

            if tg_srcid is not None:
                header.add(HeaderItem("TG_SRCID", tg_srcid.values[idx]))

            if spec_num is not None:
                header.add(HeaderItem("SPEC_NUM", spec_num.values[idx]))

            cols = [Column(col.name, values=col.values[idx].copy(),
                           desc=col.desc, unit=col.unit,
                           minval=col.minval, maxval=col.maxval)
                    for col in pha.columns
                    if col.name not in special]

            phas.append(TableBlock(pha.name, header=header, columns=cols))

    output_once = True
    datasets = []
    for p in phas:
        data = _process_pha_block(filename, p,
                                  output_once=output_once,
                                  use_errors=use_errors,
                                  use_background=use_background)
        output_once = False
        datasets.append(data)

    if len(datasets) == 1:
        return datasets[0]

    return datasets


def _empty_header(creator=False) -> Header:
    """Create an empty header."""

    out = Header([])
    if creator:
        _add_creator(out)

    return out


def _pack_header(header: HdrTypeArg) -> Header:
    """Convert a header and add a CREATOR keyword if not set."""

    out = Header([HeaderItem(name=name, value=value)
                  for name, value in header.items()])
    _add_creator(out)
    return out


def _pack_table(dataset: Data1D) -> BlockList:
    """Identify the columns in the data.

    This relies on the _fields attribute containing the data columns,
    and _extra_fields the extra information (other than the name
    column). We only return values from _fields that contain data.

    """

    cols = []
    for name in dataset._fields:
        if name == 'name':
            continue

        val = getattr(dataset, name)
        if val is None:
            continue

        cols.append(Column(name.upper(), val))

    empty = _empty_header()

    blocks: list[BlockType]
    blocks = [TableBlock(name="TABLE", header=empty, columns=cols)]

    return BlockList(header=empty, blocks=blocks)


def _pack_image(dataset: Data2D) -> ImageBlock:
    """Extract the needed data of a 2D dataset."""

    if not isinstance(dataset, Data2D):
        raise IOErr('notimage', dataset.name)

    img = dataset.get_img()

    # Data2D does not have a header or sky/eqpos
    #
    hdr = getattr(dataset, "header", {})
    sky = getattr(dataset, "sky", None)
    eqpos = getattr(dataset, "eqpos", None)

    header = Header([HeaderItem(name, value)
                     for name, value in hdr.items()])

    return ImageBlock("IMAGE", header=header, image=img,
                      sky=sky, eqpos=eqpos)


# This is to match NAXIS1 and the like, but we do not make it
# very prescriptive (ie enforce the FITS keyword standards).
#
NUMBERED_KEYWORD_PAT = re.compile(r"^([^\d]+)(\d+)$")


def _is_structural_keyword(key: str) -> bool:
    """Do we think this is a "structural" keyword?

    Parameters
    ----------
    key : str
       The keyword (need not be upper case).

    Returns
    -------
    flag : bool
       Is this considered a structural keyword.

    Notes
    -----
    The decision of what is a structural header and what isn't is
    rather ad-hoc and may need to be updated over time. The current
    code is based in part on the `HEASARC keyword list
    <https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict>`_.

    This does not strip out the OGIP keywords like HDUCLAS1 or
    HDUVERS.

    It does not remove secondary WCS keywords, like CTYPE1P.

    """

    ukey = key.upper()
    if ukey in ["BITPIX", "BLANK", "BLOCKED", "BSCALE", "BUNIT",
                "BZERO", "DATAMAX", "DATAMIN", "END", "EXTEND",
                "EXTLEVEL", "EXTNAME", "EXTVER", "GCOUNT", "GROUPS",
                "NAXIS", "PCOUNT", "SIMPLE", "TFIELDS", "THEAP",
                "XTENSION",
                # other common ones
                "CHECKSUM", "DATASUM", "CHECKVER",
                # HEASARC
                "TSORTKEY",
                # CIAO uses this as well as EXTNAME
                "HDUNAME"
                ]:
        return True

    # Keywords which end in 1, 2, ..
    match = NUMBERED_KEYWORD_PAT.match(ukey)
    if match is None:
        return False

    # The pattern is not completely robust so catch a case when
    # it does not restrict to an integer > 0.
    #
    if match[2].startswith("0"):
        return False

    return match[1] in ["CDELT", "CROTA", "CRPIX", "CRVAL", "CTYPE", "CUNIT",
                        "NAXIS", "PSCAL", "PTYPE", "PZERO", "TBCOL",
                        "TDIM", "TDISP", "TFORM", "TNULL", "TSCAL",
                        "TTYPE", "TUNIT", "TZERO",
                        "TCTYP", "TCUNI", "TCRPX", "TCRVL", "TCDLT",
                        # HEASARC
                        "TDMAX", "TDMIN", "TLMAX", "TLMIN",
                        # CIAO Data-Model keywords
                        "DSTYP", "DSVAL", "DSFORM", "DSUNIT", "TDBIN",
                        "MTYPE", "MFORM",
                        ]


def _remove_structural_items(header: Header) -> list[HeaderItem]:
    """Remove FITS keywords relating to file structure.

    The aim is to allow writing out a header that was taken from a
    file and not copy along "unwanted" values, since the data may
    no-longer be relevant. Not all FITS header keys may be passed
    along by a particular backend (e.g. crates).

    Parameters
    ----------
    header : Header
       The input header

    Returns
    -------
    nheader : list of HeaderItem
       Header with unwanted keywords returned (including those
       set to None).

    Notes
    -----
    The tricky part is knowing what is unwanted, since copying
    along a WCS transform for a column may be useful, or may
    break things.

    """

    return [item for item in header.values
            if not _is_structural_keyword(item.name)]


# By making this a package-level value, users can actually change the
# text (and in fact the name) by using code like
#
#    from sherpa.astro.io import CREATOR
#    CREATOR.value = "not sherpa"
#
# but it does not seem worth documenting this at this time.
#
CREATOR = HeaderItem(name="CREATOR",
                     value=f"sherpa {importlib.metadata.version('sherpa')}",
                     desc="Program creating this file")


def _add_creator(header: Header) -> None:
    """Add a CREATOR card if not set."""

    if header.get("CREATOR") is not None:
        return

    header.add(CREATOR)


def _pack_pha(dataset: DataPHA) -> BlockList:
    """Extract FITS column and header information.

    Notes
    -----
    The `PHA Data Extension header page
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/node6.html>`_
    lists the following keywords as either required or
    we-really-want-them:

        EXTNAME = "SPECTRUM"
        TELESCOP - the "telescope" (i.e. mission/satellite name).
        INSTRUME - the instrument/detector.
        FILTER - the instrument filter in use (if any)
        EXPOSURE - the integration time (in seconds) for the PHA data (assumed to be corrected for deadtime, data drop-outs etc. )
        BACKFILE - the name of the corresponding background file (if any)
        CORRFILE - the name of the corresponding correction file (if any)
        CORRSCAL - the correction scaling factor.
        RESPFILE - the name of the corresponding (default) redistribution matrix file (RMF; see George et al. 1992a).
        ANCRFILE - the name of the corresponding (default) ancillary response file (ARF; see George et al. 1992a).
        HDUCLASS = "OGIP"
        HDUCLAS1 = "SPECTRUM"
        HDUVERS = "1.2.1"
        POISSERR - whether Poissonian errors are appropriate to the data (see below).
        CHANTYPE - whether the channels used in the file have been corrected in any way (see below).
        DETCHANS - the total number of detector channels available.

    We also add in the following, defaulting to the first value - we
    should do better to support HDUCLAS3=RATE data!

        HDUCLAS2 - indicating the type of data stored.
          Allowed values are:
            'TOTAL' for a gross PHA Spectrum (source + bkgd)
            'NET' for a bkgd-subtracted PHA Spectrum
            'BKG' for a bkgd PHA Spectrum
        HDUCLAS3 - indicating further details of the type of data stored.
          Allowed values are:
            'COUNT' for PHA data stored as counts (rather than count/s)
            'RATE' for PHA data stored in count/s
        HDUCLAS4 - indicating whether this is a type I or II extension.
          Allowed values are:
            'TYPE:I' for type I (single spectrum) data
            'TYPE:II' for type II (multiple spectra) data

    The POISSERR keyword is not required if a STAT_ERR column is
    present however it is recommended in this case for clarity. If
    STAT_ERR is to be used for the errors then POISSERR is set to
    false.

    The TLMIN1 and TLMAX1 values are always set to match the CHANNEL
    data range. This was changed in 4.17.0 as previously they were
    only set when the channel range did not start at 1.

    """

    # The logic here repeats some of the checks that probably should
    # be done by the DataPHA class itself. However, it is likely
    # that we don't want to make the DataPHA class always reject
    # inconsistent state, as this could preclude certain workflows,
    # so we need some validation here.
    #
    if not isinstance(dataset, DataPHA):
        raise IOErr("notpha", dataset.name)

    arf, rmf = dataset.get_response()
    bkg = dataset.get_background()

    # The default keywords; these will be over-ridden by
    # anything set by the input.
    #
    default_header = {
        "EXTNAME": "SPECTRUM",
        "HDUCLASS": "OGIP",
        "HDUCLAS1": "SPECTRUM",
        "HDUCLAS2": "TOTAL",
        "HDUCLAS3": "COUNT",
        "HDUCLAS4": "TYPE:I",
        "HDUVERS": "1.2.1",
        "HDUDOC": "Arnaud et al. 1992a Legacy 2  p 65",

        # The PHA file should have an EXPOSURE value, but just in case.
        "EXPOSURE": 1.0,

        # Rely on the DataPHA class to have set up TELESCOP/INSTRUME/FILTER
        # based on any associated background or response. If the user has
        # changed them then so be it.
        #
        "TELESCOP": "none",
        "INSTRUME": "none",
        "FILTER": "none",
        "CORRFILE": "none",
        "CORRSCAL": 0,
        "CHANTYPE": "PI",
        "RESPFILE": "none",
        "ANCRFILE": "none",
        "BACKFILE": "none"
    }

    # Merge the keywords
    #
    merged_header = default_header | dataset.header
    header = _pack_header(merged_header)

    if dataset.channel is not None and header.get("DETCHANS") is None:
        header.add(HeaderItem("DETCHANS", len(dataset.channel),
                              desc="Number of channels"))

    # Over-write the header value (if set). This value should really
    # exist (OGIP standards) but may not, particularly for testing.
    #
    if dataset.exposure is not None:
        header.delete("EXPOSURE")
        header.add(HeaderItem("EXPOSURE", float(dataset.exposure),
                              unit="s", desc="Exposure time"))

    # See #1885 for a discussion about the current behaviour.
    #
    if rmf is not None:
        header.delete("RESPFILE")
        header.add(HeaderItem("RESPFILE", rmf.name.split("/")[-1],
                              desc="RMF"))

    if arf is not None:
        header.delete("ANCRFILE")
        header.add(HeaderItem("ANCRFILE", arf.name.split("/")[-1],
                              desc="ARF"))
    if bkg is not None:
        header.delete("BACKFILE")
        header.add(HeaderItem("BACKFILE", bkg.name.split("/")[-1],
                              desc="Background"))

    # Create the columns. Special case the CHANNEL and COUNTS
    # columns.
    #
    # TODO: there is currently no way to know we should write out a
    # RATE column instead of COUNTS.
    #
    cols = []
    if dataset.channel is not None:
        channel = Column("CHANNEL", dataset.channel.astype(np.int32),
                         desc="Channel values")

        # Ensure the min/max ranges are written out as integers since
        # data.channel may contain floating point values. Always
        # re-write so that the file matches the data.
        #
        channel.minval = int(dataset.channel[0])
        channel.maxval = int(dataset.channel[-1])

        cols.append(channel)

    if dataset.counts is not None:
        if np.issubdtype(dataset.counts.dtype, np.integer):
            countvals = dataset.counts.astype(np.int32)
        elif np.issubdtype(dataset.counts.dtype, np.floating):
            countvals = dataset.counts.astype(np.float32)
        else:
            raise DataErr("ogip-error", "PHA dataset", dataset.name,
                          "contains an unsupported COUNTS column")

        cols.append(Column("COUNTS", countvals, desc="Counts"))

    if dataset.staterror is not None:
        desc = "Statistical error"
        if np.ptp(dataset.staterror) == 0:
            header.add(HeaderItem("STAT_ERR", dataset.staterror[0],
                                  desc=desc))
        else:
            cols.append(Column("STAT_ERR",
                               dataset.staterror.astype(np.float32),
                               desc=desc))

    def addfield(name, value, dtype, desc):
        """Add as keyword or array of the expected type"""

        if value is None:
            header.add(HeaderItem(name, dtype(0), desc=desc))
        elif np.ptp(value) == 0:
            header.add(HeaderItem(name, dtype(value[0]), desc=desc))
        else:
            cols.append(Column(name, value.astype(dtype), desc=desc))

    # Should this convert the value/values back into a fractional error?
    addfield("SYS_ERR", dataset.syserror, np.float32,
             "Fractional systematic error")

    addfield("GROUPING", dataset.grouping, np.int16,
             "Grouping (1,-1,0)")
    addfield("QUALITY", dataset.quality, np.int16,
             "Quality (0 for okay)")

    if dataset.bin_lo is not None and dataset.bin_hi is not None:
        cols.extend([Column("BIN_LO", dataset.bin_lo),
                     Column("BIN_HI", dataset.bin_hi)])

    def addscal(value, name, label):
        """Add the scaling value."""

        if value is None:
            header.add(HeaderItem(name, 1.0, desc=label))
        elif np.isscalar(value):
            header.add(HeaderItem(name, float(value), desc=label))
        else:
            # Ensure we have no scalar value in the header
            header.delete(name)
            # OGIP standard has this as a 4-byte real column/value
            cols.append(Column(name, value.astype(np.float32),
                               desc=label))

    # Note that CORRSCAL is set in default_header. Can it also
    # be a vector?
    #
    addscal(dataset.backscal, "BACKSCAL", "Background scaling factor")
    addscal(dataset.areascal, "AREASCAL", "Area scaling factor")

    # Try to set the POISSERR keyword *IF NOT SET*, but Sherpa does
    # not use this reliably.
    #
    if header.get("POISSERR") is None:
        header.add(HeaderItem("POISSERR",
                              dataset.staterror is None))

    pheader = _empty_header()

    blocks: list[BlockType]
    blocks = [TableBlock("SPECTRUM", header=header, columns=cols)]

    return BlockList(header=pheader, blocks=blocks)


def _pack_arf(dataset: ARFType) -> BlockList:
    """Extract FITS column and header information.

    There is currently no support for Type II ARF files.

    Notes
    -----
    The `ARF Data Extension header page
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html>`_
    lists the following keywords as either required or
    we-really-want-them:

        EXTNAME = 'SPECRESP'
        TELESCOP - the "telescope" (ie mission/satellite name).
        INSTRUME - the instrument/detector.
        FILTER - the instrument filter in use (if any)
        HDUCLASS = 'OGIP' - file format is OGIP standard.
        HDUCLAS1 = 'RESPONSE' - extension contains response data.
        HDUCLAS2 = 'SPECRESP' - extension contains an ARF.
        HDUVERS = '1.1.0' - version of the file format.

    We do not add the ARFVERSN, HDUVERS1, or HDUVERS2 keywords which
    are marked as obsolete.

    """

    # The UI layer wraps up a DataARF into ARF1D quite often, so
    # support both types as input.
    #
    from sherpa.astro.instrument import ARF1D
    if not isinstance(dataset, (DataARF, ARF1D)):
        raise IOErr("data set is not an ARF")

    # Check we have data. Should we allow missing columns?
    #
    for col in ["energ_lo", "energ_hi", "specresp"]:
        if getattr(dataset, col) is None:
            raise ArgumentErr("bad", "ARF",
                              f"{col.upper()} column is missing")

    # The default keywords; these will be over-ridden by
    # anything set by the input.
    #
    default_header: HdrType = {
        "EXTNAME": "SPECRESP",
        "HDUCLASS": "OGIP",
        "HDUCLAS1": "RESPONSE",
        "HDUCLAS2": "SPECRESP",
        "HDUVERS": "1.1.0",
        "TELESCOP": "none",
        "INSTRUME": "none",
        "FILTER": "none",
    }

    # Merge the keywords
    #
    merged_header = default_header | dataset.header
    aheader = _pack_header(merged_header)

    # The exposure time is not an OGIP-mandated value but
    # may be used by CIAO, so copy it across if set.
    #
    if dataset.exposure is not None:
        aheader.add(HeaderItem("EXPOSURE", float(dataset.exposure),
                               desc="Total exposure time",
                               unit="s"))

    # The column ordering for the output file is determined by the
    # order the keys are added to the data dict. Ensure the
    # data type meets the FITS standard (Real4).
    #
    acols = [Column("ENERG_LO", unit="keV", desc="Energy",
                    values=dataset.energ_lo.astype(np.float32)),
             Column("ENERG_HI", unit="keV", desc="Energy",
                    values=dataset.energ_hi.astype(np.float32)),
             Column("SPECRESP", unit="cm^2", desc="Effective Area",
                    values=dataset.specresp.astype(np.float32))]

    # Chandra files can have BIN_LO/HI values, so copy
    # across if both set.
    #
    blo = dataset.bin_lo
    bhi = dataset.bin_hi
    if blo is not None and bhi is not None:
        acols.extend([Column("BIN_LO", blo),
                      Column("BIN_HI", bhi)])

    pheader = _empty_header(creator=True)

    blocks: list[BlockType]
    blocks = [SpecrespBlock(name="SPECRESP", header=aheader, columns=acols)]

    return BlockList(header=pheader, blocks=blocks)


def _find_int_dtype(rows: Sequence[Sequence[int]]) -> type:
    """What data type should represent the matrix of integers?

    There is no guarantee that each row has the same number of
    elements.

    """

    for row in rows:
        if len(row) == 0:
            continue

        if np.max(row) > 32767:
            return np.int32

    return np.int16


def _make_int_array(rows: Sequence[Sequence[int]],
                    ncols: int) -> np.ndarray:
    """Convert a list of rows into a 2D array of "width" ncols.

    The conversion is to a type determined by the maximum value in
    rows (it selects between 2-byte and 4-byte integers).

    Parameters
    ----------
    rows : list of 1D arrays
        The values to convert. Expected to be integers (>= 0).
    ncols : int
        The size of each row in the output. There must be no row
        with more elements.

    Returns
    -------
    out : ndarray of size (nrows, ncols)

    """

    nrows = len(rows)
    dtype = _find_int_dtype(rows)
    out: np.ndarray = np.zeros((nrows, ncols), dtype=dtype)
    for idx, row in enumerate(rows):
        if len(row) == 0:
            continue

        out[idx, 0:len(row)] = row

    return out


def _make_float32_array(rows: Sequence[Sequence[float]],
                        ncols: int) -> np.ndarray:
    """Convert a list of rows into a 2D array of "width" ncols.

    The output has type numpy.float32.

    Parameters
    ----------
    rows : list of 1D arrays
        The values to convert. Expected to have type float32.
    ncols : int
        The size of each row in the output. There must be no row
        with more elements.

    Returns
    -------
    out : ndarray of size (nrows, ncols)

    """

    nrows = len(rows)
    out = np.zeros((nrows, ncols), dtype=np.float32)
    for idx, row in enumerate(rows):
        out[idx, 0:len(row)] = row

    return out


def _make_int_vlf(rows: Sequence[Sequence[int]]) -> np.ndarray:
    """Convert a list of rows into a VLF.

    The conversion is to a type determined by the maximum value in
    rows (it selects between 2-byte and 4-byte integers).

    Parameters
    ----------
    rows : list of 1D arrays
        The values to convert. Expected to be integers (>= 0).

    Returns
    -------
    out : ndarray of dtype object

    """

    dtype = _find_int_dtype(rows)
    out: list[np.ndarray] = []
    for row in rows:
        out.append(np.asarray(row, dtype=dtype))

    return np.asarray(out, dtype=object)


def _reconstruct_rmf(rmf: RMFType) -> DataType:
    """Recreate the structure needed to write out as a FITS file.

    This does not guarantee to create byte-identical data in a round
    trip, but it should create equivalent data (e.g. the choice of
    whether a column should be a Variable Length Field may differ, as
    can the data types).

    Parameters
    ----------
    rmf : DataRMF or RMF1D instance

    Returns
    -------
    out : dict
        The re-constructed data and header values.

    Notes
    -----
    The F_CHAN, N_CHAN, and MATRIX columns are stored as either
      - 1D arrays (unlikely for MATRIX); N_GRP = 1 for all rows
      - 2D arrays when N_GRP is a constant for all rows
      - 2D arrays when N_GRP <= 4 (so some rows will have 0's in
        when N_GRP for that row is less than the max N_GRP)
      - 2D array with dtype=object indicates a VLF (each row
        should be the correct type though)

    """

    n_grp = []
    f_chan: list[list[int]] = []
    n_chan: list[list[int]] = []
    matrix: list[list[float]] = []

    # Used to reconstruct the original data
    idx = 0
    start = 0

    # Header keyword
    numelt = 0

    # How big is the matrix?
    matrix_size = set()

    for ng in rmf.n_grp:
        n_grp.append(ng)

        f_chan.append([])
        n_chan.append([])
        matrix.append([])

        # Short-cut when no data
        if ng == 0:
            # Record we have a zero-element row. Should we instead
            # remove the last element of f_chan, n_chan, matrix?
            #
            matrix_size.add(0)
            continue

        # Grab the next ng elements from rmf.f_chan/n_chan
        # and the corresponding rmf.matrix elements. Fortunately
        # we can ignore F_CHAN when extracting data from matrix.
        #
        for _ in range(ng):
            # ensure this is an integer and not a floating-point number
            # which can happen when rmf.n_chan is stored as Unt64.
            #
            end = start + rmf.n_chan[idx].astype(np.int32)
            mdata = rmf.matrix[start:end]

            f_chan[-1].append(rmf.f_chan[idx])
            n_chan[-1].append(rmf.n_chan[idx])
            matrix[-1].extend(mdata.astype(np.float32))

            idx += 1
            start = end

        matrix_size.add(len(matrix[-1]))
        numelt += np.sum(n_chan[-1])

    # N_GRP should be 2-byte integer.
    #
    n_grp_out = np.asarray(n_grp, dtype=np.int16)
    numgrp = n_grp_out.sum()

    # Can we convert F_CHAN/N_CHAN to fixed-length if either:
    #  - N_GRP is the same for all rows
    #  - max(N_GRP) < 4
    #
    # The decision to convert to the int16 or int32 types is made
    # within the _make_int_xxx routines (the maximum value is used to
    # decide what type to use).
    #
    if len(set(n_grp_out)) == 1 or n_grp_out.max() < 4:
        ny = n_grp_out.max()
        f_chan_out = _make_int_array(f_chan, ny)
        n_chan_out = _make_int_array(n_chan, ny)
    else:
        f_chan_out = _make_int_vlf(f_chan)
        n_chan_out = _make_int_vlf(n_chan)

    # We can convert the matrix to fixed size if each row in matrix
    # has the same size.
    #
    # The individual matrix elements are of type np.float32 so we
    # should not need to do any conversion, but we are explicit in
    # the fixed-length case.
    #
    if len(matrix_size) == 1:
        ny = matrix_size.pop()
        matrix_out = _make_float32_array(matrix, ny)
    else:
        # Since the elements can be lists, ensure they get converted
        # to ndarray.
        #
        matrix_out = np.asarray([np.asarray(m, dtype=np.float32)
                                 for m in matrix],
                                dtype=object)

    return {"N_GRP": n_grp_out,
            "F_CHAN": f_chan_out,
            "N_CHAN": n_chan_out,
            "MATRIX": matrix_out,
            # Ensure these are integer types
            "NUMGRP": int(numgrp),
            "NUMELT": int(numelt)}


def _pack_rmf(dataset: RMFType) -> BlockList:
    """Extract FITS column and header information.

    Unlike the other pack routines this returns data for
    multiple blocks. At present this ignores the ORDER column
    in the MATRIX block.

    Notes
    -----
    The `RMF Data Extension header page
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html>`_
    lists the following keywords as either required or
    we-really-want-them. First for the MATRIX block:

        EXTNAME = 'MATRIX' or 'SPECRESP MATRIX'
        TELESCOP - the "telescope" (ie mission/satellite name).
        INSTRUME - the instrument/detector.
        FILTER - the instrument filter in use (if any)
        CHANTYPE - PI or PHA
        DETCHANS - the total number of raw detector PHA channels
        HDUCLASS = 'OGIP' - file format is OGIP standard.
        HDUCLAS1 = 'RESPONSE' - extension contains response data.
        HDUCLAS2 = 'RSP_MATRIX' - extension contains a response matrix.
        HDUVERS = '1.3.0' - version of the file format.
        TLMIN# - the first channel in the response.

    These are optional but mechanical so can be added:

        NUMGRP - sum(N_GRP)
        NUMELT - sum(N_CHAN)

    These are optional:

        LO_THRES - minimum probability threshold
        HDUCLAS3 - 'REDIST', 'DETECTOR', 'FULL'

    For the EBOUNDS extension we have:

        EXTNAME  = 'EBOUNDS'
        TELESCOP
        INSTRUME
        FILTER
        CHANTYPE
        DETCHANS
        HDUCLASS = 'OGIP'
        HDUCLAS1 = 'RESPONSE'
        HDUCLAS2 = 'EBOUNDS'
        HDUVERS = '1.2.0'

    We do not add the RMFVERSN, HDUVERS1, or HDUVERS2 keywords which
    are marked as obsolete.

    """

    # The UI layer wraps up a DataRMF into RMF1D quite often, so
    # support both types as input.
    #
    from sherpa.astro.instrument import RMF1D
    if not isinstance(dataset, (DataRMF, RMF1D)):
        raise IOErr("data set is not a RMF")

    # We convert the input into a dict (rmfdata) and then more
    # dictionaries (both for headers and columns) before eventually
    # creating the Block objects we return. This separates out the
    # data manipulation - e.g. creating the correct column types and
    # ensuring the headers contain the needed data - from creating the
    # file contents (i.e. Block objects). The two steps could be
    # combined but leave as is for now since this is not code that
    # needs to be fast, but does need to be readable.
    #
    rmfdata = _reconstruct_rmf(dataset)

    # The default keywords; these will be over-ridden by
    # anything set by the input.
    #
    default_matrix_header = {
        "EXTNAME": "MATRIX",
        "HDUCLASS": "OGIP",
        "HDUCLAS1": "RESPONSE",
        "HDUCLAS2": "RSP_MATRIX",
        "HDUVERS": "1.3.0",
        "TELESCOP": "none",
        "INSTRUME": "none",
        "FILTER": "none",
        "CHANTYPE": "PI",
        "DETCHANS": "none",
        "NUMGRP": 0,
        "NUMELT": 0
    }

    ebounds_header: HdrType = {
        "EXTNAME": "EBOUNDS",
        "HDUCLASS": "OGIP",
        "HDUCLAS1": "RESPONSE",
        "HDUCLAS2": "EBOUNDS",
        "HDUVERS": "1.2.0",
        "TELESCOP": "none",
        "INSTRUME": "none",
        "FILTER": "none",
        "CHANTYPE": "PI",
        "DETCHANS": "none"
    }

    # Header Keys
    header = dataset.header

    # Merge the keywords (at least for the MATRIX block).
    #
    matrix_header = default_matrix_header | header

    matrix_header["NUMGRP"] = rmfdata["NUMGRP"]
    matrix_header["NUMELT"] = rmfdata["NUMELT"]
    matrix_header["DETCHANS"] = dataset.detchans

    # Copy values over.
    #
    for copykey in ["TELESCOP", "INSTRUME", "FILTER", "CHANTYPE",
                    "DETCHANS"]:
        ebounds_header[copykey] = matrix_header[copykey]

    # The column ordering for the output file is determined by the
    # order the keys are added to the data dict.
    #
    # To allow the backend convert to Variable-Length fields, we
    # have F_CHAN/N_CHAN/MATRIX either be a 2D ndarray (so fixed
    # output) OR a list of rows (use VLF). It's not ideal: we could
    # wrap up in a local VLF type just to indicate this to the
    # backend.
    #
    matrix_data = {
        "ENERG_LO": dataset.energ_lo,
        "ENERG_HI": dataset.energ_hi,
        "N_GRP": rmfdata["N_GRP"],
        "F_CHAN": rmfdata["F_CHAN"],
        "N_CHAN": rmfdata["N_CHAN"],
        "MATRIX": rmfdata["MATRIX"],
        "OFFSET": dataset.offset  # copy over this value
    }

    # TODO: is this correct?
    nchan = dataset.offset + dataset.detchans - 1
    dchan = np.int32 if nchan > 32767 else np.int16

    # Technically e_min/max can be empty, but we do not expect
    # this, and this support should probably be removed. For
    # now error out if we are sent such data.
    #
    if dataset.e_min is None or dataset.e_max is None:
        raise IOErr(f"RMF {dataset.name} has no E_MIN or E_MAX data")

    ebounds_data = {
        "CHANNEL": np.arange(dataset.offset, nchan + 1, dtype=dchan),
        "E_MIN": dataset.e_min.astype(np.float32),
        "E_MAX": dataset.e_max.astype(np.float32)
        }

    # Now convert these into TableBlock types.
    #
    # Create the MATRIX block:
    #     ENERG_LO
    #     ENERG_HI
    #     N_GRP
    #     F_CHAN
    #     N_CHAN
    #     MATRIX
    #
    mheader = _pack_header(matrix_header)

    mcols = [Column(name=name, values=value)
             for name, value in matrix_data.items()
             if name != "OFFSET"]

    # Adjust the column objects.
    #
    for col in mcols:
        if col.name.startswith("ENERG"):
            col.unit = "keV"
            continue

        if col.name == "F_CHAN":
            col.minval = matrix_data["OFFSET"]
            continue

    # Create the EBOUNDS block:
    #     CHANNEL
    #     E_MIN
    #     E_MAX
    #
    eheader = _pack_header(ebounds_header)

    ecols = [Column(name=name, values=value)
             for name, value in ebounds_data.items()]

    for col in ecols:
        if col.name == "CHANNEL":
            continue

        col.unit = "keV"

    header = _empty_header()
    blocks: list[BlockType]
    blocks = [MatrixBlock(name="MATRIX", header=mheader, columns=mcols),
              EboundsBlock(name="EBOUNDS", header=eheader, columns=ecols)]

    return BlockList(header=header, blocks=blocks)


def write_arrays(filename: str,
                 args: Sequence[np.ndarray],
                 fields: Optional[NamesType] = None,
                 ascii: bool = True,
                 clobber: bool = False) -> None:
    """Write out a collection of arrays.

    Parameters
    ----------
    filename : str
       The name of the file.
    args
       The data to write out.
    fields : None or list of str, optional
       The names to use for each column. If ``None`` then the names
       default to `col1` ... `colN`.
    ascii : bool, optional
       If `True` use an ASCII format, otherwise a binary format.
    clobber : bool, optional
       If `True` then the output file will be over-written if it
       already exists, otherwise an error is raised.

    See Also
    --------
    read_arrays

    """

    backend.set_arrays(filename, args, fields, ascii=ascii, clobber=clobber)


def write_table(filename: str,
                dataset: Data1D,
                ascii: bool = True,
                clobber: bool = False) -> None:
    """Write out a table.

    Parameters
    ----------
    filename : str
       The name of the file.
    dataset
       The data to write out.
    ascii : bool, optional
       If `True` use an ASCII format, otherwise a binary format.
    clobber : bool, optional
       If `True` then the output file will be over-written if it
       already exists, otherwise an error is raised.

    See Also
    --------
    read_table

    """

    data = _pack_table(dataset)
    backend.set_table_data(filename, data, ascii=ascii, clobber=clobber)


def write_image(filename: str,
                dataset: Data2D,
                ascii: bool = True,
                clobber: bool = False) -> None:
    """Write out an image.

    Parameters
    ----------
    filename : str
       The name of the file.
    dataset
       The data to write out.
    ascii : bool, optional
       If `True` use an ASCII format, otherwise a binary format.
    clobber : bool, optional
       If `True` then the output file will be over-written if it
       already exists, otherwise an error is raised.

    See Also
    --------
    read_image

    """

    data = _pack_image(dataset)
    backend.set_image_data(filename, data, ascii=ascii, clobber=clobber)


def write_pha(filename: str,
              dataset: DataPHA,
              ascii: bool = True,
              clobber: bool = False) -> None:
    """Write out a PHA dataset.

    Parameters
    ----------
    filename : str
       The name of the file.
    dataset
       The data to write out.
    ascii : bool, optional
       If `True` use an ASCII format, otherwise a binary format.
    clobber : bool, optional
       If `True` then the output file will be over-written if it
       already exists, otherwise an error is raised.

    See Also
    --------
    read_pha

    """

    block = _pack_pha(dataset)
    backend.set_pha_data(filename, block, ascii=ascii,
                         clobber=clobber)


def write_arf(filename: str,
              dataset: ARFType,
              ascii: bool = True,
              clobber: bool = False) -> None:
    """Write out an ARF.

    This does not handle Type II files.

    .. versionadded:: 4.16.0

    Parameters
    ----------
    filename : str
       The name of the file.
    dataset
       The data to write out.
    ascii : bool, optional
       If `True` use an ASCII format, otherwise a binary format.
    clobber : bool, optional
       If `True` then the output file will be over-written if it
       already exists, otherwise an error is raised.

    See Also
    --------
    read_arf

    """

    blocks = _pack_arf(dataset)
    backend.set_arf_data(filename, blocks,
                         ascii=ascii, clobber=clobber)


def write_rmf(filename: str,
              dataset: RMFType,
              clobber: bool = False) -> None:
    """Write out a RMF.

    .. versionadded:: 4.16.0

    Parameters
    ----------
    filename : str
       The name of the file.
    dataset
       The data to write out.
    clobber : bool, optional
       If `True` then the output file will be over-written if it
       already exists, otherwise an error is raised.

    See Also
    --------
    read_rmf

    """

    blocks = _pack_rmf(dataset)
    backend.set_rmf_data(filename, blocks, clobber=clobber)


def pack_table(dataset: Data1D) -> object:
    """Convert a Sherpa data object into an I/O item (tabular).

    Parameters
    ----------
    dataset : a sherpa.data.Data1D derived object

    Returns
    -------
    obj
       An object representing the data, the format of which depends
       on the I/O backend.

    See Also
    --------
    pack_image, pack_pha

    Examples
    --------

    >>> d = sherpa.data.Data1D('tmp', [1, 2, 3], [4, 10, 2])
    >>> tbl = pack_table(d)

    """

    data = _pack_table(dataset)
    return backend.pack_table_data(data)


def pack_image(dataset: Data2D) -> Any:
    """Convert a Sherpa data object into an I/O item (image).

    Parameters
    ----------
    dataset : a sherpa.data.Data2D derived object

    Returns
    -------
    obj
       An object representing the data, the format of which depends
       on the I/O backend.

    See Also
    --------
    pack_image, pack_pha

    Examples
    --------

    >>> y, x = np.mgrid[:10, :5]
    >>> z = (x-2)**2 + (y-2)**3
    >>> d = sherpa.data.Data2D('img', x.flatten(), y.flatten(),
                               z.flatten(), shape=z.shape)
    >>> img = pack_image(d)

    """

    data = _pack_image(dataset)
    return backend.pack_image_data(data)


def pack_pha(dataset: DataPHA) -> Any:
    """Convert a Sherpa PHA data object into an I/O item (tabular).

    Parameters
    ----------
    dataset : a sherpa.astro.data.DataPHA derived object

    Returns
    -------
    obj
       An object representing the data, the format of which depends
       on the I/O backend.

    See Also
    --------
    pack_image, pack_table

    """

    block = _pack_pha(dataset)
    return backend.pack_pha_data(block)


def read_table_blocks(arg,
                      make_copy: bool = False
                      ) -> tuple[str,
                                 dict[int, dict[str, np.ndarray]],
                                 dict[int, HdrType]]:
    """Return the HDU elements (columns and header) from a FITS table.

    Parameters
    ----------
    arg
       The data file, which can be the name of the file or an object
       representing the opened file (the type of this depends on the I/O
       backend in use).
    make_copy : bool, optional
       This argument is currently unused.

    Returns
    -------
    filename, cols, hdr : str, dict, dict
       The name of the file, the column data, and the header data
       from the HDUs in the file. The keys of the dictionaries are the
       block name and the values are dictionaries, with the keys
       being the column or header name.

    """

    def getkeys(hdr):
        if hdr is None:
            return {}

        return {key.name: key.value for key in hdr.values}

    def getcols(cols):
        if cols is None:
            return {}

        return {col.name: col.values for col in cols}

    # Desconstruct the data. It is not clear what to do if the file
    # contains any image blocks, so we explicitly ignore them at this
    # time (i.e. there is currently no requirement on the
    # backend.read_table_blocks call on what to do in this situation).
    #
    blist, filename = backend.read_table_blocks(arg, make_copy=make_copy)
    hdr = {1: getkeys(blist.header)}

    cols: dict[int, dict[str, np.ndarray]]
    cols = {1: {}}
    for idx, block in enumerate(blist.blocks, 2):
        if isinstance(block, ImageBlock):
            continue

        hdr[idx] = getkeys(block.header)
        cols[idx] = getcols(block.columns)

    return filename, cols, hdr
