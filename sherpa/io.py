#
#  Copyright (C) 2007, 2015, 2016, 2019 - 2021, 2023, 2024
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

"""I/O routines for Sherpa.

These routines are currently restricted to reading from ASCII files.
"""

import os
from typing import Optional, Sequence

import numpy as np

from sherpa.data import Data, Data1D
from sherpa.utils import is_subclass, get_num_args, is_binary_file
from sherpa.utils.err import IOErr
from sherpa.utils.numeric_types import SherpaFloat
from sherpa.utils.types import ArrayType


__all__ = ('read_data', 'write_data', 'get_ascii_data', 'read_arrays',
           'write_arrays')


NamesType = Sequence[str]


def _check_args(size: int, dstype) -> None:
    # Find the number of required args minus self, filename
    req_args = get_num_args(dstype.__init__)[1] - 2
    if size >= req_args:
        return

    raise TypeError(f"data set '{dstype.__name__}' takes at "
                    f"least {req_args} args")


def read_file_data(filename: str,
                   sep: str = ' ',
                   comment: str = '#',
                   require_floats: bool = True
                   ) -> tuple[list[str], list[np.ndarray]]:
    """Read in column data from a file."""

    bad_chars = '\t\n\r,;: |'
    raw_names = []
    rows = []

    ncols = None
    with open(filename, 'r', encoding="utf-8") as fh:
        for line in fh:
            for char in bad_chars:
                if char in line:
                    # replace any bad chars in line with sep for tokenize
                    line = line.replace(char, sep)

            line = line.strip()

            # look for last commented line before data
            if len(line) > 0 and line[0] == comment:
                # Slice off the comment
                # TODO: why is this not just `line = line[1:]`?
                line = line.replace(comment, ' ')
                raw_names = line.strip().split(sep)

            elif line == '':
                continue
            else:
                # split line at sep
                elems = line.strip().split(sep)
                row = [elem for elem in elems if elem != '']

                # make list of row elements
                rows.append(row)

                if ncols is None:
                    ncols = len(row)
                elif ncols != len(row):
                    raise IOErr('arraysnoteq')

    if ncols is None:
        raise IOErr(f"No column data found in {filename}")

    # rotate rows into list of columns
    cols = np.column_stack(rows)

    # cast columns to appropriate type
    args = []
    for col in cols:
        try:
            args.append(col.astype(SherpaFloat))
        except ValueError as ve:
            if require_floats:
                raise ValueError(f"The file {filename} could not "
                                 "be loaded, probably because it "
                                 "contained spurious data and/or "
                                 "strings") from ve
            args.append(col)

    names = [name.strip(bad_chars)
             for name in raw_names if name != '']
    nargs = len(args)
    # TODO: should this error out if nargs == 0?

    if len(names) == 0:
        names = [f'col{i}' for i in range(1, nargs + 1)]

    # TODO: This could error out if len(names) > nargs, but this might
    # break existing code, since there is such a check in
    # get_ascii_data but it's only triggered when the colkeys argument
    # is set. To avoid breaking code we leave as is for now.
    #
    return names, args


def get_column_data(*args) -> list[Optional[np.ndarray]]:
    """
    get_column_data( *NumPy_args )
    """
    if len(args) == 0:
        raise IOErr('noarrays')

    cols = []
    for arg in args:
        if arg is None or isinstance(arg, (np.ndarray, list, tuple)):
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
                   ncols: int = 1,
                   colkeys: Optional[NamesType] = None,
                   sep: str = ' ',
                   dstype: type = Data1D,
                   comment: str = '#',
                   require_floats: bool = True
                   ) -> tuple[list[str], list[np.ndarray], str]:
    r"""Read in columns from an ASCII file.

    Parameters
    ----------
    filename : str
       The name of the ASCII file to read in.
    ncols : int, optional
       The number of columns to read in (the first ``ncols`` columns
       in the file). This is ignored if ``colkeys`` is given.
    colkeys : array of str, optional
       An array of the column name to read in. The default is
       `None`.
    sep : str, optional
       The separator character. The default is ``' '``.
    dstype : data class to use, optional
       Used to check that the data file contains enough columns.
    comment : str, optional
       The comment character. The default is ``'#'``.
    require_floats : bool, optional
       If `True` (the default), non-numeric data values will
       raise a `ValueError`.

    Returns
    -------
    (colnames, coldata, filename)
       The column names read in, the data for the columns
       as an array, with each element being the data for the column
       (the order matches ``colnames``), and the name of the file.

    Raises
    ------
    sherpa.utils.err.IOErr
       Raised if a requested column is missing or the file appears
       to be a binary file.
    ValueError
       If a column value can not be converted into a numeric value
       and the `require_floats` parameter is `True`.

    See Also
    --------
    read_arrays, read_data, write_arrays, write_data

    Notes
    -----
    The file is processed by reading in each line, stripping out any
    unsupported characters (replacing them by the ``sep`` argument),
    skipping empty lines, and then identifying comment and data lines.

    The list of unsupported characters are: ``\t``, ``\n``,
    ``\r``, comma, semi-colon, colon, space, and ``|``.

    The last comment line before the data is used to define the
    column names, splitting the line by the ``sep`` argument.
    If there are no comment lines then the columns are named
    starting at ``col1``, ``col2``, up to the number of columns.

    Data lines are separated into columns - splitting by the
    ``sep`` comment - and then converted to NumPy arrays.
    If the ``require_floats`` argument is `True` then the
    column will be converted to a floating-point number
    type, with an error raised if this fails.

    An error is raised if the number of columns per row
    is not constant.

    If the ``colkeys`` argument is used then a case-sensitive
    match is used to determine what columns to return.

    Examples
    --------

    Read in the first column from the file:

    >>> (colnames, coldata, fname) = get_ascii_data('src.dat')

    Read in the first three columns from the file:

    >>> colinfo = get_ascii_data('src.dat', ncols=3)

    Read in a histogram data set, using the columns XLO, XHI,
    and Y:

    >>> cols = ['XLO', 'XHI', 'Y']
    >>> res = get_ascii_data('hist.dat', colkeys=cols,
                             dstype=sherpa.data.Data1DInt)

    Read in the first and third column from the file cols.dat,
    where the file has no header information:

    >>> res = get_ascii_data('cols.dat', colkeys=['col1', 'col3'])

    """

    if is_binary_file(filename):
        raise IOErr('notascii', filename)

    names, args = read_file_data(filename, sep=sep, comment=comment,
                                 require_floats=require_floats)

    if colkeys is None:
        kwargs = []
        if ncols != 1:
            _check_args(ncols, dstype)
        kwargs.extend(args[:ncols])
        return (names, kwargs, filename)

    kwargs = []
    colkeys = list(colkeys)

    nnames = len(names)
    nargs = len(args)
    if nnames > nargs:
        raise IOErr('wrongnumcols', nargs, nnames)

    for key in colkeys:
        if key not in names:
            raise IOErr('reqcol', key, names)
        kwargs.append(args[names.index(key)])

    _check_args(len(kwargs), dstype)
    return (colkeys, kwargs, filename)


def read_data(filename: str,
              ncols: int = 2,
              colkeys: Optional[NamesType] = None,
              sep: str = ' ',
              dstype=Data1D,
              comment: str = '#',
              require_floats: bool = True) -> Data:
    """Create a data object from an ASCII file.

    Parameters
    ----------
    filename : str
       The name of the ASCII file to read in.
    ncols : int, optional
       The number of columns to read in (the first ``ncols`` columns
       in the file). This is ignored if `colkeys` is given.
    colkeys : array of str, optional
       An array of the column name to read in. The default is
       `None`.
    sep : str, optional
       The separator character. The default is ``' '``.
    dstype : data class to use, optional
       The class of the data object to create.
    comment : str, optional
       The comment character. The default is ``'#'``.
    require_floats : bool, optional
       If `True` (the default), non-numeric data values will
       raise a `ValueError`.

    Returns
    -------
    data
       The data object (created by calling the dstype constructor
       with the filename and then the data columns from the file).

    Raises
    ------
    sherpa.utils.err.IOErr
       Raised if a requested column is missing or the file appears
       to be a binary file.
    ValueError
       If a column value can not be converted into a numeric value
       and the `require_floats` parameter is True.

    See Also
    --------
    get_ascii_data, read_arrays, write_data

    Notes
    -----

    The file format is described in `get_ascii_data`.

    Examples
    --------

    Create a 1D data object from the first two columns in the file:

    >>> dat = read_data('src.dat')

    Use the third column as the error column (statistical):

    >>> dat = read_data('src.dat', ncols=3)

    Read in a histogram data set, using the columns XLO, XHI,
    and Y:

    >>> cols = ['XLO', 'XHI', 'Y']
    >>> dat = read_data('hist.dat', colkeys=cols,
                        dstype=sherpa.data.Data1DInt)

    Use the first and third column from the file cols.dat,
    where the file has no header information:

    >>> dat = read_data('cols.dat', colkeys=['col1', 'col3'])

    """

    _, args, name = get_ascii_data(filename, ncols=ncols,
                                   colkeys=colkeys, sep=sep,
                                   dstype=dstype, comment=comment,
                                   require_floats=require_floats)
    return dstype(name, *args)


def read_arrays(*args) -> Data:
    """Create a data object from arrays.

    Parameters
    ----------
    col1, ... coln : array_like
       The data columns.
    dstype : optional
       The data type to create. It must be a subclass of
       `sherpa.data.Data` and defaults to `sherpa.data.Data1D`

    Returns
    -------
    data
       The data object (created by calling the dstype constructor
       with the filename and then the data columns from the file).

    Raises
    ------
    sherpa.utils.err.IOErr
       Raised if no arrays are sent in.

    See Also
    --------
    get_ascii_data, write_arrays

    Examples
    --------

    Create a 1D data object from the x and y arrays:

    >>> dat = read_arrays(x, y)

    Include a statistical error column:

    >>> dat = read_arrays(x, y, dy)

    Create an integrated (i.e. histogram) data set:

    >>> dat = read_arrays(xlo, xhi, y, dstype=sherpa.data.Data1DInt)

    """
    largs = list(args)
    if len(largs) == 0:
        raise IOErr('noarrays')

    if is_subclass(largs[-1], Data):
        dstype = largs.pop()
    else:
        dstype = Data1D

    dargs = get_column_data(*largs)

    # Determine max number of args for dataset constructor
    _check_args(len(dargs), dstype)

    return dstype('', *dargs)


def write_arrays(filename: str,
                 args: Sequence[ArrayType],
                 fields: Optional[NamesType] = None,
                 sep: str = ' ',
                 comment: str = '#',
                 clobber: bool = False,
                 linebreak: str = '\n',
                 format: str = '%g') -> None:
    """Write a list of arrays to an ASCII file.

    Parameters
    ----------
    filename : str
       The name of the file to write the array to.
    args : array_like
       The arrays to write out.
    fields : array_like of str
       The column names (should match the size of `args` if given).
    sep : str, optional
       The separator character. The default is ``' '``.
    comment : str, optional
       The comment character. The default is ``'#'``. This is only used
       to write out the column names when `fields` is not `None`.
    clobber : bool, optional
       If `filename` is not `None`, then this flag controls
       whether an existing file can be overwritten (``True``)
       or if it raises an exception (``False``, the default
       setting).
    linebreak : str, optional
       Indicate a new line. The default is ``'\\n'``.
    format : str, optional
       The format used to write out the numeric values. The
       default is ``'%g%'``.

    Raises
    ------
    sherpa.utils.err.IOErr
       If `filename` already exists and `clobber` is `False`
       or if there is no data to write.

    See Also
    --------
    get_ascii_data

    Examples
    --------

    Write the x and y arrays to the file 'src.dat':

    >>> write_arrays('src.dat', [x, y])

    Use the column names "r" and "surbri" for the columns:

    >>> write_arrays('prof.txt', [x, y], fields=["r", "surbri"],
                     clobber=True)

    """
    if os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    # We assume the values are numeric but we never test for this
    # explicitly, nor do we require it in the types. This can make the
    # typing code get confused about what is allowed.
    #
    # In numpy 1.24 it became an error to pass in irregularly-gridded
    # data to asarray. Prior to that it would return an ndarray with a
    # dtype of object (and generate a deprecation warning).
    #
    narg = set()
    try:
        for arg in args:
            try:
                narg.add(len(arg))
            except TypeError:
                # len(arg) fails, so assume a scalar.
                narg.add(0)

    except TypeError:
        # args is not iterable, in which case narg will be empty and
        # caught below
        pass

    # Allow args to be a sequence of non-sequences or of sequences of
    # the same size. The former is technically not in the spirit of
    # the call but users may be taking advantage of it so do not error
    # out.
    #
    if len(narg) == 0 or 0 in narg:
        raise IOErr('noarrayswrite')

    if len(narg) != 1:
        raise IOErr('arraysnoteq')

    cols = np.column_stack(np.asarray(args))
    if cols.ndim < 2:
        raise IOErr('noarrayswrite')

    lines = []
    for col in cols:
        line = [format % elem for elem in col]
        lines.append(sep.join(line))

    with open(filename, 'w', encoding="utf-8") as fh:

        if fields is not None:
            fh.write(comment + sep.join(fields) + linebreak)

        fh.write(linebreak.join(lines))

        # add a newline at end
        fh.write(linebreak)


def write_data(filename: str,
               dataset: Data,
               fields: Optional[NamesType] = None,
               sep: str = ' ',
               comment: str = '#',
               clobber: bool = False,
               linebreak: str = '\n',
               format: str = '%g'
               ) -> None:
    """Write out a dataset as an ASCII file.

    Parameters
    ----------
    filename : str
       The name of the file to write the array to.
    dataset :
       The data object to write out.
    fields : array_like of str
       The column names (should match the size of ``args`` if given).
       Any unknown columns are skipped. If not given then the field
       names from the data set will be used (for those columns which
       contain data).
    sep : str, optional
       The separator character. The default is ``' '``.
    comment : str, optional
       The comment character. The default is ``'#'``. This is used to
       write out the column names (after converting to upper case)
       before the data.
    clobber : bool, optional
       If `filename` is not `None`, then this flag controls
       whether an existing file can be overwritten (`True`)
       or if it raises an exception (`False`, the default
       setting).
    linebreak : str, optional
       Indicate a new line. The default is ``'\\n'``.
    format : str, optional
       The format used to write out the numeric values. The
       default is ``'%g%'``.

    Raises
    ------
    sherpa.utils.err.IOErr
       If `filename` already exists and `clobber` is `False`
       or if there is no data to write.

    See Also
    --------
    get_ascii_data, read_data

    Examples
    --------

    Write the x and y arrays to the file 'src.dat':

    >>> write_data('src.dat', dat)

    """

    if fields is None:
        fields = dataset._fields

    cols = []
    col_names = []

    for name in fields:
        field = getattr(dataset, name, None)
        if field is not None and name != 'name':
            col_names.append(name.upper())
            cols.append(field)

    write_arrays(filename, cols, fields=col_names, sep=sep,
                 comment=comment, clobber=clobber,
                 linebreak=linebreak, format=format)
