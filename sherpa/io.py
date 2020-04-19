#
#  Copyright (C) 2007, 2015, 2016, 2019  Smithsonian Astrophysical Observatory
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

import numpy
from sherpa.utils import SherpaFloat, get_num_args, is_binary_file
from sherpa.utils.err import IOErr
from sherpa.data import Data1D, BaseData
import os


__all__ = ('read_data', 'write_data', 'get_ascii_data', 'read_arrays',
           'write_arrays')


def _is_subclass(t1, t2):
    return isinstance(t1, type) and issubclass(t1, t2) and (t1 is not t2)


def _check_args(size, dstype):
    # Find the number of required args minus self, filename
    req_args = (get_num_args(dstype.__init__)[1] - 2)

    if size < req_args:
        # raise IOErr('badargs', dstype.__name__, req_args)
        raise TypeError("data set '%s' takes at least %s args" %
                        (dstype.__name__, req_args))


def read_file_data(filename, sep=' ', comment='#', require_floats=True):
    bad_chars = '\t\n\r,;: |'
    fp = open(filename, 'r')
    raw_names = []
    rows = []
    try:
        for line in fp:
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

    finally:
        fp.close()

    # rotate rows into list of columns
    cols = numpy.column_stack(rows)

    # cast columns to appropriate type
    args = []
    for col in cols:
        try:
            args.append(col.astype(SherpaFloat))
        except ValueError:
            if require_floats:
                raise ValueError("The file {} could not ".format(filename) +
                                 "be loaded, probably because it contained " +
                                 "spurious data and/or strings")
            args.append(col)

    names = [name.strip(bad_chars) for name in raw_names if name != '']

    if len(names) == 0:
        names = ['col%i' % (i + 1) for i in range(len(args))]

    return names, args


def get_column_data(*args):
    """
    get_column_data( *NumPy_args )
    """
    if len(args) == 0:
        raise IOErr('noarrays')

    cols = []
    for arg in args:
        if arg is None or isinstance(arg, (numpy.ndarray, list, tuple)):
            vals = arg
        else:
            raise IOErr('badarray', arg)

        if arg is not None:
            vals = numpy.asanyarray(vals)
            for col in numpy.atleast_2d(vals.T):
                cols.append(col)
        else:
            cols.append(vals)

    return cols


def get_ascii_data(filename, ncols=1, colkeys=None, sep=' ', dstype=Data1D,
                   comment='#', require_floats=True):
    """Read in columns from an ASCII file.

    Parameters
    ----------
    filename : str
       The name of the ASCII file to read in.
    ncols : int, optional
       The number of columns to read in (the first ``ncols`` columns
       in the file). This is ignored if ``colkeys`` is given.
    colkeys : array of str, optional
       An array of the column name to read in. The default is
       ``None``.
    sep : str, optional
       The separator character. The default is ``' '``.
    dstype : data class to use, optional
       Used to check that the data file contains enough columns.
    comment : str, optional
       The comment character. The default is ``'#'``.
    require_floats : bool, optional
       If ``True`` (the default), non-numeric data values will
       raise a `ValueError`.

    Returns
    -------
    (colnames, coldata, filename)
       The column names read in, the data for the columns
       as an array, with each element being the data for the column
       (the order matches ``colnames``), and the name of the file.

    Raises
    ------
    sherpa.utils.IOErr
       Raised if a requested column is missing or the file appears
       to be a binary file.
    ValueError
       If a column value can not be converted into a numeric value
       and the ``require_floats`` parameter is True.

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
    If the ``require_floats`` argument is ``True`` then the
    column will be converted to the `sherpa.utils.SherpaFloat`
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

    names, args = read_file_data(filename, sep, comment, require_floats)

    if colkeys is None:
        kwargs = []
        if ncols != 1:
            _check_args(ncols, dstype)
        kwargs.extend(args[:ncols])
        return (names, kwargs, filename)

    kwargs = []
    colkeys = list(colkeys)

    if len(names) > len(args):
        raise IOErr('toomanycols')

    assert(len(names) <= len(args))

    for key in colkeys:
        if key not in names:
            raise IOErr('reqcol', key, numpy.asarray(names, numpy.string_))
        kwargs.append(args[names.index(key)])

    _check_args(len(kwargs), dstype)
    return (colkeys, kwargs, filename)


def read_data(filename, ncols=2, colkeys=None, sep=' ', dstype=Data1D,
              comment='#', require_floats=True):
    """Create a data object from an ASCII file.

    Parameters
    ----------
    filename : str
       The name of the ASCII file to read in.
    ncols : int, optional
       The number of columns to read in (the first ``ncols`` columns
       in the file). This is ignored if ``colkeys`` is given.
    colkeys : array of str, optional
       An array of the column name to read in. The default is
       ``None``.
    sep : str, optional
       The separator character. The default is ``' '``.
    dstype : data class to use, optional
       The class of the data object to create.
    comment : str, optional
       The comment character. The default is ``'#'``.
    require_floats : bool, optional
       If ``True`` (the default), non-numeric data values will
       raise a `ValueError`.

    Returns
    -------
    data
       The data object (created by calling the dstype constructor
       with the filename and then the data columns from the file).

    Raises
    ------
    sherpa.utils.IOErr
       Raised if a requested column is missing or the file appears
       to be a binary file.
    ValueError
       If a column value can not be converted into a numeric value
       and the ``require_floats`` parameter is True.

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

    colnames, args, name = get_ascii_data(filename, ncols, colkeys,
                                          sep, dstype, comment, require_floats)
    return dstype(name, *args)


def read_arrays(*args):
    """Create a data object from arrays.

    Parameters
    ----------
    col1, ... coln : array_like
       The data columns.
    dstype : optional, default=`sherpa.data.Data1D`
       The data type to create. It must be a subclass of
       `sherpa.data.BaseData`.

    Returns
    -------
    data
       The data object (created by calling the dstype constructor
       with the filename and then the data columns from the file).

    Raises
    ------
    sherpa.utils.IOErr
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
    args = list(args)
    if len(args) == 0:
        raise IOErr('noarrays')

    dstype = Data1D
    if _is_subclass(args[-1], BaseData):
        dstype = args.pop()

    args = get_column_data(*args)

    # Determine max number of args for dataset constructor
    _check_args(len(args), dstype)

    return dstype('', *args)


def write_arrays(filename, args, fields=None, sep=' ', comment='#',
                 clobber=False, linebreak='\n', format='%g'):
    """Write a list of arrays to an ASCII file.

    Parameters
    ----------
    filename : str
       The name of the file to write the array to.
    args : array_like
       The arrays to write out.
    fields : array_like of str
       The column names (should match the size of ``args`` if given).
    sep : str, optional
       The separator character. The default is ``' '``.
    comment : str, optional
       The comment character. The default is ``'#'``. This is only used
       to write out the column names when ``fields`` is not None.
    clobber : bool, optional
       If ``filename`` is not ``None``, then this flag controls
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
       If ``filename`` already exists and ``clobber`` is ``False``
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

    if not numpy.iterable(args) or len(args) == 0:
        raise IOErr('noarrayswrite')

    if not numpy.iterable(args[0]):
        raise IOErr('noarrayswrite')

    size = len(args[0])
    for arg in args:
        if not numpy.iterable(arg):
            raise IOErr('noarrayswrite')
        elif len(arg) != size:
            raise IOErr('arraysnoteq')

    args = numpy.column_stack(numpy.asarray(args))

    f = open(filename, 'w')

    if fields is not None:
        f.write(comment + sep.join(fields) + linebreak)

    lines = []
    for arg in args:
        line = [format % elem for elem in arg]
        lines.append(sep.join(line))

    f.write(linebreak.join(lines))

    # add a newline at end
    f.write(linebreak)
    f.close()


def write_data(filename, dataset, fields=None, sep=' ', comment='#',
               clobber=False, linebreak='\n', format='%g'):
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
       If ``filename`` is not ``None``, then this flag controls
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
       If ``filename`` already exists and ``clobber`` is ``False``
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

    write_arrays(filename, cols, col_names, sep, comment, clobber,
                 linebreak, format)
