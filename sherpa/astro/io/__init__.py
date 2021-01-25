#
#  Copyright (C) 2007, 2015, 2016, 2017, 2018, 2019, 2021
#     Smithsonian Astrophysical Observatory
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
Provide Astronomy-specific I/O routines for Sherpa.

This module contains read and write routines for handling FITS [1]_
and ASCII format data. The actual support is provided by
the selected I/O backend package (currently Crates, provided
by CIAO, or the FITS support in the AstroPy package, if installed).

This module should not be imported directly if there is no
I/O backend available.

References
----------

.. [1] Flexible Image Transport System, https://en.wikipedia.org/wiki/FITS

"""

from configparser import ConfigParser
import logging
import os
import os.path
import sys
import importlib

import numpy

import sherpa.io
from sherpa.utils.err import IOErr
from sherpa.utils import SherpaFloat
from sherpa.data import Data2D, Data1D, BaseData, Data2DInt
from sherpa.astro.data import DataIMG, DataIMGInt, DataARF, DataRMF, DataPHA, DataRosatRMF
from sherpa.astro.utils import reshape_2d_arrays
from sherpa import get_config

config = ConfigParser()
config.read(get_config())

# What should we use for the default package?
io_opt = config.get('options', 'io_pkg', fallback='pyfits')
io_opt = str(io_opt).strip().lower()

ogip_emin = config.get('ogip', 'minimum_energy', fallback='1.0e-10')

if ogip_emin.upper() == 'NONE':
    ogip_emin = None
else:
    emsg = "Invalid value for [ogip] minimum_energy config value; " + \
           "it must be None or a float > 0"
    try:
        ogip_emin = float(ogip_emin)
    except ValueError:
        raise ValueError(emsg)

    if ogip_emin <= 0.0:
        raise ValueError(emsg)

if io_opt.startswith('pycrates') or io_opt.startswith('crates'):
    io_opt = 'crates_backend'

elif io_opt.startswith('pyfits'):
    io_opt = 'pyfits_backend'

try:
    importlib.import_module('.' + io_opt, package='sherpa.astro.io')
    backend = sys.modules['sherpa.astro.io.' + io_opt]
except ImportError:
    raise ImportError("""Cannot import selected FITS I/O backend {}.
    If you are using CIAO, this is most likely an error and you should contact the CIAO helpdesk.
    If you are using Standalone Sherpa, please install astropy."""
                      .format(io_opt))

warning = logging.getLogger(__name__).warning
info = logging.getLogger(__name__).info


__all__ = ('read_table', 'read_image', 'read_arf', 'read_rmf', 'read_arrays',
           'read_pha', 'write_image', 'write_pha', 'write_table',
           'pack_table', 'pack_image', 'pack_pha', 'read_table_blocks')


def _is_subclass(t1, t2):
    return isinstance(t1, type) and issubclass(t1, t2) and (t1 is not t2)


# Note: write_arrays is not included in __all__, so don't add to the
#       See Also section.
#
def read_arrays(*args):
    """Create a dataset from multiple arrays.

    The return value defaults to a `sherpa.data.Data1D` instance,
    but this can be changed by supplying the required class
    as the last argument (anything that is derived from
    `sherpa.data.BaseData`).

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
    data : a sherpa.data.BaseData derived object

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
    args = list(args)
    if len(args) == 0:
        raise IOErr('noarrays')

    dstype = Data1D
    if _is_subclass(args[-1], BaseData):
        dstype = args.pop()

    args = backend.get_column_data(*args)

    # Determine max number of args for dataset constructor
    sherpa.io._check_args(len(args), dstype)

    return dstype('', *args)


def read_table(arg, ncols=2, colkeys=None, dstype=Data1D):
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
        `sherpa.data.BaseData` interface).

    Returns
    -------
    data : a sherpa.data.BaseData derived object

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

    Create A `sherpa.data.Data1DInt` object from the first three
    columns in the file ``src.fits``:

    >>> d = read_table('src.fits', ncols=3, dstype=Data1DInt)

    Create a `sherpa.data.Data1D` data set from the specified
    columns in ``tbl.fits``, where ``WLEN`` is used for the
    independent axis, ``FLUX`` the dependent axis, and
    ``FLUXERR`` for the statistical error on the dependent axis:

    >>> d = read_table('tbl.fits', colkeys=['WLEN', 'FLUX', 'FLUXERR'])

    """
    colnames, cols, name, hdr = backend.get_table_data(arg, ncols, colkeys)

    # Determine max number of args for dataset constructor
    sherpa.io._check_args(len(cols), dstype)

    return dstype(name, *cols)


def read_ascii(filename, ncols=2, colkeys=None, dstype=Data1D, **kwargs):
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
        `sherpa.data.BaseData` interface).
    **kwargs
        The remaining arguments are passed through to the
        ``get_ascii_data`` routine of the I/O backend. It is
        expected that these include ``sep`` and ``comment``,
        which describe the column separators and indicator of
        a comment line, respectively.

    Returns
    -------
    data : a sherpa.data.BaseData derived object

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
    colnames, cols, name = backend.get_ascii_data(filename, ncols=ncols,
                                                  colkeys=colkeys,
                                                  dstype=dstype, **kwargs)

    # Determine max number of args for dataset constructor
    sherpa.io._check_args(len(cols), dstype)

    return dstype(name, *cols)


def read_image(arg, coord='logical', dstype=DataIMG):
    """Create an image dataset from a file.

    Parameters
    ----------
    arg
        The name of the file or a representation of the file
        (the type depends on the I/O backend).
    coord : {'logical', 'physical', 'world'}, optional
        The type of axis coordinates to use. An error is raised if the
        file does not contain the necessary metadata for the selected
        coordinate system.
    dstype : optional
        The data type to create (it is expected to follow the
        `sherpa.data.BaseData` interface).

    Returns
    -------
    data : a sherpa.data.BaseData derived object

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
    data, filename = backend.get_image_data(arg)
    axlens = data['y'].shape

    x0 = numpy.arange(axlens[1], dtype=SherpaFloat) + 1.
    x1 = numpy.arange(axlens[0], dtype=SherpaFloat) + 1.
    x0, x1 = reshape_2d_arrays(x0, x1)

    data['y'] = data['y'].ravel()
    data['coord'] = coord
    data['shape'] = axlens

    if issubclass(dstype, DataIMGInt):
        dataset = dstype(filename, x0 - 0.5, x1 - 0.5, x0 + 0.5, x1 + 0.5,
                         **data)
    elif issubclass(dstype, Data2DInt):
        for name in ['coord', 'eqpos', 'sky', 'header']:
            data.pop(name, None)
        dataset = dstype(filename, x0 - 0.5, x1 - 0.5, x0 + 0.5, x1 + 0.5,
                         **data)
    else:
        dataset = dstype(filename, x0, x1, **data)

    return dataset


def read_arf(arg):
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
    data, filename = backend.get_arf_data(arg)

    # It is unlikely that the backend will set this, but allow
    # it to override the config setting.
    #
    if 'emin' not in data:
        data['ethresh'] = ogip_emin

    return DataARF(filename, **data)


def read_rmf(arg):
    """Create a DataRMF object.

    Parameters
    ----------
    arg
        The name of the file or a representation of the file
        (the type depends on the I/O backend) containing the
        RMF data.

    Returns
    -------
    data : sherpa.astro.data.DataRMF

    """
    data, filename = backend.get_rmf_data(arg)

    # It is unlikely that the backend will set this, but allow
    # it to override the config setting.
    #
    if 'emin' not in data:
        data['ethresh'] = ogip_emin

    return _rmf_factory(filename, data)


def _rmf_factory(filename, data):
    response_map = {
        'ROSAT': DataRosatRMF,
        'DEFAULT': DataRMF,
    }

    rmf_telescop = data['header'].get('TELESCOP', 'DEFAULT')
    rmf_class = response_map.get(rmf_telescop, DataRMF)

    return rmf_class(filename, **data)


def _read_ancillary(data, key, label, dname,
                    read_func, output_once=True):
    """Read in a file if the keyword is set.

    Parameters
    ----------
    data
        The header information, which behaves like a dictionary.
    key
        The key to look for in data. If not set, or its value compares
        case insensitively to "none" then nothing is read in.
    label : str
        This is used to identify the file being read in. It is only
        used if output_once is set.
    dname : str
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
    data : None or a sherpa.data.BaseData derived object

    """

    if not(data[key]) or data[key].lower() == 'none':
        return None

    out = None
    try:
        if os.path.dirname(data[key]) == '':
            data[key] = os.path.join(dname, data[key])

        out = read_func(data[key])
        if output_once:
            info('read {} file {}'.format(label, data[key]))

    except:
        if output_once:
            warning(str(sys.exc_info()[1]))

    return out


def read_pha(arg, use_errors=False, use_background=False):
    """Create a DataPHA object.

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
    datasets, filename = backend.get_pha_data(arg,
                                              use_background=use_background)
    phasets = []
    output_once = True
    for data in datasets:
        if not use_errors:
            if data['staterror'] is not None or data['syserror'] is not None:
                if data['staterror'] is None:
                    msg = 'systematic'
                elif data['syserror'] is None:
                    msg = 'statistical'
                    if output_once:
                        wmsg = "systematic errors were not found in " + \
                               "file '{}'".format(filename)
                        warning(wmsg)
                else:
                    msg = 'statistical and systematic'
                if output_once:
                    imsg = msg + " errors were found in file " + \
                           "'{}' \nbut not used; ".format(filename) + \
                           "to use them, re-read with use_errors=True"
                    info(imsg)
                data['staterror'] = None
                data['syserror'] = None

        dname = os.path.dirname(filename)
        albl = 'ARF'
        rlbl = 'RMF'
        if use_background:
            albl = albl + ' (background)'
            rlbl = rlbl + ' (background)'

        arf = _read_ancillary(data, 'arffile', albl, dname, read_arf,
                              output_once)
        rmf = _read_ancillary(data, 'rmffile', rlbl, dname, read_rmf,
                              output_once)

        backgrounds = []

        if data['backfile'] and data['backfile'].lower() != 'none':
            try:
                if os.path.dirname(data['backfile']) == '':
                    data['backfile'] = os.path.join(os.path.dirname(filename),
                                                    data['backfile'])

                bkg_datasets = []
                # Do not read backgrounds of backgrounds
                if not use_background:
                    bkg_datasets = read_pha(data['backfile'], use_errors, True)

                    if output_once:
                        info('read background file {}'.format(
                            data['backfile']))

                if numpy.iterable(bkg_datasets):
                    for bkg_dataset in bkg_datasets:
                        if bkg_dataset.get_response() == (None, None) and \
                           rmf is not None:
                            bkg_dataset.set_response(arf, rmf)
                        backgrounds.append(bkg_dataset)
                else:
                    if bkg_datasets.get_response() == (None, None) and \
                       rmf is not None:
                        bkg_datasets.set_response(arf, rmf)
                    backgrounds.append(bkg_datasets)

            except:
                if output_once:
                    warning(str(sys.exc_info()[1]))

        for bkg_type, bscal_type in zip(('background_up', 'background_down'),
                                        ('backscup', 'backscdn')):
            if data[bkg_type] is not None:
                b = DataPHA(filename,
                            channel=data['channel'],
                            counts=data[bkg_type],
                            bin_lo=data['bin_lo'],
                            bin_hi=data['bin_hi'],
                            grouping=data['grouping'],
                            quality=data['quality'],
                            exposure=data['exposure'],
                            backscal=data[bscal_type],
                            header=data['header'])
                b.set_response(arf, rmf)
                if output_once:
                    info("read {} into a dataset from file {}".format(
                        bkg_type, filename))
                backgrounds.append(b)

        for k in ['backfile', 'arffile', 'rmffile', 'backscup', 'backscdn',
                  'background_up', 'background_down']:
            data.pop(k, None)

        pha = DataPHA(filename, **data)
        pha.set_response(arf, rmf)
        for id, b in enumerate(backgrounds):
            if b.grouping is None:
                b.grouping = pha.grouping
                b.grouped = (b.grouping is not None)
            if b.quality is None:
                b.quality = pha.quality
            pha.set_background(b, id + 1)

        # set units *after* bkgs have been set
        pha._set_initial_quantity()
        phasets.append(pha)
        output_once = False

    if len(phasets) == 1:
        phasets = phasets[0]

    return phasets


def _pack_table(dataset):
    names = dataset._fields
    cols = [(name, getattr(dataset, name)) for name in names]
    data = dict(cols)
    return data, names


def _pack_image(dataset):
    if not(isinstance(dataset, Data2D) or
           isinstance(dataset, DataIMG)):
        raise IOErr('notimage', dataset.name)
    data = {}

    header = {}
    if hasattr(dataset, 'header') and type(dataset.header) is dict:
        header = dataset.header.copy()

    data['pixels'] = numpy.asarray(dataset.get_img())
    data['sky'] = getattr(dataset, 'sky', None)
    data['eqpos'] = getattr(dataset, 'eqpos', None)

    return data, header


def _set_keyword(header, label, value):
    """Extract the string form of the value, and store
    the last path element in the header using the label
    key.
    """

    if value is None:
        return

    name = getattr(value, 'name', 'none')
    if name is not None and name.find('/') != -1:
        name = name.split('/')[-1]

    header[label] = name


def _pack_pha(dataset):
    if not isinstance(dataset, DataPHA):
        raise IOErr('notpha', dataset.name)

    data = {}

    arf, rmf = dataset.get_response()
    bkg = dataset.get_background()

    # Header Keys
    header = {}
    if hasattr(dataset, 'header'):  # and type(dataset.header) is dict:
        header = dataset.header.copy()

    header['EXPOSURE'] = getattr(dataset, 'exposure', 'none')

    _set_keyword(header, 'RESPFILE', rmf)
    _set_keyword(header, 'ANCRFILE', arf)
    _set_keyword(header, 'BACKFILE', bkg)

    # Columns
    col_names = ['channel', 'counts', 'stat_err', 'sys_err',
                 'bin_lo', 'bin_hi', 'grouping', 'quality']

    data['channel'] = getattr(dataset, 'channel', None)
    data['counts'] = getattr(dataset, 'counts', None)
    data['stat_err'] = getattr(dataset, 'staterror', None)
    data['sys_err'] = getattr(dataset, 'syserror', None)
    data['bin_lo'] = getattr(dataset, 'bin_lo', None)
    data['bin_hi'] = getattr(dataset, 'bin_hi', None)
    data['grouping'] = getattr(dataset, 'grouping', None)
    data['quality'] = getattr(dataset, 'quality', None)

    backscal = getattr(dataset, 'backscal', None)
    if backscal is not None:
        if numpy.isscalar(backscal):
            header['BACKSCAL'] = backscal
        else:
            data['backscal'] = backscal
            col_names.append('backscal')

    areascal = getattr(dataset, 'areascal', None)
    if areascal is not None:
        if numpy.isscalar(areascal):
            header['AREASCAL'] = areascal
        else:
            data['areascal'] = areascal
            col_names.append('areascal')

    return data, col_names, header


def write_arrays(filename, args, fields=None, ascii=True, clobber=False):
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


def write_table(filename, dataset, ascii=True, clobber=False):
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
    data, names = _pack_table(dataset)
    backend.set_table_data(filename, data, names, ascii=ascii, clobber=clobber)


def write_image(filename, dataset, ascii=True, clobber=False):
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
    data, hdr = _pack_image(dataset)
    backend.set_image_data(filename, data, hdr, ascii=ascii, clobber=clobber)


def write_pha(filename, dataset, ascii=True, clobber=False):
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
    data, col_names, hdr = _pack_pha(dataset)
    backend.set_pha_data(filename, data, col_names, hdr, ascii=ascii,
                         clobber=clobber)


def pack_table(dataset):
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
    data, names = _pack_table(dataset)
    return backend.set_table_data('', data, names, packup=True)


def pack_image(dataset):
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
    data, hdr = _pack_image(dataset)
    return backend.set_image_data('', data, hdr, packup=True)


def pack_pha(dataset):
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
    data, col_names, hdr = _pack_pha(dataset)
    return backend.set_pha_data('', data, col_names, hdr, packup=True)


def read_table_blocks(arg, make_copy=False):
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

    return backend.read_table_blocks(arg, make_copy=make_copy)
