#
#  Copyright (C) 2024
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
"""Defines the routines a backend must provide.

Since this uses the ABCMeta metaclass, it does not provide a usable
backend. For this, one of the *_backend modules must be used.

The actual backend-specific code should

- be placed in a file called xxx_backend.py (where xxx is the label
  that is used in the options.io_pkg setting by the user, such as
  crates),

- and this module should contain a class called Backend, which derives
  from BaseBackend, over-rides all the abstract methods of the class,
  and has a name attribute set to the string matching xxx.

This approach is simpler than that used for the plotting backends,
which uses a registry for the backends, and does not require a fixed
class name.

"""

from abc import ABCMeta, abstractmethod
import os
from typing import Any, Mapping, Optional, Sequence, Union

import numpy as np

from ...utils.err import IOErr
from ..data import Data1D

from .xstable import TableHDU


__all__ = ('BaseBackend', )


# Some variants are named xxx and xxxArg, where the former is the
# return value (an invariant type, like dict) and the latter is
# covariant (such as Mapping) as it's used as an argument to a
# function.
#
KeyType = Union[str, bool, int, float]
NamesType = Sequence[str]
HdrTypeArg = Mapping[str, KeyType]
HdrType = dict[str, KeyType]
ColumnsType = Mapping[str, Union[np.ndarray, list, tuple]]

# It's hard to type the arguments to the Data constructors
DataTypeArg = Mapping[str, Any]
DataType = dict[str, Any]


class BaseBackend(metaclass=ABCMeta):
    """Define the backend API."""

    name: str = "base"
    """The name of the backend."""

    def check_clobber(self, filename: str, clobber: bool) -> None:
        """Error out if the file exists and clobber is not set."""

        if clobber or not os.path.isfile(filename):
            return

        raise IOErr("filefound", filename)

    @abstractmethod
    def get_table_data(self, arg,
                       ncols: int = 1,
                       colkeys: Optional[NamesType] = None,
                       make_copy: bool = True,
                       fix_type: bool = True,
                       blockname: Optional[str] = None,
                       hdrkeys: Optional[NamesType] = None
                       ) -> tuple[list[str], list[np.ndarray], str, HdrType]:
        """Read columns from a file or object.

        The columns to select depend on the ncols and colkeys arguments,
        as well as the backend.

        Parameters
        ----------
        arg
            The data to read the columns from. This depends on the
            backend but is expected to be a file name or a tabular
            data structure supported by the backend.
        ncols: int, optional
            The number of columns to read in when colkeys is not
            set (the first ncols columns are chosen).
        colkeys: sequence of str or None, optional
            If given, what columns from the table should be selected,
            otherwise the backend selects. The default is `None`.
            Names are compared using a case insensitive match.
        make_copy: bool, optional
            If set then the returned NumPy arrays are expictly copied,
            rather than using a reference from the data structure
            created by the backend. Backends are not required to
            honor this setting. The default is `True`.
        fix_type: bool, optional
            Should the returned arrays be converted to the
            `sherpa.utils.numeric_types.SherpaFloat` type. The default
            is `True`.
        blockname: str or None, optional
            The name of the "block" (HDU) to read the column data (useful
            for data structures containing multiple blocks/HDUs) or None.
            Names are compared using a case insensitive match.
        hdrkeys: sequence of str or None, optional
            If set, the table structure must contain these keys, and the
            values are returned.  Names are compared using a case
            insensitive match.

        Returns
        -------
        names, data, filename, hdr
            The column names, as a list of strings, and the data as
            a list of NumPy arrays (matching the order and length of
            the names array). The filename is the name of the file (a
            string) and hdr is a dictionary with the requested keywords
            (when hdrkeys is `None` this dictionary will be empty).

        Raises
        ------
        sherpa.utils.err.IOErr
            The arg argument is invalid, or a required column or keyword
            is missing.

        """

    @abstractmethod
    def get_header_data(self, arg,
                        blockname: Optional[str] = None,
                        hdrkeys: Optional[NamesType] = None
                        ) -> HdrType:
        """Read the metadata.

        Parameters
        ----------
        arg
            The data to read the header values from. This depends on the
            backend but is expected to be a file name or a data structure
            supported by the backend.
        blockname: str or None, optional
            The name of the "block" (HDU) to read the data (useful
            for data structures containing multiple blocks/HDUs) or None.
            Names are compared using a case insensitive match.
        hdrkeys: sequence of str or None, optional
            If set, the table structure must contain these keys, and the
            values are returned. Names are compared using a case
            insensitive match.

        Returns
        -------
        hdr: dict
            A dictionary with the keyword data (only the name and values
            are returned, any sort of metadata, such as comments or units,
            are not returned).

        Raises
        ------
        sherpa.utils.err.IOErr
            The arg argument is invalid or a keyword is missing.

        """

    @abstractmethod
    def get_image_data(self, arg,
                       make_copy: bool = True,
                       fix_type: bool = True
                       ) -> tuple[DataType, str]:
        """Read image data.

        Parameters
        ----------
        arg
            The data to read the header values from. This depends on the
            backend but is expected to be a file name or a data structure
            supported by the backend.
        make_copy: bool, optional
            If set then the returned NumPy arrays are expictly copied,
            rather than using a reference from the data structure
            created by the backend. Backends are not required to
            honor this setting. The default is `True`.
        fix_type: bool, optional
            Should the returned arrays be converted to the
            `sherpa.utils.numeric_types.SherpaFloat` type. The default
            is `True`.

        Returns
        -------
        data, filename
            The data, as a dictionary, and the filename. The keys of the
            dictionary match the arguments when creating a
            sherpa.astro.data.DataIMG object.

        Raises
        ------
        sherpa.utils.err.IOErr
            The arg argument is invalid or not an image.

        """

    @abstractmethod
    def get_column_data(self, *args) -> list[np.ndarray]:
        """Extract the column data.

        Parameters
        ----------
        *args
            Extract column information from each argument. It can be an
            ndarray, list, or tuple, or a data structure from the backend,
            with each argument representing a column. 2D arguments are
            separated by column.

        Returns
        -------
        data: list of ndarray
            The column data.

        Raises
        ------
        sherpa.utils.err.IOErr
            There are no arguments or an argument is not supported.

        Notes
        -----

        An argument can be None, which is just passed back to the caller.
        This means the typing rules for the function are not quite
        correct.

        """

    # Follow sherpa.io.get_ascii_data API.
    #
    @abstractmethod
    def get_ascii_data(self,
                       filename: str,
                       ncols: int = 1,
                       colkeys: Optional[NamesType] = None,
                       sep: str = ' ',
                       dstype: type = Data1D,
                       comment: str = '#',
                       require_floats: bool = True
                       ) -> tuple[list[str], list[np.ndarray], str]:
        """Read columns from an ASCII file.

        The `sep`, `dstype`, `comment`, and `require_floats` arguments
        may be ignored by the backend.q

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
        colnames, coldata, filename
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


        """

    @abstractmethod
    def get_arf_data(self, arg,
                     make_copy: bool = True
                     ) -> tuple[DataType, str]:
        """Read in the ARF.

        Parameters
        ----------
        arg
            The data to read the ARF from. This depends on the backend but
            is expected to be a file name or a data structure supported by
            the backend.
        make_copy: bool, optional
            If set then the returned NumPy arrays are expictly copied,
            rather than using a reference from the data structure
            created by the backend. Backends are not required to
            honor this setting. The default is `True`.

        Returns
        -------
        data, filename
            The data, as a dictionary, and the filename. The keys of the
            dictionary match the arguments when creating a
            sherpa.astro.data.DataARF object.

        Raises
        ------
        sherpa.utils.err.IOErr
            The arg argument is invalid or not an ARF.

        """

    @abstractmethod
    def get_rmf_data(self, arg,
                     make_copy: bool = True
                     ) -> tuple[DataType, str]:
        """Read in the RMF.

        Parameters
        ----------
        arg
            The data to read the RMF from. This depends on the backend but
            is expected to be a file name or a data structure supported by
            the backend.
        make_copy: bool, optional
            If set then the returned NumPy arrays are expictly copied,
            rather than using a reference from the data structure
            created by the backend. Backends are not required to
            honor this setting. The default is `True`.

        Returns
        -------
        data, filename
            The data, as a dictionary, and the filename. The keys of the
            dictionary match the arguments when creating a
            sherpa.astro.data.DataRMF object.

        Raises
        ------
        sherpa.utils.err.IOErr
            The arg argument is invalid or not a RMF.

        """

    @abstractmethod
    def get_pha_data(self, arg,
                     make_copy: bool = True,
                     use_background: bool = False
                     ) -> tuple[list[DataType], str]:
        """Read in the PHA.

        Parameters
        ----------
        arg
            The data to read the PHA from. This depends on the backend but
            is expected to be a file name or a data structure supported by
            the backend.
        make_copy: bool, optional
            If set then the returned NumPy arrays are expictly copied,
            rather than using a reference from the data structure
            created by the backend. Backends are not required to
            honor this setting. The default is `True`.
        use_background: bool, optional
            Should the data be read in as a background file (only relevant
            for files that contain the background data in a separate block
            of the same file, such as Chandra Level 3 PHA files, as used
            by the Chandra Source Catalog). The default is `False`.

        Returns
        -------
        datas, filename
            A list of dictionaries, containing the PHA data (since there
            can be multiple datasets with a PHA-II file) and the filename.
            The keys of the dictionary match the arguments when creating
            a sherpa.astro.data.DataPHA object.

        Raises
        ------
        sherpa.utils.err.IOErr
            The arg argument is invalid or not a PHA.

        """

    @abstractmethod
    def pack_table_data(self,
                        data: ColumnsType,
                        col_names: NamesType,
                        header: Optional[HdrTypeArg] = None) -> Any:
        """Create the tabular data.

        Parameters
        ----------
        data : dict
            The table data, where the key is the column name and the value
            the data.
        col_names : sequence of str
            The column names from data to use (this also sets the order).
        header : dict or None, optional
            Any header information to include.

        Returns
        -------
        table
            A data structure used by the backend to represent tabular data.

        """

    @abstractmethod
    def pack_image_data(self,
                        data: DataTypeArg,
                        header: HdrTypeArg) -> Any:
        """Create the image data.

        Parameters
        ----------
        data : dict
            The image data, where the keys are arguments used to create a
            sherpa.astro.data.DataIMG object.
        header : dict
            The header information to include.

        Returns
        -------
        image
            A data structure used by the backend to represent image data.

        """

    @abstractmethod
    def pack_pha_data(self,
                      data: ColumnsType,
                      col_names: NamesType,
                      header: Optional[HdrTypeArg] = None) -> Any:
        """Create the PHA.

        Parameters
        ----------
        data : dict
            The table data, where the key is the column name and the value
            the data.
        col_names : sequence of str
            The column names from data to use (this also sets the order).
        header : dict or None, optional
            Any header information to include.

        Returns
        -------
        pha
            A data structure used by the backend to represent PHA data.

        """

    @abstractmethod
    def pack_arf_data(self,
                      data: ColumnsType,
                      col_names: NamesType,
                      header: Optional[HdrTypeArg] = None) -> Any:
        """Create the ARF.

        Parameters
        ----------
        data : dict
            The table data, where the key is the column name and the value
            the data.
        col_names : sequence of str
            The column names from data to use (this also sets the order).
        header : dict or None, optional
            Any header information to include.

        Returns
        -------
        arf
            A data structure used by the backend to represent ARF data.

        """

    @abstractmethod
    def pack_rmf_data(self, blocks) -> Any:
        """Create the RMF.

        Parameters
        ----------
        blocks : sequence of pairs
            The RMF data, stored as pairs of (data, header), where data is
            a dictionary of column name (keys) and values, and header is a
            dictionary of key and values. The first element is the MATRIX
            block and the second is for the EBOUNDS block.

        Returns
        -------
        rmf
            A data structure used by the backend to represent RMF data.

        """

    @abstractmethod
    def pack_hdus(self,
                  blocks: Sequence[TableHDU]) -> Any:
        """Create a dataset.

        Parameters
        ----------
        blocks : sequence of TableHDU
            The blocks (HDUs) to store.

        Returns
        -------
        hdus
            A data structure used by the backend to represent the data.

        """

    @abstractmethod
    def set_table_data(self,
                       filename: str,
                       data: ColumnsType,
                       col_names: NamesType,
                       header: Optional[HdrTypeArg] = None,
                       ascii: bool = False,
                       clobber: bool = False) -> None:
        """Write out the tabular data.

        Parameters
        ----------
        filename : str
            The name of the file to create.
        data : dict
            The table data, where the key is the column name and the value
            the data.
        col_names : sequence of str
            The column names from data to use (this also sets the order).
        header : dict or None, optional
            Any header information to include.
        ascii : bool, optional
            Is the file to be written out as a text file (`True`) or a
            binary file? The default is `False`.
        clobber : bool, optional
            If the file already exists can it be over-written (`True`) or
            will a sherpa.utils.err.IOErr error be raised? The default is
            `False`.

        """

    @abstractmethod
    def set_image_data(self,
                       filename: str,
                       data: DataTypeArg,
                       header: HdrTypeArg,
                       ascii: bool = False,
                       clobber: bool = False) -> None:
        """Write out the image data.

        Parameters
        ----------
        filename : str
            The name of the file to create.
        data : dict
            The image data, where the keys are arguments used to create a
            sherpa.astro.data.DataIMG object.
        header : dict
            The header information to include.
        ascii : bool, optional
            Is the file to be written out as a text file (`True`) or a
            binary file? The default is `False`.
        clobber : bool, optional
            If the file already exists can it be over-written (`True`) or
            will a sherpa.utils.err.IOErr error be raised? The default is
            `False`.

        """

    @abstractmethod
    def set_pha_data(self,
                     filename: str,
                     data: ColumnsType,
                     col_names: NamesType,
                     header: Optional[HdrTypeArg] = None,
                     ascii: bool = False,
                     clobber: bool = False) -> None:
        """Write out the PHA.

        Parameters
        ----------
        filename : str
            The name of the file to create.
        data : dict
            The table data, where the key is the column name and the value
            the data.
        col_names : sequence of str
            The column names from data to use (this also sets the order).
        header : dict or None, optional
            Any header information to include.
        ascii : bool, optional
            Is the file to be written out as a text file (`True`) or a
            binary file? The default is `False`.
        clobber : bool, optional
            If the file already exists can it be over-written (`True`) or
            will a sherpa.utils.err.IOErr error be raised? The default is
            `False`.

        """

    @abstractmethod
    def set_arf_data(self,
                     filename: str,
                     data: ColumnsType,
                     col_names: NamesType,
                     header: Optional[HdrTypeArg] = None,
                     ascii: bool = False,
                     clobber: bool = False) -> None:
        """Write out the ARF.

        Parameters
        ----------
        filename : str
            The name of the file to create.
        data : dict
            The table data, where the key is the column name and the value
            the data.
        col_names : sequence of str
            The column names from data to use (this also sets the order).
        header : dict or None, optional
            Any header information to include.
        ascii : bool, optional
            Is the file to be written out as a text file (`True`) or a
            binary file? The default is `False`.
        clobber : bool, optional
            If the file already exists can it be over-written (`True`) or
            will a sherpa.utils.err.IOErr error be raised? The default is
            `False`.

        """

    @abstractmethod
    def set_rmf_data(self,
                     filename: str,
                     blocks,
                     clobber: bool = False) -> None:
        """Write out the RMF.

        Parameters
        ----------
        filename : str
            The name of the file to create.
        blocks : sequence of pairs
            The RMF data, stored as pairs of (data, header), where data is
            a dictionary of column name (keys) and values, and header is a
            dictionary of key and values. The first element is the MATRIX
            block and the second is for the EBOUNDS block.
        clobber : bool, optional
            If the file already exists can it be over-written (`True`) or
            will a sherpa.utils.err.IOErr error be raised? The default is
            `False`.

        Notes
        -----
        There is currently no support for writing out a RMF as an ASCII
        file.

        """

    @abstractmethod
    def set_hdus(self,
                 filename: str,
                 blocks: Sequence[TableHDU],
                 clobber: bool = False) -> None:
        """Write out (possibly multiple) blocks.

        Parameters
        ----------
        filename : str
            The name of the file to create.
        blocks : sequence of TableHDU
            The blocks (HDUs) to store.
        clobber : bool, optional
            If the file already exists can it be over-written (`True`) or
            will a sherpa.utils.err.IOErr error be raised? The default is
            `False`.

        """

    @abstractmethod
    def read_table_blocks(self, arg,
                          make_copy: bool = False
                          ) -> tuple[str,
                                     dict[int, dict[str, np.ndarray]],
                                     dict[int, HdrType]]:
        """Read in tabular data with no restrictions on the columns.

        Parameters
        ----------
        arg
            The data to read the columns from. This depends on the
            backend but is expected to be a file name or a tabular
            data structure supported by the backend.
        make_copy: bool, optional
            If set then the returned NumPy arrays are expictly copied,
            rather than using a reference from the data structure
            created by the backend. Backends are not required to
            honor this setting. The default is `True`.

        Returns
        -------
        filename, blockdata, hdrdata
            The filename as a string. The blockdata and hdrdata values are
            dictionaries where the key is an integer representing the
            block (or HDU) number (where the first block is numbered 1 and
            represents the first tabular block, that is it does not
            include the primary HDU) and the values are dictionaries
            representing the column data or header data for each block.

        Raises
        ------
        sherpa.utils.err.IOErr
            The arg argument is invalid, or a required column or keyword
            is missing.

        """

    def set_arrays(self,
                   filename: str,
                   args: Sequence[np.ndarray],
                   fields: Optional[NamesType] = None,
                   ascii: bool = True,
                   clobber: bool = False) -> None:
        """Write out columns.

        Parameters
        ----------
        filename : str
            The file name.
        args : sequence of ndarray
            The column data.
        fields : sequence of str or None, optional
            The column names to use. If set to `None` then the columns are
            named ``col1``, ``col2``, ...
        ascii : bool, optional
            Is the file to be written out as a text file (`True`) or a
            binary file? The default is `True`.
        clobber : bool, optional
            If the file already exists can it be over-written (`True`) or
            will a sherpa.utils.err.IOErr error be raised? The default is
            `False`.

        """

        self.check_clobber(filename, clobber)

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

        # This is the backend-specific functionality.
        #
        self.write_arrays(filename, args, fieldnames, ascii=ascii)

    @abstractmethod
    def write_arrays(self,
                     filename: str,
                     args: Sequence[np.ndarray],
                     fields: Optional[NamesType] = None,
                     ascii: bool = True) -> None:
        """Write out columns.

        The clobber check has been made in save_arrays.

        Parameters
        ----------
        filename : str
            The file name.
        args : sequence of ndarray
            The column data. This has been checked to have the same
            size for each column.
        fields : sequence of str
            The column names to use. This matches the number of columns
            in args.
        ascii : bool, optional
            Is the file to be written out as a text file (`True`) or a
            binary file? The default is `True`.

        """
