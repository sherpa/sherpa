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

"""Corner-case handling of the I/O code.

This is named test_sherpa_io.py as there is already a test_io.py
in the astro tree.

"""

import pytest

from sherpa.data import Data1D
from sherpa.io import get_ascii_data, get_column_data, read_arrays, \
    write_arrays, write_data
from sherpa.utils.err import IOErr


def test_get_column_data_no_data():
    """Regression test"""

    with pytest.raises(IOErr, match="^no arrays found to be loaded$"):
        get_column_data()


def test_get_column_data_not_a_sequence():
    """Regression test"""

    with pytest.raises(IOErr, match="^'1' must be a Numpy array, list, or tuple$"):
        get_column_data(1, 2, 3.4)


def test_read_arrays_no_data():
    """Regression test"""

    with pytest.raises(IOErr, match="^no arrays found to be loaded$"):
        read_arrays()


def test_write_arrays_no_data(tmp_path):
    """Regression test"""

    tfile = tmp_path / 'x.dat'
    with pytest.raises(IOErr,
                       match=r"^please supply array\(s\) to write to file$"):
        write_arrays(str(tfile), 23, fields=["x"])


def test_get_ascii_data_irregular_data(tmp_path):
    """What happens when the number of columns is not constant?"""

    tfile = tmp_path / "col.dat"
    tfile.write_text("23 1\n24 2\n25 3 4\n")

    with pytest.raises(IOErr,
                       match="^not all arrays are of equal length$"):
        get_ascii_data(str(tfile))


def test_write_arrays_irregular_data(tmp_path):
    """What happens when the number of columns is not constant?"""

    tfile = tmp_path / "col.dat"
    with pytest.raises(IOErr,
                       match="^not all arrays are of equal length$"):
        write_arrays(str(tfile), [[1, 2, 3], [4, 5], [6, 7, 8]],
                     clobber=True)


def test_get_ascii_data_not_numeric_data(tmp_path):
    """What happens when given non-numeric data?"""

    tfile = tmp_path / "col.dat"
    tfile.write_text("23 1 foo\n24 2 bar\n25 3 baz\n")
    with pytest.raises(ValueError,
                       match=" could not be loaded, probably "
                       "because it contained spurious data and/or "
                       "strings"):

        get_ascii_data(str(tfile))


def test_write_arrays_not_numeric_data(tmp_path):
    """What happens when given non-numeric data?"""

    tfile = tmp_path / "col.dat"
    with pytest.raises(TypeError,
                       match="must be real number"):
        write_arrays(str(tfile), [[23, 1, "foo"],
                                  [24, 2, "bar"],
                                  [25, 3, "baz"]],
                     clobber=True)


def test_get_ascii_data_get_one_col(tmp_path):
    """This just checks that the test data is correct"""

    tfile = tmp_path / "1col.dat"
    tfile.write_text("# col2\n23\n24\n")

    (colnames, coldata, fname) = get_ascii_data(str(tfile))
    assert colnames == ["col2"]
    assert len(coldata) == 1
    assert coldata[0] == pytest.approx([23, 24])
    assert fname.endswith("1col.dat")


def test_get_ascii_data_more_columns_than_colnames_no_colkeys_ncols1(tmp_path):
    """What happens if more columns than column names?"""

    tfile = tmp_path / "col.dat"
    tfile.write_text("# col2\n-0.2 23\n-0.4 24\n")

    # Although ncols=1 is the default, be explicit in this test
    (colnames, coldata, fname) = get_ascii_data(str(tfile), ncols=1)
    assert colnames == ["col2"]
    assert len(coldata) == 1
    assert coldata[0] == pytest.approx([-0.2, -0.4])


def test_get_ascii_data_more_columns_than_colnames_no_colkeys_ncols2(tmp_path):
    """What happens if more columns than column names?"""

    tfile = tmp_path / "col.dat"
    tfile.write_text("# col2 col7\n-0.2 23 1\n-0.4 24 0\n")

    (colnames, coldata, fname) = get_ascii_data(str(tfile), ncols=2)
    assert colnames == ["col2", "col7"]
    assert len(coldata) == 2
    assert coldata[0] == pytest.approx([-0.2, -0.4])
    assert coldata[1] == pytest.approx([23, 24])


def test_get_ascii_data_more_columns_than_colnames_colkeys_ncols1(tmp_path):
    """What happens if more columns than column names?

    This also tests a mix-match between ncols and len(colkeys).
    """

    tfile = tmp_path / "col.dat"
    tfile.write_text("# col2 col7\n-0.2 23 1\n-0.4 24 0\n")

    # Although ncols=1 is the default, be explicit in this test
    (colnames, coldata, fname) = get_ascii_data(str(tfile), ncols=1,
                                                colkeys=["col7", "col2"])
    assert colnames == ["col7", "col2"]
    assert len(coldata) == 2
    assert coldata[0] == pytest.approx([23, 24])
    assert coldata[1] == pytest.approx([-0.2, -0.4])


def test_get_ascii_data_more_columns_than_colnames_colkeys_ncols2(tmp_path):
    """What happens if more columns than column names?

    We need to set ncols=2 to pass the check for creating a
    dstype=Data1D object (for some reason, when colkeys is not set,
    then there is no check if ncols=1).

    """

    tfile = tmp_path / "col.dat"
    tfile.write_text("# col2 col7\n-0.2 23 1\n-0.4 24 0\n")

    (colnames, coldata, fname) = get_ascii_data(str(tfile), ncols=2,
                                                colkeys=["col7", "col2"])
    assert colnames == ["col7", "col2"]
    assert len(coldata) == 2
    assert coldata[0] == pytest.approx([23, 24])
    assert coldata[1] == pytest.approx([-0.2, -0.4])


def test_get_ascii_data_more_colnames_than_columns_no_colkeys_ncols1(tmp_path):
    """What happens if more column names than columns: ncols=1?"""

    tfile = tmp_path / "col.dat"
    tfile.write_text("# col2 xcol ycol\n-0.2 23\n-0.4 24\n")

    # Although ncols=1 is the default be explicit in this test
    (colnames, coldata, fname) = get_ascii_data(str(tfile), ncols=1)
    assert colnames == ["col2", "xcol", "ycol"]
    assert len(coldata) == 1
    assert coldata[0] == pytest.approx([-0.2, -0.4])


def test_get_ascii_data_more_colnames_than_columns_no_colkeys_ncols2(tmp_path):
    """What happens if more column names than columns: ncols=2?"""

    tfile = tmp_path / "col.dat"
    tfile.write_text("# col2 xcol ycol\n-0.2 23\n-0.4 24\n")

    (colnames, coldata, fname) = get_ascii_data(str(tfile), ncols=2)
    assert colnames == ["col2", "xcol", "ycol"]
    assert len(coldata) == 2
    assert coldata[0] == pytest.approx([-0.2, -0.4])
    assert coldata[1] == pytest.approx([23, 24])


def test_get_ascii_data_more_colnames_than_columns_colkeys_ncols1(tmp_path):
    """What happens if more column names than columns? colkeys/ncols=1

    There is also a check for specifying an invalid column name but
    that check is not reached.
    """

    tfile = tmp_path / "col.dat"
    tfile.write_text("# col2 xcol ycol\n-0.2 23\n-0.4 24\n")

    # Although ncols=1 is the default be explicit in this test
    with pytest.raises(IOErr,
                       match="^Expected 2 columns but found 3$"):
        get_ascii_data(str(tfile), ncols=1, colkeys=["col2", "col3"])


def test_get_ascii_data_more_colnames_than_columns_colkeys_ncols2(tmp_path):
    """What happens if more column names than columns? colkebys/ncols=2

    There is also a check for specifying an invalid column name but
    that check is not reached.
    """

    tfile = tmp_path / "col.dat"
    tfile.write_text("# col2 xcol ycol\n-0.2 23\n-0.4 24\n")

    with pytest.raises(IOErr,
                       match="^Expected 2 columns but found 3$"):
        get_ascii_data(str(tfile), ncols=2, colkeys=["col2", "col3"])


def test_get_ascii_data_too_many_colkeys(tmp_path):
    """Regression test"""

    tfile = tmp_path / "1col.dat"
    tfile.write_text("23 400\n24 -0.2\n")

    with pytest.raises(IOErr,
                       match=r"Required column 'col3' not found in \['col1'. 'col2'\]"):
        get_ascii_data(str(tfile), ncols=3, colkeys=["col1", "col2", "col3"])


def test_get_ascii_data_unknown_colkey(tmp_path):
    """Regression test"""

    tfile = tmp_path / "1col.dat"
    tfile.write_text("# col2\n23\n24\n")

    with pytest.raises(IOErr,
                       match=r"Required column 'col1' not found in \['col2'\]"):
        get_ascii_data(str(tfile), colkeys=["col1", "col2", "col3"])


def test_write_arrays_single_column(tmp_path):
    """What happens if do not give a list of lists?

    Should this error out? Users may be relying on this behavior, so
    we do not error out, but that seems unlikely.

    """

    tfile = tmp_path / "not-changed.dat"
    with pytest.raises(IOErr,
                       match=r"^please supply array\(s\) to write to file$"):
        write_arrays(str(tfile), [1, 2, 3], format="<%02d>", clobber=True)


def test_write_data_no_fields(tmp_path):
    """Basic check of write_data with no explicit fields"""

    tfile = tmp_path / "test.dat"

    d = Data1D("x", [1, 2], [4, 5], syserror=[0.2, 0.3])
    write_data(str(tfile), d,
               format="%.1f", clobber=True)
    assert tfile.read_text() == "#X Y SYSERROR\n1.0 4.0 0.2\n2.0 5.0 0.3\n"


def test_write_data_with_fields(tmp_path):
    """Basic check of write_data with fields set"""

    tfile = tmp_path / "test.dat"

    d = Data1D("x", [1, 2], [4, 5], syserror=[0.2, 0.3])
    write_data(str(tfile), d, fields=["y", "syserror"],
               format="%.1f", clobber=True)
    assert tfile.read_text() == "#Y SYSERROR\n4.0 0.2\n5.0 0.3\n"
