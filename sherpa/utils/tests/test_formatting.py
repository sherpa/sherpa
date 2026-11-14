#
#  Copyright (C) 2026, MIT
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
import numpy as np
import pytest

from sherpa.utils.formatting import smallest_significant_order, latex_val, latex_val_errors


@pytest.mark.parametrize(
    "val, expected",
    [
        (0.01, -2),
        (0.03, -2),
        (0.06, -2),
        (0.1, -1),
        (0.3, -1),
        (0.6, -1),
        (1, 0),
        (1000, 3),
    ]
)
def test_smallest_significant_order(val, expected):
    '''Test with a range of values that are powers of 10 and some in between.'''
    assert smallest_significant_order(val) == expected


def test_smallest_significant_negative():
    """Test that the function works for negative values as well."""
    for n in np.random.uniform(0.001, 10, 100):
        # There is nothing special about 2.34, it is just a number that is not a power of 10
        val = 2.34**n
        assert smallest_significant_order(-val) == smallest_significant_order(val)


def test_smallest_significant_order_array():
    """Test that the function works for arrays of values."""
    vals = np.array([0.01, -0.003, 12, 1e12])
    assert smallest_significant_order(vals) == -3
    # And make a more random test with a range of values
    vals = np.random.uniform(0, 10, 10) ** np.random.uniform(0.001, 4, 10)
    assert smallest_significant_order(vals) == smallest_significant_order(np.min(np.abs(vals)))


def test_latex_val():
    """Test that the latex_val function works as expected."""
    assert latex_val(0.000123) == "$0.000123$"
    assert latex_val(0.0000123) == '$1.23 \\times 10^{-5}$'


def test_latex_val_errors():
    """Test that the latex_val_errors function works as expected."""
    assert latex_val_errors(-12345, 100, 200) == r'$-12300^{+100}_{-200}$'
    val = 0.000123
    up = 0.0000134
    down = 0.00002123
    assert latex_val_errors(val, up, down) == '$12^{+1}_{-2} \\times 10^{-5}$'
    assert latex_val_errors(val, up, -down) == '$12^{+1}_{-2} \\times 10^{-5}$'
    # Only one number small enough to be in scientific notation
    up = 0.0134
    down = -0.00003
    assert latex_val_errors(val, up, down) == '$12^{+1340}_{-3} \\times 10^{-5}$'
    assert latex_val_errors(val, up, -down) == '$12^{+1340}_{-3} \\times 10^{-5}$'


def test_latex_val_errors_trailing_zeros():
    """We want to keep trailing zeros so that the printed precision
    matches the precision of the error bars."""
    assert latex_val_errors(1., 0.001, -0.32451) == '$1.000^{+0.001}_{-0.325}$'


def test_latex_val_errors_none():
    """Test that the latex_val_errors function works as expected with None values."""
    val = 0.000123
    up = None
    down = 0.00002123
    assert latex_val_errors(val, up, down) == '$12^{+\\infty}_{-2} \\times 10^{-5}$'