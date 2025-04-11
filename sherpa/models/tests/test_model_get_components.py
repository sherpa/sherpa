#
#  Copyright (C) 2025
#  MIT
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
"""Test several methods to get the components of a more complex model.

These tests could be in the test_model.py file, but that file is already pretty
long, so we move them out to a separate file for easier readability.
"""
from typing import Literal
import pytest

from sherpa.models.basic import Gauss1D, Polynom1D
from sherpa.models.model import ArithmeticConstantModel, BinaryOpModel, ArithmeticModel

class Gauss1DSubClass(Gauss1D):
    """A Gauss1D model with a different name to test get_components_by_name."""
    pass

l1 = Gauss1D('l1')
l2 = Gauss1D('l_other_name')
l_same_name = Gauss1D('l1')
l1_subclass = Gauss1DSubClass('l1_subclass')
l2_subclass_same_name = Gauss1DSubClass('l1')
b = Polynom1D('b')


@pytest.mark.parametrize("func", ['get_parts', '_get_parts'])
@pytest.mark.parametrize("remove_duplicates", [True, False])
def test_get_parts_composite(func: Literal['get_parts'] | Literal['_get_parts'], remove_duplicates: bool):
    """Test the get_components method of a composite model."""
    mdl = l1 + (0.5 * l2) + b + b
    parts = getattr(mdl, func)(remove_duplicates=remove_duplicates)
    if remove_duplicates:
        assert len(parts) == 8
    else:
        assert len(parts) == 9
        assert parts == [p for p in mdl]
    for i in [0, 1, 2, 4]:
        assert isinstance(parts[i], BinaryOpModel)
    assert parts[3] == l1
    assert isinstance(parts[5], ArithmeticConstantModel)
    assert parts[6] == l2
    assert parts[7] == b
    if not remove_duplicates:
        assert parts[8] == b


@pytest.mark.parametrize("remove_duplicates", [True, False])
def test_get_parts(remove_duplicates: bool):
    """Test the get_components method of a non-composite model."""
    parts = l1.get_parts(remove_duplicates=remove_duplicates)
    assert len(parts) == 1
    assert parts[0] == l1
    assert parts == [p for p in l1]


@pytest.mark.parametrize("func", ['get_parts', '_get_parts'])
@pytest.mark.parametrize("remove_duplicates", [True, False])
def test_get_parts_composite_no_composites(func: Literal['get_parts'] | Literal['_get_parts'], remove_duplicates: bool):
    """Test the get_components method of a model."""
    # Test the get_components method
    mdl = l1 + (0.5 * l2) + b + b
    parts = getattr(mdl, func)(include_composites=False, remove_duplicates=remove_duplicates)
    if remove_duplicates:
        assert len(parts) == 4
    else:
        assert len(parts) == 5

    assert parts[0] == l1
    assert isinstance(parts[1], ArithmeticConstantModel)
    assert parts[2] == l2
    assert parts[3] == b
    if not remove_duplicates:
        assert parts[4] == b


def test_get_components_by_name():
    """Test the get_components method of a model."""
    # Simple case: One component with the name
    assert [l2] == l2.get_components_by_name('l_other_name')

    with pytest.raises(KeyError, match="No components found with name 'ls'"):
        parts = l2.get_components_by_name('ls')


def test_get_components_by_name_composite():
    """Test the get_components method of a composite model."""
    # Simple case: One component with the name
    mdl = l1 + (0.5 * l2) + b + b + l_same_name
    assert [l2] == mdl.get_components_by_name('l_other_name')

    # One component appears twice
    assert [b] == mdl.get_components_by_name('b')

    # More than one component with the same string name
    assert [l1, l_same_name] == mdl.get_components_by_name('l1')


    with pytest.raises(KeyError, match="No components found with name 'ls'"):
        parts = mdl.get_components_by_name('ls')


def test_components_by_class():
    """Test the get_components method of a model."""
    assert [l1] == l1.get_components_by_class(Gauss1D)


def test_components_by_class_composite():
    """Test the get_components method of a model."""
    mdl = l1 + (0.5 * l2) + b + b + l_same_name

    # Simple case: One component with the class
    parts = mdl.get_components_by_class(ArithmeticConstantModel)
    assert len(parts) == 1
    assert isinstance(parts[0], ArithmeticConstantModel)
    assert parts[0].name == '0.5'

    # Simple case: Each component appears only once
    mdl = l1 + (0.5 * l2) + b + b + l_same_name
    assert [l1, l2, l_same_name] == mdl.get_components_by_class(Gauss1D)

    # One component appears twice
    assert [b] == mdl.get_components_by_class(Polynom1D)


def test_components_by_class_subclass_composite():
    """Test the get_components method of a model in the presence of subclasses."""
    mdl = l1 * l1_subclass - abs(l_same_name)
    assert [l1_subclass] == mdl.get_components_by_class(Gauss1DSubClass)

    assert [l1, l1_subclass, l_same_name] == mdl.get_components_by_class(Gauss1D)

    assert [l1, l_same_name] == mdl.get_components_by_class(Gauss1D, subclass_ok=False)

    with pytest.raises(KeyError, match="No components found with type '<class 'sherpa.models.basic.Polynom1D'>'"):
        parts = mdl.get_components_by_class(Polynom1D)

    parts = mdl.get_components_by_class(ArithmeticModel, subclass_ok=True)
    assert len(parts) == 6

    with pytest.raises(KeyError, match="No components found with type '<class 'sherpa.models.model.ArithmeticModel'>'"):
        parts = mdl.get_components_by_class(ArithmeticModel, subclass_ok=False)
