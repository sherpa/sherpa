#
# Copyright (C) 2020  Smithsonian Astrophysical Observatory
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

"""Very-basic tests of the HTML representation of models and parameters.

"""

import numpy as np

from sherpa.astro.optical import AbsorptionEdge
from sherpa.models.parameter import Parameter
from sherpa.models.basic import Const1D, Gauss1D, Gauss2D


def test_parameter():
    p = Parameter('foo', 'x1', 1.2)
    r = p._repr_html_()

    assert r is not None

    assert '<summary>Parameter</summary>' in r
    assert '<table class="model">' in r
    assert '<tr><th class="model-odd">foo</th><td>x1</td><td><input disabled type="checkbox" checked></input></td><td>1.2</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r


def test_parameter_frozen():
    p = Parameter('foo', 'x1', 1.2)
    p.freeze()
    r = p._repr_html_()

    assert r is not None

    assert '<summary>Parameter</summary>' in r
    assert '<table class="model">' in r
    assert '<tr><th class="model-odd">foo</th><td>x1</td><td><input disabled type="checkbox"></input></td><td>1.2</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r


def test_parameter_linked():
    p = Parameter('foo', 'x1', 1.2)
    q = Parameter('bar', 'x2', 2.2)
    p.val = 2 + q
    r = p._repr_html_()

    assert r is not None

    assert '<summary>Parameter</summary>' in r
    assert '<table class="model">' in r
    assert '<th class="model-odd">foo</th><td>x1</td><td>linked</td><td>4.2</td><td colspan="2">&#8656; 2 + bar.x2</td><td></td></tr>' in r


def test_model():
    m = Gauss1D('ff')
    r = m._repr_html_()

    assert r is not None

    assert '<summary>Model</summary>' in r
    assert '<table class="model">' in r

    assert '<tr><th class="model-odd" scope="rowgroup" rowspan=3>ff</th><td>fwhm</td><td><input disabled type="checkbox" checked></input></td><td>10.0</td><td>TINY</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>pos</td><td><input disabled type="checkbox" checked></input></td><td>0.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>ampl</td><td><input disabled type="checkbox" checked></input></td><td>1.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r


def test_model_tau():
    """Check special-case of 2pi"""
    m = Gauss2D('g2')
    r = m._repr_html_()

    assert r is not None

    assert '<tr><th class="model-odd" scope="rowgroup" rowspan=6>g2</th><td>fwhm</td><td><input disabled type="checkbox" checked></input></td><td>10.0</td><td>TINY</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>xpos</td><td><input disabled type="checkbox" checked></input></td><td>0.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>ellip</td><td><input disabled type="checkbox"></input></td><td>0.0</td><td>0.0</td><td>0.999</td><td></td></tr>'
    assert '<tr><td>theta</td><td><input disabled type="checkbox"></input></td><td>0.0</td><td>-2&#960;</td><td>2&#960;</td><td>radians</td></tr>' in r


def test_model_pi():
    """Check special-case of pi"""
    m = Gauss2D('g2')
    m.theta.min = -np.pi
    m.theta.max = np.pi
    r = m._repr_html_()

    assert r is not None

    assert '<tr><th class="model-odd" scope="rowgroup" rowspan=6>g2</th><td>fwhm</td><td><input disabled type="checkbox" checked></input></td><td>10.0</td><td>TINY</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>xpos</td><td><input disabled type="checkbox" checked></input></td><td>0.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>ellip</td><td><input disabled type="checkbox"></input></td><td>0.0</td><td>0.0</td><td>0.999</td><td></td></tr>'
    assert '<tr><td>theta</td><td><input disabled type="checkbox"></input></td><td>0.0</td><td>-&#960;</td><td>&#960;</td><td>radians</td></tr>' in r


def test_model_linked():
    """Check linking of models"""
    m = Gauss1D('g1')
    c = Const1D('c1')
    m.fwhm = 8 * c.c0
    r = m._repr_html_()

    assert r is not None

    assert '<tr><th class="model-odd" scope="rowgroup" rowspan=3>g1</th><td>fwhm</td><td>linked</td><td>8.0</td><td colspan=2>&#8656; 8 * c1.c0</td><td></td></tr>' in r
    assert '<tr><td>pos</td><td><input disabled type="checkbox" checked></input></td><td>0.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>ampl</td><td><input disabled type="checkbox" checked></input></td><td>1.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r


def test_model_combined():
    """We can show a binary op"""
    m1 = Gauss1D('g1')
    m2 = Const1D('c1')
    m = m1 + m2
    r = m._repr_html_()

    assert r is not None

    assert '<th class="model-odd" scope="rowgroup" rowspan=3>g1</th>' in r
    assert '<th class="model-even" scope="rowgroup" rowspan=1>c1</th>' in r

    assert '<tr><th class="model-odd" scope="rowgroup" rowspan=3>g1</th><td>fwhm</td><td><input disabled type="checkbox" checked></input></td><td>10.0</td><td>TINY</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>pos</td><td><input disabled type="checkbox" checked></input></td><td>0.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r
    assert '<tr><td>ampl</td><td><input disabled type="checkbox" checked></input></td><td>1.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r
    assert '<tr class="block"><th class="model-even" scope="rowgroup" rowspan=1>c1</th><td>c0</td><td><input disabled type="checkbox" checked></input></td><td>1.0</td><td>-MAX</td><td>MAX</td><td></td></tr>' in r


def test_model_hidden():
    """We can show a hidden parameter"""
    m = AbsorptionEdge('mdl')
    r = m._repr_html_()

    assert r is not None

    assert '<tbody><tr><th class="model-odd" scope="rowgroup" rowspan=2>mdl</th><td>edgew</td><td><input disabled type="checkbox"></input></td><td>5000.0</td><td>TINY</td><td>MAX</td><td>angstroms</td></tr><tr><td>tau</td><td><input disabled type="checkbox" checked></input></td><td>0.5</td><td>-MAX</td><td>MAX</td><td></td></tr></tbody>' in r

def test_model_combined_samename():
    """We can show a binary op"""
    m1 = Gauss1D('name')
    m2 = Gauss1D('name')
    m = m1 + m2
    r = m._repr_html_()

    assert r is not None

    assert '<th class="model-odd" scope="rowgroup" rowspan=3>name</th>' in r
    assert '<th class="model-even" scope="rowgroup" rowspan=3>name</th>' in r
