#
#  Copyright (C) 2025
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
"""Tests related to caching model evaluations

More tests with 1D data are in test_model.py and might be moved here at a later point.
"""
import logging

import numpy as np
import pytest

from sherpa.data import Data1D
from sherpa.models import Polynom2D, Gauss2D
from sherpa.astro.models import Beta2D, Lorentz2D
from sherpa.models.basic import Polynom1D, Gauss1D, Sin
from sherpa.models.model import Model, UnaryOpModel, BinaryOpModel, RegridWrappedModel, \
    hashfunc, modelCacher1d, ArithmeticConstantModel
from sherpa.fit import Fit

from sherpa.models.tests.test_model import ReportKeywordsModel


@pytest.mark.parametrize("modelclass", [Polynom2D, Gauss2D, Beta2D, Lorentz2D])
def test_default_cache_size(modelclass):
    """Check the default cache size for some models."""

    mdl = modelclass()
    assert mdl.cache == 0


def test_evaluate_no_cache2d():
    """Check we can turn off caching: 2d"""

    x0grid = np.array([2, 10, 1.5])
    x1grid = np.array([-2, 0, 3.])

    mdl = Polynom2D()
    mdl.integrate = False
    mdl.cache = 0
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(3)
    assert mdl(x0grid, x1grid) == pytest.approx(expected)
    assert len(mdl._cache) == 0

    mdl.cx1 = 5
    mdl.cy1 = 2

    expected = 1 + 2 * x0grid + 5 * x1grid
    assert mdl(x1grid, x0grid) == pytest.approx(expected)
    assert len(mdl._cache) == 0


def test_evaluate_cache2d():
    """Check we run with caching on: 1d"""

    x0grid = np.array([2, 10, 1.5])
    x1grid = np.array([-2, 0, 3.])

    mdl = Polynom2D()
    mdl.integrate = False
    mdl.cache = 5
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(3)
    assert mdl(x1grid, x0grid) == pytest.approx(expected)
    # Model has been run once, so there is one entry in the cache
    assert len(mdl._cache) == 1
    assert mdl._cache[list(mdl._cache.keys())[0]] == pytest.approx(expected)

    mdl.cx1 = 5
    mdl.cy1 = 2

    expected2 = 1 + 2 * x0grid + 5 * x1grid
    assert mdl(x1grid, x0grid) == pytest.approx(expected2)
    assert len(mdl._cache) == 2

    # Now, the model has run a second time, so there are two entries in the cache
    assert mdl._cache[list(mdl._cache.keys())[0]] == pytest.approx(expected)
    assert mdl._cache[list(mdl._cache.keys())[1]] == pytest.approx(expected2)


def test_cache_is_actually_used_2():
    """Most other tests check that the cache has values in it,
    but not that those values are actually returned.
    Here, we manipulate the cached value and then call the model
    to check that the cached value is used.
    """
    x0grid = [2, 10, 1.5]
    x1grid = [-2, 0, 3.]

    mdl = Polynom2D()
    mdl.integrate = False
    mdl.cache = 5

    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(3)
    assert mdl(x0grid, x1grid) == pytest.approx(expected)

    # Manipulate the values in the cache
    mdl._cache[list(mdl._cache.keys())[0]] = 2 * expected
    assert mdl(x0grid, x1grid) == pytest.approx(2 * expected)


def test_evaluate_cache_many_parameters():
    """Check that we cache correctly with > 2 parameters, e.g.
    for a Data2DInt (with x0lo, x1lo, x0hi, x1hi) or an ND model.
    Since Data2DInt is used more often, there are models defined for
    that in Sherpa already and that's an easier approach to write this test.
    """
    x0lo = np.arange(-20, -10, 1)
    x0hi = x0lo + 1
    x1lo = np.arange(-3.2, 8.4, 1.2)
    x1hi = x1lo + 1

    mdl = Polynom2D()
    mdl.integrate = True
    mdl.cache = 5

    assert len(mdl._cache) == 0
    out = mdl(x0lo, x1lo, x0hi, x1hi)
    assert len(mdl._cache) == 1
    # Manipulate the values in the cache to check if it's used
    mdl._cache[list(mdl._cache.keys())[0]] = 2 * out

    # second call should return value from cache
    out1 = mdl(x0lo, x1lo, x0hi, x1hi)
    assert len(mdl._cache) == 1
    assert pytest.approx(out1) == 2 * out
    # Call with third array different, so should not return cached value
    x0hi2 = x1hi + 0.1
    out2 = mdl(x0lo, x1lo, x0hi2, x1hi)
    assert len(mdl._cache) == 2
    assert pytest.approx(out2) != out1
    assert pytest.approx(out2) != 2 * out
    # And just for good measure, try the first call again
    out3 = mdl(x0lo, x1lo, x0hi, x1hi)
    assert len(mdl._cache) == 2
    assert pytest.approx(out3) == 2 * out


def test_evaluate_cache_regrid2d():
    """How about a regridded model?"""

    mdl = Lorentz2D()

    x0 = np.arange(0, 10, 0.1)
    x1 = np.arange(2, 20, 0.1)
    rmdl = mdl.regrid(x0, x1)

    assert isinstance(rmdl, RegridWrappedModel)
    assert not hasattr(rmdl, 'cache')


def test_evaluate_cache2dint():
    """Check we run with caching on: 2dint"""

    x0lo = np.array([-20, -1, 10])
    x1lo = np.array([-25, -1, -10])
    x0hi = np.array([-10, 1, 15])
    x1hi = np.array([-10, 1, -8])

    mdl = Polynom2D()
    mdl.cache = 12
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.array([150, 4, 10])
    assert mdl(x0lo, x1lo, x0hi, x1hi) == pytest.approx(expected)
    assert len(mdl._cache) == 1
    assert mdl._cache[list(mdl._cache.keys())[0]] == pytest.approx(expected)

    # Manipulate the values in the cache to check if it's used
    mdl._cache[list(mdl._cache.keys())[0]] = 2 * expected
    assert mdl(x0lo, x1lo, x0hi, x1hi) == pytest.approx(2 * expected)



def check_cache(mdl, expected, x, xhi=None, cache_size=1):
    """Check the cache contents.

    The code matches that in sherpa.models.model.modelCacher1d,
    so all it does is check we are using this method.
    """

    cache = mdl._cache
    assert len(cache) == cache_size

    pars = [p.val for p in mdl.pars]
    data = [np.asarray(pars).tobytes(),
            b'1' if mdl.integrate else b'0',
            x.tobytes()]
    if xhi is not None:
        data.append(xhi.tobytes())

    token = b''.join(data)
    digest = hashfunc(token).digest()
    assert digest in cache
    assert cache[digest] == pytest.approx(expected)


def test_deprecated_use_caching():
    """_use_caching is rundundant with `cache`, but kept for
    backwards compatibility. In this test, we check that it
    still works.
    """
    mdl = Polynom1D()
    assert mdl._use_caching is True
    with pytest.warns(DeprecationWarning):
        mdl._use_caching = False
    assert mdl._use_caching is False
    assert mdl.cache == 0
    with pytest.warns(DeprecationWarning):
        mdl._use_caching = True
    assert mdl._use_caching is True
    assert mdl.cache == 5


def test_evaluate_no_cache1d_use_caching():
    """Check we can turn off caching: 1d"""

    xgrid = np.arange(2, 10, 1.5)

    mdl = Polynom1D()
    mdl.integrate = False
    with pytest.warns(DeprecationWarning):
        mdl._use_caching = False
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(6)
    assert mdl(xgrid) == pytest.approx(expected)
    assert len(mdl._cache) == 0

    mdl.c0 = 5
    mdl.c1 = 2

    expected = 5 + 2 * xgrid
    assert mdl(xgrid) == pytest.approx(expected)
    assert len(mdl._cache) == 0


def test_evaluate_no_cache1d():
    """Check we can turn off caching: 1d"""

    xgrid = np.arange(2, 10, 1.5)

    mdl = Polynom1D()
    mdl.integrate = False
    mdl.cache = 0
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(6)
    assert mdl(xgrid) == pytest.approx(expected)
    assert len(mdl._cache) == 0

    mdl.c0 = 5
    mdl.c1 = 2

    expected = 5 + 2 * xgrid
    assert mdl(xgrid) == pytest.approx(expected)
    assert len(mdl._cache) == 0


def test_evaluate_cache1d():
    """Check we run with caching on: 1d"""

    xgrid = np.arange(2, 10, 1.5)

    mdl = Polynom1D()
    mdl.cache = 1
    mdl.integrate = False
    mdl.cache = 5
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(6)
    assert mdl(xgrid) == pytest.approx(expected)
    # Model has been run once, so there is one entry in the cache
    check_cache(mdl, expected, xgrid, cache_size=1)

    mdl.c0 = 5
    mdl.c1 = 2

    expected = 5 + 2 * xgrid
    assert mdl(xgrid) == pytest.approx(expected)
    # Now, the model has run a second time, so there are two entries in the cache
    check_cache(mdl, expected, xgrid, cache_size=2)


def test_cache_is_actually_used():
    """Most other tests check that the cache has values in it,
    but not that those values are actually returned.
    Here, we manipulte the cached value and then call the model
    to check that the cached value is used.
    """
    xgrid = np.arange(2, 10, 1.5)

    mdl = Polynom1D()
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(6)
    assert mdl(xgrid) == pytest.approx(expected)

    # Manipulate the values in the cache
    mdl._cache[list(mdl._cache.keys())[0]] = 2 * expected
    assert mdl(xgrid) == pytest.approx(2 * expected)


def test_evaluate_no_cache1dint():
    """Check we can turn off caching: 1dint"""

    xgrid = np.arange(2, 10, 1.5)
    xlo, xhi = xgrid[:-1], xgrid[1:]

    mdl = Polynom1D()
    mdl.cache = 0
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(5) * 1.5
    assert mdl(xlo, xhi) == pytest.approx(expected)
    assert len(mdl._cache) == 0

    mdl.c0 = 5
    mdl.c1 = 2

    def xval(x):
        return 5 + 2 * x

    dx = xhi - xlo
    ylo = xval(xlo)
    yhi = xval(xhi)

    expected = dx * (yhi + ylo) / 2
    assert mdl(xlo, xhi) == pytest.approx(expected)
    assert len(mdl._cache) == 0


def test_evaluate_cache1dint():
    """Check we run with caching on: 1dint"""

    xgrid = np.arange(2, 10, 1.5)
    xlo, xhi = xgrid[:-1], xgrid[1:]

    mdl = Polynom1D()
    mdl.cache = 1
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(5) * 1.5
    assert mdl(xlo, xhi) == pytest.approx(expected)
    check_cache(mdl, expected, xlo, xhi)

    mdl.c0 = 5
    mdl.c1 = 2

    def xval(x):
        return 5 + 2 * x

    dx = xhi - xlo
    ylo = xval(xlo)
    yhi = xval(xhi)

    expected = dx * (yhi + ylo) / 2
    assert mdl(xlo, xhi) == pytest.approx(expected)
    check_cache(mdl, expected, xlo, xhi)


def test_evaluate_cache_swap():
    """Check issue #959 when swapping integrate flag caused problems.

    Note that the problem causing #959 is actually tested in
    test_evaluate_cache1dint but it's nice to have this
    separate check in case things change.
    """

    xgrid = np.arange(2, 10, 1.5)
    xlo, xhi = xgrid[:-1], xgrid[1:]

    mdl = Polynom1D()
    mdl.cache = 1

    mdl.c0 = 5
    mdl.c1 = 2

    mdl.integrate = False
    expected = 5 + 2 * xlo

    y1 = mdl(xlo, xhi)
    assert y1 == pytest.approx(expected)
    check_cache(mdl, expected, xlo, xhi)

    mdl.integrate = True

    def xval(x):
        return 5 + 2 * x

    dx = xhi - xlo
    ylo = xval(xlo)
    yhi = xval(xhi)

    expected = dx * (yhi + ylo) / 2

    y2 = mdl(xlo, xhi)
    assert y2 == pytest.approx(expected)
    check_cache(mdl, expected, xlo, xhi)


def test_evaluate_cache_arithmeticconstant():
    """Check we run with caching: ArihmeticConstant"""

    mdl = ArithmeticConstantModel(2.3)
    assert not hasattr(mdl, 'cache')


def test_evaluate_cache_unaryop():
    """UnaryOp has no cache"""

    mdl = Polynom1D()
    assert hasattr(mdl, 'cache')

    fmdl = -mdl
    assert isinstance(fmdl, UnaryOpModel)
    assert not hasattr(fmdl, 'cache')


def test_evaluate_cache_binaryop():
    """BinaryOp has no cache"""

    mdl = Polynom1D()
    assert hasattr(mdl, 'cache')

    fmdl = mdl + 2
    assert isinstance(fmdl, BinaryOpModel)
    assert not hasattr(fmdl, 'cache')


def test_evaluate_cache_regrid1d():
    """How about a regridded model?"""

    mdl = Polynom1D()

    x = np.arange(2, 20, 0.5)
    rmdl = mdl.regrid(x)

    assert isinstance(rmdl, RegridWrappedModel)
    assert not hasattr(rmdl, 'cache')


class DoNotUseModel(Model):
    """This is only used to get a integrate-free model.

    ArithmeticModel sets a integrate field, so this is derived
    from Model. This requires adding some fields from ArithmeticModel
    to support use with modelCacher1d.
    """

    def __init__(self, *args, **kwargs) -> None:
        # Model caching ability
        self.cache_clear()
        Model.__init__(self, *args, **kwargs)

    def cache_clear(self) -> None:
        """Clear the cache."""

        self._cache: dict[bytes, np.ndarray] = {}
        self._cache_ctr: dict[str, int] = {'hits': 0, 'misses': 0, 'check': 0}
        self.cache: int = 2


    @modelCacher1d
    def calc(self, p, *args, **kwargs):
        """p is ignored."""

        return np.ones(args[0].size)


def get_cache_classes():
    """This is a function because we want to conditionally
    include an XSPEC model. Within a function, we can simply
    pass that if XSPEC is not available.
    """
    cls_list = [DoNotUseModel, Polynom1D]

    try:
        from sherpa.astro.xspec import XSphabs
        cls_list.append(XSphabs)
    except ImportError:
        pass

    return cls_list


@pytest.mark.parametrize('cls', get_cache_classes())
def test_cache_uses_instance_attributes(cls):
    """Check that the cache uses the instance attributes.

    This tests both the real code (e.g. ArithmetricModels)
    but also the DoNotUseModel class which we defined above
    just for testing, because it previously failed just in that
    test class in a hard-to-debug way.
    """
    mdl = cls("some-name")

    # cache is an int, not a mutable object, so it's OK if it's the same
    # as the class attribute, unless we set it - then it ought to be different.
    mdl.cache = 2234

    for attr in ["cache", "_cache", "_cache_ctr"]:
        assert hasattr(mdl, attr)
        if hasattr(cls, attr):
            assert id(getattr(mdl, attr)) != id(getattr(cls, attr))


def test_cache_reset_when_size_changes():
    """Check the cache is reset when the cache size changes."""

    mdl = Polynom1D()
    mdl.cache = 2

    x = np.arange(2, 10, 1.5)
    mdl(x)

    assert len(mdl._cache) == 1
    assert mdl._cache_ctr['check'] == 1

    mdl.cache = 3
    assert len(mdl._cache) == 0
    assert mdl._cache_ctr['check'] == 0

    mdl(x)
    assert len(mdl._cache) == 1
    assert mdl._cache_ctr['check'] == 1


def test_caching_not_used_when_set_to_zero():
    """Check the cache is not used when the cache size is zero."""

    mdl = Polynom1D()
    mdl.cache = 0

    x = np.arange(2, 10, 1.5)
    mdl(x)

    assert len(mdl._cache) == 0
    # We checked, but found nothing (because the cache is empty)
    assert mdl._cache_ctr['check'] == 1


def test_cache_not_used_in_fit():
    """Check the cache is not used in a fit when `fit=False`.

    How do we do that without relying too much on the internal
    implementation?
    """
    mdl = Gauss1D('con1')
    mdl.ampl = 1
    mdl.pos = 0
    mdl.fwhm = 1
    dat = Data1D('data', [-2, -1, 0, 1, 2,], [0.05, 1, 2, 1, 0.05], np.ones(5))
    fit = Fit(dat, mdl)
    res = fit.fit()
    assert mdl.ampl.val == pytest.approx(2.02, rel=1e-2)
    assert mdl.pos.val == pytest.approx(0, abs=1e-6)
    assert len(mdl._cache) == 5

    mdl.cache_clear()
    mdl.ampl = 1
    mdl.pos = 0
    mdl.fwhm = 1
    res = fit.fit(cache=False)
    assert mdl.ampl.val == pytest.approx(2.02, rel=1e-2)
    assert mdl.pos.val == pytest.approx(0, abs=1e-6)
    assert len(mdl._cache) == 0
    # We cleared the cache and then called the fit with `cache=False`.
    # After the fit, the cache is still empty, so presumably it was never used.


def test_cache_integrate_fall_through_no_integrate():
    """Try and test the fall-through of the integrate setting.

    This is a bit contrived.
    """

    mdl = DoNotUseModel('notme')
    x = np.asarray([2, 3, 7, 100])
    y = mdl(x)

    expected = [1, 1, 1, 1]
    assert y == pytest.approx(expected)

    # A simplified version of check_cache (as we have no
    # integrate setting).
    #
    cache = mdl._cache
    assert len(cache) == 1

    pars = []
    data = [np.asarray(pars).tobytes(),
            b'0', # not integrated
            x.tobytes()]

    token = b''.join(data)
    digest = hashfunc(token).digest()
    assert digest in cache
    assert cache[digest] == pytest.approx(expected)


def test_cache_integrate_fall_through_integrate_true():
    """See also test_cache_integrate_fall_through_no_integrate."""

    mdl = DoNotUseModel('notme')

    cache = mdl._cache
    assert len(cache) == 0

    x = np.asarray([2, 3, 7, 100])
    y = mdl(x, integrate=True)

    expected = [1, 1, 1, 1]
    assert y == pytest.approx(expected)

    # A simplified version of check_cache (as we have no
    # integrate setting).
    #
    cache = mdl._cache
    assert len(cache) == 1

    pars = []
    data = [np.asarray(pars).tobytes(),
            b'1', # integrated
            x.tobytes(),
            # The integrate setting is included twice because we can
            # not guarantee it has been sent in with a keyword
            # argument.
            b'integrate',
            np.asarray(True).tobytes()]

    token = b''.join(data)
    digest = hashfunc(token).digest()
    assert digest in cache
    assert cache[digest] == pytest.approx(expected)


def test_cache_integrate_fall_through_integrate_false():
    """See also test_cache_integrate_fall_through_no_integrate."""

    mdl = DoNotUseModel('notme')
    x = np.asarray([2, 3, 7, 100])
    y = mdl(x, integrate=False)

    expected = [1, 1, 1, 1]
    assert y == pytest.approx(expected)

    # A simplified version of check_cache (as we have no
    # integrate setting).
    #
    cache = mdl._cache
    assert len(cache) == 1

    pars = []
    data = [np.asarray(pars).tobytes(),
            b'0', # not integrated
            x.tobytes(),
            # The integrate setting is included twice because we can
            # not guarantee it has been sent in with a keyword
            # argument.
            b'integrate',
            np.asarray(False).tobytes()]

    token = b''.join(data)
    digest = hashfunc(token).digest()
    assert digest in cache
    assert cache[digest] == pytest.approx(expected)


def test_cache_status_single(caplog):
    """Check cache_status for a single model."""

    p = Polynom1D()
    with caplog.at_level(logging.INFO, logger='sherpa'):
        p.cache_status()

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == 'sherpa.models.model'
    assert lvl == logging.INFO
    toks = msg.split()
    assert toks[0] == 'polynom1d'
    assert toks[1] == 'size:'
    assert toks[2] == '0'
    assert toks[3] == 'hits:'
    assert toks[4] == '0'
    assert toks[5] == 'misses:'
    assert toks[6] == '0'
    assert toks[7] == 'check:'
    assert toks[8] == '0'
    assert len(toks) == 9


def test_cache_status_multiple(caplog):
    """Check cache_status for a multi-component model.

    Unlike test_cache_syayus_single we also have evaluated the model
    so we can check that the cache status has changed.
    """

    # The model expression includes an ArithmeticConstant model (the
    # term 2) which does not have a cache and so is ignored by
    # cache_status.
    #
    p = Polynom1D()
    b = Sin()
    c = Gauss1D()
    mdl = c * (2 * p + b)

    # One model is not cached
    b.cache = 0

    mdl([0.1, 0.2, 0.3])
    mdl([0.1, 0.2, 0.3])
    mdl([0.1, 0.2, 0.3, 0.4])

    with caplog.at_level(logging.INFO, logger='sherpa'):
        mdl.cache_status()

    assert len(caplog.records) == 3

    tokens = []
    for lname, lvl, msg in caplog.record_tuples:
        assert lname == 'sherpa.models.model'
        assert lvl == logging.INFO
        toks = msg.split()
        assert len(toks) == 9
        assert toks[1] == 'size:'
        assert toks[3] == 'hits:'
        assert toks[5] == 'misses:'
        assert toks[7] == 'check:'
        assert toks[8] == '3'

        tokens.append(toks)

    toks = tokens[0]
    assert toks[0] == 'gauss1d'
    assert toks[2] == '2'
    assert toks[4] == '1'
    assert toks[6] == '2'

    toks = tokens[1]
    assert toks[0] == 'polynom1d'
    assert toks[2] == '2'
    assert toks[4] == '1'
    assert toks[6] == '2'

    toks = tokens[2]
    assert toks[0] == 'sin'
    assert toks[2] == '0'
    assert toks[4] == '0'
    assert toks[6] == '0'


def test_cache_clear_single():
    """Check cache_clear for a single model."""

    p = Polynom1D()

    # There's no official API for accessing the cache data,
    # so do it directly.
    #
    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0

    p([1, 2, 3])
    p([1, 2, 3])
    p([1, 2, 3, 4])

    assert len(p._cache) == 2
    assert p._cache_ctr['check'] == 3
    assert p._cache_ctr['hits'] == 1
    assert p._cache_ctr['misses'] == 2

    p.cache_clear()

    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0


def test_cache_clear_multiple():
    """Check cache_clear for a combined model."""

    p = Polynom1D()
    b = Sin()
    c = Gauss1D()
    mdl = c * (p + 2 * b)

    # Ensure one component doesn't use the cache
    c.cache = 0

    # There's no official API for accessing the cache data,
    # so do it directly.
    #
    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0

    assert len(b._cache) == 0
    assert b._cache_ctr['check'] == 0
    assert b._cache_ctr['hits'] == 0
    assert b._cache_ctr['misses'] == 0

    assert len(c._cache) == 0
    assert c._cache_ctr['check'] == 0
    assert c._cache_ctr['hits'] == 0
    assert c._cache_ctr['misses'] == 0

    mdl([1, 2, 3])
    mdl([1, 2, 3])
    mdl([1, 2, 3, 4])

    assert len(p._cache) == 2
    assert p._cache_ctr['check'] == 3
    assert p._cache_ctr['hits'] == 1
    assert p._cache_ctr['misses'] == 2

    assert len(b._cache) == 2
    assert b._cache_ctr['check'] == 3
    assert b._cache_ctr['hits'] == 1
    assert b._cache_ctr['misses'] == 2

    assert len(c._cache) == 0
    assert c._cache_ctr['check'] == 3
    assert c._cache_ctr['hits'] == 0
    assert c._cache_ctr['misses'] == 0

    mdl.cache_clear()

    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0

    assert len(b._cache) == 0
    assert b._cache_ctr['check'] == 0
    assert b._cache_ctr['hits'] == 0
    assert b._cache_ctr['misses'] == 0

    assert len(c._cache) == 0
    assert c._cache_ctr['check'] == 0
    assert c._cache_ctr['hits'] == 0
    assert c._cache_ctr['misses'] == 0


def test_cache_keeps_limited_size():
    """Check cache_clear for a single model."""

    p = Polynom1D()
    p.cache = 2

    # There's no official API for accessing the cache data,
    # so do it directly.
    #
    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0

    p([1, 2, 3])
    p([1, 2, 3])
    p([1, 2, 3, 4])

    assert len(p._cache) == 2
    assert p._cache_ctr['check'] == 3
    assert p._cache_ctr['hits'] == 1
    assert p._cache_ctr['misses'] == 2

    p([1.2, 2, 3, 4])

    assert len(p._cache) == 2
    assert p._cache_ctr['check'] == 4
    assert p._cache_ctr['hits'] == 1
    assert p._cache_ctr['misses'] == 3

    # Check what numbers are in the cache.
    # In order to not have to reconstruct the binary representation,
    # we just check the length - the first call had three elements and
    # and that one should have been dropped.
    for val in p._cache.values():
        assert len(val) == 4

    p.cache_clear()

    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0



def test_model_keyword_cache():
    """Check what happens with the cache and keywords"""

    mdl = ReportKeywordsModel()

    store = mdl._keyword_store
    assert len(store) == 0

    x = [1, 3, 4]

    # We do not care about the returned value, we just
    # want to see what happens with the model.
    #
    mdl(x, user_arg=1)
    assert len(store) == 1
    assert store[0][1] == {"user_arg": 1}

    mdl(x, user_arg=1)
    assert len(store) == 1
    assert store[0][1] == {"user_arg": 1}

    mdl(x, user_arg=2)
    assert len(store) == 2
    assert store[1][1] == {"user_arg": 2}

    # Explicit check what happens when we turn off the cache
    mdl.cache = 0
    store.clear()
    mdl(x, user_arg=3)
    assert len(store) == 1
    assert store[0][1] == {"user_arg": 3}

    mdl(x, user_arg=4)
    assert len(store) == 2
    assert store[1][1] == {"user_arg": 4}


def test_model1d_cache_xhi_positional():
    """Do we cache the model when xhi argument is not named?"""

    xlo = [1, 2, 3]
    xhi = [2, 2.5, 5]
    mdl = ReportKeywordsModel()
    y1 = mdl(xlo, xhi)

    store = mdl._keyword_store
    assert len(store) == 1

    y2 = mdl(xlo, xhi)
    assert len(store) == 1
    assert y2 == pytest.approx(y1)
    assert y1 == pytest.approx(10, 5, 20)

    # Now change xhi so we can see if it is used in the cache?
    xhi = [x + 1 for x in xlo]
    y3 = mdl(xlo, xhi)
    assert len(store) == 2
    assert y3 == pytest.approx(10, 10, 10)


def test_model1d_cache_xhi_named():
    """Do we cache the model when xhi argument is named?"""

    xlo = [1, 2, 3]
    xhi = [2, 2.5, 5]
    mdl = ReportKeywordsModel()
    y1 = mdl(xlo, xhi=xhi)

    store = mdl._keyword_store
    assert len(store) == 1

    y2 = mdl(xlo, xhi=xhi)
    assert len(store) == 1
    assert y2 == pytest.approx(y1)
    assert y1 == pytest.approx(10, 5, 20)

    # Now change xhi so we can see if it is used in the cache?
    xhi = [x + 1 for x in xlo]
    y3 = mdl(xlo, xhi=xhi)
    assert len(store) == 2
    assert y3 == pytest.approx(10, 10, 10)
