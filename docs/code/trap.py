import numpy as np

from sherpa.models import model

__all__ = ('Trap1D', )


def _trap1d(pars, x):
    """Evaluate the Trapezoid.

    Parameters
    ----------
    pars: sequence of 4 numbers
        The order is amplitude, center, width, and slope.
        These numbers are assumed to be valid (e.g. width
        is 0 or greater).
    x: sequence of numbers
        The grid on which to evaluate the model. It is expected
        to be a floating-point type.

    Returns
    -------
    y: sequence of numbers
        The model evaluated on the input grid.

    Notes
    -----
    This is based on the interface described at
    https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Trapezoid1D.html
    but implemented without looking at the code, so any errors
    are not due to AstroPy.
    """

    (amplitude, center, width, slope) = pars

    # There are five segments:
    #    xlo = center - width/2
    #    xhi = center + width/2
    #    x0  = xlo - amplitude/slope
    #    x1  = xhi + amplitude/slope
    #
    #    flat   xlo <= x < xhi
    #    slope  x0 <= x < xlo
    #           xhi <= x < x1
    #    zero   x < x0
    #           x >= x1
    #
    hwidth = width / 2.0
    dx = amplitude / slope
    xlo = center - hwidth
    xhi = center + hwidth
    x0 = xlo - dx
    x1 = xhi + dx

    out = np.zeros(x.size)
    out[(x >= xlo) & (x < xhi)] = amplitude

    idx = np.where((x >= x0) & (x < xlo))
    out[idx] = slope * x[idx] - slope * x0

    idx = np.where((x >= xhi) & (x < x1))
    out[idx] = - slope * x[idx] + slope * x1

    return out


class Trap1D(model.RegriddableModel1D):
    """A one-dimensional trapezoid.

    The model parameters are:

    ampl
        The amplitude of the central (flat) segment (zero or greater).
    center
        The center of the central segment.
    width
        The width of the central segment (zero or greater).
    slope
        The gradient of the slopes (zero or greater).

    """

    def __init__(self, name='trap1d'):
        self.ampl = model.Parameter(name, 'ampl', 1, min=0, hard_min=0)
        self.center = model.Parameter(name, 'center', 1)
        self.width = model.Parameter(name, 'width', 1, min=0, hard_min=0)
        self.slope = model.Parameter(name, 'slope', 1, min=0, hard_min=0)

        model.RegriddableModel1D.__init__(self, name,
                                          (self.ampl, self.center, self.width,
                                           self.slope))

    def calc(self, pars, x, *args, **kwargs):
        """Evaluate the model"""

        # If given an integrated data set, use the center of the bin
        if len(args) == 1:
            x = (x + args[0]) / 2

        return _trap1d(pars, x)
