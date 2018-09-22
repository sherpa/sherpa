from sherpa.models import model
from astropy.modeling.polynomial import Polynomial2D

__all__ = ('WrapPoly2D', )


class WrapPoly2D(model.ArithmeticModel):
    """A two-dimensional polynomial from AstroPy, restricted to degree=2.

    The model parameters (with the same meaning as the underlying
    AstroPy model) are:

    c0_0
    c1_0
    c2_0
    c0_1
    c0_2
    c1_1

    """

    def __init__(self, name='wrappoly2d'):
        self._actual = Polynomial2D(degree=2)
        self.c0_0 = model.Parameter(name, 'c0_0', 0)
        self.c1_0 = model.Parameter(name, 'c1_0', 0)
        self.c2_0 = model.Parameter(name, 'c2_0', 0)
        self.c0_1 = model.Parameter(name, 'c0_1', 0)
        self.c0_2 = model.Parameter(name, 'c0_2', 0)
        self.c1_1 = model.Parameter(name, 'c1_1', 0)

        model.ArithmeticModel.__init__(self, name,
                                       (self.c0_0, self.c1_0, self.c2_0,
                                        self.c0_1, self.c0_2, self.c1_1))

    def calc(self, pars, x0, x1, *args, **kwargs):
        """Evaluate the model"""

        # This does not support 2D integrated data sets
        mdl = self._actual
        for n in ['c0_0', 'c1_0', 'c2_0', 'c0_1', 'c0_2', 'c1_1']:
            pval = getattr(self, n).val
            getattr(mdl, n).value = pval

        return mdl(x0, x1)
