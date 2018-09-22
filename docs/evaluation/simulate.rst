***************
Simulating data
***************

.. todo::

   Is this section worth it? There is an example which the user will
   already have read which talks about adding noise to a model. I
   do not think there is much more to say, unless we want to through
   some domain-specific examples here (e.g. fake a grouped X-ray
   spectrum, simulate an image with different noise/response components).

Simulating a data set normally involves:

 1. evaluate the model
 2. add in noise

This may need to be repeated several times for complex models, such
as when different components have different noise models or the noise
needs to be added before evaluation by a component.

The model evaluation would be performed using the techniques
described in this section, and then the noise term can be
handled with :py:func:`sherpa.utils.poisson_noise` or routines from
NumPy or SciPy to evaluate noise, such as ``numpy.random.standard_normal``.

::

   >>> import numpy as np
   >>> from sherpa.models.basic import Polynom1D

   >>> np.random.seed(235)
   
   >>> x = np.arange(10, 100, 12)
   >>> mdl = Polynom1D('mdl')
   >>> mdl.offset = 35
   >>> mdl.c1 = 0.5
   >>> mdl.c2 = 0.12

   >>> ymdl = mdl(x)

   >>> from sherpa.utils import poisson_noise
   >>> ypoisson = poisson_noise(ymdl)

   >>> from numpy.random import standard_normal, normal
   >>> yconst = ymdl + standard_normal(ymdl.shape) * 10
   >>> ydata = ymdl + normal(scale=np.sqrt(ymdl))
