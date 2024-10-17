***************
Simulating data
***************

.. plot::
   :context:
   :nofigs:

   # This is set up for us with doctest but not for sphinxext_plot_directive,
   # so we need to redo it. This means that the RTD builds for PRs
   # end up installing the latest version of sherpa/sherpa-test-data
   # to try and pick up the correct files, but if a change to that
   # repository has been made but not merged in then this documentation
   # may be out of data.
   #
   from sherpa.utils.testing import get_datadir
   data_dir = get_datadir()
   if data_dir is None:
       assert False, "data directory is not set"  # simple error log
   if not data_dir.endswith('/'):
       data_dir += '/'


A simple case
=============

Simulating a data set normally involves:

 1. evaluate the model
 2. add in noise

This may need to be repeated several times for complex models, such
as when different components have different noise models or the noise
needs to be added before evaluation by a component.

The model evaluation would be performed using the techniques
described in this section, and then the noise term can be
handled with :py:func:`sherpa.utils.random.poisson_noise` or routines from
`NumPy <https://numpy.org/doc/stable/reference/random/index.html>`_
or SciPy to evaluate noise.

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> from matplotlib import pyplot as plt
   >>> import numpy as np
   >>> rng = np.random.default_rng(235)  # a "repeatable" set of random values
   >>> from sherpa.utils.random import poisson_noise

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> from sherpa.models.basic import Polynom1D
   >>> mdl = Polynom1D('mdl')
   >>> mdl.offset = 35
   >>> mdl.c1 = 0.5
   >>> mdl.c2 = 0.12
   >>> print(mdl)
   mdl
      Param        Type          Value          Min          Max      Units
      -----        ----          -----          ---          ---      -----
      mdl.c0       thawed            1 -3.40282e+38  3.40282e+38
      mdl.c1       frozen          0.5 -3.40282e+38  3.40282e+38
      mdl.c2       frozen         0.12 -3.40282e+38  3.40282e+38
      mdl.c3       frozen            0 -3.40282e+38  3.40282e+38
      mdl.c4       frozen            0 -3.40282e+38  3.40282e+38
      mdl.c5       frozen            0 -3.40282e+38  3.40282e+38
      mdl.c6       frozen            0 -3.40282e+38  3.40282e+38
      mdl.c7       frozen            0 -3.40282e+38  3.40282e+38
      mdl.c8       frozen            0 -3.40282e+38  3.40282e+38
      mdl.offset   frozen           35 -3.40282e+38  3.40282e+38

.. plot::
   :context:
   :include-source:

   >>> x = np.arange(10, 100, 12)
   >>> ymdl = mdl(x)
   >>> plt.plot(x, ymdl, label="model")
   [<matplotlib.lines.Line2D object at ...>]
   >>> plt.legend()
   <matplotlib.legend.Legend object at ...>

.. plot::
   :context:
   :include-source:

   >>> ypoisson = poisson_noise(ymdl, rng=rng)
   >>> plt.clf()
   >>> plt.plot(x, ymdl, label="model")
   [<matplotlib.lines.Line2D object at ...>]
   >>> plt.plot(x, ypoisson, ".", label="with noise")
   [<matplotlib.lines.Line2D object at ...>]
   >>> plt.legend()
   <matplotlib.legend.Legend object at ...>

.. _data_pha_fake:

X-ray data (:py:mod:`~sherpa.astro.data.DataPHA`)
=================================================

.. note::
   The interface to :py:func:`~sherpa.astro.fake.fake_pha` changed
   significantly in the Sherpa 4.16.1 release and a number of options
   previously supported now generate a warning message.

In principle, the same steps apply when simulating :term:`PHA` data
(:py:class:`~sherpa.astro.data.DataPHA` objects), however, the mechanics are a
little more complicated because we need to account for the
instrumental response (:term:`ARF` and :term:`RMF`) and possibly also
the background, which may contribute to the source spectrum that we
want to simulate.  Readers not interested in X-ray data analysis may
want to skip this section.

Sherpa offers a dedicated function :py:func:`sherpa.astro.fake.fake_pha` for
simulations of :term:`PHA` data. A :py:class:`~sherpa.astro.data.DataPHA` object
needs to be set up with the responses for the source and an exposure
time. So, we first create a :py:class:`~sherpa.astro.data.DataPHA` object with
the correct exposure time, but we can leave the settings for channel
and counts empty, because these will be filled in by the simulation:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> from sherpa.astro.data import DataPHA
   >>> from sherpa.astro.io import read_arf, read_rmf, read_pha
   >>> data = DataPHA(name='simulation', channel=None, counts=None, exposure=10000.)
   >>> data.set_arf(read_arf(data_dir + '9774.arf'))
   >>> data.set_rmf(read_rmf(data_dir + '9774.rmf'))
   >>> data.channel = np.arange(1, 1025, dtype=np.int16)

Alternatively, one could read in a PHA file (``data =
read_pha('9774_bg.pi')``). In this case, the response and backgrounds
will be automatically loaded, if the relevant header keywords are
set. Also, the exposure time, background scaling etc. will be taken
from the header of the file. When new data is simulated later, the
counts in ``data`` will be overwritten, but all other information
stays the same.

Next, we set up a model. In this case, we start with a powerlaw source
where the slope and normalization of that powerlaw are already known,
e.g. from the literature. We then add a weak emission line. Our
simulation will show us if this emission line would be detectable in
a real observation:

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> from sherpa.models.basic import PowLaw1D, Gauss1D
   >>> pl = PowLaw1D()
   >>> line = Gauss1D()
   >>> pl.gamma = 1.8
   >>> pl.ampl = 2e-05
   >>> line.pos = 6.7
   >>> line.ampl = .0003
   >>> line.fwhm = .1
   >>> srcmdl = pl + line
   >>> print(srcmdl)
   powlaw1d + gauss1d
      Param        Type          Value          Min          Max      Units
      -----        ----          -----          ---          ---      -----
      powlaw1d.gamma thawed          1.8          -10           10
      powlaw1d.ref frozen            1 -3.40282e+38  3.40282e+38
      powlaw1d.ampl thawed        2e-05            0  3.40282e+38
      gauss1d.fwhm thawed          0.1  1.17549e-38  3.40282e+38
      gauss1d.pos  thawed          6.7 -3.40282e+38  3.40282e+38
      gauss1d.ampl thawed       0.0003 -3.40282e+38  3.40282e+38

The model we use with :py:func:`~sherpa.astro.fake.fake_pha` **must**
include the instrument response, which we can access with the
:py:meth:`~sherpa.astro.data.DataPHA.get_full_response` method
call.

.. plot::
   :context:
   :nofigs:
   :include-source:

   >>> resp = data.get_full_response()
   >>> instmdl = resp(srcmdl)
   >>> print(instmdl)
   apply_rmf(apply_arf(10000.0 * (powlaw1d + gauss1d)))
      Param        Type          Value          Min          Max      Units
      -----        ----          -----          ---          ---      -----
      powlaw1d.gamma thawed          1.8          -10           10
      powlaw1d.ref frozen            1 -3.40282e+38  3.40282e+38
      powlaw1d.ampl thawed        2e-05            0  3.40282e+38
      gauss1d.fwhm thawed          0.1  1.17549e-38  3.40282e+38
      gauss1d.pos  thawed          6.7 -3.40282e+38  3.40282e+38
      gauss1d.ampl thawed       0.0003 -3.40282e+38  3.40282e+38


The simplest case: Simulate the source spectrum only
----------------------------------------------------

With this model, it is now easy to run the simulation, which will
calculate the expected number of counts in each spectral channel
(where the number and width of the channels is given by the responses)
and then draw from a Poisson distribution with this expected
number. Thus, the simulated number of counts in each channel is always
an integer and includes Poisson noise - running
:py:func:`~sherpa.astro.fake.fake_pha` twice with identical settings will give
slightly different answers. With default settings, the input model is
convolved with the RMF and multiplied by the ARF, and properly scaled
for the exposure time.

.. plot::
   :context:
   :include-source:

   >>> from sherpa.astro.fake import fake_pha
   >>> from sherpa.astro.plot import DataPHAPlot
   >>> fake_pha(data, instmdl)  # NOTE: include the instrument response
   >>> dplot = DataPHAPlot()
   >>> dplot.prepare(data)
   >>> dplot.plot()

We bin the counts into bins of at least 5 counts per bin and display
an image of the simulated spectrum (see :ref:`model_evaluate_example_pha`
for details):

.. plot::
   :context:
   :include-source:

   >>> data.set_analysis('energy')
   >>> data.notice(0.3, 8)
   >>> data.group_counts(5)
   >>> dplot.prepare(data)
   >>> dplot.plot(xlog=True, ylog=True)

We can overplot the model (both grouped and ungrouped):

.. plot::
   :context:
   :include-source:

   >>> from sherpa.astro.plot import ModelPHAHistogram, ModelHistogram
   >>> dplot.plot(xlog=True, ylog=True)
   >>> m1plot = ModelPHAHistogram()
   >>> m1plot.prepare(data, instmdl)
   >>> m1plot.overplot(alpha=0.5, label="Grouped model")
   >>> m2plot = ModelHistogram()
   >>> m2plot.prepare(data, instmdl)
   >>> m2plot.overplot(alpha=0.5, label="Ungrouped model")

Sometimes, more than one response is needed, e.g. in Chandra LETG/HRC-S
different orders of the grating overlap on the detector, so they all
contribute to the observed source
spectrum. :py:func:`~sherpa.astro.fake.fake_pha` works if one or more responses
are set for the input :py:class:`~sherpa.astro.data.DataPHA` object.

It is also possible that the input model already includes the
appropriate responses and scalings. In this case, ``is_source=False``
can be set to indicate to sherpa that the model is not a naked source
model, but includes all required instrumental effects. In this
way, :py:func:`~sherpa.astro.fake.fake_pha` can be used with arbitrarily
complex models which may include such components as the instrumental
background (which should not be folded through the ARF) or arbitrary
other components::

   >>> fake_pha(data, model=inst_bkg + my_arf(my_rmf(srcmodel)), is_source=False)  # doctest: +SKIP

.. image:: ../_static/evaluation/fake_pha.png

.. todo::
   Can we create the model "in line" so we can create the image rather
   than hard code it.


Adding background
-----------------

Weak spectral features can be hidden in a large (instrumental or
astrophysical) background, thus it can be important to include
background in the simulation.

Sample background from a PHA file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One way to include background is to sample it from a
:py:class:`~sherpa.astro.data.DataPHA` object. To do so, a background need to be
loaded into the dataset before running the simulation and, if not done
before, the scale of the background scaling has to be set:

.. plot::
   :context:
   :include-source:

   >>> bgdata = read_pha(data_dir + '9774_bg.pi')
   >>> data.set_background(bgdata)
   >>> data.backscal = 9.6e-06
   >>> fake_pha(data, instmdl, add_bkgs=True)
   >>> bplot = DataPHAPlot()
   >>> bplot.prepare(data)
   >>> dplot.plot(linestyle="solid", ylog=True, label="data")
   >>> bplot.overplot(linestyle="solid", alpha=0.5, label="with background")

The `fake_pha` function simulates the source spectrum as above, but
then it samples from the background PHA. For each bin, it treats the
background count number as the expected value and performs a Poisson
draw. The background drawn from the Poisson distribution is than added
to the simulated source spectrum (the sum of two Poisson distributions
is a Poisson distribution again). This works best if the background is
well exposed and has a large number or counts.

Why do we need to set the ``add_bkgs=True`` argument and do not simply
use all available backgrounds? Often is is
useful to read in a file with ``data = read_pha('9774.pi')``, which
might automatically read in the background as well. Using the
`add_bkgs` to switch the background on or off in the simulation
makes it easy to compare the same simulation with and without a
background.

Background models
^^^^^^^^^^^^^^^^^

When the number of counts in the background is low, the above
procedure amplifies the noise. The experimental background already
contains noise just from Poisson statistics. Using the observed count
numbers as background value and then performing a Poisson draw on them
again gives a higher level of noise than a real observation would
have. To avoid this problem, a background model may be used instead::

   >>> from sherpa.models.basic import Const1D
   >>> bkgmdl = Const1D('bmdl')
   >>> bkgmdl.c0 = 2
   >>> fake_pha(data, mdl, add_bkgs=True, bkg_models={'1': bkgmdl})  # doctest: +SKIP

The keys in the `bkg_models` dictionary are the identifiers of the
backgrounds. Above, we loaded a background with the default identifier
(which is ``1``).

.. todo::
   Looks like we need `instmdl` here, but does this mean that
   bkg_models needs to be sent `resp(bkgmdl)` ?


More than one background
^^^^^^^^^^^^^^^^^^^^^^^^

If more than one background is set up, then the expected backgrounds
will be averaged before the Poisson draw. For some backgrounds that
expected value might be the value of the counts in a PHA file,
while for others the expected value might be calculated from a model.


Reference/API
=============

.. currentmodule:: sherpa.astro.fake

.. automodule:: sherpa.astro.fake

   .. rubric:: Functions

   .. autosummary::
      :toctree: api

      fake_pha
