
****************
Utility routines
****************

.. todo::

   This section needs looking at.
   Should they be documented here or elsewhere?

There are a number of utility routines provided by Sherpa that may
be useful. Unfortunately it is not always obvious whether a routine is for use
with the Object-Oriented API or the Session API.

.. _resource-file:

The Sherpa resource file
========================

The `sherpa.get_config` routine returns the location of the default
resource file used to configure Sherpa. The location of this file is
taken from (in order):

  1. the ``SHERPARC`` environment file,

  2. the file ``.sherpa.rc`` in the user's home directory,

  3. and the file returned by `~sherpa.get_config`.

.. _multi-core:

Taking advantage of multiple CPUs
=================================

The `sherpa.utils.parallel` module provides code to run code on
multiple CPUs. There are two settings in the Sherpa resource file
that can be changed:

- ``parallel.numcores``
- ``multiprocessing.multiprocessing_start_method``

If ``parallel.numcores`` is set to 1 then the code will run in serial
and not parallel.

The ``multiprocessing.multiprocessing_start_method`` field controls
how the parallel code is run. Using methods other than "fork" are
unlikely to provide any benefit for Sherpa, and it is probably better to set
``parallel.numcores`` to 1 if the "fork" method can not be used.
The `Python documentation on multiprocessing <https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods>`_
explains the methods that Python offers and discusses some pitfalls.

Sherpa relies on mutable state, in particular for handling parameter
values, and so is unlikely to benefit from using threads.

.. _verbosity:

Controlling the verbosity of Sherpa
===================================

Sherpa uses `Python logging
<https://docs.python.org/3/library/logging.html>`_ for most
messages. This allows the user to redirected the output to a file or
suppress it by setting the logging level. The following example will
globally change the level for all sherpa modules, such that debug and
informational messages are no longer displayed:

  >>> import logging
  >>> sherpalog = logging.getLogger('sherpa')
  >>> sherpalog.setLevel('WARNING')

Sherpa also provides a context manager -
:py:class:`~sherpa.utils.logging.SherpaVerbosity` - to change the
logging level only for a specific portion of the code. This can be
used, e.g., to hide the long default output printed after fitting a
model:

  >>> from sherpa.utils.logging import SherpaVerbosity
  >>> import numpy as np
  >>> from sherpa.astro import ui
  >>> ui.load_arrays("mydata", np.arange(5), np.ones(5))
  >>> ui.set_model("mydata", "polynom1d.poly")
  >>> with SherpaVerbosity('WARNING'):
  ...    ui.fit("mydata")


Random Numbers
==============

.. todo::

   Is this the best place for this information?

Sherpa uses random numbers

- in an optimisation routine (`~sherpa.optmethods.optfcts.montecarlo`
  used by `~sherpa.optmethods.MonCar`);
- when evaluating a MCMC chain such as `~sherpa.sim.mh.MetropolisMH`;
- estimating the likelihood ratio of two models
  (`sherpa.sim.simulate`);
- when sampling distributions using routines from
  `sherpa.sim.sample`;
- and  when simulating data, such as `~sherpa.ui.fake` and
  `~sherpa.astro.ui.fake_pha`.

The primary way to control the random numbers used by Sherpa is to
create a `NumPy Random Generator
<https://numpy.org/doc/stable/reference/random/generator.html>`_ and
pass it to routines either with the ``rng`` parameter or, for users of
the Session class, the `~sherpa.ui.utils.Session.set_rng` method.
However, for the optimiser code, which uses a C++ implementation
of the Mersenne Twister random number generator,
the important value is the ``seed`` argument (although this can be
derived from the random number generator if the ``seed`` value is
set to `None`).

.. note::

   Prior to Sherpa 4.16.0 the random numbers were controlled by a
   combination of the legacy NumPy random-number API, that is,
   calling `numpy.random.seed()` , and - in a few places, the
   `Python random module <https://docs.python.org/3/library/random.html>`_,
   as well as the ``seed`` argument for the optimiser code.


Reference/API
=============

.. toctree::
   :maxdepth: 2

   sherpa
   err
   logging
   parallel
   random
   guess
   utils
   testing
   io
   astro
   astro_io
   astro_io_types
   astro_io_wcs
   astro_io_xstable
   astro_utils
   astro_utils_xspec
