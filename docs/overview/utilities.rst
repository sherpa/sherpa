
****************
Utility routines
****************

.. todo::

   This section needs looking at.
   Should they be documented here or elsewhere?

There are a number of utility routines provided by Sherpa that may
be useful. Unfortunately it is not always obvious whether a routine is for use
with the Object-Oriented API or the Session API.

Contolling the verbosity of Sherpa
==================================

Sherpa uses `Python logging
<https://docs.python.org/3/library/logging.htm>`_ for most
messages. This allows the user to redirected the output to a file or
suppress it by setting the logging level. The following example will
globally change the level for all sherpa moduls, such that debug and
informational messages are no longer displayed:

  >>> import logging
  >>> sherpalog = logging.getLogger('sherpa')
  >>> sherpalog.setLevel('WARNING')

Sherpa also provides a context manager to change the logging level
only for a specific portion of the code. This can be used, e.g., to
hide the long default output printed after fitting a model:

  >>> from sherpa.utils.logging import SherpaVerbosity
  >>> import numpy as np
  >>> from sherpa.astro import ui
  >>> ui.load_arrays("mydata", np.arange(5), np.ones(5))
  >>> ui.set_model("mydata", "polynom1d.poly")
  >>> with SherpaVerbosity('WARNING'):
  ...    ui.fit("mydata")

 
Reference/API
=============

.. toctree::
   :maxdepth: 2

   sherpa
   err
   utils
   testing
   io
   astro_io
