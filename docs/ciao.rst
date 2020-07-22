
***************
Sherpa and CIAO
***************

The Sherpa package was developed by the
Chandra X-ray Center (:term:`CXC`)
as a general purpose fitting and modeling tool, with specializations
for handling X-ray Astronomy data. It is provided as part of the
:term:`CIAO` analysis package,
where the code is the same as that available
from the
`Sherpa GitHub page <https://github.com/sherpa/sherpa>`_,
with the following modifications:

* the I/O backend uses the CIAO library :term:`Crates` rather than
  :term:`Astropy`;

* a set of customized IPython routines are provided as part of
  CIAO that automatically loads Sherpa and adjusts the appearance
  of IPython (mainly changes to the prompt);

* and the CIAO version of Sherpa includes the optional XSPEC model
  library (:py:mod:`sherpa.astro.xspec`).

The online documentation provided for Sherpa as part of CIAO,
namely https://cxc.harvard.edu/sherpa/, can be used with the
standalone version of Sherpa, but note that the focus of this
documentation is the
:doc:`session-based API <ui/index>`
provided by the
:py:mod:`sherpa.astro.ui` and :py:mod:`sherpa.ui` modules.
These are wrappers around the Object-Oriented
interface described in this document, and  data management
and utility routines.
