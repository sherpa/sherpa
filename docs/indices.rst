******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Glossary
--------

.. glossary::

   ARF
       The Ancillary Response Function used to describe the effective
       area curve of an X-ray telescope; that is, the area of the telescope
       and detector tabulated as a function of energy. The
       :term:`FITS` format used to represent ARFs is defined in
       the :term:`OGIP` Calibration Memo
       `CAL/GEN/02-002 <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html>`_.

   Astropy
       A community Python library for Astronomy:
       https://www.astropy.org/.

   CIAO
       The data reduction and analysis provided by the :term:`CXC`
       for users of the Chandra X-ray telescope. Sherpa is provided
       as part of CIAO, but can also be used separately. The
       CIAO system is available from https://cxc.harvard.edu/ciao/.

   Crates
       The Input/Output library provided as part of :term:`CIAO`.
       It provides read and write access to FITS data files, with
       speciality support for X-ray data formats that follow
       :term:`OGIP` standards (such as :term:`ARF` and :term:`RMF`
       files).

   CXC
       The `Chandra X-ray Center <https://cxc.harvard.edu/>`_.

   DS9
       An external image viewer designed to allow users to interact with
       gridded data sets (2D and 3D). Is is used by Sherpa to display
       image data, and is available from https://ds9.si.edu/. It uses
       the :term:`XPA` messaging system to communicate with external
       processes.

   FITS
       The `Flexible Image Transport System
       <https://en.wikipedia.org/wiki/FITS>`_ is a common data format in
       Astronomy, originally defined to represent imaging data from radio
       telescopes, but has since been extended to contain a mixture of
       imaging and tabular data. Information on the various standards
       related to the format are available from the
       `FITS documentation page <https://fits.gsfc.nasa.gov/fits_documentation.html>`_
       at :term:`HEASARC`.

   HEASARC
       NASA's High Energy Astrophysics Science Archive Research
       Center at Goddard Space Flight Center:
       https://heasarc.gsfc.nasa.gov/.

   matplotlib
       The matplotlib plotting package, which is documented at
       https://matplotlib.org/, is used to provide the plotting support
       in Sherpa.

   OGIP
       The Office for Guest Investigator Programs (OGIP) was a division of
       the Laboratory for High Energy Astrophysics at Goddard Space Flight
       Center. The activities of that group have now become the responsibility
       of the `HEASARC FITS Working Group (HFWG)
       <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/ofwg_intro.html>`_,
       and supports the use of high-energy astrophysics data through
       multimission standards and archive access. Of particular note for
       users of Sherpa are the standard documents produced by this group
       that define the data formats and standards used by high-energy
       Astrophysics missions.

   PHA
       The standard file format used to store astronomical X-ray
       spectral data. The format is defined as part of the
       :term:`OGIP` set of standards, in particular OGIP memos
       `OGIP/92-007
       <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html>`_
       and
       `OGIP/92-007a
       <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007a/ogip_92_007a.html>`_.
       Confusingly, PHA can also refer to the Pulse Height Amplitude (the
       amount of charge detected) of an event, which is one of the
       two channel types that can be found in a PHA format file.

   PSF
       The Point Spread Function. This represents the response of an
       imaging system to a delta function: e.g. what is the shape that
       a point source would produce when observed by a system. It is
       dependent on the optical design of the system but can also be
       influenced by other factors (e.g. for ground-based observatories
       the atmosphere can add additional blurring).

   RMF
       The Redistribution Matrix Function used to describe the response
       of an Astronomical X-ray detector. It is a matrix containing the
       probability of detecting a photon of a given energy at a
       given detector channel.  The :term:`FITS` format used to
       represent RMFs is defined in the
       :term:`OGIP` Calibration Memo
       `CAL/GEN/02-002 <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html>`_.

   WCS
       The phrase World Coordinate System for an Astronomical data set
       represents the mapping between the measured position on the detector
       and a "celestial" coordinate. The most common case is in providing
       a location on the sky (e.g. in
       `Equatorial
       <https://en.wikipedia.org/wiki/Equatorial_coordinate_system>`_
       or `Galactic <https://en.wikipedia.org/wiki/Galactic_coordinate_system>`_
       coordinates)
       for a given image pixel, but it can also be used to map between
       row on a spectrograph and the corresponding wavelength of light.

   XPA
       The `XPA messaging system
       <https://hea-www.harvard.edu/saord/xpa/>`_
       is used by :term:`DS9` to communicate
       with external programs. Sherpa uses this functionality to
       control DS9 - by sending it images to display and retriving
       any regions a used may have created on the image data.
       The command-line tools used for this commiunication may be
       available via the package manager for a particular
       operating system, such as
       `xpa-tools for Ubuntu
       <https://packages.ubuntu.com/xenial/xpa-tools>`_,
       or they can be
       `built from source <https://github.com/ericmandel/xpa>`_.

   XSPEC
       This can refer to either the X-ray Spectral fitting package,
       or the models from this package. XSPEC is distributed by
       :term:`HEASARC` and its home page is
       https://heasarc.gsfc.nasa.gov/xanadu/xspec/. Sherpa can be
       built with support for the
       `models from XSPEC
       <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixExternal.html>`_.

       Sherpa can be built to use XSPEC versions 12.12.1, 12.12.0, 12.11.1, 12.11.0,
       12.10.1 (patch level `a` or later), 12.10.0, 12.9.1, or 12.9.0.
