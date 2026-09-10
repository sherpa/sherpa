.. _Image-backends:

*********************
Image backend details
*********************

Image display in Sherpa is done through *image backends*, although
unlike the :ref:`plotting case <Plotting-backends>` there is currently
only one usable backend, that uses :term:`XPA` to connect to the :term:`DS9`
application.

Which backend is used?
======================

When the `sherpa.image` module is first imported, Sherpa tries to
import the `sherpa.image.ds9_backend` module, falling back to
`sherpa.image.dummy_backend` if DS9 and the XPA tools can not be
found. The logic for sending data to and from DS9, as well as
starting and stopping the DS9 application, it provided by the
``sherpa.image.DS9`` module.
