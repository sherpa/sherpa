************************
Plotting backend details
************************

.. warning::
   The plotting backend API is currently being re-designed and should be considered
   experimental.

Plotting in Sherpa is done through *plotting backends*.
The term *backend* is used in two related ways in Sherpa. First, a *backend* is an
external package that performs the plotting such as :term:`matplotlib`, i.e. it creates a
canvas and puts color on paper or on the screen. The Sherpa plotting classes are
written in a backend-independent way; they work with any plotting backend that
is supported by Sherpa. We also use the term *backend* to refer to Sherpa
classes that connect the Sherpa plot objects with the plotting backend by
translating Sherpa options into backend-specific commands. The details of that
are explained below, but, in short, a backend class would translate the sherpa
backend-independent "plot a line in red from here to there" into backend
specific ``plt.plot(x,y, linecolor='r')`` (if the backend is :term:`matplotlib`).

Instead of utilizing Sherpa classes and commands to plot data, users can
alternatively just access the data in the Sherpa data objects (e.g. the x and y
values of a dataset) and perform plotting operations directly with any plotting
backend that is available to them. This method may be less convenient, but it
works for plotting backends not (yet) supported by Sherpa.

Which backend is used?
======================

When the `sherpa.plot` module is first imported, Sherpa tries to import the
backends installed with Sherpa in the order listed in the
``options.plot_pkg`` setting from the ``sherpa.rc`` startup file.
The first module that imports successfully is set as the active
backend. The following command prints the name and the location
on disk of that module::

   >>> from sherpa import plot
   >>> print(plot.backend.name)
   pylab

.. warning::
    Of course, you could be tempted to write::

      >>> from sherpa.plot import backend
      >>> print(backend.name)
      pylab

    However, ``backend`` is now a reference to the backend that
    was active when the import was done. If the backend is changed
    later, ``backend`` will still refer to the old backend. Thus,
    always use ``plot.backend`` to access the active backend.

Change the backend
------------------

.. warning::
   When the plotting backend is changed, new plotting objects must be
   initialized because objects like e.g. `sherpa.plot.DataPlot` keep a
   copy of all default parameters of the backend that was active when they
   are initialized. A different backend might not understand all those
   parameters.

   In particular, Sherpa's UI layer (which contains functions such as
   `sherpa.astro.ui.plot_data`) keeps a reference to specific plot objects.
   Changing the backend by itself will not change those references, and thus
   those function might not work with a new backend.

After the initial import, the backend can be changed by loading one of
the plotting backends shipped with sherpa (or any other module that
provides the same interface)::

  >>> from sherpa.plot import set_backend
  >>> set_backend('pylab')

Sherpa also provides a context manager to change the backend for just one plot::

  >>> from sherpa.plot import TemporaryPlottingBackend
  >>> with TemporaryPlottingBackend('pylab'):
  ...     x = [1,2,3,4]
  ...     # do some plotting with x

The list of available plotting backend names and the classes that implement them is::

  >>> from sherpa.plot.backends import PLOT_BACKENDS
  >>> print(PLOT_BACKENDS.keys())
  dict_keys(['BaseBackend', 'BasicBackend', 'IndepOnlyBackend', 'pylab', 'PylabErrorArea', 'BokehBackend'])

Which backends are available depends on which packages are installed in your Python
environment, e.g. the ``"pylab"`` backend requires :term:`matplotlib`.

.. _backend-independent-plotting-options:

Backend-independent plotting options
====================================

Sherpa defines a number of plotting options that can be used with any backend,
thus code that limits itself to those options can run independent of which
backend is set up in Sherpa. All default settings and most examples use only
backend-independent plotting options. The resulting plots will not look
identical in each backend, e.g. thickness of a line or the font type of an
annotation might differ, but they will convey the same information. For example,
"a blue dotted line" will generally appear blue and dotted, even though the
shade of blue or the size of the dots might differ between plotting backends.

In some cases, a plotting backend might not support all Sherpa plot options (for
example, a plotting backend might not have a command to change the line style).
In those rare cases, a setting might be ignored; but a setting from Sherpa's
list of backend independent values will never raise an error.

The names of plotting options and values that all Sherpa plotting backends
accept are chosen to match common :term:`matplotlib` settings which are familiar
to many scientific Python users, and to maintain backwards compatibility with
previous version of Sherpa. They offer limited choices, but those are sufficient
for most plots.

The following settings are accepted for all Sherpa plotting backends:

Colors
------

Colors are used for the color of lines, symbols, error bars etc.

- ``'b'`` (blue)
- ``'r'`` (red)
- ``'g'`` (green)
- ``'k'`` (black)
- ``'w'`` (white),
- ``'c'`` (cyan)
- ``'y'`` (yellow)
- ``'m``` (magenta)
- `None` (plotting backend default)

Line styles
------------

Because of legacy from both Chips and :term:`matplotlib` backends, line styles
can be specified in more than one form:

- ``'noline'`` or  ``''`` (empty string) or ``'None'`` (as a string)
  means that no line shall be plotted.
- ``'solid'`` or ``'-'``  for solid lines,
- `None`  for the default style - usually a solid line, too,
- ``'dot'`` or ``':'`` for dotted lines,
- ``'dash'`` or ``'--'`` (dashed) for dashed lines, and
- ``'dashdot'`` or ``'-.'`` for a dot-dashed line.

Markers
--------
- ``"None"`` as a string or ``""`` (empty string) means that no marker will be
  shown,
- ``"."`` shows dots,
-  ``"o"`` shows circles (filled or unfilled depending on the backend),
- ``"+"`` shows plus signs,
- ``"s"`` shows squares.

Additional backend-specific settings
====================================

Most plotting backends accept more than just the backend-independent options
listed above. For example, :term:`matplotlib` allows `many different ways to
specify colors to cover the entire RGB range
<https://matplotlib.org/stable/tutorials/colors/colors.html>`_, such as
``color=(.3, .4, .5, .2)`` or ``color='xkcd:eggshell'``. The Sherpa plotting
methods will pass any value to the underlying backend::

  >>> from sherpa.data import Data1D
  >>> from sherpa.plot import DataPlot
  >>> d = Data1D('example data', [1, 2, 3], [3, 2, 5])
  >>> dplot = DataPlot()
  >>> dplot.prepare(d)
  >>> dplot.plot(markerfacecolor='xkcd:eggshell')

will succeed with the `~sherpa.plot.pylab_backend.PylabBackend`
(activated with ``sherpa.plot.set_backend("pylab")``), but raise an
error if the active backend is `~sherpa.plot.bokeh_backend.BokehBackend`. In
contrast, ``color='k'`` is also not understood by bokeh natively, but because it is on
the backend-independent list, Sherpa will translate ``'k'`` to a form that bokeh
does understand (``'black'`` in this case).

Backends may also accept additional keywords to specify more plotting properties
such as the transparency of an element or a URL that is opened when clicking on
an element. Those can simply be passed to the Sherpa plotting command, which
will pass them through to the plotting backend:

  >>> from sherpa.data import Data1D
  >>> from sherpa.plot import DataPlot
  >>> d = Data1D('example data', [1, 2, 3], [3, 2, 5])
  >>> dplot = DataPlot()
  >>> dplot.prepare(d)
  >>> dplot.plot(marker='*')

Since Sherpa does not process those options itself, but just passes them on to
the underlying backend module, they are not documented here - see the
documentation of the specific plotting module for details. Also, they will fail and
raise an error if the plotting backend in use does not understand the ``url`` keyword.

In some cases, the Sherpa plotting commands create several visualization
elements at the same time (lines, symbols, error bars, axes, labels). This makes
using Sherpa classes convenient, but it also means that the plotting functions
do not offer options to customize each and every part. In general, the plotting
functions pass color, line style etc. to the elements that describes the data
(line, marker) and generate labels or axes grids using default settings. Backend
specific code can be used to change the properties of the current figure after
the Sherpa plotting.

Backend interface
=================

.. note::

   This section is mostly relevant for developers or advanced users who write new
   Sherpa plot classes or new backends.

This section describes the API that all Sherpa backends offer to explain how to
use it and why it was designed this way. See `sherpa.plot.backends.BaseBackend`
for a complete listing of the calling signature for each function.
The `~sherpa.plot.backends.BasicBackend` backend extends
`~sherpa.plot.backends.BaseBackend` by raising a warning message for
plotting functions (plot, image, histogram etc.) that are not implemented.
It is the base for any real functional backend, which will override those
methods, but offer useful user feedback for any method not provided.
This future-proofs any backend derived from this class: When sherpa adds new
functions to its backend definition, they will be added here with a warning
message. Thus, any backend derived from this class will always provide the
interface that sherpa requires from a plotting backend.


Plotting functions
------------------

Each backend shall support the plotting functions listed below, where "support"
means "has to provide these functions and accept a standard list of arguments
without crashing or raising an exception". We explicitly allow for backends that
implement some of these as a no-op, e.g. because the underlying plotting library
does not support 2D data. In that case, the backend would typically issue a
warning.

The plotting functions are not separated by "how things look on paper" (thus "plot" is
a long method that is responsible for points, lines, and error bars), but
by "what is the input data type":

- `~sherpa.plot.backend.BaseBackend.plot` (for scatter plots with marker style
  set, for line plots with line style set, and for error bars with ``xerr`` or
  ``yerr`` set to `True`); accepts (x, y) data with optional error bars in each
  dimension. Data can be scalar (for a single marker), or array-like.
  Note that ``x`` and ``y`` can also be `None`, which should create an empty plot.
- `~sherpa.plot.backend.BaseBackend.histo` (similar to plot, but with
  "histogram-style" lines); accepts (xlo, xhi, y) data with optional xerr, yerr.
- `~sherpa.plot.backend.BaseBackend.contour` for (x0, x1, z) data.

If called with empty data, a plotting function shall at a minimum create an
empty plot.

Annotations
-----------

Backends should also implement the following annotation functions. They do not
depend on the data plotted, but just annotate the plot, e.g. a
`~sherpa.plot.RatioPlot` shows the ratio between data and model and can use an
annotation to mark the ``ratio=1`` line.

- `~sherpa.plot.backend.BaseBackend.hline` (horizontal across the entire axes)
- `~sherpa.plot.backend.BaseBackend.vline` (vertical across the entire axes)

Other annotations (e.g. text labels) might be added to the API in the future.
For this reason new backends should inherit from
`~sherpa.plot.backend.BasicBackend`. Any function added to the API will be
implemented in `~sherpa.plot.backend.BasicBackend` as a no-op with a
warning to the user like "Feature XYZ is not available in your backend". That
way, all Sherpa plots can immediately make use of newly added functions without
breaking existing plotting backends; the worst that happens is that not all
annotation will be visible in every backend.

Return values
-------------

Sherpa does not expect a specific return argument from any plotting function,
but they are allowed to have return values if that is helpful for their internal
implementation, e.g. in the `~sherpa.plot.pylab_backend.PylabBackend` backend, plotting a line might return a
line object so that error bars plotted later can use ``line.color`` to match the
color of that line.

Creating plots and panels, clearing and overplotting
----------------------------------------------------
.. todo::
   Add details

Each of the plotting functions above accepts
the following arguments: title, xlabel, ylabel, xlog, ylog, overplot, clearwindow

Multi-panels plot can be set with `~sherpa.plot.backends.BaseBackend.set_subplot``
and `~sherpa.plot.backends.BaseBacken.set_jointplot` (from `sherpa.plot.backends.BaseBackend`).


Interaction with interactive plots in the UI
--------------------------------------------
Each backend also acts as a context manager: All plotting commands
in the UI are wrapped in a ``with`` statement like this::

    >>> with sherpa.plot.backend():  # doctest: +SKIP
    ...     plotobj.plot()  # doctest: +SKIP

That allows the backend to execute any finishing code after the plots are done,
e.g. to save the plot to disk or to display it in a window.

Other methods
--------------

Backends need to have a few more methods:

- ``as_html_XXX`` (where XXX is a plot type) that are used for interactive
  display in the notebook with ``_repr_html_``.  These functions take a plot
  object and return a html representation as a string.
- ``get_XXX_plot/hist_prefs`` (where XXX is a plot type) which returns a
  dictionary of preferences that is used for displaying this plot.
- `~sherpa.plot.backend.BasicBackend.get_latex_for_string` to format latex in strings.
- ``colorlist(n)`` generates a list of n distinct colors. In backends that only have a
  limited number of colors available, the list might repeat.


How do I implement a new backend?
=================================

Example for a modified backend
------------------------------
Since Sherpa backends are defined using inheritance, new backends can be
created to change the appearance of plots in ways that are beyond the
options offered as keyword arguments for the existing backends. As an example,
the `~sherpa.plot.pylab_area_backend.PylabErrorArea` backend changes the way that
x-errors are visualized in plots from simple error bars to a shaded area.

.. literalinclude:: ../../sherpa/plot/pylab_area_backend.py
   :language: python
   :lines: 38-99

.. todo::

   Put code example in here that uses it?

Testing backends and plotting code
----------------------------------
Currently, Sherpa does not employ pixel-level tests that compare a generated
image pixel-by-pixel to a reference image. While testing like that guarantees
that any and all changes are found, they are susceptible to failing for
reasons unrelated to Sherpa, such as minor changes in the default of
:term:`matplotlib`.

Instead, the plotting tests in Sherpa fall into the following categories:

  - Tests that do not depend on the output of a plotting package. Those tests
    could check e.g. the ``prepare`` stage of a plotting object.
    They work with any backend (since Sherpa has `sherpa.plot.backends.BasicBackend`
    there is always a backend). For those tests, no special
    setup is needed and they will just be run with whatever backend is loaded
    based on the installed packages and the ``.sherpa.rc`` file.
  - Tests that should work with *all* backends. In particular, the
    `sherpa.plot.backends.IndepOnlyBackend` does not perform any plotting,
    but it does raise a warning when plotting options that are not on the list of
    backend-independent parameters are used. Users will never use this backend,
    but running our tests with it ensures that we don't hardcode
    defaults specific to matplotlib into any of our objects.
    In principle, most plotting tests fall into this
    category, since most code needs to work for all backends. However, running every
    plotting test with every possible backend can make the runtime for the tests long.
    In practice, we thus do not run every test that we could for all backends.
    To run a test for all backends, add the `sherpa.conftest.all_plot_backends`
    `pytest fixture <https://docs.pytest.org/en/stable/explanation/fixtures.html>`_::

        >>> def test_this_plotting_feature(all_plot_backends):
        ...     x = [1,2,3,4]
        ...     y = [2,3,4,5]
        ...     # plot x and y in some way in your test.

    This fixture will execute the test several times with different active backends.
    Alternatively, sherpa also has a fixture that runs all backends *except* for the
    `sherpa.plot.backends.IndepOnlyBackend`: `plot_backends`. This is useful if a test
    makes use of options that are not on the Backend-independent list and needs to avoid the
    extra warning that `sherpa.plot.backends.IndepOnlyBackend` would emit.
  - Tests for a specific backend using code such as ``assert plt.gca().xlabel == "text")``.
    These tests can either use the "require" decorators from `sherpa.utils.testing` or, if parts
    of the test are useful for other backends as well, but skip specific statements if the wrong
    backend is active::

      >>> import numpy as np
      >>> from sherpa import plot
      >>> from sherpa.data import Data1D
      >>> def test_something_that_also_uses_matplotlib(all_plot_backends):
      ...     d = Data1D('x', np.asarray([2, 4, 10]), np.asarray([2, 4, 0]))
      ...     r = d._repr_html_()
      ...     if plot.backend.name == 'pylab':
      ...         assert f'<summary>{summary} Plot</summary>' in r
      ...     else:
      ...         assert f'<summary>{summary} Data (' in r
