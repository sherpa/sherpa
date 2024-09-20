************************************
The sherpa.plot.bokeh_backend module
************************************

:term:`bokeh` is a plotting library that can be used to create interactive
plots, which are rendered in a webbrowser. Those plots allow various options to zoom,
pan, and scroll out of the box. Many additional interactions can be added by the
user later, and we also plan to make Sherpa use more features of this library
in the future.

.. note::
    The :term:`bokeh` backend is new in Sherpa 4.16 and should be considered
    experimental.
    We appreciate `feedback <https://github.com/sherpa/sherpa/issues/new>`_
    on how it works for you, or what features are missing; and we anticipate
    changes to how this backend works in the future.


How to display and save a plot
------------------------------
Bokeh plots are rendered by a webbrowser and, by default, are displayed in a
new tab. The plot can be saved to a file by clicking on the save icon in the
toolbar at the top of the plot or by using `bokeh.io.save`.

Once displayed, bokeh has no direct control over the browser tab anymore,
so the plot cannot be updated automatically. When doing plots with bokeh,
that means that the user needs to call `bokeh.plotting.show` when all elements
have been added.


Jupyter notebooks
-----------------
Similar to e.g. the ``%matplotlib inline`` command in a notebook that specifies
how to display :term:`matplotlib` plots, some setup is required to display
:term:`bokeh` plots inline in a notebook.

Before a plot is created, the following commands need to be run::

    >>> from bokeh.io import output_notebook
    >>> output_notebook()

See `the bokeh documentation <https://docs.bokeh.org/en/latest/docs/user_guide/output/jupyter.html>`_
for details.

Python scripts and (I)Python terminals
--------------------------------------
Bokeh plots can either be displayed in a new tab in the default webbrowser or be
saved to a file. In both cases, a plot is not continuously updated.
Users must explicitly call `bokeh.plotting.show` to display the plot or
`bokeh.plotting.save` to save the plot to a html file. The reference to
the current plot is stored in the ``current_plot`` attribute of the backend object::

    >>> from bokeh.plotting import show
    >>> from sherpa import plot, data
    >>> from sherpa.data import Data1D
    >>> plot.set_backend('BokehBackend')
    >>> x1 = [100, 200, 600, 1200]
    >>> y1 = [2000, 2100, 1400, 3050]
    >>> d1 = data.Data1D('oned', x1, y1)
    >>> plot1 = plot.DataPlot()
    >>> plot1.prepare(d1)
    >>> plot1.plot()
    >>> show(plot.backend.current_fig)

The same is true for the UI interface::

    >>> from bokeh.plotting import show, save
    >>> from sherpa.astro import ui
    >>> from sherpa import plot
    >>> from sherpa.astro.data import DataPHA
    >>> ui.load_arrays(1, [1, 2, 3], [1, 2, 3], DataPHA)
    >>> ui.plot_data()
    >>> show(plot.backend.current_fig)
    >>> # or
    >>> save(plot.backend.current_fig, filename='plot.html')  # doctest: +SKIP

Suggestions on how to integrate that better into the Sherpa UI are welcome.


Differences to matplotlib
-------------------------
The :term:`bokeh` backend is intrinsically different from the :term:`matplotlib` backend,
but Sherpa will automatically translate some settings.
For example, if the color `'k'` (black) is passed into a sherpa plotting routine, it
will be internally translated to bokeh's representation of black (`'black'`).
However, this is limited to a subset of the available options that are commonly
used, see :ref:`backend-independent-plotting-options`.
A user can pass additional backend specific options to the plotting
routines, such as `marker="*"` in matplotlib. This will cause an error in
bokeh, as bokeh does not have a `"*"` marker; a similar look can be achieved
with `marker="star"`.

Similarly, in bokeh, text for labels can be formatted with LaTeX by
enclosing it in `'$$'` (e.g. `label=r"$$\alpha$$"`), while in matplotlib
it is done with `label=r"$\alpha$"`. All labels that sherpa itself sets will be set
correctly, but if a user wants to manually set a label, they need adjust it
for the backend that is currently active.

Limitations
-----------
:term:`bokeh` itself is quite different from :term:`matplotlib` in the way plots
are displayed in an external program (a webbrowser). Thus, the workflow is
different from what we might be used to.
See the `bokeh documentation on output options <https://docs.bokeh.org/en/latest/docs/user_guide/output.html>`_
for details. In particular:

  - Bokeh plots cannot easily be exported to png or pdf. They are always
    rendered in javascript. Bokeh does
    `offer a mechanism to export plots to static images <https://docs.bokeh.org/en/latest/docs/user_guide/output/export.html>`_,
    but it requires additional software and drivers.
  - For display, bokeh saves temporary files that are then automatically opened in a
    webbrowser. In some linux systems, browsers are not allowed to open
    files in the local ``\temp`` directory for security reasons, so the user
    has to manually change the location the bokeh plot is written to with
    `bokeh.io.output_file`.
  - Once displayed, bokeh no longer has control over the webbrowser.
    For plotting in scripts or when using virtual frame buffers, that might
    lead to those browser processes staying open and consuming resources.
    Use `bokeh.io.output_file` together with `bokeh.io.save` instead of
    `bokeh.io.show`.



.. currentmodule:: sherpa.plot.bokeh_backend

.. automodule:: sherpa.plot.bokeh_backend

   .. rubric:: Classes

   .. autosummary::
      :toctree: api

      BokehBackend

Class Inheritance Diagram
=========================

.. inheritance-diagram:: BokehBackend
   :parts: 1