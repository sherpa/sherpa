************************************
The sherpa.plot.bokeh_backend module
************************************

The :term:`bokeh` plotting module is not integrated as tightly into
IPython and Jupyter notebooks as :term:`matplotlib` is, so there is a little
more setup required to make this work well.

.. note::
    The :term:`bokeh` backend is new in Sherpa 4.16 and should be considered
    experimental.
    We appreciate `feedback <https://github.com/sherpa/sherpa/issues/new>`_
    on how it works for you, or what features are missing; and we anticipte
    changes to how this backend works in the future.


How to display and save a plot
------------------------------
Bokeh plots are rendered by a webbrowser and, by default, are displayed in a
new tab. The plot can be saved to a file by clicking on the save icon in the
toolbar at the top of the plot or by using `bokeh.io.save`.

Once displayed, bokeh has no direct control over the browser tab any more,
so the plot cannot be updated automatically. In bokeh, that means that the
user needs to call `bokeh.plotting.show` when all elements have been added.
In Sherpa, this is done automatically in the UI, but users of the OO interface
need to call `bokeh.plotting.show` themselves.

In the UI, a new plot is shown after every plot command. That's pretty
automatic and works fine with jupyter notebook or on the command line.
However, it is annoying if you use ``overplot=True`` repeatedly,
because you'll get a new plot each time - and in the notebook, you end up with a
list of plots in the cell where the top one shows only what you plotted first,
and with each overplot, another plot is appended at the bottom that is
like the previous one, just with the overplotted data added.

To avoid it, in the OO interface, you need to call `bokeh.plotting.show` yourself
or activate the `bokeh.io.output_notebook` mode.


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
In the UI, plotting is wrapped into a context manager that displays the plot after
plotting is finished; users of Sherpa's OO interface
must explicitly call `bokeh.plotting.show` to display the plot. The reference to
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