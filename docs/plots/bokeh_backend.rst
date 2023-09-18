************************************
The sherpa.plot.bokeh_backend module
************************************

The :term:`bokeh` plotting module is not integrated as tightly into
IPython and Jupyter notebooks as :term:`matplotlib` is, so there is a little
more setup required to make this work well.

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
    >>> show(plot.backend.current_fig['all_axes'])

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