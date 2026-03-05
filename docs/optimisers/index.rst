.. _optimizers:

*******************************************************
Optimisers: How to improve the current parameter values
*******************************************************

The optimiser varies the model parameters in an attempt to find
the solution which minimises the chosen
:doc:`statistic <../statistics/index>`.

In general, it is expected that the optimiser will be used by
a :py:class:`~sherpa.fit.Fit` object to
:doc:`perform the fit <../fit/index>`, but
it can be used directly using the
:py:meth:`~sherpa.optmethods.OptMethod.fit` method. The optimiser
object allows configuration values to be changed which can
tweak the behavior; for instance the tolerance to determine whether
the fit has converged, the maximum number of iterations to use,
or how much information to display whilst optimising a model.

As an example, the default parameter values for the
:py:class:`Levenberg-Marquardt <sherpa.optmethods.LevMar>`
optimiser are::

    >>> from sherpa.optmethods import LevMar
    >>> lm = LevMar()
    >>> print(lm)
    name     = levmar
    ftol     = 1.1920928955078125e-07
    xtol     = 1.1920928955078125e-07
    gtol     = 1.1920928955078125e-07
    maxfev   = None
    epsfcn   = 1.1920928955078125e-07
    factor   = 100.0
    numcores = 1
    verbose  = 0

These settings are available both as fields of the object and via
the :py:attr:`~sherpa.optmethods.OptMethod.config` dictionary
field.

In addition to the built-in optimisers in `sherpa.optmethods.buildin`,
Sherpa includes interfaces to several of the optimizers in
`scipy.optimize` or :term:`optimagic` packages, which will be
available when :term:`SciPy` or :term:`optimagic`, respectively, are
installed in your Python environment.  Available optimisers are listed
in the `Reference/API` section below.

More optimisers can be built with the
:py:class:`sherpa.optmethods.OptMethod` class; a detailed example is given
in the documentation of that class. This can be used
to provide access to external packages such as
`CERN's MINUIT optimisation library <https://iminuit.readthedocs.io>`_.


Choosing an optimiser
=====================

Different optimisers have different strengths and weaknesses, and
there is no single best choice for all problems - that is why Sherpa
provides a range of them. Since optimisers are usually used with
`~sherpa.fit.Fit`, we discuss the strategies for choosing an
optimiser in practice in :doc:`the section on fitting <../fit/index>`.

Below is some background information on the mathematics of optimisation
and the optimisers implemented in Sherpa. A reader more interested in
running Sherpa can skip this section and go straight to the
:doc:`section on fitting <../fit/index>`.

For further reading, we also encourage to look at the very
accessible `Scipy User Guide <https://docs.scipy.org/doc/scipy/tutorial/optimize.html>`_ and
the `"explanation" chapters of the optimagic package <https://docs.scipy.org/doc/scipy/tutorial/optimize.html>`_.
For the optimisers in Sherpa, references to the original papers are given in the
documentation for the individual optimisers.


In general, the best assurance that the correct minimum has been found in a
particular calculation is careful examination of the nature of the solution
(e.g., by plotting a fitted function over data), and some confidence that the
full region that the minimum may lie in has been well searched by the algorithm
used, we will show examples in the :doc:`section on fitting <../fit/index>`.


Some points to take away from the discussions in the rest of this document.

  1. Never accept the result of a minimization using a single optimisation run;
     always test the minimum using a different method.

  2. Check that the result of the minimization does not have parameter values
     at the edges of the parameter space. If this happens, then the fit must be
     disregarded since the minimum lies outside the space that has been
     searched, or the minimization missed the minimum.

  3. Get a feel for the range of values of the target function (in Sherpa this
     is the fit statistic), and the stability of the solution, by starting the
     minimization from several different parameter values.

  4. Always check that the minimum "looks right" by visualizing the
     model and the data.

.. _optimizers-localglobal:

Background information on different types of optimisers
=======================================================

The following sections are adapted from a memo written by
Mark Birkinshaw (1998).

The minimization of mathematical functions is a difficult operation. A general
function :math:`f({\bf x})` of the vector argument :math:`\bf x` may
have many isolated local minima, non-isolated minimum hypersurfaces, or
even more complicated topologies. No finite minimization routine can
guarantee to locate the unique, global, minimum of :math:`f({\bf x})`
without being fed intimate knowledge about the function by the user.

This does not mean that minimization is a hopeless task. For many problems
there are techniques which will locate a local minimum which may be "close
enough" to the global minimum, and there are techniques which will find the
global minimum a large fraction of the time (in a probabilistic
sense). However, the reader should be aware of my philosophy is that there is
no "best" algorithm for finding the minimum of a general function. Instead,
Sherpa provides tools which will allow the user to look at the overall behavior
of the function and find plausible local minima, which will often contain
the physically-meaningful minimum in the types of problem with which Sherpa
deals.

Sherpa contains two types of routine for minimizing a fit statistic. I will
call them the "single-shot" routines, which start from a guessed set of
parameters, and then try to improve the parameters in a continuous fashion, and
the "scatter-shot" routines, which try to look at parameters over the entire
permitted hypervolume to see if there are better minima than near the starting
guessed set of parameters.

Single-shot techniques (or "local" techniques)
----------------------------------------------

As the reader might expect, the single-shot routines are relatively quick, but
depend critically on the guessed initial parameter values :math:`{\bf x}_0`
being near (in some sense) to the minimum :math:`{\bf x}_{\rm min}`. All the
single-shot routines investigate the local behaviour of the function near
:math:`{\bf x}_0`, and then make a guess at the best direction and distance to
move to find a better minimum. After testing at the new point, they accept that
point as the next guess, :math:`{\bf x}_1`, if the fit statistic is smaller
than at the first point, and modify the search procedure if it isn't
smaller. The routines continue to run until one of the following occurs:

  1. all search directions result in an increased value of the fit statistic
     (i.e., the search has converged to a local or global minimum);
  2. an excessive number of steps have been taken; or
  3. something strange happens to the fit statistic (e.g., it turns out to be
     discontinuous in some horrible way).

This description indicates that for the single-shot routines, there is a
considerable emphasis on the initial search position, :math:`{\bf x}_0`, being
reasonable. It may also be apparent that the values of these parameters should
be moderate; neither too small (:math:`10^{-12}`, say), nor too large
(:math:`10^{12}`, say). This is because the initial choice of step size in
moving from :math:`{\bf x}_0` towards the next improved set of parameters,
:math:`{\bf x}_1`, is based on the change in the fit statistic, :math:`f({\bf
x})` as components of :math:`{\bf x}` are varied by amounts :math:`{\cal
O}(1)`.  If :math:`f` varies little as :math:`{\bf x}` is varied by this
amount, then the calculation of the distance to move to reach the next root may
be inaccurate. On the other hand, if :math:`f` has a lot of structure (several
maxima and minima) as :math:`{\bf x}` is varied by the initial step size, then
these single-shot minimizers may mistakenly jump entirely over the
"interesting" region of parameter space.

These considerations suggest that the user should arrange that the search
vector is scaled so that the range of parameter space to be searched is neither
too large nor too small. To take a concrete example, it would not be a good
idea to have :math:`x_7` parameterize the Hydrogen column density
(:math:`N_{\rm H}`) in a spectral fit, with an initial guess of
:math:`10^{20}\ {\rm cm}^{-2}`, and a search range
(in units of :math:`{\rm cm}^{-2}`) of
:math:`10^{16}` to :math:`10^{24}`. The minimizers will look for variations in
the fit statistic as :math:`N_{\rm H}` is varied by
:math:`1\ {\rm cm}^{-2}`, and
finding none (to the rounding accuracy likely for the code), will conclude that
:math:`x_7` is close to being a null parameter and can be ignored in the
fitting. It would be much better to have :math:`x_7 = \log_{10}(N_{\rm H})`,
with a search range of 16 to 24. Significant variations in the fit statistic
will occur as :math:`x_7` is varied by :math:`\pm 1`, and the code has a
reasonable chance of finding a useful solution.

Bearing this in mind, the single-shot minimizers in Sherpa are listed below:

:py:class:`~sherpa.optmethods.NelderMead`
  This technique - also known as Simplex - creates a polyhedral search element
  around the initial position, :math:`{\bf x}_0`, and then grows or shrinks in
  particular directions while crawling around parameter space, to try to place
  a minimum within the final search polyhedron. This technique has some
  hilarious ways of getting stuck in high-dimension parameter spaces (where the
  polyhedron can become a strange shape), but is very good at finding minima in
  regions where the fit statistic has a moderately well-defined topology. Since
  it works in a different way than Levenberg-Marquardt minimization, a good
  strategy is to combine both minimization to test whether an apparent minimum
  found by one technique is stable when searched by the other. I regard
  NelderMead searching as good in smooth and simple parameter spaces,
  particularly when looking at regions where the fit statistic depends on a
  parameter in a linear or parabolic fashion, and bad where surfaces of equal
  value of the fit statistic are complicated. In either case, it is essential
  that the initial size of the polyhedron (with sides of length 1 unit) is a
  smallish fraction of the search space.

:py:class:`Levenberg-Marquardt <sherpa.optmethods.LevMar>`
  This can be considered to be a censored maximum-gradients technique which,
  starting from a first guess, moves towards a minimum by finding a good
  direction in which to move, and calculating a sensible distance to go.
  Its principal drawback is that to calculate the distance to move it has to
  make some assumptions about how large a step size to take, and hence there is
  an implicit assumption that the search space is reasonably well scaled (to
  :math:`\pm 10` units in each of the search directions, say). It is also
  important that in finding these gradients, the steps do not miss a lot of
  important structure; i.e. there should not be too many subsidiary minima.
  The search directions and distances to move are based on the shape of the
  target function near the initial guessed minimum, :math:`{\bf x}_0`,
  with progressive movement towards the dominant local minimum. Since
  this technique uses information about the local curvature of the fit
  statistic as well as its local gradients, the approach tends to stabilize the
  result in some cases. I regard the techniques implemented in Sherpa as being
  good minimum-refiners for simple local topologies, since more assumptions
  about topology are made than in the NelderMead approach, but bad at
  finding global minima for target functions with complicated topologies.

Scatter-shot techniques (or "global" techniques)
------------------------------------------------

Although a bit ad hoc, these techniques attempt to locate a decent minimum over
the entire range of the search parameter space. Because they involve searching
a lot of the parameter space, they involve many function evaluations, and are
somewhere between quite slow and incredibly-tediously slow.

The routines are listed below:

:py:class:`~sherpa.optmethods.GridSearch`
  This routine simply searches a grid in each of the search parameters,
  where the spacing is uniform between the minimum and maximum
  value of each parameter.
  This often requires many function evaluations and takes a long time
  to run. Even if there are only seven parameters (e.g. two Gaussians with
  position, amplitude, and width plus a constant background) and we
  want to search each parameter in 10 steps, we already need 100000
  function evaluations, while a Levenberg-Marquardt optimiser might find
  that same minimum in a few hundred evaluations.

  The coarseness of the grid sets how precise a root will be found,
  and if the fit statistic has significant structure on a
  smaller scale, then the grid-searcher will miss it completely.

  On the other hand, GridSearch is very useful together with
  :py:class:`template models <sherpa.models.template.TemplateModel>`.

:py:class:`Monte Carlo <sherpa.optmethods.MonCar>`
  This is a simple population based, stochastic function minimizer. At
  each iteration it combines population vectors - each containing a set of
  parameter values - using a weighted difference. This optimiser can
  be used to find solutions to complex search spaces but is not guaranteed
  to find a global minimum. It is over-kill for relatively simple problems.


Reference/API
=============

.. toctree::
   :maxdepth: 2

   optmethods
   optscipy
   optoptimagic
   optfcts
