Fitting Poisson distributed data
================================

Data that originates in a `Poisson process <https://en.wikipedia.org/wiki/Poisson_point_process>`_
can be encountered in many areas of science where the data consists of
discrete countable events. Examples include the number of photons detected
in a given time interval, the number of particles detected in a given
area, or the number of stars in a given mass range.

In principle, fitting Poisson distributed data follows the same steps as fitting
continuously distributed data. However, there are a few pitfalls to keep in mind.

Do not subtract the background
------------------------------
While the sum of two Poisson distributions is also Poisson distributed, the difference
is not. Thus, we cannot simply subtract the background from the data and then fit the result.
Instead, the background should be modeled separately and included in the fit as an additional
component. This ensures that the Poisson nature of the data is preserved.



Pick an appropriate statistic
-----------------------------

While very few distributions in science are truly normally distributed, in practice
they are often close enough that a :math:`\chi^2` statistic can be used and plots can
display symmetric error bars. Neither is true for Poisson distributed data, in particular
when there are fewer than 25 counts or so in any bin.

Sherpa implements `sherpa.stats.Chi2Gehrels` which applies a
weighting scheme to allow a :math:`\chi^2`-like statistic to be used (details and
references are given in the documentation for that statistic).
In general, though, Poisson distributed data should be fit by minimizing the
likelihood. Two options are available in Sherpa:
`sherpa.stats.Cash` and `sherpa.stats.CStat` differ slightly in the equations used (details
in the documentation for each class).


Keep your models positive
-------------------------
Both `sherpa.stats.Cash` and `sherpa.stats.CStat` calculate the logarithm
of number of counts in each bin and the logarithm of the model prediction.
To prevent the code from crashing when the model prediction is zero or negative, Sherpa
replaces such values with a small positive number before taking the logarithm.
This way, there is a step in the statistic when one or more model bins are zero or below
and the optimizer will hopefully move back into the allowed parameter space. However, if
the optimizer starts in the wrong place, it might be stuck.

Example
-------
In this example, we will look at the distributions of count rates from
young stars in the star forming region IRAS 20050+2720. The data is from
[Günther et al. (2012)](https://ui.adsabs.harvard.edu/abs/2012AJ....144..101G/abstract).
For this example, we downloaded the data and binned it.
The numbers are hardcoded below for simplicity.
We ignore the fainter half of the objects because the sample is probably
incomplete in this range. The cumulative distribution of observed X-ray fluxes
is often modelled as a lognormal distribution, but here we want to try to model
the non-cumulative distribution using a Gaussian distribution.

.. plot::
    :include-source:
    :context:
    :nofigs:

    >>> import numpy as np
    >>> from sherpa.data import Data1DInt
    >>> from sherpa.models import Gauss1D
    >>> from sherpa.stats import CStat
    >>> from sherpa.optmethods import LevMar
    >>> from sherpa.fit import Fit
    >>>
    >>> hist = np.array([45, 50, 43, 30, 27, 32,  5,  8,  0,  4,  1,  0])
    >>> edges = np.array([-6.5, -6.3, -6.1, -5.9, -5.7, -5.5, -5.3, -5.1, -4.9, -4.7, -4.5,
    ...    -4.3, -4.1])
    >>> xrayflux = Data1DInt('xrayflux', edges[:-1], edges[1:], hist)

The data is clearly in the Poisson regime. The histogram has "number of stars"
in each bin, where the X-axis is the logarithm of the count rate. There are
even bins with zero counts.
Thus, we pick `~sherpa.stats.CStat` as the statistic and use the default
`~sherpa.optmethods.LevMar` optimizer which usually works well for count data.

.. plot::
    :include-source:
    :context:

    >>> from sherpa.plot import DataPlot, ModelPlot
    >>> xraydistribution = Gauss1D('xraydistribution')
    >>> xfit = Fit(data=xrayflux, model=xraydistribution, stat=CStat(), method=LevMar())
    >>> xfit.guess()
    >>> print(xfit.fit())
    datasets       = None
    itermethodname = none
    methodname     = levmar
    statname       = cstat
    succeeded      = True
    parnames       = ('xraydistribution.fwhm', 'xraydistribution.pos', 'xraydistribution.ampl')
    parvals        = (1.502103404881948, -6.254572554858249, 235.94188464394244)
    statval        = 24.122296228170086
    istatval       = 60.989781628740936
    dstatval       = 36.86748540057085
    numpoints      = 12
    dof            = 9
    qval           = 0.004112093257497959
    rstat          = 2.6802551364633427
    message        = successful termination
    nfev           = 25

    >>> dplot = DataPlot()
    >>> dplot.prepare(xfit.data)
    >>> mplot = ModelPlot()
    >>> mplot.prepare(xfit.data, xfit.model)
    >>> dplot.plot()
    >>> mplot.overplot()



