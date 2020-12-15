# mp-interp #

## Mean-preserving interpolation ## 

This repository contains Fortran 90/95 code that implements four approaches for mean-preserving interpolation, for example, the interpolation of pseudo-daily values from monthly means, that when averaged, reproduce the monthly means.  (The often-used linear interpolation does not do so.)  The four methods compared here include those by:

- Epstein, E.S. (1991), On obtaining daily climatological values from monthly means, *J. Climate* 4:365-368;  
- Harzallah, A. (1995) The interpolation of data series using a constrained iterating technique, *Monthly Weather Review* 123:2251-2254;
- Killworth, P.D. (1996) Time interpolation of forcing fields in ocean models, *J. Physical Oceanography* 26:136-143; and
- Rymes, M.D. and D.R. Meyers (2001) Mean preserving alogrithm for smoothing interpolating averaged data, *Solar Energy* 71:225-231.

In the discussion here, the example cases are considered to be single- or multiple year time series of monthly data, from we seek mean-preserving daily values.  The approaches can be generalized to other data sets; for example, the Harzallah (1991) approach was illustrated using a time series of annual values, in which interpolated monthly values were sought.

The Epstein (1991) approach is based on a harmonic analysis (or harmonic regression) of, for example, a year's worth of monthly mean values, where the basic harmonic analysis is recast as the integration of the values over the year, which assures preservation of the interpolated values.  It is periodic, in the sense that the curve defined by the interpolated daily values "wraps around".  It calculated directly; no iterative fitting or improvement is involved.  Like all of these methods, it is prone to "overshooting" (which, in practice, is a necessary feature (as opposed to problem) for reproducing the means).  Because of this, when applied to a time series of many years of monthly data (a year at a time), small discontinuities in the interpolated daily data occur between years.  In the implementation here, these are "fixed" by smoothing over the end of one year and beginning of the following year.

The Harzallah (1995) approach involves fitting, for example, a cubic spline to, e.g. the monthly data.  Interpolated daily values using the original fit are not necessarily mean preserving.  The method involves finding the "residuals" between the original (e.g. monthly) means and the means of the interpolated (e.g. daily) values, and fitting a second spline to those values.  This procedure is iterated until the residuals are acceptably small, and then the interpolated values for a particular "target" day are taken as the sum of the values across iterations.  The method is not periodic, but can be made so by padding the beginning and end of a year  with data from the end and beginning of the year, respectively.  Year-to-year continuity can be had by padding the beginning of any particular year with data from the end of the previous year, and padding the end with data from the next year.

The Killworth (1996) approach is a direct, as opposed to iterative approach, and is non periodic.  It involves...
