#! /usr/bin/env python3
# Last Modification: $Id$
#=======================================================================
# find.py
#
# Facilities for analyzing the Pencil Code data.
#
# Chao-Chin Yang, 2013-10-21
#=======================================================================
def avgt1d(tmin=None, **kwarg):
    """Finds the time average of the 1D averages.

    Keyword Arguments:
        tmin
            Minimum time to be included in the average.  If None, the
            whole time series is used.
        **kwarg
            Keyword arguments to be passed to read.avg1d().
    """
    # Chao-Chin Yang, 2014-10-27
    from . import read
    import numpy as np
    from scipy import integrate
    # Read the 1D averages.
    t, avg1d = read.avg1d(**kwarg)
    t, indices = np.unique(t, return_index=True)
    avg1d = avg1d[indices,:]
    # Check the time span.
    tmax = t.max()
    if tmin is None:
        tmin = t.min()
    if tmax <= tmin:
        print("The minimum time has not been reached. ")
        print("  tmin, tmax = ", tmin, tmax)
        tmin = t.min()
    indices = t >= tmin
    dtinv = 1 / (tmax - tmin)
    # Find the time average at each location for each variable.
    var = avg1d.dtype.names
    nvar = len(var)
    nz = avg1d.shape[1]
    avg = np.core.records.array(nvar * [np.zeros(nz,)], names=var)
    sd = np.core.records.array(nvar * [np.zeros(nz,)], names=var)
    for v in var:
        for j in range(nz):
            avg[v][j] = integrate.simps(avg1d[v][indices,j], t[indices])
            sd[v][j] = integrate.simps(avg1d[v][indices,j]**2, t[indices])
        avg[v] *= dtinv
        sd[v] = np.sqrt(dtinv * sd[v] - avg[v]**2)
    return avg, sd
#=======================================================================
def sigma0(datadir='./data'):
    """Returns the background gas column density inside the
    computational domain.

    Keyword Arguments:
        datadir
            Name of the data directory.
    """
    # Chao-Chin Yang, 2014-11-04
    from . import read
    from math import pi, sqrt, erf
    # Read the parameters.
    par = read.parameters(datadir=datadir)
    # Find the column density.
    if par.gztype in {'zero', 'none'}:
        s = par.rho0 * par.Lxyz[2]
    elif par.gztype == 'linear':
        h = par.cs0 / par.gz_coeff
        z0 = par.xyz0[2] / (sqrt(2) * h)
        z1 = par.xyz1[2] / (sqrt(2) * h)
        s = sqrt(0.5 * pi) * h * par.rho0 * (erf(z1) - erf(z0))
    return s
#=======================================================================
def stratz(datadir='./data', trim=False):
    """Finds the vertical background stratification.

    Returned Values:
        The z coordinates and the corresponding density stratification.

    Keyword Arguments:
        datadir
            Name of the data directory.
        trim
            Whether or not to trim ghost cells.
    """
    # Chao-Chin Yang, 2014-11-03
    from . import read
    import numpy as np
    # Read the dimensions and the parameters.
    dim = read.dimensions(datadir=datadir)
    par = read.parameters(datadir=datadir)
    z = read.grid(datadir=datadir).z
    if trim:
        z = z[dim.nghost:-dim.nghost]
    # Find the density stratification.
    if par.gztype in {'zero', 'none'}:
        rho = par.rho0 * np.ones(dim.nzgrid,)
    elif par.gztype == 'linear':
        h = par.cs0 / par.gz_coeff
        rho = par.rho0 * np.exp(-0.5 * (z / h)**2)
    return z, rho
#=======================================================================
def time_average(datadir='./data', diagnostics=None, tmin=0, verbose=True):
    """Finds the time average of each given diagnostic variable.

    Returned Values:
        The mean and standard deviation of each diagnostics for
        t >= tmin.

    Keyword Arguments:
        datadir
            Name of the data directory.
        diagnostics
            (A list of) diagnostic variable(s).  If None, all
            diagnostics but it and t are processed.
        tmin
            Starting time of the time average.
        verbose
            Directly print out the averages when True.
    """
    # Chao-Chin Yang, 2014-07-31
    from . import read
    from collections import namedtuple
    from numpy import sqrt
    from scipy import integrate
    # Read the time series.
    ts = read.time_series(datadir=datadir)
    # Check the time span.
    tmax = max(ts.t)
    if tmax <= tmin:
        print("The minimum time has not been reached. ")
        print("  tmax = ", tmax)
        print("  tmin = ", tmin)
        return None
    indices = ts.t >= tmin
    dtinv = 1 / (tmax - tmin)
    # Local function for conducting statistics.
    def stats(diag):
        mean = dtinv * integrate.simps(ts[diag][indices], ts.t[indices])
        stddev = sqrt(dtinv * integrate.simps((ts[diag][indices] - mean)**2, ts.t[indices]))
        if verbose:
            print("<", diag, "> = ", mean, "+/-", stddev)
        return mean, stddev
    # Default diagnostics is all except it and t.
    if diagnostics is None:
        diagnostics = list(ts.dtype.names)
        for name in ['it', 't']:
            try:
                diagnostics.remove(name)
            except ValueError:
                pass
    # Find and return the time average and its standard deviation of each diagnostics.
    if type(diagnostics) is str:
        return stats(diagnostics)
    else:  # assuming diagnostics is iterable.
        Mean = namedtuple('Mean', diagnostics)
        StdDev = namedtuple('StandardDeviation', diagnostics)
        mean, stddev = zip(*(stats(diag) for diag in diagnostics))
        return Mean(*mean), StdDev(*stddev)
