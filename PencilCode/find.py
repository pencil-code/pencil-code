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
    # Chao-Chin Yang, 2014-10-29
    from . import read
    import numpy as np
    from scipy import integrate
    # Read the 1D averages.
    t, avg1d = read.avg1d(**kwarg)
    t, indices = np.unique(t, return_index=True)
    avg1d = avg1d[indices,:]
    # Check the time span.
    indices, tmin, tmax = _get_indices(t, tmin=tmin)
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
def midplane(**kwarg):
    """Finds horizontally-averaged mid-plane properties as a function of
    time.

    Returned Values:
        Time and mid-plane properties.

    Keyword Arguments:
        **kwarg
            Keyword arguments to be passed to read.avg1d(), except
            keyword plane.
    """
    # Author: Chao-Chin Yang
    # Created: 2014-11-27
    # Last Modified: 2017-06-08
    from . import read
    import numpy as np
    from scipy.interpolate import interp1d

    # Read the horizontal averages.
    dir = {"datadir": kwarg["datadir"]} if "datadir" in kwarg else {}
    t, avg = read.avg1d(plane='xy', **kwarg)
    g = read.grid(trim=True, **dir)

    # Find the mid-plane properties.
    print("Finding mid-plane properties...")
    var = avg.dtype.names
    nvar = len(var)
    nt = avg.shape[0]
    avg0 = np.rec.array(nvar * [np.zeros(nt,)], names=var)
    for v in var:
        for i in range(nt):
            avg0[v][i] = interp1d(g.z, avg[v][i,:])(0)

    return t, avg0
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
def slice(axis=2, var=None, z0=0, **kwarg):
    """Finds a slice normal to one of the axis direction in one
    snapshot.

    Keyword Arguments:
        axis
            Axis that is normal to the slice.
        z0
            Coordinate of the slice on the axis.
        var
            List of variables; all variables are retrieved if None.
        **kwarg
            Keyword arguments to be passed to read.var().

    Returned Values:
        Time of and fields on the slice.
    """
    # Author: Chao-Chin Yang
    # Created: 2017-04-30
    # Last Modified: 2017-06-08
    from . import read
    import numpy as np
    from scipy.interpolate import interp1d

    # Read the snapshot.
    f = read.var(compact=False, **kwarg)

    # Check the normal of the slice.
    if axis == 0:
        dim = f.y.size, f.z.size
        z = f.x
    elif axis == 1:
        dim = f.x.size, f.z.size
        z = f.y
    elif axis == 2:
        dim = f.x.size, f.y.size
        z = f.z
    else:
        raise ValueError("axis is out of range. ")

    # Get the names of the variables.
    if var is None:
        datadir = {"datadir": kwarg["datadir"]} if "datadir" in kwarg else {}
        var = read.varname(**datadir)
    else:
        if type(var) is str: var = [var]
        for v in var:
            if v not in dir(f):
                print("Variable " + v + " does not exist. ")
                var.remove(v)
        if len(var) == 0: return f.t, None

    # Interpolate onto the slice.
    s = np.rec.array(len(var) * [np.zeros(dim)], names=var)
    for v in var:
        s[v] = interp1d(z, getattr(f,v), axis=axis, kind="cubic")(z0)

    return f.t, s
#=======================================================================
def stratz(datadir='./data', par=None, trim=False):
    """Finds the vertical background stratification.

    Returned Values:
        The z coordinates and the corresponding density stratification.

    Keyword Arguments:
        datadir
            Name of the data directory.
        par
            If not None, parameters given by read.parameters().
        trim
            Whether or not to trim ghost cells.
    """
    # Author: Chao-Chin Yang
    # Created: 2014-10-08
    # Last Modified: 2017-03-02
    from . import read
    import numpy as np

    # Read the dimensions and the parameters.
    dim = read.dimensions(datadir=datadir)
    if par is None:
        par = read.parameters(datadir=datadir)
    z = read.grid(datadir=datadir).z
    if trim:
        z = z[dim.nghost:-dim.nghost]

    # Find the density stratification.
    if not par.lstratz or par.gztype in {'zero', 'none'}:
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
    # Author: Chao-Chin Yang
    # Created: 2013-10-21
    # Last Modified: 2017-01-20
    from . import read
    from collections import namedtuple
    from numpy import sqrt, unique
    from scipy import integrate

    # Read the time series.
    ts = read.time_series(datadir=datadir)
    t, indices = unique(ts.t, return_index=True)
    ts = ts[indices]

    # Check the time span.
    indices, tmin, tmax = _get_indices(t, tmin=tmin)
    dtinv = 1 / (tmax - tmin)

    # Define function for conducting statistics.
    def stats(diag):
        v = ts[diag][indices]

        if diag[0:2] == "dt" or diag[0:4] == "npar" or diag[0:4] == "nmig":
            mean = v.mean()
            stddev = v.std()
        else:
            t = ts.t[indices]
            mean = dtinv * integrate.simps(v, t)
            stddev = sqrt(dtinv * integrate.simps((v - mean)**2, t))

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

    # Find and return the time average and its standard deviation of each
    # diagnostics.
    if type(diagnostics) is str:
        return stats(diagnostics)
    else:  # assuming diagnostics is iterable.
        Mean = namedtuple('Mean', diagnostics)
        StdDev = namedtuple('StandardDeviation', diagnostics)
        mean, stddev = zip(*(stats(diag) for diag in diagnostics))
        return Mean(*mean), StdDev(*stddev)
#=======================================================================
def _get_indices(t, tmin=None, tmax=None):
    """Gets the indices in t that is in the range [tmin, tmax].

    Positional Arguments:
        t
            A numpy array.
        tmin
            Minimum value; no minimum if None.
        tmax
            Maximum value; no maximum if None.

    Returned Values:
        indices
            Indices in t that is in the range [tmin, tmax].
        tmin
            Real minimum in t[indices].
        tmax
            Real maximum in t[indices].
    """
    # Created: 2015-10-29
    # Author: Chao-Chin Yang
    if tmin is None: tmin = t.min()
    if tmax is None: tmax = t.max()
    if tmin >= tmax:
        print("tmin = ", tmin)
        print("tmax = ", tmax)
        raise ValueError("tmin >= tmax")
    indices = (tmin <= t) & (t <= tmax)
    tmin = t[indices].min()
    tmax = t[indices].max()
    return indices, tmin, tmax
