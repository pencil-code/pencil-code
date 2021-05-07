#!/usr/bin/env python3
#=======================================================================
# find.py
#
# Facilities for analyzing the Pencil Code data.
#
# Author: Chao-Chin Yang
# Created: 2013-10-21
# Last Modified: 2021-05-07
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
def par_disp(datadir="./data", save_to=None):
    """Finds the displacement of each particle as a function of time.

    Keyword Arguments:
        datadir
            Path to the data directory.
        save_to
            If not None, a string of the filename (without .npz
            extension) to save the results (t and dxp) to a numpy data
            file under datadir.

    Returned Values:
        t
            A numpy array of times, starting from zero.
        dxp, dyp, dzp
            A tuple of three 2D numpy arrays, where dxp[i,j], dyp[i,j],
            and dzp[i,j] are the components of the displacement of
            particle j at t[i].  A component is returned None if the
            dimension is neither active nor periodic.
    """
    # Author: Chao-Chin Yang
    # Created: 2018-01-16
    # Last Modified: 2021-05-07
    import numpy as np
    from . import read

    # Safe factor for detecting jumps.
    SAFE = 1.1

    # Determine numbers of particles and snapshots.
    pdim = read.pardim(datadir=datadir)
    pvarfiles, t = read.varlist(datadir=datadir, listname="pvarN.list")
    npar = pdim.npar
    nt = len(pvarfiles)
    print(f"Total number of particles: {npar}")
    print(f"Total number of snapshots: {nt}")

    # Determine which dimension(s) should be considered.
    par = read.parameters(datadir=datadir)
    dim = read.dimensions(datadir=datadir)
    active = (dim.nxgrid > 1 and par.lperi[0],
              dim.nygrid > 1 and par.lperi[1],
              dim.nzgrid > 1 and par.lperi[2])
    print(f"Active dimensions: {active}")

    # Allocate memory.
    t = np.empty(nt,)
    dxp = np.empty((nt,npar)) if active[0] else None
    dyp = np.empty((nt,npar)) if active[1] else None
    dzp = np.empty((nt,npar)) if active[2] else None
    vpxmin, vpxmax = float("inf"), -float("inf")
    vpymin, vpymax = float("inf"), -float("inf")
    vpzmin, vpzmax = float("inf"), -float("inf")

    for i, pvarfile in enumerate(pvarfiles):
        # Read in the particle data.
        print("\rReading PVAR files ({:6.1%})......".format((i+1)/nt),
              end='', flush=True)
        fp = read.pvar(datadir=datadir, pvarfile=pvarfile, verbose=False)

        # Record time, positions, and monitor velocity extrema.
        t[i] = fp.t
        if active[0]:
            dxp[i] = fp.xp
            vpxmin, vpxmax = min(vpxmin, min(fp.vpx)), max(vpxmax, max(fp.vpx))
        if active[1]:
            dyp[i] = fp.yp
            vpymin, vpymax = min(vpymin, min(fp.vpy)), max(vpymax, max(fp.vpy))
        if active[2]:
            dzp[i] = fp.zp
            vpzmin, vpzmax = min(vpzmin, min(fp.vpz)), max(vpzmax, max(fp.vpz))
    print("Done. ")

    # Check the dimensions.
    dtmax = max(t[1:] - t[:-1])
    if ((active[0] and SAFE * (vpxmax - vpxmin) * dtmax > par.lxyz[0]) or
        (active[1] and SAFE * (vpymax - vpymin) * dtmax > par.lxyz[1]) or
        (active[2] and SAFE * (vpzmax - vpzmin) * dtmax > par.lxyz[2])):
        print(f"dtmax = {dtmax}")
        if active[0]: print("vpxmin, vpxmax = ", vpxmin, vpxmax)
        if active[1]: print("vpymin, vpymax = ", vpymin, vpymax)
        if active[2]: print("vpzmin, vpzmax = ", vpzmin, vpzmax)
        print("Warning: Boundary jumping may not be properly detected. ")

    # Remove duplicate data.
    print("Removing duplicate data......", end='', flush=True)
    t, indices = np.unique(t, return_index=True)
    if active[0]: dxp = dxp[indices]
    if active[1]: dyp = dyp[indices]
    if active[2]: dzp = dzp[indices]
    nt = len(t)
    print("Done. ")

    # Find the displacement.
    print("Computing the displacement......", end='', flush=True)
    if active[0]:
        dx = dxp[1:] - dxp[:-1]
        dxp -= dxp[0]
    if active[1]:
        dy = dyp[1:] - dyp[:-1]
        dyp -= dyp[0]
    if active[2]:
        dz = dzp[1:] - dzp[:-1]
        dzp -= dzp[0]

    # Detect boundary jumps.
    for i in range(nt-1):
        dt = t[i+1] - t[i]
        if active[0]:
            if vpxmin < 0:
                dxp[i+1:, dx[i] > vpxmax * dt] -= par.lxyz[0]
            if vpxmax > 0:
                dxp[i+1:, dx[i] < vpxmin * dt] += par.lxyz[0]
        if active[1]:
            if vpymin < 0:
                dyp[i+1:, dy[i] > vpymax * dt] -= par.lxyz[1]
            if vpymax > 0:
                dyp[i+1:, dy[i] < vpymin * dt] += par.lxyz[1]
        if active[2]:
            if vpzmin < 0:
                dzp[i+1:, dz[i] > vpzmax * dt] -= par.lxyz[2]
            if vpzmax > 0:
                dzp[i+1:, dz[i] < vpzmin * dt] += par.lxyz[2]

    # Reset the starting time to zeo.
    t -= t[0]
    print("Done. ")

    # Save the results if requested.
    if save_to is not None:
        file = datadir + '/' + save_to
        print("Saving the results to " + file + ".npz......",
              end='', flush=True)
        kw = dict(t=t)
        if active[0]: kw["dxp"] = dxp
        if active[1]: kw["dyp"] = dyp
        if active[2]: kw["dzp"] = dzp
        np.savez(file, **kw)
        print("Done. ")

    return t, (dxp, dyp, dzp)
#=======================================================================
def eta_z(z, datadir="./data"):
    """Finds the vertical profile of the resistivity, if any.

    Positional Arguments:
        z
            Numpy array of vertical coordinates.

    Keyword Arguments:
        datadir
            Data directory.

    Returned values:
        Magnetic diffusivity at the corresponding z.
    """
    # Author: Chao-Chin Yang
    # Created: 2017-08-13
    # Last Modified: 2017-08-13
    from . import read
    from math import sqrt
    import numpy as np
    from scipy.special import erfc

    # Read the parameters.
    par = read.parameters(datadir=datadir, par2=True)

    # Initialize the profile.
    eta = np.zeros_like(z)

    # Process the resistivities.
    print("iresistivity = ", par.iresistivity)
    for res in par.iresistivity:
        if res == "zdep":
            if par.zdep_profile == "fs":
                # Fleming & Stone (2003, ApJ, 585, 908)
                if par.cs0 > 0 and par.omega > 0:
                    zoh = z / (sqrt(2) * par.cs0 / par.omega)
                else:
                    zoh = z
                eta += par.eta * np.exp(-0.5 * zoh**2 +
                                        0.25 * par.sigma_ratio * erfc(abs(zoh)))
            else:
                raise NotImplementedError("zdep_profile = " + par.zdep_profile)
            print("Processed 'zdep'. ")

    return eta
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
