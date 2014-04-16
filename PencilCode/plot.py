#! /usr/bin/env python3
# Last Modification: $Id$
#=======================================================================
# plot.py
#
# Facilities for plotting the Pencil Code data.
#
# Chao-Chin Yang, 2013-10-22
#=======================================================================
def avg1d(datadir='./data', plane='xy', tsize=1024, var=None, **kwargs):
    """Plots the space-time diagram of a 1D average.

    Keyword Arguments:
        datadir
            Name of the data directory.
        plane
            Plane of the average.
        tsize
            Number of regular time intervals.
        var
            Name of the variable; if None, first variable is used.
        **kwargs
            Sent to matplotlib.pyplot.imshow.
    """
    # Chao-Chin Yang, 2013-10-29

    # Check the plane of the average.
    if plane == 'xy':
        xlabel = '$z$'
        xdir = 2
    elif plane == 'xz':
        xlabel = '$y$'
        xdir = 1
    elif plane == 'yz':
        xlabel = '$x$'
        xdir = 0
    else:
        raise ValueError("Keyword plane only accepts 'xy', 'xz', or 'yz'. ")

    # Read the data.
    print("Reading 1D averages...")
    from . import read
    time, avg = read.avg1d(datadir=datadir, plane=plane, verbose=False)
    par = read.parameters(datadir=datadir)
    xmin, xmax = par.xyz0[xdir], par.xyz1[xdir]

    # Default variable name.
    if var is None:
        var = avg.dtype.names[0]

    # Interpolate the time series.
    print("Interpolating", var, "...")
    import numpy as np
    from scipy.interpolate import interp1d
    tmin, tmax = np.min(time), np.max(time)
    ns = avg.shape[1]
    t = np.linspace(tmin, tmax, tsize)
    a = np.empty((tsize, ns))
    for j in range(ns):
        a[:,j] = interp1d(time, avg[var][:,j])(t)

    # Plot the space-time diagram.
    print("Plotting...")
    import matplotlib.pyplot as plt
    img = plt.imshow(a, origin='bottom', extent=[xmin,xmax,tmin,tmax], aspect='auto', **kwargs)
    ax = plt.gca()
    ax.set_ylabel('$t$')
    ax.set_xlabel(xlabel)
    cb = plt.colorbar(img)
    cb.set_label(var)
    plt.show()

#=======================================================================
def time_series(datadir='./data', diagnostics='dt', trange=None, xlog=False, ylog=False):
    """Plots diagnostic variable(s) as a function of time.

    Keyword Arguments:
        datadir
            Name of the data directory.
        diagnostics
            (A list of) diagnostic variable(s).
        trange
            A tuple of (tmin, tmax) for the time range to be shown; if
            None, all time is shown.
        xlog
            A boolean value for turning on or off logarithmic scale in
            x axis.
        ylog
            A boolean value for turning on or off logarithmic scale in
            y axis.
    """
    # Chao-Chin Yang, 2014-04-07
    from . import read
    import matplotlib.pyplot as plt

    # Read the time series.
    ts = read.time_series(datadir=datadir)

    # Determine the time range.
    if trange is None:
        it = list(range(len(ts.t)))
    else:
        it = (trange[0] <= ts.t) & (ts.t <= trange[1])
        if not it.any():
            print("No data lies in the time range. ")
            it = list(range(len(ts.t)))

    # Determine the axis scales.
    if xlog and ylog:
        p = plt.loglog
    elif xlog:
        p = plt.semilogx
    elif ylog:
        p = plt.semilogy
    else:
        p = plt.plot

    # Plot the diagnostics.
    if type(diagnostics) is list:
        for diag in diagnostics:
            p(ts.t[it], ts[diag][it], label=diag)
    else:
        p(ts.t[it], ts[diagnostics][it], label=diagnostics)
    plt.minorticks_on()
    plt.xlabel('t')
    plt.legend(loc='best')
    plt.show()
