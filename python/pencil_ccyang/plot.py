#! /usr/bin/env python3
# Last Modification: $Id$
#=======================================================================
# plot.py
#
# Facilities for plotting the Pencil Code data.
#
# Chao-Chin Yang, 2013-10-22
#=======================================================================
def avg1d(name, datadir='./data', logscale=False, plane='xy', tsize=None,
          **kwargs):
    """Plots the space-time diagram of a 1D average.

    Positional Arguments:
        name
            Name of the average to be plotted.

    Keyword Arguments:
        datadir
            Name of the data directory.
        logscale
            Take logarithm or not.
        plane
            Plane of the average.
        tsize
            Number of regular time intervals; if None, it is
            automatically determined.
        **kwargs
            Sent to matplotlib.pyplot.imshow.
    """
    # Author: Chao-Chin Yang
    # Created: 2013-10-28
    # Last Modified: 2022-11-07
    from . import read
    import math
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import interp1d

    # Check the plane of the average.
    par = read.parameters(datadir=datadir)
    g = read.grid(datadir=datadir, par=par)
    if plane == 'xy':
        xlabel = '$z$'
        xdir = 2
        x = g.z
    elif plane == 'xz':
        xlabel = '$y$'
        xdir = 1
        x = g.y
    elif plane == 'yz':
        xlabel = '$x$'
        xdir = 0
        x = g.x
    else:
        raise ValueError("Keyword plane only accepts 'xy', 'xz', or 'yz'. ")

    # Read the data.
    print("Reading 1D averages...")
    time, avg = read.avg1d(datadir=datadir, plane=plane, verbose=False)
    xmin, xmax = par.xyz0[xdir], par.xyz1[xdir]

    # Set colorbar label.
    if logscale:
        cblabel = r'$\log(\tt{' + name + '})$'
    else:
        cblabel = name

    # Interpolate the time series.
    print("Interpolating", name, "...")
    tmin, tmax = np.min(time), np.max(time)
    if tsize is None:
        dt = (time[1:] - time[:-1]).min()
        tsize = int(math.ceil((tmax - tmin) / dt)) + 1
    ns = avg.shape[1]
    t = np.linspace(tmin, tmax, tsize)
    a = np.empty((tsize, ns))
    for j in range(ns):
        if logscale:
            a[:,j] = interp1d(time, np.log10(avg[name][:,j]))(t)
        else:
            a[:,j] = interp1d(time, avg[name][:,j])(t)

    # Plot the space-time diagram.
    print("Plotting...")
    img = plt.pcolormesh(x, t, a, **kwargs)
    ax = plt.gca()
    ax.set_ylabel('$t$')
    ax.set_xlabel(xlabel)
    cb = plt.colorbar(img)
    cb.set_label(cblabel)
    plt.show()
#=======================================================================
def time_series(diagnostics, datadir='./data', trange=None, xlog=False, ylog=False):
    """Plots diagnostic variable(s) as a function of time.

    Positional Arguments:
        diagnostics
            (A list of) diagnostic variable(s).

    Keyword Arguments:
        datadir
            Name of the data directory.
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
    # Chao-Chin Yang, 2014-10-11
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
