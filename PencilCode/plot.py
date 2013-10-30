#=======================================================================
# plot.py
#
# Facilities for plotting the Pencil Code data.
#
# Chao-Chin Yang, 2013-10-22
# Last Modification: $Id$
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
    elif plane == 'xz':
        xlabel = '$y$'
    elif plane == 'yz':
        xlabel = '$x$'
    else:
        raise ValueError("Keyword plane only accepts 'xy', 'xz', or 'yz'. ")

    # Read the data.
    from . import read
    time, avg = read.avg1d(datadir=datadir, plane=plane)

    # Default variable name.
    if var is None:
        var = avg.dtype.names[0]

    # Interpolate the time series.
    import numpy as np
    from scipy.interpolate import interp1d
    tmin, tmax = np.min(time), np.max(time)
    ns = avg.shape[1]
    t = tmin + (np.arange(tsize) + 0.5) * (tmax - tmin) / tsize
    a = np.empty((tsize, ns))
    for j in range(ns):
        a[:,j] = interp1d(time, avg[var][:,j])(t)

    # Plot the space-time diagram.
    import matplotlib.pyplot as plt
    img = plt.imshow(a, origin='bottom', **kwargs)
    ax = plt.gca()
    ax.set_ylabel('$t$')
    ax.set_xlabel(xlabel)
    cb = plt.colorbar(img)
    cb.set_label(var)
    plt.show()

#=======================================================================
def time_series(datadir='./data', diagnostics='dt'):
    """Plots diagnostic variable(s) as a function of time.

    Keyword Arguments:
        datadir
            Name of the data directory.
        diagnostics
            (A list of) diagnostic variable(s).
    """
    # Chao-Chin Yang, 2013-10-22

    from . import read
    import matplotlib.pyplot as plt

    # Read the time series.
    ts = read.time_series(datadir=datadir)

    # Plot the diagnostics.
    if type(diagnostics) is list:
        for diag in diagnostics:
            plt.plot(ts.t, ts[diag])
    else:
        plt.plot(ts.t, ts[diagnostics])
    plt.xlabel('t')
    plt.ylabel(diagnostics)
    plt.show()

