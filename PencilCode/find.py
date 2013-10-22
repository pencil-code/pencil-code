#
# find.py
#
# Facilities for analyzing the Pencil Code data.
#
# Chao-Chin Yang, 2013-10-21
# Last Modification: $Id$
#
def time_average(datadir='./data', diagnostics=['urms'], tmin=0):
    """Finds the time average of each given diagnostic variable.

    Keyword Arguments:
        datadir:  Name of the data directory.
        diagnostics:  A list of the diagnostic variables.
        tmin:  Starting time of the time average.

    A dictionary is returned with the mean and standard deviation of
    each diagnostics.
    """
    #
    # Chao-Chin Yang, 2013-10-21
    #
    from . import read

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
        from scipy import integrate
        from numpy import sqrt
        mean = dtinv * integrate.simps(ts[diag][indices], ts.t[indices])
        stddev = sqrt(dtinv * integrate.simps((ts[diag][indices] - mean)**2, ts.t[indices]))
        print("<", diag, "> = ", mean, "+/-", stddev)
        return mean, stddev

    # Find and return the time average and its standard deviation of each diagnostics.
    return dict((diag, stats(diag)) for diag in diagnostics)

