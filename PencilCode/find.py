#
# find.py
#
# Facilities for analyzing the Pencil Code data.
#
# Chao-Chin Yang, 2013-10-21
# Last Modification: $Id$
#
def time_average(datadir='./data', diagnostics=['urms'], tmin=0):
    """Find the time average of each given diagnostic variable.
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
        mean = dtinv * integrate.simps(ts[diag], ts.t)
        stddev = sqrt(dtinv * integrate.simps((ts[diag] - mean)**2, ts.t))
        print("<", diag, "> = ", mean, " +/- ", stddev)
        return (mean, stddev)

    return dict((diag, stats(diag)) for diag in diagnostics)
