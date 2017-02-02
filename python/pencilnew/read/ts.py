# ts.py
#
# Read the time_series.dat and return a TimeSeries class of 1D numpy
# arrrays
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read the time series file.
"""


def ts(*args, **kwargs):
    """
    Read Pencil Code time series data.

    call signature:

    ts(file_name='time_series.dat', datadir='data',
       quiet=False, comment_char='#')

    Keyword arguments:

    *file_name*:
      Name of the time series file.

    *datadir*:
      Directory where the data is stored.

    *quiet*
      Flag for switching off output.

    *comment_char*
      Comment character in the time series file.
    """

    ts_tmp = TimeSeries()
    ts_tmp.read(*args, **kwargs)
    return ts_tmp


class TimeSeries(object):
    """
    TimeSeries -- holds Pencil Code time series data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.keys = []


    def read(self, file_name='time_series.dat', datadir='data',
             quiet=False, comment_char='#', sim=False):
        """
        Read Pencil Code time series data.

        call signature:

        read(self, file_name='time_series.dat', datadir='data',
             double=0, quiet=0, comment_char='#')

        Keyword arguments:

        *SIM*:
          return time series for this simulation object

        *file_name*:
          Name of the time series file.

        *datadir*:
          Directory where the data is stored.

        *quiet*
          Flag for switching of output.

        *comment_char*
          Comment character in the time series file.
        """

        import numpy as np
        import os.path
        import re

        if str(sim.__class__) == "<class 'pencilnew.sim.simulation.__Simulation__'>":
            datadir = sim.datadir

        datadir = os.path.expanduser(datadir)
        infile = open(os.path.join(datadir, file_name), "r")
        lines = infile.readlines()
        infile.close()

        # Need to handle cases where restart AND print.in changes,
        # but not right away.
        nlines_init = len(lines)
        data = np.zeros((nlines_init, len(self.keys)))
        nlines = 0
        for line in lines:
            if re.search("^%s--" % comment_char, line):
                # Read header and create keys for dictionary.
                line = line.strip("%s-\n" % comment_char)
                keys_new = re.split("-+", line)
                if keys_new != self.keys:
                    n_newrows = abs(len(keys_new) - len(self.keys))
                    data = np.append(data, np.zeros((nlines_init, n_newrows)),
                                     axis=1)
                    self.keys = keys_new
            else:
                try:
                    row = np.array(map(float, re.split(" +", line.strip(" \n"))))
                    data[nlines, :] = row
                    nlines += 1
                except ValueError:
                    print("Invalid data on line {0}. Skipping.".format(nlines))

        # Clean up data.
        data = np.resize(data, (nlines, len(self.keys)))

        if not quiet:
            print("Read {0} lines".format(nlines))

        # Assemble into a TimeSeries class.
        for i in range(0, len(self.keys)):
            setattr(self, self.keys[i], data[:, i])
