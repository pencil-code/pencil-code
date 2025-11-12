# timeseries.py
#
# Read the time_series.dat and return a TimeSeries object of 1D numpy
# arrrays
"""
Contains the classes and methods to read the time series file.
"""

from pencil.util import copy_docstring


class TimeSeries(object):
    """
    TimeSeries -- holds Pencil Code time series data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = []
        self.keys = []

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(
        self,
        file_name="time_series.dat",
        datadir="data",
        quiet=False,
        comment_char="#",
        sim=None,
        unique_clean=False,
        time_range=None,
        precision="f",
    ):
        """
        read(file_name='time_series.dat', datadir='data',
             quiet=False, comment_char='#', sim=None, unique_clean=False)

        Read Pencil Code time series data.

        Parameters
        ----------
        file_name : string
            Name of the time series file.
            For supernova data change file_name to 'sn_series.dat'.

        datadir : string
            Directory where the data is stored.

        quiet : bool
            Flag for switching off output.

        comment_char : string
            Comment character in the time series file.

        sim : obj
          Simulation object from which to take the datadir.

        unique_clean : bool
          Set True, np.unique is used to clean up the ts, e.g. remove errors
          at the end of crashed runs.

        time_range : bool
          List of length 2, start and end time, of float with end time.

        precision : float
          "f" (single,default) or "d" (double) or "h" (half).
        """

        import numpy as np
        import os.path
        import re

        if precision == "h":
            precision = "half"
        if sim:
            from pencil.sim import __Simulation__

            if isinstance(sim, __Simulation__):
                datadir = sim.datadir

        datadir = os.path.expanduser(datadir)
        with open(os.path.join(datadir, file_name), "r") as infile:
            lines = infile.readlines()

        nlines_init = len(lines)
        data = np.zeros((nlines_init, len(self.keys)),dtype=precision)
        nlines = 0
        for i, line in enumerate(lines):
            if re.search("^%s--" % comment_char, line):
                # Read header and create keys for dictionary.
                line = line.strip("{0}-\n".format(comment_char))
                keys_new = re.split("-+", line)
                if keys_new != self.keys:
                    n_newrows = abs(len(keys_new) - len(self.keys))
                    data = np.append(data, np.zeros((nlines_init, n_newrows), dtype=precision), axis=1)
                    self.keys = keys_new
            else:
                try:
                    row = np.array(list(map(float, re.split(" +", line.strip(" \n")))))
                    data[nlines, :] = row
                    nlines += 1
                except ValueError:
                    print(f"Invalid data on line {i}. Skipping.")

        # Clean up data.
        data = np.resize(data, (nlines, len(self.keys)))

        if not quiet:
            print("Read {0} lines.".format(nlines))

        # Assemble into a TimeSeries class.
        for i in range(0, len(self.keys)):
            setattr(self, self.keys[i], data[:, i])

        # Do unique clean up.
        if unique_clean:
            clean_t, unique_indices = np.unique(self.t, return_index=True)

            if np.size(clean_t) != np.size(self.t):
                for key in self.keys:
                    setattr(self, key, getattr(self, key)[unique_indices])
        if time_range:
            if isinstance(time_range, list):
                time_range = time_range
            else:
                time_range = [time_range]
            if len(time_range) == 1:
                start_time = 0.
                end_time = time_range[0]
            elif len(time_range) == 2:
                start_time = time_range[0]
                end_time = time_range[1]
            ilist = list()
            for i, time in zip(range(self.t.size),self.t):
                if time >= start_time:
                    if time <= end_time:
                        ilist.append(i)
            for key in self.keys:
                tmp = self.__getattribute__(key)[ilist]
                self.__delattr__(key)
                setattr(self, key, tmp)

@copy_docstring(TimeSeries.read)
def ts(*args, **kwargs):
    """
    Wrapper for :py:meth:`TimeSeries.read`
    """
    ts_tmp = TimeSeries()
    ts_tmp.read(*args, **kwargs)
    return ts_tmp
