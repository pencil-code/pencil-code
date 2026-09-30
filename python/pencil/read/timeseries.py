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
        return list(self.__dict__.keys())

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

        precision : str
          "f" (single,default) or "d" (double) or "h" (half).
        """

        import numpy as np
        import os.path
        import re

        if precision == "h":
            precision = "half"
        if sim:
            from pencil.sim import Simulation

            if isinstance(sim, Simulation):
                datadir = sim.datadir

        datadir = os.path.expanduser(datadir)
        with open(os.path.join(datadir, file_name), "r") as infile:
            lines = infile.readlines()

        # The pencil code writes a new "#--key1--key2--..." header line
        # whenever the set of diagnostic variables changes (variables
        # added, removed or reordered). Split the file into segments
        # between headers, each with its own fixed column layout, so
        # that a later merge can match columns by name rather than by
        # position.
        segments = []
        keys = None
        rows = None
        for i, line in enumerate(lines):
            if re.search("^%s--" % comment_char, line):
                keys_new = re.split("-+", line.strip("{0}-\n".format(comment_char)))
                if keys_new != keys:
                    if keys is not None and rows:
                        segments.append((keys, rows))
                    keys = keys_new
                    rows = []
                continue
            if keys is None:
                # Data lines before the first header cannot be interpreted.
                continue
            fields = re.split(" +", line.strip(" \n"))
            if len(fields) != len(keys):
                print(f"Invalid data on line {i}. Skipping.")
                continue
            try:
                rows.append(np.array(fields, dtype=precision))
            except ValueError:
                print(f"Invalid data on line {i}. Skipping.")
        if keys is not None and rows:
            segments.append((keys, rows))

        # Union of all keys across segments, in order of first appearance.
        self.keys = []
        for seg_keys, _ in segments:
            for key in seg_keys:
                if key not in self.keys:
                    self.keys.append(key)

        nlines = sum(len(seg_rows) for _, seg_rows in segments)
        data = np.full((nlines, len(self.keys)), np.nan, dtype=precision)

        irow = 0
        for seg_keys, seg_rows in segments:
            icols = [self.keys.index(key) for key in seg_keys]
            for row in seg_rows:
                data[irow, icols] = row
                irow += 1

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
