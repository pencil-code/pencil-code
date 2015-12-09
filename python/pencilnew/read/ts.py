# $Id$
#
# read time_series.dat and return a TimeSeries class of 1D numpy
# arrrays
#
#
import os.path
import re
import numpy as N


def ts(*args, **kwargs):
    """Read Pencil Code time series data.
    params:
     string: filename  ='time_series.dat'
     string: datadir   = 'data'
     logical: double    = 0
     logical: print_std = 0
     logical: quiet     = 0
    """
    return TimeSeries(*args, **kwargs)


class TimeSeries(object):
    """
    TimeSeries -- holds pencil code time series data. each variable is
    represented by a data member of the class.
    """

    def __init__(self, filename='time_series.dat', datadir='data',
                 double=0, print_std=0, quiet=0, comment_char='#'):
        """
        constructor:
        -----------
        __init__ -- reads Pencil Code time series data.
        Modeled after idl function of same name.

        params:
        ______
         string: filename  ='time_series.dat'
         string: datadir   = 'data'
         logical: double    = 0
         logical: print_std = 0
         logical: quiet     = 0
        """

        datadir = os.path.expanduser(datadir)
        infile = open(os.path.join(datadir, filename), "r")
        lines = infile.readlines()
        infile.close()

        # need to handle cases where restart AND print.in changes,
        # but not right away
        # idl version uses input_table function with a STOP_AT
        # and FILEPOSITION keywords
        nlines_init = len(lines)
        self.keys = []
        data = N.zeros((nlines_init, len(self.keys)))
        nlines = 0
        for line in lines:
            if re.search("^%s--" % comment_char, line):
                # read header and create keys for dictionary
                line = line.strip("%s-\n" % comment_char)
                keys_new = re.split("-+", line)
                if keys_new != self.keys:
                    n_newrows = abs(len(keys_new) - len(self.keys))
                    data = N.append(data, N.zeros((nlines_init, n_newrows)),
                                 axis=1)
                    self.keys = keys_new
            else:
                try:
                    row = N.array(list(map(float, re.split(" +", line.strip(" \n")))))
                    data[nlines, :] = row
                    nlines += 1
                except ValueError:
                    print(("Invalid data on line %i. Skipping." % nlines))
        #clean up data
        data = N.resize(data, (nlines, len(self.keys)))

        if (not quiet):
            print(("Read",nlines,"lines."))

        #assemble into a TimeSeries class
        for i in range(0, len(self.keys)):
            setattr(self, self.keys[i], data[:,i])

    def __repr__(self):
        count = len(self.t)
        t_min = self.t.min()
        t_max = self.t.max()
        return "TimeSeries(t=%g..%g, rows=%d): %s" \
            % (t_min, t_max, count, str(self.keys), )


if __name__=='__main__':
    print((TimeSeries.__doc__))
