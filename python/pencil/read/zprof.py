# ts.py
#
# Read vertical profiles written in data/proc*/zprof_varname.dat.
"""
Contains the classes and methods to read the time series file.
"""


def zprof(*args, **kwargs):
    """
    zprof(var_name, datadir='data', dim=None, nfield=1)

    Read vertical profiles written in data/proc*/zprof_varname.dat.
    Returns a ZProfile object with z and profiles(z).

    Parameters
    ----------
    var_name : string
      Name of the zprof var file.

    datadir : string
      Directory where the data is stored.

    dim : obj
      Dimension object.

    nfield : int
      Number of fields to be read.
    """

    zprof_tmp = ZProfile()
    zprof_tmp.read(*args, **kwargs)
    return zprof_tmp


class ZProfile(object):
    """
    ZProfile -- holds the zprofile data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.prof = 0
        self.z = 0

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(self, var_name, datadir="data", dim=None, nfield=1):
        """
        read(var_name, datadir='data', dim=None, nfield=1)

        Read vertical profiles written in data/proc*/zprof_varname.dat.
        Returns a ZProfile object with z and profiles(z).

        Parameters
        ----------
        var_name : string
          Name of the zprof var file.

        datadir : string
          Directory where the data is stored.

        dim : obj
          Dimension object.

        nfield : int
          Number of fields to be read.
        """

        import os as os
        import numpy as np
        from pencil import read

        if not dim:
            dim = read.dim()

        nz = int(dim.nzgrid / dim.nprocz)
        self.z = np.zeros(nz * dim.nprocz, dtype=np.float32)
        if nfield > 1:
            self.prof = np.zeros((nfield, dim.nzgrid), dtype=np.float32)
        else:
            self.prof = np.zeros(dim.nzgrid, dtype=np.float32)

        # Loop over all processors and records in file.
        izcount = 0
        for iprocz in range(0, dim.nprocz):
            proc_name = "proc{0}".format(iprocz)
            file_name = os.path.join(datadir, proc_name, "zprof_", var_name, ".dat")
            fd = open(file_name, "r")

            #  When reading a zprof_once_X file, the first dim.nghostz gridpoints are
            #  not saved.
            if var_name.find("once") != -1:
                for i in range(dim.nghostz):
                    line = fd.readline()
            for i in range(nz):
                line = fd.readline()
                data = np.asarray(line.split()).astype(np.float32)
                self.z[izcount] = data[0]
                if nfield > 1:
                    for j in range(nfield):
                        self.prof[j, izcount] = data[j + 1]
                else:
                    self.prof[izcount] = data[1]
                izcount = izcount + 1
        fd.close()
