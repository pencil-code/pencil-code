# pdims.py
"""
Contains the perticle dimension class and its reading routine.
"""

from pencil.util import copy_docstring

class PDim(object):
    """
    Class holding the data from pdim.dat and its methods.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.npar = 0
        self.mpvar = 0
        self.npar_stalk = 0
        self.mpaux = 0

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(self, datadir="data"):
        """
        read(sim=None, datadir='data')

        Read the pdim.dat file.

        Parameters
        ----------
        sim : obj
            Specifies the simulation object from which to obtain the datadir.
        
        datadir : string
            Directory where the data is stored.

        Returns
        -------
        Object containing the particle dimensions.
        """

        import os
        import numpy as np

        # Contains the global particle properties.
        filename = os.path.join(datadir, "pdim.dat")

        try:
            filename = os.path.expanduser(filename)
            file = open(filename, "r")
        except IOError:
            print("File " + filename + " could not be opened.")
        else:
            lines = file.readlines()[0].split()
            file.close()
            if np.size(lines) == 3:
                npar, mpvar, mpaux = tuple(map(int, lines))
                npar_stalk = 0
            if np.size(lines) == 4:
                npar, mpvar, npar_stalk, mpaux = tuple(map(int, lines))

        self.npar = npar
        self.mpvar = mpvar
        self.npar_stalk = npar_stalk
        self.mpaux = mpaux

@copy_docstring(PDim.read)
def pdim(*args, **kwargs):
    """
    Wrapper for :py:meth:`PDim.read`
    """
    pdim_tmp = PDim()
    pdim_tmp.read(*args, **kwargs)
    return pdim_tmp
