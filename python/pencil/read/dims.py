# dims.py
#
# Read the dimensions of the simulation.
"""
Contains the classes and methods to read the simulation dimensions.
"""
import numpy as np


def dim(*args, **kwargs):
    """
    dim(datadir='data', proc=-1)

    Read the dim.dat file.

    Parameters
    ----------
    datadir : string
      Directory where the data is stored.

    proc : int
      Processor to be read. If proc is -1, then read the 'global'
      dimensions. If proc is >=0, then read the dim.dat in the
      corresponding processor directory.

    Returns
    -------
    Class containing the domain dimension information.
    """

    dim_tmp = Dim()
    dim_tmp.read(*args, **kwargs)
    return dim_tmp


class Dim(object):
    """
    Dim -- holds pencil code dimension data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.mx = self.my = self.mz = 0
        self.mvar = 0
        self.maux = 0
        self.mglobal = 0

        self.precision = "S"
        self.nghostx = self.nghosty = self.nghostz = 0

        self.nprocx = self.nprocy = self.nprocz = 0
        self.iprocz_slowest = 0
        self.ipx = self.ipy = self.ipz = 0

        # Add derived quantities to the dim object.
        self.nx = self.ny = self.nz = 0
        self.mw = 0
        self.l1 = self.l2 = 0
        self.m1 = self.m2 = 0
        self.n1 = self.n2 = 0

        self.nxgrid = self.nygrid = self.nzgrid = 0
        self.mxgrid = self.mygrid = self.mzgrid = 0

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(self, datadir="data", proc=-1, ogrid=False, down=False):
        """
        dim(datadir='data', proc=-1)

        Read the dim.dat file.

        Parameters
        ----------
        datadir : string
          Directory where the data is stored.

        proc : int
          Processor to be read. If proc is -1, then read the 'global'
          dimensions. If proc is >=0, then read the dim.dat in the
          corresponding processor directory.

        Returns
        -------
        Class containing the domain dimension information.
        """

        import os

        if os.path.exists(os.path.join(datadir, "grid.h5")):
            import h5py

            with h5py.File(os.path.join(datadir, "grid.h5"), "r") as tmp:
                self.mx = np.array(tmp["settings"]["mx"][()]).item()
                self.my = np.array(tmp["settings"]["my"][()]).item()
                self.mz = np.array(tmp["settings"]["mz"][()]).item()
                self.mvar = np.array(tmp["settings"]["mvar"][()]).item()
                self.maux = np.array(tmp["settings"]["maux"][()]).item()
                self.mglobal = np.array(tmp["settings"]["mglobal"][()])
                self.precision = np.array(tmp["settings"]["precision"][()]).item()
                self.nghostx = np.array(tmp["settings"]["nghost"][()]).item()
                self.nghosty = np.array(tmp["settings"]["nghost"][()]).item()
                self.nghostz = np.array(tmp["settings"]["nghost"][()]).item()
                self.nprocx = np.array(tmp["settings"]["nprocx"][()]).item()
                self.nprocy = np.array(tmp["settings"]["nprocy"][()]).item()
                self.nprocz = np.array(tmp["settings"]["nprocz"][()]).item()
                self.nx = np.array(tmp["settings"]["nx"][()]).item()
                self.ny = np.array(tmp["settings"]["ny"][()]).item()
                self.nz = np.array(tmp["settings"]["nz"][()]).item()
                self.l1 = np.array(tmp["settings"]["l1"][()]).item()
                self.l2 = np.array(tmp["settings"]["l2"][()]).item()
                self.m1 = np.array(tmp["settings"]["m1"][()]).item()
                self.m2 = np.array(tmp["settings"]["m2"][()]).item()
                self.n1 = np.array(tmp["settings"]["n1"][()]).item()
                self.n2 = np.array(tmp["settings"]["n2"][()]).item()
                self.iprocz_slowest = 0
                self.ipx = self.ipy = self.ipz = 0
                self.nxgrid = np.array(tmp["settings"]["nx"][()]).item()
                self.nygrid = np.array(tmp["settings"]["ny"][()]).item()
                self.nzgrid = np.array(tmp["settings"]["nz"][()]).item()
                self.mxgrid = np.array(tmp["settings"]["mx"][()]).item()
                self.mygrid = np.array(tmp["settings"]["my"][()]).item()
                self.mzgrid = np.array(tmp["settings"]["mz"][()]).item()
                self.mw = self.mx * self.my * self.mz
        else:
            if not ogrid:
                if down:
                    file_name = "dim_down.dat"
                else:
                    file_name = "dim.dat"
            else:
                file_name = "ogdim.dat"

            if proc < 0:
                file_name = os.path.join(datadir, file_name)
            else:
                file_name = os.path.join(datadir, "proc{0}".format(proc), file_name)

            try:
                file_name = os.path.expanduser(file_name)
                dim_file = open(file_name, "r")
            except IOError:
                print("File {0} could not be opened.".format(file_name))
                return -1
            else:
                lines = dim_file.readlines()
                dim_file.close()

            if len(lines[0].split()) == 6:
                self.mx, self.my, self.mz, self.mvar, self.maux, self.mglobal = tuple(
                    map(int, lines[0].split())
                )
            else:
                self.mx, self.my, self.mz, self.mvar, self.maux = tuple(
                    map(int, lines[0].split())
                )
                self.mglobal = 0

            self.precision = lines[1].strip("\n")
            self.nghostx, self.nghosty, self.nghostz = tuple(map(int, lines[2].split()))
            if proc < 0:
                # Set global parameters.
                self.nprocx, self.nprocy, self.nprocz, self.iprocz_slowest = tuple(
                    map(int, lines[3].split())
                )
                self.ipx = self.ipy = self.ipz = -1
            else:
                # Set local parameters to this proc.
                self.ipx, self.ipy, self.ipz = tuple(map(int, lines[3].split()))
                self.nprocx = self.nprocy = self.nprocz = self.iprocz_slowest = -1

            # Add derived quantities to the dim object.
            self.nx = self.mx - (2 * self.nghostx)
            self.ny = self.my - (2 * self.nghosty)
            self.nz = self.mz - (2 * self.nghostz)
            self.mw = self.mx * self.my * self.mz
            self.l1 = self.nghostx
            self.l2 = self.mx - self.nghostx - 1
            self.m1 = self.nghosty
            self.m2 = self.my - self.nghosty - 1
            self.n1 = self.nghostz
            self.n2 = self.mz - self.nghostz - 1
            if self.ipx == self.ipy == self.ipz == -1:
                # Set global parameters.
                self.nxgrid = self.nx
                self.nygrid = self.ny
                self.nzgrid = self.nz
                self.mxgrid = self.nxgrid + (2 * self.nghostx)
                self.mygrid = self.nygrid + (2 * self.nghosty)
                self.mzgrid = self.nzgrid + (2 * self.nghostz)
            else:
                # Set local parameters to this proc.
                self.nxgrid = self.nygrid = self.nzgrid = 0
                self.mxgrid = self.mygrid = self.mzgrid = 0

        return 0
