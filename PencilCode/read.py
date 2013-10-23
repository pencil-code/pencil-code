#
# read.py
#
# Facilities for reading the Pencil Code data.
#
# Chao-Chin Yang, 2013-05-06
# Last Modification: $Id$
#
def dimensions(datadir='./data'):
    """Returns the dimensions of the Pencil Code data from datadir. """
    #
    # Chao-Chin Yang, 2013-10-23
    #
    # Read dim.dat.
    f = open(datadir.strip() + '/dim.dat')
    a = f.read().rsplit()
    f.close()
    # Sanity check
    if a[6] == '?':
        print('Warning: unknown data precision. ')
    if not a[7] == a[8] == a[9]:
        raise Exception('unequal number of ghost zones in different directions. ')
    # Define and return a dimension object.
    class Dimension:
        """Contains the dimensions of the Pencil Code data."""
        @property
        def nghost(self):
            """number of ghost cells"""
            return int(a[7])
        @property
        def mxgrid(self):
            """number of cells in x direction, including ghost cells"""
            return int(a[0])
        @property
        def mygrid(self):
            """number of cells in y direction, including ghost cells"""
            return int(a[1])
        @property
        def mzgrid(self):
            """number of cells in z direction, including ghost cells"""
            return int(a[2])
        @property
        def nxgrid(self):
            """number of cells in x direction, excluding ghost cells"""
            return int(a[0]) - 2 * int(a[7])
        @property
        def nygrid(self):
            """number of cells in y direction, excluding ghost cells"""
            return int(a[1]) - 2 * int(a[7])
        @property
        def nzgrid(self):
            """number of cells in z direction, excluding ghost cells"""
            return int(a[2]) - 2 * int(a[7])
        @property
        def mvar(self):
            """number of state variables"""
            return int(a[3])
        @property
        def maux(self):
            """number of auxiliary variables"""
            return int(a[4])
        @property
        def mglobal(self):
            """number of global variables"""
            return int(a[5])
        @property
        def double_precision(self):
            """True if the data is in double precision; False otherwise"""
            return a[6] == 'D'
        @property
        def nprocx(self):
            """number of processors in x direction"""
            return int(a[10])
        @property
        def nprocy(self):
            """number of processors in x direction"""
            return int(a[11])
        @property
        def nprocz(self):
            """number of processors in x direction"""
            return int(a[12])
        @property
        def procz_last(self):
            """True if the z direction is the last; False otherwise"""
            return int(a[13]) == 1
    return Dimension()


class DimProc:
    """object holding the dimensions of the data held by one processor.

       Data attributes:

         double:     True if the data is in double precision.
         maux:       number of auxiliary variables.
         mglobal:    number of global variables.
         mvar:       number of state variables.
         mx, my, mz: number of cells including ghost cells in each direction.
         nghost:     number of ghost cells from boundary.
         iprocx, iprocy, iprocz:
                     processor ID in each direction.
    """

    def __init__(self, proc=0, datadir='./data'):
        """reads the dimensions of the data from processor of ID proc under datadir. """
        #
        # Chao-Chin Yang, 2013-05-06
        #
        f = open(path_proc(proc, datadir.strip()) + '/dim.dat')
        a = f.read().rsplit()
        self.mx, self.my, self.mz, self.mvar, self.maux, self.mglobal = (int(b) for b in a[0:6])
        self.double = a[6] == 'D'
        if a[6] == '?':
            print('Warning: unknown data precision. ')
        if not a[7] == a[8] == a[9]:
            raise Exception('unequal number of ghost zones in every direction. ')
        self.nghost = int(a[7])
        self.iprocx, self.iprocy, self.iprocz = (int(b) for b in a[10:13])
        f.close()


def time_series(datadir='./data'):
    """returns a NumPy recarray from the time series data under datadir. """
    #
    # Chao-Chin Yang, 2013-05-13
    #
    from numpy import recfromtxt
    from io import BytesIO

    path = datadir.strip()

    f = open(path + '/legend.dat')
    names = f.read().replace('-', ' ').split()
    f.close()

    f = open(path + '/time_series.dat')
    ts = recfromtxt(BytesIO(f.read().encode()), names=names)
    f.close()

    return ts


def path_proc(proc=0, datadir='./data'):
    """returns the path to the data held by the processor of ID proc. """
    #
    # Chao-Chin Yang, 2013-05-13
    #
    return datadir.strip() + '/proc' + str(proc)


def varname(datadir='./data'):
    """returns the names of the farray variables. """
    #
    # Chao-Chin Yang, 2013-05-13
    #
    f = open(datadir.strip() + '/varname.dat')
    var = []
    for line in f:
        var.append(line.split()[1])
    f.close()
    return var

