#
# read.py
#
# Facilities for reading the Pencil Code data.
#
# Chao-Chin Yang, 2013-05-06
# Last Modification: $Id$
#
class Dim:
    """object holding the dimensions of the data.

       Data attributes:

         double:     True if the data is in double precision.
         procz_last: True if the z direction runs last.
         maux:       number of auxiliary variables.
         mglobal:    number of global variables.
         mvar:       number of state variables.
         mxgrid, mygrid, mzgrid:
                     number of cells including ghost cells in each direction.
         nghost:     number of ghost cells from boundary.
         nprocx, nprocy, nprocz:
                     number of processors allocated for each dimension.
    """

    def __init__(self, datadir='./data'):
        """reads the dimensions of the data from datadir. """
        #
        # Chao-Chin Yang, 2013-05-13
        #
        f = open(datadir.strip() + '/dim.dat')
        a = f.read().rsplit()
        self.mxgrid, self.mygrid, self.mzgrid, self.mvar, self.maux, self.mglobal = (int(b) for b in a[0:6])
        self.double = a[6] == 'D'
        if a[6] == '?':
            print('Warning: unknown data precision. ')
        if not a[7] == a[8] == a[9]:
            raise Exception('unequal number of ghost zones in every direction. ')
        self.nghost = int(a[7])
        self.nprocx, self.nprocy, self.nprocz = (int(b) for b in a[10:13])
        self.procz_last = a[13] == 1
        f.close()


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


