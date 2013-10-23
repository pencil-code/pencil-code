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

    # Extract the dimensions.
    mxgrid, mygrid, mzgrid, mvar, maux, mglobal = (int(b) for b in a[0:6])
    double_precision = a[6] == 'D'
    nghost = int(a[7])
    nxgrid, nygrid, nzgrid = (int(b) - 2 * nghost for b in a[0:3])
    nprocx, nprocy, nprocz = (int(b) for b in a[10:13])
    procz_last = int(a[13]) == 1

    # Define and return a named tuple.
    from collections import namedtuple
    Dimensions = namedtuple('Dimensions', ['nxgrid', 'nygrid', 'nzgrid', 'nghost', 'mxgrid', 'mygrid', 'mzgrid',
                                           'mvar', 'maux', 'mglobal', 'double_precision',
                                           'nprocx', 'nprocy', 'nprocz', 'procz_last'])
    return Dimensions(nxgrid=nxgrid, nygrid=nygrid, nzgrid=nzgrid, nghost=nghost, mxgrid=mxgrid, mygrid=mygrid, mzgrid=mzgrid,
                      mvar=mvar, maux=maux, mglobal=mglobal, double_precision=double_precision,
                      nprocx=nprocx, nprocy=nprocy, nprocz=nprocz, procz_last=procz_last)


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

