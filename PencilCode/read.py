#! /usr/bin/env python3
# Last Modification: $Id$
#=======================================================================
# read.py
#
# Facilities for reading the Pencil Code data.
#
# Chao-Chin Yang, 2013-05-06
#=======================================================================
def avg1d(datadir='./data', plane='xy', verbose=True):
    """Returns the time series of 1D averages.

    Keyword Arguments:
        datadir
            Name of the data directory.
        plane
            Plane of average: 'xy', 'xz', or 'yz'.
        verbose
            Whether or not to print information.
    """
    # Chao-Chin Yang, 2013-10-31

    # Read the dimensions and check the plane of average.
    dim = dimensions(datadir=datadir)
    if plane == 'xy':
        nc = dim.nzgrid
    elif plane == 'xz':
        nc = dim.nygrid
    elif plane == 'yz':
        nc = dim.nxgrid
    else:
        raise ValueError("Keyword plane only accepts 'xy', 'xz', or 'yz'. ")

    # Read the names of the variables.
    var = varname(datadir=datadir+'/..', filename=plane.strip()+'aver.in')
    nvar = len(var)

    # Open file and define data stream.
    f = open(datadir.strip() + '/' + plane.strip() + 'averages.dat')
    import numpy as np
    def fetch(nval):
        return np.fromfile(f, count=nval, sep=' ')

    # Check the data size.
    if verbose:
        print("Checking the data size...")
    nt = 0
    nblock = nvar * nc
    while fetch(1).size:
        if fetch(nblock).size != nblock:
            raise EOFError("incompatible data file")
        nt += 1
    f.seek(0)

    # Read the data.
    if verbose:
        print("Reading 1D averages", var, "...")
    t = np.zeros(nt)
    avg = np.core.records.array(len(var) * [np.zeros((nt,nc))], names=var)
    for i in range(nt):
        t[i] = fetch(1)
        for v in var:
            avg[v][i,:] = fetch(nc)

    # Close file.
    f.close()

    return t, avg.view(np.recarray)

#=======================================================================
def dimensions(datadir='./data'):
    """Returns the dimensions of the Pencil Code data from datadir.

    Keyword Arguments:
        datadir
            Name of the data directory
    """
    # Chao-Chin Yang, 2013-10-23

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

#=======================================================================
def parameters(datadir='./data'):
    """Returns runtime parameters.

    Keyword Arguments:
        datadir
            Name of the data directory
    """
    # Chao-Chin Yang, 2013-10-31

    # Function to convert a string to the correct type.
    def convert(v):
        if v == 'T' or v == ".TRUE.":
            return True
        elif v == 'F' or v == ".FALSE.":
            return False
        else:
            try:
                return int(v)
            except ValueError:
                try:
                    return float(v)
                except ValueError:
                    return v.strip("' ")

    # Function to parse the values.
    def parse(v):
        v = v.split(',')
        u = []
        for w in v:
            w = w.strip()
            if '*' in w:
                nrep, sep, x = w.partition('*')
                u += int(nrep) * [convert(x)]
            else:
                u += [convert(w)]
        if len(u) == 1:
            return u[0]
        else:
            return u

    # Read the parameter file.
    f = open(datadir.strip() + "/param2.nml")
    keys, values = [], []
    for line in f:
        if '=' in line:
            k, s, v = line.partition('=')
            k = k.strip().lower()
            if k not in keys:
                keys.append(k)
                values.append(parse(v.strip(" ,\n")))
            else:
                print("Duplicate parameter:", k, '=', v.strip(" ,\n"))
    f.close()

    # Define a container class and return the parameters in it.
    class Parameter:
        """A container class to hold the parameters of the Pencil Code data. """
        def __init__(self, d):
            for key, value in d.items():
                setattr(self, key, value)
    return Parameter(dict(zip(keys, values)))

#=======================================================================
def proc_dim(datadir='./data', proc=0):
    """Returns the dimensions of the data from one process.

    Keyword Arguments:
        datadir
            Name of the data directory
        proc
            Process ID
    """
    # Chao-Chin Yang, 2013-10-23

    # Read dim.dat.
    f = open(datadir.strip() + '/proc' + str(proc) + '/dim.dat')
    a = f.read().rsplit()
    f.close()

    # Sanity Check
    if a[6] == '?':
        print('Warning: unknown data precision. ')
    if not a[7] == a[8] == a[9]:
        raise Exception('unequal number of ghost zones in every direction. ')

    # Extract the dimensions.
    mx, my, mz, mvar, maux, mglobal = (int(b) for b in a[0:6])
    double_precision = a[6] == 'D'
    nghost = int(a[7])
    nx, ny, nz = (int(b) - 2 * nghost for b in a[0:3])
    iprocx, iprocy, iprocz = (int(b) for b in a[10:13])

    # Define and return a named tuple.
    from collections import namedtuple
    Dimensions = namedtuple('Dimensions', ['nx', 'ny', 'nz', 'nghost', 'mx', 'my', 'mz', 'mvar', 'maux', 'mglobal',
                                           'double_precision', 'iprocx', 'iprocy', 'iprocz'])
    return Dimensions(nx=nx, ny=ny, nz=nz, nghost=nghost, mx=mx, my=my, mz=mz, mvar=mvar, maux=maux, mglobal=mglobal,
                      double_precision=double_precision, iprocx=iprocx, iprocy=iprocy, iprocz=iprocz)

#=======================================================================
def time_series(datadir='./data'):
    """Returns a NumPy recarray from the time series.

    Keyword Arguments:
        datadir
            Name of the data directory
    """
    # Chao-Chin Yang, 2013-05-13

    # Import necessary modules.
    from numpy import recfromtxt
    from io import BytesIO

    # Save the data path.
    path = datadir.strip()

    # Read legend.dat.
    f = open(path + '/legend.dat')
    names = f.read().replace('-', ' ').split()
    f.close()

    # Read time_series.dat.
    f = open(path + '/time_series.dat')
    ts = recfromtxt(BytesIO(f.read().encode()), names=names)
    f.close()

    return ts

#=======================================================================
def varname(datadir='./data', filename='varname.dat'):
    """Returns the names of variables.

    Keyword Arguments:
        datadir
            Name of the data directory
        filename
            Name of the file containing variable names
    """
    # Chao-Chin Yang, 2013-05-13

    # Read varname.dat.
    f = open(datadir.strip() + '/' + filename.strip())
    var = []
    for line in f:
        try:
            var.append(line.split()[1])
        except:
            var.append(line.rstrip('\n'))
    f.close()

    return var

