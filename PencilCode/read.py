#! /usr/bin/env python3
# Last Modification: $Id$
#=======================================================================
# read.py
#
# Facilities for reading the Pencil Code data.
#
# Chao-Chin Yang, 2013-05-06
#=======================================================================
hsize = 4    # Heading and trailing bytes in Fortran binary.
#=======================================================================
def avg1d(datadir='./data', plane='xy', tsize=None, verbose=True):
    """Returns the time series of 1D averages.

    Keyword Arguments:
        datadir
            Name of the data directory.
        plane
            Plane of average: 'xy', 'xz', or 'yz'.
        tsize
            If not None, the time series is interpolated onto regular
            time intervals with tsize elements and endpoints included.
        verbose
            Whether or not to print information.
    """
    # Chao-Chin Yang, 2014-08-29
    import numpy as np
    from scipy.interpolate import interp1d
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
    # Interpolate the time series if requested.
    if tsize is not None:
        if verbose:
            print("Interpolating...")
        tmin, tmax = t.min(), t.max()
        ti = tmin + (tmax - tmin) / (tsize - 1) * np.arange(tsize)
        ti[0], ti[-1] = tmin, tmax
        avgi = np.core.records.array(len(var) * [np.zeros((tsize,nc))], names=var)
        for v in var:
            for k in range(nc):
                avgi[v][:,k] = interp1d(t, avg[v][:,k])(ti)
        t, avg = ti, avgi
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
def parameters(datadir='./data', par2=False):
    """Returns runtime parameters.

    Keyword Arguments:
        datadir
            Name of the data directory.
        par2
            Read param2.nml if True; read param.nml otherwise.
    """
    # Chao-Chin Yang, 2014-06-18

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
    if par2:
        f = open(datadir.strip() + "/param2.nml")
    else:
        f = open(datadir.strip() + "/param.nml")
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
def pdim(datadir='./data'):
    """Returns the numbers associated with particles.

    Keyword Arguments:
        datadir
            Name of the data directory
    """
    # Chao-Chin Yang, 2014-02-11

    # Read pdim.dat.
    f = open(datadir.strip() + '/pdim.dat')
    a = f.read().rsplit()
    f.close()

    # Extract the numbers.
    npar, mpvar, npar_stalk = (int(b) for b in a)

    # Define and return a named tuple.
    from collections import namedtuple
    ParticleNumbers = namedtuple('ParticleNumbers', ['npar', 'mpvar', 'npar_stalk'])
    return ParticleNumbers(npar=npar, mpvar=mpvar, npar_stalk=npar_stalk)
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
    from collections import namedtuple
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
    Dimensions = namedtuple('Dimensions', ['nx', 'ny', 'nz', 'nghost', 'mx', 'my', 'mz', 'mvar', 'maux', 'mglobal',
                                           'double_precision', 'iprocx', 'iprocy', 'iprocz'])
    return Dimensions(nx=nx, ny=ny, nz=nz, nghost=nghost, mx=mx, my=my, mz=mz, mvar=mvar, maux=maux, mglobal=mglobal,
                      double_precision=double_precision, iprocx=iprocx, iprocy=iprocy, iprocz=iprocz)
#=======================================================================
def proc_grid(datadir='./data', dim=None, proc=0):
    """Returns the grid controlled by one process.

    Keyword Arguments:
        datadir
            Name of the data directory.
        dim
            Dimensions supplied by proc_dim().  If None, proc_dim() will
            be called.
        proc
            Process ID.
    """
    # Chao-Chin Yang, 2014-10-29
    from collections import namedtuple
    import numpy as np
    from struct import unpack, calcsize
    # Check the dimensions and precision.
    if dim is None:
        dim = proc_dim(datadir=datadir, proc=proc)
    if dim.double_precision:
        dtype = np.float64
        fmt = 'd'
    else:
        dtype = np.float32
        fmt = 'f'
    nb = calcsize(fmt)
    # Read grid.dat.
    f = open(datadir.strip() + '/proc' + str(proc) + '/grid.dat', 'rb')
    f.read(hsize)
    t = unpack(fmt, f.read(nb))[0]
    x = np.frombuffer(f.read(nb*dim.mx), dtype=dtype)
    y = np.frombuffer(f.read(nb*dim.my), dtype=dtype)
    z = np.frombuffer(f.read(nb*dim.mz), dtype=dtype)
    dx, dy, dz = unpack(3*fmt, f.read(3*nb))
    f.read(2*hsize)
    dx, dy, dz = unpack(3*fmt, f.read(3*nb))
    f.read(2*hsize)
    Lx, Ly, Lz = unpack(3*fmt, f.read(3*nb))
    f.read(2*hsize)
    dx_1 = np.frombuffer(f.read(nb*dim.mx), dtype=dtype)
    dy_1 = np.frombuffer(f.read(nb*dim.my), dtype=dtype)
    dz_1 = np.frombuffer(f.read(nb*dim.mz), dtype=dtype)
    f.read(2*hsize)
    dx_tilde = np.frombuffer(f.read(nb*dim.mx), dtype=dtype)
    dy_tilde = np.frombuffer(f.read(nb*dim.my), dtype=dtype)
    dz_tilde = np.frombuffer(f.read(nb*dim.mz), dtype=dtype)
    f.read(hsize)
    f.close()
    # Define and return a named tuple.
    Grid = namedtuple('Grid', ['x', 'y', 'z', 'dx', 'dy', 'dz', 'Lx', 'Ly', 'Lz', 'dx_1', 'dy_1', 'dz_1',
                               'dx_tilde', 'dy_tilde', 'dz_tilde'])
    return Grid(x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, Lx=Lx, Ly=Ly, Lz=Lz, dx_1=dx_1, dy_1=dy_1, dz_1=dz_1,
                dx_tilde=dx_tilde, dy_tilde=dy_tilde, dz_tilde=dz_tilde)
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

