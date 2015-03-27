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
    # Chao-Chin Yang, 2014-11-27
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
    # Discard duplicate entries.
    t, indices = np.unique(t, return_index=True)
    avg = avg[indices]
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
    from collections import namedtuple
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
    Dimensions = namedtuple('Dimensions', ['nxgrid', 'nygrid', 'nzgrid', 'nghost', 'mxgrid', 'mygrid', 'mzgrid',
                                           'mvar', 'maux', 'mglobal', 'double_precision',
                                           'nprocx', 'nprocy', 'nprocz', 'procz_last'])
    return Dimensions(nxgrid=nxgrid, nygrid=nygrid, nzgrid=nzgrid, nghost=nghost, mxgrid=mxgrid, mygrid=mygrid, mzgrid=mzgrid,
                      mvar=mvar, maux=maux, mglobal=mglobal, double_precision=double_precision,
                      nprocx=nprocx, nprocy=nprocy, nprocz=nprocz, procz_last=procz_last)
#=======================================================================
def grid(datadir='./data', trim=False):
    """Returns the coordinates and their derivatives of the grid.

    Keyword Arguments:
        datadir
            Name of the data directory.
        trim
            Whether or not to trim the ghost cells.
    """
    # Chao-Chin Yang, 2014-11-02
    from collections import namedtuple
    import numpy as np
    # Get the dimensions.
    dim = dimensions(datadir=datadir)
    # Allocate arrays.
    x = np.zeros(dim.mxgrid,)
    y = np.zeros(dim.mygrid,)
    z = np.zeros(dim.mzgrid,)
    dx, dy, dz, Lx, Ly, Lz = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    dx_1 = np.zeros(dim.mxgrid,)
    dy_1 = np.zeros(dim.mygrid,)
    dz_1 = np.zeros(dim.mzgrid,)
    dx_tilde = np.zeros(dim.mxgrid,)
    dy_tilde = np.zeros(dim.mygrid,)
    dz_tilde = np.zeros(dim.mzgrid,)
    # Define functions for assigning local coordinates to global.
    def assign(l1, l2, xg, xl, label, proc):
        indices = xg[l1:l2] != 0.0
        if any(xg[l1:l2][indices] != xl[indices]):
            print("Warning: inconsistent ", label, " for process", proc)
        xg[l1:l2][~indices] = xl[~indices]
        return xg
    def assign1(old, new, label, proc):
        if old != 0.0 and old != new:
            print("Warning: inconsistent ", label, " for process", proc)
        return new
    # Read the grid from each process and construct the whole grid.
    for proc in range(dim.nprocx * dim.nprocy * dim.nprocz):
        dim1 = proc_dim(datadir=datadir, proc=proc)
        grid = proc_grid(datadir=datadir, dim=dim1, proc=proc)
        # x direction
        l1 = dim1.iprocx * dim1.nx
        l2 = l1 + dim1.mx
        x = assign(l1, l2, x, grid.x, 'x', proc)
        dx = assign1(dx, grid.dx, 'dx', proc)
        Lx = assign1(Lx, grid.Lx, 'Lx', proc)
        dx_1 = assign(l1, l2, dx_1, grid.dx_1, 'dx_1', proc)
        dx_tilde = assign(l1, l2, dx_tilde, grid.dx_tilde, 'dx_tilde', proc)
        # y direction
        m1 = dim1.iprocy * dim1.ny
        m2 = m1 + dim1.my
        y = assign(m1, m2, y, grid.y, 'y', proc)
        dy = assign1(dy, grid.dy, 'dy', proc)
        Ly = assign1(Ly, grid.Ly, 'Ly', proc)
        dy_1 = assign(m1, m2, dy_1, grid.dy_1, 'dy_1', proc)
        dy_tilde = assign(m1, m2, dy_tilde, grid.dy_tilde, 'dy_tilde', proc)
        # z direction
        n1 = dim1.iprocz * dim1.nz
        n2 = n1 + dim1.mz
        z = assign(n1, n2, z, grid.z, 'z', proc)
        dz = assign1(dz, grid.dz, 'dz', proc)
        Lz = assign1(Lz, grid.Lz, 'Lz', proc)
        dz_1 = assign(n1, n2, dz_1, grid.dz_1, 'dz_1', proc)
        dz_tilde = assign(n1, n2, dz_tilde, grid.dz_tilde, 'dz_tilde', proc)
    # Trim the ghost cells if requested.
    if trim:
        x = x[dim.nghost:-dim.nghost]
        y = y[dim.nghost:-dim.nghost]
        z = z[dim.nghost:-dim.nghost]
        dx_1 = dx_1[dim.nghost:-dim.nghost]
        dy_1 = dy_1[dim.nghost:-dim.nghost]
        dz_1 = dz_1[dim.nghost:-dim.nghost]
        dx_tilde = dx_tilde[dim.nghost:-dim.nghost]
        dy_tilde = dy_tilde[dim.nghost:-dim.nghost]
        dz_tilde = dz_tilde[dim.nghost:-dim.nghost]
    # Define and return a named tuple.
    Grid = namedtuple('Grid', ['x', 'y', 'z', 'dx', 'dy', 'dz', 'Lx', 'Ly', 'Lz', 'dx_1', 'dy_1', 'dz_1',
                               'dx_tilde', 'dy_tilde', 'dz_tilde'])
    return Grid(x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, Lx=Lx, Ly=Ly, Lz=Lz, dx_1=dx_1, dy_1=dy_1, dz_1=dz_1,
                dx_tilde=dx_tilde, dy_tilde=dy_tilde, dz_tilde=dz_tilde)
#=======================================================================
def parameters(datadir='./data', par2=False, warning=True):
    """Returns runtime parameters.

    Keyword Arguments:
        datadir
            Name of the data directory.
        par2
            Read param2.nml if True; read param.nml otherwise.
        warning
            Give warning messages or not.
    """
    # Chao-Chin Yang, 2015-02-23
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
            elif warning:
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
    # Chao-Chin Yang, 2014-12-31
    from collections import namedtuple
    # Read pdim.dat.
    f = open(datadir.strip() + '/pdim.dat')
    a = f.read().rsplit()
    f.close()
    # Extract the numbers.
    npar, mpvar, npar_stalk, mpaux = (int(b) for b in a)
    # Define and return a named tuple.
    ParticleNumbers = namedtuple('ParticleNumbers', ['npar', 'mpvar', 'npar_stalk', 'mpaux'])
    return ParticleNumbers(npar=npar, mpvar=mpvar, npar_stalk=npar_stalk, mpaux=mpaux)
#=======================================================================
def proc_avg2d(datadir='./data', direction='z', proc=0):
    """Returns the time series of one chunk of the 2D averages from one
    process.

    Keyword Arguments:
        datadir
            Name of the data directory.
        direction
            Direction of the average: 'x', 'y', or 'z'.
        proc
            Process ID.
    """
    # Chao-Chin Yang, 2015-03-27
    import numpy as np
    import os.path
    from struct import unpack, calcsize
    # Find the dimensions and the precision.
    dim = proc_dim(datadir=datadir, proc=proc)
    if dim.double_precision:
        dtype = np.float64
        fmt = 'd'
    else:
        dtype = np.float32
        fmt = 'f'
    nb = calcsize(fmt)
    # Check the direction of average.
    if direction == 'x':
        n1 = dim.ny
        n2 = dim.nz
    elif direction == 'y':
        n1 = dim.nx
        n2 = dim.nz
    elif direction == 'z':
        n1 = dim.nx
        n2 = dim.ny
    else:
        raise ValueError("Keyword direction only accepts 'x', 'y', or 'z'. ")
    # Read the names of the averages.
    var = varname(datadir=os.path.dirname(datadir.rstrip('/')), filename=direction.strip()+'aver.in')
    nvar = len(var)
    adim = np.array((n1, n2, nvar))
    nbavg = nb * adim.prod()
    # Open file and define data stream.
    f = open(datadir.strip() + '/proc' + str(proc) + '/' + direction.strip() + 'averages.dat', 'rb')
    def get_time():
        if len(f.read(hsize)) > 0:
            buf = f.read(nb)
            if len(buf) > 0:
                f.read(hsize)
                return unpack(fmt, buf)[0]
            else:
                raise EOFError("incompatible data file")
    def get_avg():
        f.read(hsize)
        try:
            a = np.frombuffer(f.read(nbavg), dtype=dtype).reshape(adim, order='F')
        except:
            raise EOFError("incompatible data file")
        f.read(hsize)
        return a
    # Check the data size.
    nt = 0
    while get_time() is not None:
        get_avg()
        nt += 1
    f.seek(0)
    # Read the data.
    t = np.zeros(nt)
    avg = np.core.records.array(nvar * [np.zeros((n1,n2,nt))], names=var)
    for i in range(nt):
        t[i] = get_time()
        a = get_avg()
        for j in range(nvar):
            avg[var[j]][:,:,i] = a[:,:,j]
    # Close file.
    f.close()
    return t, avg.view(np.recarray)
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
def proc_var(datadir='./data', dim=None, par=None, proc=0, varfile='var.dat'):
    """Returns the patch of one snapshot saved by one process.

    Keyword Arguments:
        datadir
            Name of the data directory.
        dim
            Dimensions supplied by proc_dim().  If None, proc_dim() will
            be called.
        par
            Parameters supplied by parameters().  If None, parameters()
            will be called.
        proc
            Process ID.
        varfile
            Name of the snapshot file.
    """
    # Chao-Chin Yang, 2015-02-20
    from collections import namedtuple
    import numpy as np
    from struct import unpack, calcsize
    # Check the dimensions and precision.
    if par is None:
        par = parameters(datadir=datadir)
    if dim is None:
        dim = proc_dim(datadir=datadir, proc=proc)
    if dim.double_precision:
        dtype = np.float64
        fmt = 'd'
    else:
        dtype = np.float32
        fmt = 'f'
    if par.lwrite_aux:
        adim = np.array((dim.mx, dim.my, dim.mz, dim.mvar + dim.maux))
    else:
        adim = np.array((dim.mx, dim.my, dim.mz, dim.mvar))
    nb = calcsize(fmt)
    # Read the snapshot.
    f = open(datadir.strip() + '/proc' + str(proc) + '/' + varfile.strip(), 'rb')
    f.read(hsize)
    a = np.frombuffer(f.read(nb*adim.prod()), dtype=dtype).reshape(adim, order='F')
    f.read(2*hsize)
    t = unpack(fmt, f.read(nb))[0]
    x = np.frombuffer(f.read(nb*dim.mx), dtype=dtype)
    y = np.frombuffer(f.read(nb*dim.my), dtype=dtype)
    z = np.frombuffer(f.read(nb*dim.mz), dtype=dtype)
    dx, dy, dz = unpack(3*fmt, f.read(3*nb))
    if par.lshear:
        deltay = unpack(fmt, f.read(nb))[0]
    else:
        deltay = None
    f.read(hsize)
    f.close()
    # Define and return a named tuple.
    Var = namedtuple('Var', ['f', 't', 'x', 'y', 'z', 'dx', 'dy', 'dz', 'deltay'])
    return Var(f=a, t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, deltay=deltay)
#=======================================================================
def time_series(datadir='./data'):
    """Returns a NumPy recarray from the time series.

    Keyword Arguments:
        datadir
            Name of the data directory
    """
    # Chao-Chin Yang, 2013-05-13
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
def var(datadir='./data', varfile='var.dat', verbose=True):
    """Returns one snapshot.

    Keyword Arguments:
        datadir
            Name of the data directory.
        varfile
            Name of the snapshot file.
        verbose
            Verbose output or not.
    """
    # Chao-Chin Yang, 2015-02-23
    from collections import namedtuple
    import numpy as np
    # Get the dimensions.
    par = parameters(datadir=datadir, warning=verbose)
    dim = dimensions(datadir=datadir)
    var = varname(datadir=datadir)
    if par.lwrite_aux:
        mvar = dim.mvar + dim.maux
        if len(var) == dim.mvar:
            for i in range(dim.maux):
                var.append("aux" + str(i+1))
    else:
        mvar = dim.mvar
    # Check the precision.
    if dim.double_precision:
        dtype = np.float64
    else:
        dtype = np.float32
    # Allocate arrays.
    fdim = [dim.nxgrid, dim.nygrid, dim.nzgrid, mvar]
    f = np.zeros(fdim, dtype=dtype)
    x = np.zeros((dim.nxgrid,), dtype=dtype)
    y = np.zeros((dim.nygrid,), dtype=dtype)
    z = np.zeros((dim.nzgrid,), dtype=dtype)
    t, dx, dy, dz = 0.0, 0.0, 0.0, 0.0
    if par.lshear:
        deltay = 0.0
    else:
        deltay = None
    # Define functions for assigning local coordinates to global.
    def assign(l1, l2, xg, xl, label, proc):
        indices = xg[l1:l2] != 0.0
        if any(xg[l1:l2][indices] != xl[indices]):
            print("Warning: inconsistent ", label, " for process", proc)
        xg[l1:l2][~indices] = xl[~indices]
        return xg
    def assign1(old, new, label, proc):
        if old != 0.0 and old != new:
            print("Warning: inconsistent ", label, " for process", proc)
        return new
    # Loop over each process.
    for proc in range(dim.nprocx * dim.nprocy * dim.nprocz):
        # Read data.
        if verbose:
            print("Reading", datadir + "/proc" + str(proc) + "/" + varfile, "...")
        dim1 = proc_dim(datadir=datadir, proc=proc)
        snap1 = proc_var(datadir=datadir, par=par, proc=proc, varfile=varfile) 
        # Assign the data to the corresponding block.
        l1 = dim1.iprocx * dim1.nx
        l2 = l1 + dim1.nx
        m1 = dim1.iprocy * dim1.ny
        m2 = m1 + dim1.ny
        n1 = dim1.iprocz * dim1.nz
        n2 = n1 + dim1.nz
        f[l1:l2,m1:m2,n1:n2,:] = snap1.f[dim1.nghost:-dim1.nghost,dim1.nghost:-dim1.nghost,dim1.nghost:-dim1.nghost,:]
        # Assign the coordinates and other scalars.
        t = assign1(t, snap1.t, 't', proc)
        x = assign(l1, l2, x, snap1.x[dim1.nghost:-dim1.nghost], 'x', proc)
        y = assign(m1, m2, y, snap1.y[dim1.nghost:-dim1.nghost], 'y', proc)
        z = assign(n1, n2, z, snap1.z[dim1.nghost:-dim1.nghost], 'z', proc)
        dx = assign1(dx, snap1.dx, 'dx', proc)
        dy = assign1(dy, snap1.dy, 'dy', proc)
        dz = assign1(dz, snap1.dz, 'dz', proc)
        if par.lshear:
            deltay = assign1(deltay, snap1.deltay, 'deltay', proc)
    # Define and return a named tuple.
    fdim.pop()
    for i in range(fdim.count(1)):
        fdim.remove(1)
    keys = ['t', 'x', 'y', 'z', 'dx', 'dy', 'dz', 'deltay'] + var
    values = [t, x, y, z, dx, dy, dz, deltay]
    for i in range(len(var)):
        values.append(f[:,:,:,i].reshape(fdim))
    Var = namedtuple('Var', keys)
    return Var(**dict(zip(keys, values)))
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
