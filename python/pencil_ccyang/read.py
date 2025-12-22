#!/usr/bin/env python3
#=======================================================================
# read.py
#
# Facilities for reading the Pencil Code data.
#
# Author: Chao-Chin Yang
# Created: 2013-05-06
# Last Modified: 2021-04-19
#=======================================================================
hsize = 4    # Heading and trailing bytes in Fortran binary.
#=======================================================================
class Dict:
    """A class that emulates NamedTuple but is picklable.

    Instance Variables
        All the keyword arguments when constructed.
    """
    # Author: Chao-Chin Yang
    # Created: 2021-03-11
    # Last Modified: 2021-03-16

    def __init__(self, **kw):
        """Initiate the object with all the keyword arguments.

        **kw
            Keyword values to be recorded in the object.
        """
        # Author: Chao-Chin Yang
        # Created: 2021-03-11
        # Last Modified: 2021-03-16

        self._keys = list(kw.keys())
        for key, value in kw.items():
            setattr(self, key, value)

    def keys(self):
        """Returns all the keywords recorded in the object. """
        # Author: Chao-Chin Yang
        # Created: 2021-03-16
        # Last Modified: 2021-03-16

        return self._keys
#=======================================================================
def allprocs_avg2d(datadir='./data', direction='z'):
    """Returns the time series of the 2D averages read from allprocs/.

    Keyword Arguments:
        datadir
            Name of the data directory.
        direction
            Direction of the average: 'x', 'y', or 'z'.
    """
    # Author: Chao-Chin Yang
    # Created: 2022-12-14
    # Last Modified: 2022-12-14
    import numpy as np
    from pathlib import Path
    from struct import unpack

    # Find the dimensions and the precision.
    dim = dimensions(datadir=datadir)
    fmt, dtype, nb = _get_precision(dim)

    # Check the direction of average.
    if direction == 'x':
        n1, n2 = dim.nygrid, dim.nzgrid
    elif direction == 'y':
        n1, n2 = dim.nxgrid, dim.nzgrid
    elif direction == 'z':
        n1, n2 = dim.nxgrid, dim.nygrid
    else:
        raise ValueError("Keyword direction only accepts 'x', 'y', or 'z'. ")

    # Read the names of the averages.
    workdir = Path(datadir).parent.as_posix()
    var = varname(datadir=workdir, filename=direction.strip()+'aver.in')
    nvar = len(var)
    adim = np.array((n1, n2, nvar))
    nbavg = nb * adim.prod()

    # Open file.
    path = Path(datadir) / "allprocs" / (direction.strip() + "averages.dat")
    f = open(path, 'rb')

    # Define data stream.
    def get_time():
        buf = f.read(nb)
        if len(buf) > 0:
            return unpack(fmt, buf)[0]
        else:
            return None

    def get_avg():
        try:
            a = np.frombuffer(f.read(nbavg), dtype=dtype)
            a = a.reshape(adim, order='F')
        except:
            raise EOFError("incompatible data file")
        return a

    # Read the data.
    times = []
    avg = [[] for i in range(nvar)]
    while True:
        t1 = get_time()
        if t1 is None: break
        times.append(t1)
        a = get_avg()
        for i in range(nvar):
            avg[i].append(a[:,:,i])

    # Close file.
    f.close()

    # Return the data.
    times = np.array(times)
    avg = [np.stack(a, axis=0) for a in avg]
    avg = np.rec.array(avg, names=var)
    return times, avg
#=======================================================================
def allprocs_grid(datadir='./data', dim=None):
    """Returns the grid under allprocs/.

    Keyword Arguments:
        datadir
            Name of the data directory.
        dim
            Dimensions supplied by dimensions().  If None, they will be
            read in at runtime.
    """
    # Author: Chao-Chin Yang
    # Created: 2020-11-15
    # Last Modified: 2021-03-11
    import numpy as np
    from struct import unpack, calcsize

    # Check the dimensions and precision.
    if dim is None: dim = dimensions(datadir=datadir)
    fmt, dtype, nb = _get_precision(dim)

    # Read grid.dat.
    f = open(datadir.strip() + "/allprocs/grid.dat", "rb")
    f.read(hsize)
    t = unpack(fmt, f.read(nb))[0]
    x = np.frombuffer(f.read(nb*dim.mxgrid), dtype=dtype)
    y = np.frombuffer(f.read(nb*dim.mygrid), dtype=dtype)
    z = np.frombuffer(f.read(nb*dim.mzgrid), dtype=dtype)
    dx, dy, dz = unpack(3*fmt, f.read(3*nb))
    f.read(2*hsize)
    dx, dy, dz = unpack(3*fmt, f.read(3*nb))
    f.read(2*hsize)
    Lx, Ly, Lz = unpack(3*fmt, f.read(3*nb))
    f.read(2*hsize)
    dx_1 = np.frombuffer(f.read(nb*dim.mxgrid), dtype=dtype)
    dy_1 = np.frombuffer(f.read(nb*dim.mygrid), dtype=dtype)
    dz_1 = np.frombuffer(f.read(nb*dim.mzgrid), dtype=dtype)
    f.read(2*hsize)
    dx_tilde = np.frombuffer(f.read(nb*dim.mxgrid), dtype=dtype)
    dy_tilde = np.frombuffer(f.read(nb*dim.mygrid), dtype=dtype)
    dz_tilde = np.frombuffer(f.read(nb*dim.mzgrid), dtype=dtype)
    f.read(hsize)
    f.close()

    # Return the data.
    return Dict(x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, Lx=Lx, Ly=Ly, Lz=Lz,
                dx_1=dx_1, dy_1=dy_1, dz_1=dz_1,
                dx_tilde=dx_tilde, dy_tilde=dy_tilde, dz_tilde=dz_tilde)
#=======================================================================
def allprocs_pvar(datadir='./data', pdim=None, pvarfile='pvar.dat'):
    """Returns the particles in one snapshot under allprocs/.

    Keyword Arguments:
        datadir
            Name of the data directory.
        pdim
            Particle dimensions returned by pardim().  If None, they
            will be read in at runtime.
        pvarfile
            Name of the snapshot file.
    """
    # Author: Chao-Chin Yang
    # Created: 2020-11-08
    # Last Modified: 2021-03-11
    import numpy as np
    from struct import calcsize, unpack

    # Check the dimensions and precision.
    fmt, dtype, nb = _get_precision(dimensions(datadir=datadir))
    fmti, nbi = 'i', calcsize('i')
    if pdim is None:
        pdim = pardim(datadir=datadir)
    mparray = pdim.mpvar + pdim.mpaux

    # Read number of particles.
    f = open(datadir.strip() + "/allprocs/" + pvarfile.strip(), 'rb')
    npar = unpack(fmti, f.read(nbi))[0]

    # Read particle data.
    if npar > 0:
        fpdim = np.array((npar, mparray))
        ipar = np.frombuffer(f.read(nbi*npar), dtype='i4')
        fp = np.frombuffer(f.read(nb*fpdim.prod()),
                           dtype=dtype).reshape(fpdim, order='F')
    else:
        fp = None
        ipar = None

    # Read time.
    t = unpack(fmt, f.read(nb))[0]
    f.close()

    # Return the data.
    return Dict(npar=npar, fp=fp, ipar=ipar, t=t)
#=======================================================================
def allprocs_var(datadir='./data', dim=None, par=None, varfile='var.dat'):
    """Returns one snapshot under allprocs/.

    Keyword Arguments:
        datadir
            Name of the data directory.
        dim
            Dimensions supplied by dimensions().  If None, the
            dimensions will be read in at runtime.
        par
            Parameters supplied by parameters().  If None, the
            parameters will be read in at runtime.
        varfile
            Name of the snapshot file.
    """
    # Author: Chao-Chin Yang
    # Created: 2020-11-03
    # Last Modified: 2021-03-11
    import numpy as np
    from struct import unpack

    # Get the dimensions and precision.
    if dim is None:
        dim = dimensions(datadir=datadir)
    fmt, dtype, nb = _get_precision(dim)

    # Get the total number of variables.
    nvar = dim.mvar
    if par is None:
        par = parameters(datadir=datadir)
    if par.lwrite_aux:
        nvar += dim.maux

    # Allocate memory space.
    adim = np.array((dim.mxgrid, dim.mygrid, dim.mzgrid, nvar))

    # Read the snapshot.
    f = open(datadir.strip() + '/allprocs/' + varfile.strip(), 'rb')
    a = np.frombuffer(f.read(nb*adim.prod()), dtype=dtype)
    a = a.reshape(adim, order='F')
    f.read(hsize)
    t = unpack(fmt, f.read(nb))[0]
    x = np.frombuffer(f.read(nb*dim.mxgrid), dtype=dtype)
    y = np.frombuffer(f.read(nb*dim.mygrid), dtype=dtype)
    z = np.frombuffer(f.read(nb*dim.mzgrid), dtype=dtype)
    dx, dy, dz = unpack(3*fmt, f.read(3*nb))
    if par.lshear:
        deltay = unpack(fmt, f.read(nb))[0]
    else:
        deltay = None
    f.read(hsize)
    f.close()

    # Return the data.
    return Dict(f=a, t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, deltay=deltay)
#=======================================================================
def avg1d(datadir='./data', plane='xy', tsize=None, unformatted=True,
          verbose=True):
    """Returns the time series of 1D averages.

    Keyword Arguments:
        datadir
            Name of the data directory.
        plane
            Plane of average: 'xy', 'xz', or 'yz'.
        tsize
            If not None, the time series is interpolated onto regular
            time intervals with tsize elements and endpoints included.
        unformatted
            Set True if the data file is written in binary format.
        verbose
            Whether or not to print information.
    """
    # Author: Chao-Chin Yang
    # Created: 2013-10-28
    # Last Modified: 2018-02-09
    import numpy as np
    import os.path
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

    # Check the file format.
    mode = 'r'
    if unformatted:
        fmt, dtype, nb = _get_precision(dim)
        mode += 'b'

    # Read the names of the variables.
    var = varname(datadir=os.path.dirname(datadir),
                  filename=plane.strip()+'aver.in')
    nvar = len(var)

    # Open file and define data stream.
    f = open(datadir.strip() + '/' + plane.strip() + 'averages.dat', mode)
    def fetch(nval):
        if unformatted:
            f.read(hsize)
            a = np.fromfile(f, dtype=dtype, count=nval)
            f.read(hsize)
        else:
            a = np.fromfile(f, count=nval, sep=' ')
        return a

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
    avg = np.rec.array(len(var) * [np.zeros((nt,nc))], names=var)
    for i in range(nt):
        t[i] = fetch(1)
        block = fetch(nblock).reshape((nc, nvar), order='F')
        for j, v in enumerate(var):
            avg[v][i,:] = block[:,j]

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
        avgi = np.rec.array(len(var) * [np.zeros((tsize,nc))], names=var)
        for v in var:
            for k in range(nc):
                avgi[v][:,k] = interp1d(t, avg[v][:,k])(ti)
        t, avg = ti, avgi

    return t, avg
#=======================================================================
def avg2d(datadir='./data', direction='z', par=None):
    """Returns the time series of the 2D averages.

    Keyword Arguments:
        datadir
            Name of the data directory.
        direction
            Direction of average: 'x', 'y', or 'z'.
        par
            Parameters supplied by parameters().  If None, the
            parameters will be read in at runtime.
    """
    # Author: Chao-Chin Yang
    # Created: 2015-03-27
    # Last Modified: 2022-12-14
    import numpy as np

    # Delegate to allprocs_avg2d() if using MPI IO.
    if par is None:
        par = parameters(datadir=datadir)
    if par.io_strategy == "MPI-IO":
        return allprocs_avg2d(datadir=datadir, direction=direction)

    # Get the dimensions.
    dim = dimensions(datadir=datadir)
    if direction == 'x':
        n1 = dim.nygrid
        n2 = dim.nzgrid
    elif direction == 'y':
        n1 = dim.nxgrid
        n2 = dim.nzgrid
    elif direction == 'z':
        n1 = dim.nxgrid
        n2 = dim.nygrid
    else:
        raise ValueError("Keyword direction only accepts 'x', 'y', or 'z'. ")

    # Allocate data arrays.
    t, avg1 = proc_avg2d(datadir=datadir, direction=direction)
    nt = len(t)
    var = avg1.dtype.names
    avg = np.rec.array(len(var) * [np.zeros((nt,n1,n2))], names=var)

    # Read the averages from each process and assemble the data.
    for proc in range(dim.nprocx * dim.nprocy * dim.nprocz):
        dim1 = proc_dim(datadir=datadir, proc=proc)
        if direction == 'x' and dim1.iprocx == 0:
            iproc1 = dim1.iprocy
            iproc2 = dim1.iprocz
            n1 = dim1.ny
            n2 = dim1.nz
        elif direction == 'y' and dim1.iprocy == 0:
            iproc1 = dim1.iprocx
            iproc2 = dim1.iprocz
            n1 = dim1.nx
            n2 = dim1.nz
        elif direction == 'z' and dim1.iprocz == 0:
            iproc1 = dim1.iprocx
            iproc2 = dim1.iprocy
            n1 = dim1.nx
            n2 = dim1.ny
        else:
            continue
        t1, avg1 = proc_avg2d(datadir=datadir, direction=direction, proc=proc)
        if any(t1 != t):
            print("Warning: inconsistent time series. ")
        avg[:,iproc1*n1:(iproc1+1)*n1,iproc2*n2:(iproc2+1)*n2] = avg1

    # Return the time and the averages.
    return t, avg
#=======================================================================
def dimensions(datadir="./data"):
    """Returns the dimensions of the Pencil Code data from datadir.

    Keyword Arguments:
        datadir
            Name of the data directory.
    """
    # Author: Chao-Chin Yang
    # Created: 2013-10-23
    # Last Modified: 2021-03-11

    # Read dim.dat.
    f = open(datadir.strip() + "/dim.dat")
    a = f.read().rsplit()
    f.close()

    # Sanity check
    if a[6] == '?':
        print("Warning: unknown data precision. ")
    if not a[7] == a[8] == a[9]:
        raise Exception(
                "unequal number of ghost zones in different directions. ")

    # Extract the dimensions.
    mxgrid, mygrid, mzgrid, mvar, maux, mglobal = (int(b) for b in a[0:6])
    double_precision = a[6] == 'D'
    nghost = int(a[7])
    nxgrid, nygrid, nzgrid = (int(b) - 2 * nghost for b in a[0:3])
    nprocx, nprocy, nprocz = (int(b) for b in a[10:13])
    procz_last = int(a[13]) == 1

    # Return the data.
    return Dict(nxgrid=nxgrid, nygrid=nygrid, nzgrid=nzgrid, nghost=nghost,
               mxgrid=mxgrid, mygrid=mygrid, mzgrid=mzgrid,
               mvar=mvar, maux=maux, mglobal=mglobal,
               double_precision=double_precision,
               nprocx=nprocx, nprocy=nprocy, nprocz=nprocz,
               procz_last=procz_last)
#=======================================================================
def grid(datadir='./data', interface=False, par=None, trim=True):
    """Returns the coordinates and their derivatives of the grid.

    Keyword Arguments:
        datadir
            Name of the data directory.
        interface
            If True, the coordinates of the interfaces between cells are
            returned; otherwise, those of the cell centers are returned.
        par
            Parameters supplied by parameters() needed when the keyword
            interface is True.  If None, parameters() will be called.
        trim
            Whether or not to trim the ghost cells.  If the keyword
            interface is True, ghost cells are automatically trimmed and
            this keyword has no effect.
    """
    # Author: Chao-Chin Yang
    # Created: 2014-11-02
    # Last Modified: 2021-03-11
    import numpy as np

    # Get the dimensions and parameters.
    dim = dimensions(datadir=datadir)
    if par is None: par = parameters(datadir=datadir)
    allprocs = par.io_strategy == "MPI-IO"

    if allprocs:
        # Read grid from under allprocs/ and unpack.
        grid = allprocs_grid(datadir=datadir, dim=dim)
        x, y, z = grid.x, grid.y, grid.z
        dx, dy, dz = grid.dx, grid.dy, grid.dz
        Lx, Ly, Lz = grid.Lx, grid.Ly, grid.Lz
        dx_1, dy_1, dz_1 = grid.dx_1, grid.dy_1, grid.dz_1
        dx_tilde = grid.dx_tilde
        dy_tilde = grid.dy_tilde
        dz_tilde = grid.dz_tilde

    else:
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
    if interface or trim:
        x = x[dim.nghost:-dim.nghost]
        y = y[dim.nghost:-dim.nghost]
        z = z[dim.nghost:-dim.nghost]
        dx_1 = dx_1[dim.nghost:-dim.nghost]
        dy_1 = dy_1[dim.nghost:-dim.nghost]
        dz_1 = dz_1[dim.nghost:-dim.nghost]
        dx_tilde = dx_tilde[dim.nghost:-dim.nghost]
        dy_tilde = dy_tilde[dim.nghost:-dim.nghost]
        dz_tilde = dz_tilde[dim.nghost:-dim.nghost]

    # Find the interfaces if requested.
    if interface:
        f = (lambda x, dx_1, x0, x1:
                np.concatenate((
                    (x0,),
                    (x[:-1] * dx_1[:-1] + x[1:] * dx_1[1:]) /
                        (dx_1[:-1] + dx_1[1:]),
                    (x1,)))
                if len(x) > 1 else x)
        x = f(x, dx_1, par.xyz0[0], par.xyz1[0])
        y = f(y, dy_1, par.xyz0[1], par.xyz1[1])
        z = f(z, dz_1, par.xyz0[2], par.xyz1[2])

    # Return the data.
    return Dict(x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, Lx=Lx, Ly=Ly, Lz=Lz,
                dx_1=dx_1, dy_1=dy_1, dz_1=dz_1,
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
    # Author: Chao-Chin Yang
    # Created: 2013-10-31
    # Last Modified: 2021-03-11

    # Function to convert a string to the correct type.
    def convert(v):
        if v == 'T' or v == ".TRUE.":
            return True
        elif v == 'F' or v == ".FALSE.":
            return False
        elif ',' in v:
            re, im = v.split(',')
            return complex(float(re), float(im))
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
        if '(' in v:
            v = (v + ',').split("),")[:-1]
        else:
            v = v.split(',')
        u = []
        for w in v:
            w = w.strip()
            if '*' in w:
                nrep, sep, x = w.partition('*')
                u += int(nrep) * [convert(x.lstrip('('))]
            else:
                u += [convert(w.lstrip('('))]
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
    line = ""
    for next_line in f:
        lead = next_line.lstrip()
        if len(lead) > 0: lead = lead[0]
        if '=' in next_line or lead == '&' or lead == '/':
            line = line.strip().replace('\n', '')
            if len(line) > 0:
                k, s, v = line.partition('=')
                k = k.strip().lower()
                if "'" not in v: v = v.replace(' ', '')
                if k not in keys:
                    keys.append(k)
                    values.append(parse(v.strip(" ,\n")))
                elif warning:
                    print("Duplicate parameter:", k, '=', v.strip(" ,\n"))
            line = next_line if '=' in next_line else ""
        else:
            line += next_line
    f.close()

    # Return the parameters.
    return Dict(**dict(zip(keys, values)))
#=======================================================================
def pardim(datadir='./data'):
    """Returns the numbers associated with particles.

    Keyword Arguments:
        datadir
            Name of the data directory
    """
    # Author: Chao-Chin Yang
    # Created: 2014-02-11
    # Last Modified: 2021-03-11

    # Read pdim.dat.
    f = open(datadir.strip() + '/pdim.dat')
    a = f.read().rsplit()
    f.close()

    # Extract the numbers.
    if len(a) == 4:
        npar, mpvar, npar_stalk, mpaux = (int(b) for b in a)
    else:  # for backward compatibility.
        npar, mpvar, npar_stalk = (int(b) for b in a)
        mpaux = 0

    # Return the data.
    return Dict(npar=npar, mpvar=mpvar, npar_stalk=npar_stalk, mpaux=mpaux)
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
    # Chao-Chin Yang, 2015-04-20
    import numpy as np
    import os.path
    from struct import unpack
    # Find the dimensions and the precision.
    dim = proc_dim(datadir=datadir, proc=proc)
    fmt, dtype, nb = _get_precision(dim)
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
    avg = np.rec.array(nvar * [np.zeros((nt,n1,n2))], names=var)
    for i in range(nt):
        t[i] = get_time()
        a = get_avg()
        for j in range(nvar):
            avg[var[j]][i,:,:] = a[:,:,j]
    # Close file.
    f.close()
    return t, avg
#=======================================================================
def proc_dim(datadir='./data', proc=0):
    """Returns the dimensions of the data from one process.

    Keyword Arguments:
        datadir
            Name of the data directory
        proc
            Process ID
    """
    # Author: Chao-Chin Yang
    # Created: 2013-10-16
    # Last Modified: 2021-03-11

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

    # Return the data.
    return Dict(nx=nx, ny=ny, nz=nz, nghost=nghost, mx=mx, my=my, mz=mz,
            mvar=mvar, maux=maux, mglobal=mglobal,
            double_precision=double_precision,
            iprocx=iprocx, iprocy=iprocy, iprocz=iprocz)
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
    # Author: Chao-Chin Yang
    # Created: 2014-10-29
    # Last Modified: 2021-03-11
    import numpy as np
    from struct import unpack, calcsize

    # Check the dimensions and precision.
    if dim is None:
        dim = proc_dim(datadir=datadir, proc=proc)
    fmt, dtype, nb = _get_precision(dim)

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

    # Return the data.
    return Dict(x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, Lx=Lx, Ly=Ly, Lz=Lz,
                dx_1=dx_1, dy_1=dy_1, dz_1=dz_1,
                dx_tilde=dx_tilde, dy_tilde=dy_tilde, dz_tilde=dz_tilde)
#=======================================================================
def proc_pvar(datadir='./data', pdim=None, proc=0, pvarfile='pvar.dat'):
    """Returns the particles in one snapshot held by one process.

    Keyword Arguments:
        datadir
            Name of the data directory.
        pdim
            Particle dimensions returned by pardim().  If None, they
            will be read in at runtime.
        proc
            Process ID.
        pvarfile
            Name of the snapshot file.
    """
    # Author: Chao-Chin Yang
    # Created: 2015-07-12
    # Last Modified: 2021-03-11
    import numpy as np
    from struct import calcsize, unpack

    # Check the dimensions and precision.
    dim = proc_dim(datadir=datadir, proc=proc)
    fmt, dtype, nb = _get_precision(dim)
    fmti, nbi = 'i', calcsize('i')
    if pdim is None:
        pdim = pardim(datadir=datadir)
    mparray = pdim.mpvar + pdim.mpaux

    # Read number of particles.
    f = open(datadir.strip()+'/proc'+str(proc)+'/'+pvarfile.strip(), 'rb')
    f.read(hsize)
    npar_loc = unpack(fmti, f.read(nbi))[0]

    # Read particle data.
    if npar_loc > 0:
        fpdim = np.array((npar_loc, mparray))
        f.read(2*hsize)
        ipar = np.frombuffer(f.read(nbi*npar_loc), dtype='i4')
        f.read(2*hsize)
        fp = np.frombuffer(f.read(nb*fpdim.prod()),
                           dtype=dtype).reshape(fpdim, order='F')
    else:
        fp = None
        ipar = None

    # Read time and grid.
    f.read(2*hsize)
    t = unpack(fmt, f.read(nb))[0]
    x = np.frombuffer(f.read(nb*dim.mx), dtype=dtype)
    y = np.frombuffer(f.read(nb*dim.my), dtype=dtype)
    z = np.frombuffer(f.read(nb*dim.mz), dtype=dtype)
    dx, dy, dz = unpack(3*fmt, f.read(3*nb))
    f.read(hsize)
    f.close()

    # Return the data.
    return Dict(npar_loc=npar_loc, fp=fp, ipar=ipar,
                t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz)
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
    # Author: Chao-Chin Yang
    # Created: 2014-12-03
    # Last Modified: 2021-03-11
    import numpy as np
    from struct import unpack

    # Check the dimensions and precision.
    if dim is None:
        dim = proc_dim(datadir=datadir, proc=proc)
    fmt, dtype, nb = _get_precision(dim)
    if par is None:
        par = parameters(datadir=datadir)
    if par.lwrite_aux:
        adim = np.array((dim.mx, dim.my, dim.mz, dim.mvar + dim.maux))
    else:
        adim = np.array((dim.mx, dim.my, dim.mz, dim.mvar))

    # Read the snapshot.
    f = open(datadir.strip()+'/proc'+str(proc)+'/'+varfile.strip(), 'rb')
    f.read(hsize)
    a = np.frombuffer(f.read(nb*adim.prod()), dtype=dtype)
    a = a.reshape(adim, order='F')
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

    # Return the data.
    return Dict(f=a, t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, deltay=deltay)
#=======================================================================
def pvar(allprocs=True, datadir='./data', ivar=None, pvarfile='pvar.dat',
         verbose=True):
    """Returns particles in one snapshot.

    Keyword Arguments:
        allprocs
            If True, the snapshot is read under allprocs/.
        datadir
            Name of the data directory.
        ivar
            If not None, an integer specifying the snapshot number and
            the argument pvarfile takes no effect.
        pvarfile
            Name of the snapshot file.
        verbose
            Verbose output or not.
    """
    # Author: Chao-Chin Yang
    # Created: 2015-07-12
    # Last Modified: 2021-03-11
    import numpy as np

    # Get the dimensions.
    dim = dimensions(datadir=datadir)
    pdim = pardim(datadir=datadir)
    var = varname(datadir=datadir, filename="pvarname.dat")
    nvar = len(var)
    if nvar != pdim.mpvar + pdim.mpaux:
        raise RuntimeError("Inconsistent number of variables. ")

    # Check the precision.
    fmt, dtype, nb = _get_precision(dim)

    # Get the file name of the snapshot.
    if ivar is not None:
        pvarfile = "PVAR{}".format(ivar)

    # Allocate arrays.
    fp = np.zeros((pdim.npar, nvar), dtype=dtype)
    exist = np.full((pdim.npar,), False)
    t = None

    if allprocs:
        # Read data under allprocs/.
        pvar = allprocs_pvar(datadir=datadir, pdim=pdim, pvarfile=pvarfile)
        t = pvar.t

        # Check numbers of particles.
        if pvar.npar != pdim.npar:
            print("pvar.npar, pdim.npar = {}, {}".format(pvar.npar, pdim.npar))
            raise RuntimeError("Inconsistent numbers of particles. ")

        # Check the integrity of the particle IDs.
        if len(np.unique(pvar.ipar)) != len(pvar.ipar):
            raise RuntimeError("Duplicate particles IDs. ")
        ipar = pvar.ipar - 1

        # Unpack the data.
        fp[ipar,:] = pvar.fp
        exist[ipar] = True

    else:
        # Loop over each process.
        for proc in range(dim.nprocx * dim.nprocy * dim.nprocz):

            # Read data.
            if verbose:
                print("Reading", datadir.strip() + "/proc" + str(proc) + "/" +
                                 pvarfile.strip(), "...")
            p = proc_pvar(datadir=datadir, pdim=pdim, proc=proc,
                          pvarfile=pvarfile) 

            # Check the consistency of time.
            if t is None:
                t = p.t
            elif p.t != t:
                raise RuntimeError(
                        "Inconsistent time stamps in multiple processes. ")

            # Return if no particles exist.
            if p.npar_loc <= 0: continue
            ipar = p.ipar - 1

            # Check the integrity of the particle IDs.
            if any(exist[ipar]):
                raise RuntimeError(
                        "Duplicate particles in multiple processes. ")

            # Assign the particle data.
            fp[ipar,:] = p.fp
            exist[ipar] = True

    # Final integrity check.
    if not all(exist):
        raise RuntimeError("Missing some particles. ")

    # Return the data.
    pvar = dict(npar=pdim.npar, t=t)
    for i, v in enumerate(var):
        pvar[v.lstrip('i')] = fp[:,i]
    return Dict(**pvar)
#=======================================================================
def slices(field, datadir="./data", return_pos=False):
    """Reads the video slices.

    Returned Values:
        1. Array of time points.
        2. Slice planes of field at each time.
        3. (optional) Slice positions of each plane at each time.

    Positional Argument:
        field
            Field name.

    Keyword Argument:
        datadir
            Name of the data directory.
        return_pos
            Also return slice positions if True.
    """
    # Author: Chao-Chin Yang
    # Created: 2015-04-21
    # Last Modified: 2021-03-11
    from glob import glob
    import numpy as np
    import os
    from struct import unpack

    # Determine where the slice files should be searched.
    par = parameters(datadir=datadir, warning=False)
    allprocs = par.io_strategy == "MPI-IO"
    if allprocs:
        prefix = datadir + "/allprocs/slice_" + field
    else:
        prefix = datadir + "/slice_" + field

    # Find the slice planes.
    datadir = datadir.strip()
    field = field.strip()
    planes = []
    for path in glob(prefix + ".*"):
        planes.append(os.path.basename(path).rsplit('.')[1])
    if len(planes) == 0:
        raise IOError("Found no slices. ")

    # Get the dimensions and check the data precision.
    dim = dimensions(datadir=datadir)
    fmt, dtype, nb = _get_precision(dim)

    # Function to return the slice dimensions.
    def slice_dim(plane):
        if plane[0:2] == 'xy':
            sdim = dim.nxgrid, dim.nygrid
        elif plane[0:2] == 'xz':
            sdim = dim.nxgrid, dim.nzgrid
        elif plane[0:2] == 'yz':
            sdim = dim.nygrid, dim.nzgrid
        else:
            raise ValueError("Unknown plane type " + plane)
        return sdim

    # Define function to read file and check size.
    def readf(f, size):
        d = f.read(size)
        if len(d) != size:
            raise EOFError("Incompatible slice file. ")
        return d

    # Process each plane.
    s, t, pos = [], None, []
    for p in planes:
        path = prefix + '.' + p
        print("Reading " + path + "......")

        # read slice data.
        f = open(path, 'rb')
        s1, t1, pos1 = [], [], []
        sdim = slice_dim(p)
        nbs = nb * np.array(sdim).prod()
        while len(f.peek(1)) > 0:
            if not allprocs: readf(f, hsize)
            s1.append(np.frombuffer(readf(f, nbs),
                    dtype=dtype).reshape(sdim, order='F'))
            t1.append(unpack(fmt, readf(f, nb))[0])
            pos1.append(unpack(fmt, readf(f, nb))[0])
            if not allprocs: readf(f, hsize)
        f.close()

        # Check time stamps.
        if t is None:
            t = np.array(t1, dtype=dtype)
        elif np.any(np.array(t1, dtype=dtype) != t):
            raise RuntimeError("Inconsistent time points. ")

        # Record the slice data.
        s.append(np.array(s1))
        pos.append(np.array(pos1))

    # Return the data.
    dicts, dictp = {}, {}
    for k, v in zip(planes, s): dicts[k] = v
    for k, v in zip(planes, p): dictp[k] = v
    s, pos = Dict(**dicts), Dict(**dictp)
    return (t, s, pos) if return_pos else (t, s)
#=======================================================================
def time_series(datadir='./data', unique=False):
    """Returns a NumPy recarray from the time series.

    Keyword Arguments:
        datadir
            Name of the data directory
        unique
            If True, discard duplicate entries.
    """
    # Author: Chao-Chin Yang
    # Created: 2013-05-13
    # Last Modified: 2016-07-16
    import numpy as np
    from io import BytesIO

    # Save the data path.
    path = datadir.strip()

    # Read legend.dat.
    f = open(path + '/legend.dat')
    names = f.read().replace('-', ' ').split()
    f.close()

    # Read time_series.dat.
    f = open(path + '/time_series.dat')
    ts = np.recfromtxt(BytesIO(f.read().encode()), names=names)
    f.close()

    # Discard duplicate entries if requested.
    if unique:
        t, indices = np.unique(ts.t, return_index=True)
        ts = ts[indices]

    return ts
#=======================================================================
def var(allprocs=True, compact=True, datadir='./data', ivar=None, par=None,
        trim=True, unshear=False, varfile='var.dat', verbose=True):
    """Returns one snapshot.

    Keyword Arguments:
        allprocs
            If True, the snapshot is read under allprocs/.
        compact
            If True, dimensions of one are removed.
        datadir
            Name of the data directory.
        ivar
            If not None, an integer specifying the snapshot number and
            the argument varfile takes no effect.
        par
            Parameters supplied by parameters().  If None, parameters()
            will be called.
        trim
            Whether or not to trim the ghost cells.
        unshear
            Whether or not to make x periodic
        varfile
            Name of the snapshot file.
        verbose
            Verbose output or not.
    """
    # Author: Chao-Chin Yang
    # Created: 2014-12-03
    # Last Modified: 2021-03-11
    import numpy as np

    # Get the parameters.
    if par is None:
        par = parameters(datadir=datadir, warning=verbose)
    dim = dimensions(datadir=datadir)
    var = varname(datadir=datadir)
    if par.lwrite_aux:
        mvar = dim.mvar + dim.maux
        for i in range(mvar - len(var)):
            var.append("aux" + str(i+1))
    else:
        mvar = dim.mvar

    # Check the precision.
    fmt, dtype, nb = _get_precision(dim)

    # Get the file name of the snapshot.
    if ivar is not None:
        varfile = "VAR{}".format(ivar)

    # Allocate arrays.
    fdim = [dim.mxgrid, dim.mygrid, dim.mzgrid, mvar]
    f = np.zeros(fdim, dtype=dtype)
    x = np.zeros((dim.mxgrid,), dtype=dtype)
    y = np.zeros((dim.mygrid,), dtype=dtype)
    z = np.zeros((dim.mzgrid,), dtype=dtype)
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

    if allprocs:
        # Read data under allprocs/.
        data = allprocs_var(datadir=datadir, dim=dim, par=par, varfile=varfile)
        f, t = data.f, data.t
        x, y, z = data.x, data.y, data.z
        dx, dy, dz = data.dx, data.dy, data.dz
        if par.lshear:
            deltay = data.deltay

    else:
        # Loop over each process.
        for proc in range(dim.nprocx * dim.nprocy * dim.nprocz):
            # Read data.
            if verbose:
                print("Reading", datadir + "/proc" + str(proc) + "/" + varfile,
                      "...")
            dim1 = proc_dim(datadir=datadir, proc=proc)
            snap1 = proc_var(datadir=datadir, par=par, proc=proc,
                             varfile=varfile) 

            # Assign the data to the corresponding block.
            l1 = dim1.iprocx * dim1.nx
            l2 = l1 + dim1.mx
            m1 = dim1.iprocy * dim1.ny
            m2 = m1 + dim1.my
            n1 = dim1.iprocz * dim1.nz
            n2 = n1 + dim1.mz
            dl1 = 0 if dim1.iprocx == 0 else dim.nghost
            dl2 = 0 if dim1.iprocx == dim.nprocx - 1 else -dim.nghost
            dm1 = 0 if dim1.iprocy == 0 else dim.nghost
            dm2 = 0 if dim1.iprocy == dim.nprocy - 1 else -dim.nghost
            dn1 = 0 if dim1.iprocz == 0 else dim.nghost
            dn2 = 0 if dim1.iprocz == dim.nprocz - 1 else -dim.nghost
            f[l1+dl1:l2+dl2,
              m1+dm1:m2+dm2,
              n1+dn1:n2+dn2,:] = snap1.f[dl1:dim1.mx+dl2,
                                         dm1:dim1.my+dm2,
                                         dn1:dim1.mz+dn2,:]

            # Assign the coordinates and other scalars.
            t = assign1(t, snap1.t, 't', proc)
            x = assign(l1, l2, x, snap1.x, 'x', proc)
            y = assign(m1, m2, y, snap1.y, 'y', proc)
            z = assign(n1, n2, z, snap1.z, 'z', proc)
            dx = assign1(dx, snap1.dx, 'dx', proc)
            dy = assign1(dy, snap1.dy, 'dy', proc)
            dz = assign1(dz, snap1.dz, 'dz', proc)
            if par.lshear:
                deltay = assign1(deltay, snap1.deltay, 'deltay', proc)

    # Trim the ghost cells if requested.
    if trim:
        f = f[dim.nghost:-dim.nghost,
              dim.nghost:-dim.nghost,
              dim.nghost:-dim.nghost,:]
        fdim = [dim.nxgrid, dim.nygrid, dim.nzgrid, mvar]
        x = x[dim.nghost:-dim.nghost]
        y = y[dim.nghost:-dim.nghost]
        z = z[dim.nghost:-dim.nghost]

    if unshear:
        f  = _unshear(f,dim,xax=x,param=par,t=t)

    # Define the returned dimensions.
    fdim.pop()
    if compact:
        for i in range(fdim.count(1)):
            fdim.remove(1)

    # Return the data.
    data = dict(t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, deltay=deltay)
    for i, v in enumerate(var):
        data[v] = f[:,:,:,i].ravel().reshape(fdim)
    return Dict(**data)
#=======================================================================
def varlist(datadir="./data", listname="varN.list"):
    """Reads and returns the names and times of the snapshots.

    Keyword Arguments:
        datadir
            Name of the data directory.
        listname
            Name of the file containing the list of snapshots.

    Returned Values
        varfile
            A list of file names.
        time
            A list of the corresponding times.
    """
    # Author: Chao-Chin Yang
    # Created: 2021-04-19
    # Last Modified: 2021-04-19
    from pathlib import Path

    # Determine where the file list should be located.
    p = Path(datadir.strip()) / "allprocs" / listname.strip()
    if not p.exists():
        p = Path(datadir.strip()) / "proc0" / listname.strip()
        if not p.exists():
            raise RuntimeError("cannot find the list of snapshots. ")

    # Parse the file.
    varfile, time = [], []
    with open(p) as f:
        for line in f:
            v, t = line.split()
            varfile.append(v)
            time.append(float(t))

    return varfile, time
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
#=======================================================================
######  LOCAL FUNCTIONS  ######
#=======================================================================
def _get_precision(dim):
    """Checks the data precision for reading.

    Returned Values:
        fmt
            'd' for double precision and 'f' for single precision.
        dtype
            The corresponding numpy dtype.
        nb
            Number of bytes per floating-point number.

    Positional Argument:
        dim
            Dimensions supplied by dim() or proc_dim().
    """
    # Chao-Chin Yang, 2015-04-20
    import numpy as np
    from struct import calcsize
    if dim.double_precision:
        fmt = 'd'
        dtype = np.float64
    else:
        fmt = 'f'
        dtype = np.float32
    nb = calcsize(fmt)
    return fmt, dtype, nb
#=======================================================================
def _unshear(arr, dim, xax=None, x0=0.0, param=None, t=None, nowrap=False):
    # Wladimir Lyra, 2025-12-21
    from scipy.fft import fft, ifft
    import numpy as np

    if xax is None:
        raise ValueError("_unshear: must provide 1-D array of x coordinates")
    if param is None:
        raise ValueError("param must be provided")
    if t is None:
        raise ValueError("_unshear: must provide t")

    Lx = param.lxyz[0]
    Ly = param.lxyz[1]  # assuming second dimension is y
    deltay = -param.sshear * Lx * t

    # Check dimensions
    if arr.ndim != 4:
        raise ValueError("_unshear only supports 4D arrays (nxgrid, nygrid, nzgrid, mvar)")

    if len(xax) != dim.nxgrid:
        raise ValueError(f"_unshear: length of xax ({len(xax)}) must match nx ({dim.nxgrid})")

    # FFT wavenumbers along y
    ky = 2 * np.pi / Ly * np.concatenate([np.arange(dim.nygrid // 2 + 1), -np.arange(1, dim.nygrid // 2)[::-1]])

    arr_unsheared = np.empty_like(arr)

    for ix in range(dim.nxgrid):
        # Compute shift along y
        if nowrap:
            deltay_x = deltay * (xax[ix] - x0) / Lx
        else:
            deltay_x = (deltay % Ly) * (xax[ix] - x0) / Lx

        # Extract plane at this x (shape: ny, nz, mvar)
        plane = arr[ix, :, :, :]  # shape: (ny, nz, mvar)
        
        # FFT along y (axis=0 in row-major)
        plane_ky = fft(plane, axis=0)

        # Broadcast shift along all remaining axes
        shape = [dim.nygrid, 1, 1]  # broadcast along nz and mvar
        shift_array = np.exp(-1j * ky.reshape(shape) * deltay_x)

        # Apply shift
        plane_ky *= shift_array

        # Inverse FFT
        arr_unsheared[ix, :, :, :] = ifft(plane_ky, axis=0).real

    return arr_unsheared
#=======================================================================
