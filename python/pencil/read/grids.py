# grids.py
#
# Read the grid information from grid.dat.
#
# Authors:
# J. Oishi (joishi@amnh.org)
# S. Candelaresi (iomsn1@gmail.com)
#
# 27-jun-19: F. Gent added hdf5
"""
Contains classes and methods to read the grid data.
"""


def grid(*args, **kwargs):
    """
    Read the grid data from the pencil code simulation.
    If proc < 0, then load all data and assemble.
    Otherwise, load grid from specified processor.

    call signature:

    grid(datadir='data', proc=-1, quiet=False, trim=False)

    Keyword arguments:

    *datadir*:
      Directory where the data is stored.

    *proc*
      Processor to be read. If proc is -1, then read the 'global'
      grid. If proc is >=0, then read the grid.dat in the
      corresponding processor directory.

    *quiet*
      Flag for switching of output.

    *trim*
      Cuts off the ghost points.
    """

    grid_tmp = Grid()
    grid_tmp.read(*args, **kwargs)
    return grid_tmp


class Grid(object):
    """
    Grid -- holds pencil code time grid data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.x = self.y = self.z = 0
        self.dx_1 = self.dy_1 = self.dz_1 = 0
        self.dx_tilde = self.dy_tilde = self.dz_tilde = 0

        self.t = 0
        self.dx = self.dy = self.dz = 0
        self.Lx = self.Ly = self.Lz = 0


    def keys(self):
        for i in self.__dict__.keys():
            print(i)


    def read(self, datadir='data', proc=-1, quiet=False, precision='f',
             trim=False):
        """
        Read the grid data from the pencil code simulation.
        If proc < 0, then load all data and assemble.
        Otherwise, load grid from specified processor.

        call signature:

        grid(datadir='data', proc=-1, quiet=False, trim=False)

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *proc*
          Processor to be read. If proc is -1, then read the 'global'
          grid. If proc is >=0, then read the grid.dat in the
          corresponding processor directory.

        *quiet*
          Flag for switching of output.

        *precision*
          Float (f), double (d) or half (half).

        *trim*
          Cuts off the ghost points.
        """

        import numpy as np
        import os
        from scipy.io import FortranFile
        from pencil import read

        if precision == 'f':
            dtype = np.float32
        elif precision == 'd':
            dtype = np.float64
        elif precision == 'half':
            dtype = np.float16
        else:
            print('read grid: {} precision not set, using "f"'.format(
                  precision))
            dtype = np.float32

        if os.path.exists(os.path.join(datadir, 'grid.h5')):
            dim = read.dim(datadir, proc)
            import h5py

            with h5py.File(os.path.join(datadir, 'grid.h5'), 'r') as tmp:
                x = dtype(tmp['grid']['x'][()])
                y = dtype(tmp['grid']['y'][()])
                z = dtype(tmp['grid']['z'][()])
                dx_1 = dtype(tmp['grid']['dx_1'][()])
                dy_1 = dtype(tmp['grid']['dy_1'][()])
                dz_1 = dtype(tmp['grid']['dz_1'][()])
                dx_tilde = dtype(tmp['grid']['dx_tilde'][()])
                dy_tilde = dtype(tmp['grid']['dy_tilde'][()])
                dz_tilde = dtype(tmp['grid']['dz_tilde'][()])
                dx = dtype(tmp['grid']['dx'][()])
                dy = dtype(tmp['grid']['dy'][()])
                dz = dtype(tmp['grid']['dz'][()])
                Lx = dtype(tmp['grid']['Lx'][()])
                Ly = dtype(tmp['grid']['Ly'][()])
                Lz = dtype(tmp['grid']['Lz'][()])
                t = dtype(0.0)
        else:
            datadir = os.path.expanduser(datadir)
            dim = read.dim(datadir, proc)
            param = read.param(datadir=datadir, quiet=True,
                               conflicts_quiet=True)
            if dim.precision == 'D':
                read_precision = 'd'
            else:
                read_precision = 'f'

            if proc < 0:
                proc_dirs = list(filter(lambda string: string.startswith(
                    'proc'), os.listdir(datadir)))
                if (proc_dirs.count("proc_bounds.dat") > 0):
                    proc_dirs.remove("proc_bounds.dat")
                if param.lcollective_io:
                    # A collective IO strategy is being used
                    proc_dirs = ['allprocs']
            else:
                proc_dirs = ['proc' + str(proc)]

            # Define the global arrays.
            x = np.zeros(dim.mx, dtype=precision)
            y = np.zeros(dim.my, dtype=precision)
            z = np.zeros(dim.mz, dtype=precision)
            dx_1 = np.zeros(dim.mx, dtype=precision)
            dy_1 = np.zeros(dim.my, dtype=precision)
            dz_1 = np.zeros(dim.mz, dtype=precision)
            dx_tilde = np.zeros(dim.mx, dtype=precision)
            dy_tilde = np.zeros(dim.my, dtype=precision)
            dz_tilde = np.zeros(dim.mz, dtype=precision)

            for directory in proc_dirs:
                if not param.lcollective_io:
                    proc = int(directory[4:])
                    procdim = read.dim(datadir, proc)
                    if not quiet:
                        print("reading grid data from processor"+
                              " {0} of {1} ...".format(proc, len(proc_dirs)))
                else:
                    procdim = dim
                mxloc = procdim.mx
                myloc = procdim.my
                mzloc = procdim.mz

                # Read the grid data.
                file_name = os.path.join(datadir, directory, 'grid.dat')
                infile = FortranFile(file_name, 'r')
                grid_raw = infile.read_record(dtype=read_precision)
                dx, dy, dz = tuple(infile.read_record(dtype=read_precision))
                Lx, Ly, Lz = tuple(infile.read_record(dtype=read_precision))
                dx_1_raw = infile.read_record(dtype=read_precision)
                dx_tilde_raw = infile.read_record(dtype=read_precision)
                infile.close()

                # Reshape the arrays.
                t = dtype(grid_raw[0])
                x_loc = grid_raw[1:mxloc+1]
                y_loc = grid_raw[mxloc+1:mxloc+myloc+1]
                z_loc = grid_raw[mxloc+myloc+1:mxloc+myloc+mzloc+1]
                dx_1_loc = dx_1_raw[0:mxloc]
                dy_1_loc = dx_1_raw[mxloc:mxloc+myloc]
                dz_1_loc = dx_1_raw[mxloc+myloc:mxloc+myloc+mzloc]
                dx_tilde_loc = dx_tilde_raw[0:mxloc]
                dy_tilde_loc = dx_tilde_raw[mxloc:mxloc+myloc]
                dz_tilde_loc = dx_tilde_raw[mxloc+myloc:mxloc+myloc+mzloc]

                if len(proc_dirs) > 1:
                    if procdim.ipx == 0:
                        i0x = 0
                        i1x = i0x + procdim.mx
                        i0x_loc = 0
                        i1x_loc = procdim.mx
                    else:
                        i0x = procdim.ipx*procdim.nx + procdim.nghostx
                        i1x = i0x + procdim.mx - procdim.nghostx
                        i0x_loc = procdim.nghostx
                        i1x_loc = procdim.mx

                    if procdim.ipy == 0:
                        i0y = 0
                        i1y = i0y + procdim.my
                        i0y_loc = 0
                        i1y_loc = procdim.my
                    else:
                        i0y = procdim.ipy*procdim.ny + procdim.nghosty
                        i1y = i0y + procdim.my - procdim.nghosty
                        i0y_loc = procdim.nghosty
                        i1y_loc = procdim.my

                    if procdim.ipz == 0:
                        i0z = 0
                        i1z = i0z + procdim.mz
                        i0z_loc = 0
                        i1z_loc = procdim.mz
                    else:
                        i0z = procdim.ipz*procdim.nz + procdim.nghostz
                        i1z = i0z + procdim.mz - procdim.nghostz
                        i0z_loc = procdim.nghostz
                        i1z_loc = procdim.mz

                    x[i0x:i1x] = x_loc[i0x_loc:i1x_loc]
                    y[i0y:i1y] = y_loc[i0y_loc:i1y_loc]
                    z[i0z:i1z] = z_loc[i0z_loc:i1z_loc]
                    dx_1[i0x:i1x] = dx_1_loc[i0x_loc:i1x_loc]
                    dy_1[i0y:i1y] = dy_1_loc[i0y_loc:i1y_loc]
                    dz_1[i0z:i1z] = dz_1_loc[i0z_loc:i1z_loc]
                    dx_tilde[i0x:i1x] = dx_tilde_loc[i0x_loc:i1x_loc]
                    dy_tilde[i0y:i1y] = dy_tilde_loc[i0y_loc:i1y_loc]
                    dz_tilde[i0z:i1z] = dz_tilde_loc[i0z_loc:i1z_loc]

                else:
                    #x = dtype(x_loc.astype)
                    x = dtype(x_loc)
                    y = dtype(y_loc)
                    z = dtype(z_loc)
                    dx_1 = dtype(dx_1_loc)
                    dy_1 = dtype(dy_1_loc)
                    dz_1 = dtype(dz_1_loc)
                    dx_tilde = dtype(dx_tilde_loc)
                    dy_tilde = dtype(dy_tilde_loc)
                    dz_tilde = dtype(dz_tilde_loc)

        if trim:
            self.x = x[dim.l1:dim.l2+1]
            self.y = y[dim.m1:dim.m2+1]
            self.z = z[dim.n1:dim.n2+1]
            self.dx_1 = dx_1[dim.l1:dim.l2+1]
            self.dy_1 = dy_1[dim.m1:dim.m2+1]
            self.dz_1 = dz_1[dim.n1:dim.n2+1]
            self.dx_tilde = dx_tilde[dim.l1:dim.l2+1]
            self.dy_tilde = dy_tilde[dim.m1:dim.m2+1]
            self.dz_tilde = dz_tilde[dim.n1:dim.n2+1]
        else:
            self.x = x
            self.y = y
            self.z = z
            self.dx_1 = dx_1
            self.dy_1 = dy_1
            self.dz_1 = dz_1
            self.dx_tilde = dx_tilde
            self.dy_tilde = dy_tilde
            self.dz_tilde = dz_tilde

        self.t = t
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
