# tracers.py
#
# Perform a streamline tracer scan.
#
# Written by Simon Candelaresi (iomsn1@gmail.com)
"""
Reads the tracer files, composes a color map.
"""


class Tracers(object):
    """
    Tracers -- Holds the traced tracer object with the field line integrated
    quantities and the mapping.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.params = TracersParameterClass()
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.z1 = None
        self.l = None
        self.mapping = None
        self.t = None
        self.aa = None
        self.ee = None
        self.curly_A = None


    def find_tracers(self, var_file='VAR0', datadir='data', trace_field='bb',
                     ti=-1, tf=-1):
        """
        Trace streamlines of the vectofield 'field' from z = z0 to z = z1
        and integrate quantities 'int_q' along the lines. Creates a 2d
        mapping as in 'streamlines.f90'.

        call signature:

        find_tracers(var_file='VAR0', datadir='data', trace_field='bb',
                     ti=-1, tf=-1)

        Keyword arguments:

        *var_file*:
          Varfile to be read.

        *datadir*:
          Directory where the data is stored.

        *trace_field*:
          Vector field used for the streamline tracing.

        *ti*:
          Initial VAR file index for tracer time sequences. Overrides 'var_file'.

        *tf*:
          Final VAR file index for tracer time sequences. Overrides 'var_file'.
        """

        import numpy as np
        import multiprocessing as mp
        from pencil import read
        from pencil import math

        # Write the tracing parameters.
        self.params.trace_field = trace_field
        self.params.datadir = datadir

        # Multi core setup.
        if not(np.isscalar(self.params.n_proc)) or (self.params.n_proc%1 != 0):
            print("error: invalid processor number")
            return -1
        queue = mp.Queue()

        # Read the data.
        magic = []
        if trace_field == 'bb':
            magic.append('bb')
        if trace_field == 'jj':
            magic.append('jj')
        if trace_field == 'vort':
            magic.append('vort')
        if self.params.int_q == 'ee':
            magic.append('bb')
            magic.append('jj')
        dim = read.dim(datadir=datadir)
        self.params.var_file = var_file

        # Check if user wants a tracer time series.
        if (ti%1 == 0) and (tf%1 == 0) and (ti >= 0) and (tf >= ti):
            series = True
            nTimes = tf-ti+1
        else:
            series = False
            nTimes = 1

        # Initialize the arrays.
        self.x0 = np.zeros([int(self.params.trace_sub*dim.nx),
                            int(self.params.trace_sub*dim.ny), nTimes])
        self.y0 = np.zeros([int(self.params.trace_sub*dim.nx),
                            int(self.params.trace_sub*dim.ny), nTimes])
        self.x1 = np.zeros([int(self.params.trace_sub*dim.nx),
                            int(self.params.trace_sub*dim.ny), nTimes])
        self.y1 = np.zeros([int(self.params.trace_sub*dim.nx),
                            int(self.params.trace_sub*dim.ny), nTimes])
        self.z1 = np.zeros([int(self.params.trace_sub*dim.nx),
                            int(self.params.trace_sub*dim.ny), nTimes])
        self.l = np.zeros([int(self.params.trace_sub*dim.nx),
                           int(self.params.trace_sub*dim.ny), nTimes])
        if self.params.int_q == 'curly_A':
            self.curly_A = np.zeros([int(self.params.trace_sub*dim.nx),
                                     int(self.params.trace_sub*dim.ny), nTimes])
        if self.params.int_q == 'ee':
            self.ee = np.zeros([int(self.params.trace_sub*dim.nx),
                                int(self.params.trace_sub*dim.ny), nTimes])
        self.mapping = np.zeros([int(self.params.trace_sub*dim.nx),
                                 int(self.params.trace_sub*dim.ny),
                                 nTimes, 3])
        self.t = np.zeros(nTimes)

        for t_idx in range(ti, tf+1):
            if series:
                var_file = 'VAR' + str(t_idx)

            # Read the data.
            var = read.var(var_file=var_file, datadir=datadir, magic=magic,
                           quiet=True, trimall=True)
            grid = read.grid(datadir=datadir, quiet=True, trim=True)
            param2 = read.param(datadir=datadir, quiet=True)
            self.t[t_idx] = var.t

            # Extract the requested vector trace_field.
            field = getattr(var, trace_field)
            if self.params.int_q == 'curly_A':
                self.aa = var.aa
            if self.params.int_q == 'ee':
                self.ee = var.jj*param2.eta - math.cross(var.uu, var.bb)

            # Get the simulation parameters.
            self.params.dx = var.dx
            self.params.dy = var.dy
            self.params.dz = var.dz
            self.params.Ox = var.x[0]
            self.params.Oy = var.y[0]
            self.params.Oz = var.z[0]
            self.params.Lx = grid.Lx
            self.params.Ly = grid.Ly
            self.params.Lz = grid.Lz
            self.params.nx = dim.nx
            self.params.ny = dim.ny
            self.params.nz = dim.nz

            # Initialize the tracers.
            for ix in range(int(self.params.trace_sub*dim.nx)):
                for iy in range(int(self.params.trace_sub*dim.ny)):
                    self.x0[ix, iy, t_idx] = grid.x[0] + grid.dx/self.params.trace_sub*ix
                    self.x1[ix, iy, t_idx] = self.x0[ix, iy, t_idx].copy()
                    self.y0[ix, iy, t_idx] = grid.y[0] + grid.dy/self.params.trace_sub*iy
                    self.y1[ix, iy, t_idx] = self.y0[ix, iy, t_idx].copy()
                    self.z1[ix, iy, t_idx] = grid.z[0]

            proc = []
            sub_data = []
            for i_proc in range(self.params.n_proc):
                proc.append(mp.Process(target=self.__sub_tracers,
                                       args=(queue, field, t_idx, i_proc, self.params.n_proc)))
            for i_proc in range(self.params.n_proc):
                proc[i_proc].start()
            for i_proc in range(self.params.n_proc):
                sub_data.append(queue.get())
            for i_proc in range(self.params.n_proc):
                proc[i_proc].join()
            for i_proc in range(self.params.n_proc):
                # Extract the data from the single cores. Mind the order.
                sub_proc = sub_data[i_proc][0]
                self.x1[sub_proc::self.params.n_proc, :, t_idx] = sub_data[i_proc][1]
                self.y1[sub_proc::self.params.n_proc, :, t_idx] = sub_data[i_proc][2]
                self.z1[sub_proc::self.params.n_proc, :, t_idx] = sub_data[i_proc][3]
                self.l[sub_proc::self.params.n_proc, :, t_idx] = sub_data[i_proc][4]
                self.mapping[sub_proc::self.params.n_proc, :, t_idx, :] = sub_data[i_proc][5]
                if self.params.int_q == 'curly_A':
                    self.curly_A[sub_proc::self.params.n_proc, :, t_idx] = sub_data[i_proc][6]
                if self.params.int_q == 'ee':
                    self.ee[sub_proc::self.params.n_proc, :, t_idx] = sub_data[i_proc][7]
            for i_proc in range(self.params.n_proc):
                proc[i_proc].terminate()

            return 0


    # Return the tracers for the specified starting locations.
    def __sub_tracers(self, queue, field, t_idx, i_proc, n_proc):
        import numpy as np
        from pencil.calc.streamlines import Stream
        from pencil.math.interpolation import vec_int

        # Prepare the splines for the tricubis interpolation.
        if self.params.interpolation == 'tricubic':
            try:
                from eqtools.trispline import Spline
                x = np.linspace(self.params.Ox, self.params.Ox+self.params.Lx, self.params.nx)
                y = np.linspace(self.params.Oy, self.params.Oy+self.params.Ly, self.params.ny)
                z = np.linspace(self.params.Oz, self.params.Oz+self.params.Lz, self.params.nz)
                field_x = Spline(z, y, x, field[0, ...])
                field_y = Spline(z, y, x, field[1, ...])
                field_z = Spline(z, y, x, field[2, ...])
                splines = np.array([field_x, field_y, field_z])
            except:
                splines = None
        else:
            splines = None

        xx = np.zeros([(self.x0.shape[0]+n_proc-1-i_proc)//n_proc,
                       self.x0.shape[1], 3])
        xx[:, :, 0] = self.x0[i_proc:self.x0.shape[0]:n_proc, :, t_idx].copy()
        xx[:, :, 1] = self.y0[i_proc:self.x0.shape[0]:n_proc, :, t_idx].copy()
        xx[:, :, 2] = self.z1[i_proc:self.x0.shape[0]:n_proc, :, t_idx].copy()

        # Initialize the local arrays for this core.
        sub_x1 = np.zeros(xx[:, :, 0].shape)
        sub_y1 = np.zeros(xx[:, :, 0].shape)
        sub_z1 = np.zeros(xx[:, :, 0].shape)
        sub_l = np.zeros(xx[:, :, 0].shape)
        sub_curly_A = np.zeros(xx[:, :, 0].shape)
        sub_ee = np.zeros(xx[:, :, 0].shape)
        sub_mapping = np.zeros([xx[:, :, 0].shape[0], xx[:, :, 0].shape[1], 3])
        for ix in range(i_proc, self.x0.shape[0], n_proc):
            for iy in range(self.x0.shape[1]):
                time = np.linspace(0, 20*self.params.Lz/field[2, 0, iy, ix], 1000)
                stream = Stream(field, self.params, xx=xx[int(ix/n_proc), iy, :],
                                time=time, splines=splines)
                sub_x1[int(ix/n_proc), iy] = stream.tracers[-1, 0]
                sub_y1[int(ix/n_proc), iy] = stream.tracers[-1, 1]
                sub_z1[int(ix/n_proc), iy] = stream.tracers[-1, 2]
                sub_l[int(ix/n_proc), iy] = stream.total_l
                if self.params.int_q == 'curly_A':
                    for l in range(stream.total_l):
                        aaInt = vec_int((stream.tracers[l+1] + stream.tracers[l])/2, self.aa,
                                        [self.params.dx, self.params.dy, self.params.dz],
                                        [self.params.Ox, self.params.Oy, self.params.Oz],
                                        [self.params.nx, self.params.ny, self.params.nz],
                                        interpolation=self.params.interpolation)
                        sub_curly_A[int(ix/n_proc), iy] += \
                            np.dot(aaInt, (stream.tracers[l+1] - stream.tracers[l]))
                if self.params.int_q == 'ee':
                    for l in range(stream.total_l):
                        eeInt = vec_int((stream.tracers[l+1] + stream.tracers[l])/2, self.ee,
                                        [self.params.dx, self.params.dy, self.params.dz],
                                        [self.params.Ox, self.params.Oy, self.params.Oz],
                                        [self.params.nx, self.params.ny, self.params.nz],
                                        interpolation=self.params.interpolation)
                        sub_ee[int(ix/n_proc), iy] += \
                            np.dot(eeInt, (stream.tracers[l+1] - stream.tracers[l]))

                # Create the color mapping.
#                if (sub_z1[int(ix/n_proc), iy] > self.params.Oz+self.params.Lz-self.params.dz*12):
                if (self.x0[ix, iy, t_idx] - sub_x1[int(ix/n_proc), iy]) > 0:
                    if (self.y0[ix, iy, t_idx] - sub_y1[int(ix/n_proc), iy]) > 0:
                        sub_mapping[int(ix/n_proc), iy, :] = [0, 1, 0]
                    else:
                        sub_mapping[int(ix/n_proc), iy, :] = [1, 1, 0]
                else:
                    if (self.y0[ix, iy, t_idx] - sub_y1[int(ix/n_proc), iy]) > 0:
                        sub_mapping[int(ix/n_proc), iy, :] = [0, 0, 1]
                    else:
                        sub_mapping[int(ix/n_proc), iy, :] = [1, 0, 0]
#                else:
#                    sub_mapping[int(ix/n_proc), iy, :] = [1, 1, 1]

        queue.put((i_proc, sub_x1, sub_y1, sub_z1, sub_l, sub_mapping,
                   sub_curly_A, sub_ee))


    def write(self, datadir='data', destination='tracers.hdf5'):
        """
        Write the tracers into a file.

        call signature::

        write(self, datadir='data', destination='tracers.hdf5')

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *destination*:
          Destination file.
        """

        import os
        try:
            import h5py
        except:
            print("Warning: no h5py library found.")

        self.params.destination = destination

        # Write the results into hdf5 file.
        if destination != '':
            f = h5py.File(os.path.join(datadir, destination), 'w')
            # Write main data arrays.
            set_x0 = f.create_dataset("x0", self.x0.shape, dtype=self.x0.dtype)
            set_y0 = f.create_dataset("y0", self.y0.shape, dtype=self.y0.dtype)
            set_x1 = f.create_dataset("x1", self.x1.shape, dtype=self.x1.dtype)
            set_y1 = f.create_dataset("y1", self.y1.shape, dtype=self.y1.dtype)
            set_z1 = f.create_dataset("z1", self.z1.shape, dtype=self.z1.dtype)
            set_l = f.create_dataset("l", self.l.shape, dtype=self.l.dtype)
            set_x0[...] = self.x0[...]
            set_y0[...] = self.y0[...]
            set_x1[...] = self.x1[...]
            set_y1[...] = self.y1[...]
            set_z1[...] = self.z1[...]
            set_l[...] = self.l[...]
#            set_q = []
#            if not self.params.int_q == '':
#                set_q.append(f.create_dataset(self.params.int_q, getattr(self, self.params.int_q).shape,
#                                              dtype=getattr(self, self.params.int_q).dtype))
#                set_q[-1][...] = getattr(self, self.params.int_q)[...]
            set_t = f.create_dataset("t", self.t.shape, dtype=self.l.dtype)
            set_m = f.create_dataset("mapping", self.mapping.shape,
                                     dtype=self.mapping.dtype)
            set_t[...] = self.t[...]
            set_m[...] = self.mapping[...]
            # Write the parameters into their own group.
            group_params = f.create_group('params')
            for key in dir(self.params):
                if not key.startswith('_'):
                    if not key == 'int_q':
                        value = getattr(self.params, key)
                        group_params.attrs[key] = value
            f.close()
        else:
            print("error: empty destination file")


    def read(self, datadir='data', file_name='tracers.hdf5'):
        """
        Read the tracers from a file.

        call signature:

        read(self, datadir='data', file_name='tracers.hdf5')

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *file_name*:
          File with the tracer data.
        """

        import os
        try:
            import h5py
        except:
            print("Warning: no h5py library found.")

        # Open the file.
        f = h5py.File(os.path.join(datadir, file_name), 'r')

        # Extract arrays.
        self.t = f['t'].value
        self.x0 = f['x0'].value
        self.y0 = f['y0'].value
        self.x1 = f['x1'].value
        self.y1 = f['y1'].value
        self.z1 = f['z1'].value
        self.l = f['l'].value
        self.mapping = f['mapping'].value
        print(f.keys())
        if 'curly_A' in list(f.keys()):
            self.curly_A = f['curly_A'].value
        if 'ee' in list(f.keys()):
            self.ee = f['ee'].value

        # Extract parameters.
        params = f['params']
        self.params = TracersParameterClass()
        for param in params.attrs.keys():
            setattr(self.params, param, params.attrs[param])

        f.close()


# Class containing simulation and tracing parameters.
class TracersParameterClass(object):
    """
    __TracersParameterClass -- Holds the simulation and tracing parameters.
    """
    def __init__(self):
        """
        Initialize the parameters.
        """
        self.dx = 0
        self.dy = 0
        self.dz = 0
        self.Ox = 0
        self.Oy = 0
        self.Oz = 0
        self.Lx = 0
        self.Ly = 0
        self.Lz = 0
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.trace_field = ''
        self.rtol = 1e-8
        self.atol = 1e-8
        self.periodic_x = False
        self.periodic_y = False
        self.periodic_z = False
        self.interpolation = 'trilinear'
        self.method = 'RK45'
        self.trace_sub = 1
        self.int_q = ''
        self.var_file = 'VAR0'
        self.ti = -1
        self.tf = -1
        self.datadir = 'data'
        self.destination = 'tracers.hdf5'
        self.n_proc = 1
