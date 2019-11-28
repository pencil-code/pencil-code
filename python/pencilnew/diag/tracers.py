# tracers.py
# Written by Simon Candelaresi (iomsn1@gmail.com)
"""
Reads the tracer files, composes a color map.
"""

#import numpy as np
import os
try:
    import h5py
except:
    print("Warning: no h5py library found.")
import multiprocessing as mp
from ..calc.streamlines import Stream
from ..math.interpolation import vec_int


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
        self.x0 = []
        self.y0 = []
        self.x1 = []
        self.y1 = []
        self.z1 = []
        self.l = []
        self.mapping = []
        self.t = []


    def find_tracers(self, trace_field='bb', h_min=2e-3, h_max=2e4, len_max=500,
                     tol=1e-2, iter_max=1e3, interpolation='trilinear',
                     trace_sub=1, int_q=[''], varfile='VAR0', ti=-1, tf=-1,
                     integration='simple', datadir='./data', n_proc=1):
        """
        Trace streamlines from the VAR files and integrate quantity 'int_q'
        along them.

        call signature:

        find_tracers(self, trace_field='bb', h_min=2e-3, h_max=2e4, len_max=500,
                     tol=1e-2, iter_max=1e3, interpolation='trilinear',
                     trace_sub=1, int_q=[''], varfile='VAR0', ti=-1, tf=-1,
                     integration='simple', datadir='data/', n_proc=1)

        Trace streamlines of the vectofield 'field' from z = z0 to z = z1
        and integrate quantities 'int_q' along the lines. Creates a 2d
        mapping as in 'streamlines.f90'.

        Keyword arguments:

        *trace_field*:
          Vector field used for the streamline tracing.

        *h_min*:
          Minimum step length for and underflow to occur.

        *h_max*:
          Parameter for the initial step length.

        *len_max*:
          Maximum length of the streamline. Integration will stop if
          len >= len_max.

        *tol*:
          Tolerance for each integration step. Reduces the step length if
          error >= tol.

        *iter_max*:
          Maximum number of iterations.

        *interpolation*:
          Interpolation of the vector field.
          'mean': takes the mean of the adjacent grid point.
          'trilinear': weights the adjacent grid points according to
          their distance.

        *trace_sub*:
          Number of sub-grid cells for the seeds.

        *int_q*:
          Quantities to be integrated along the streamlines.

        *varfile*:
          Varfile to be read.

        *integration*:
          Integration method.
          'simple': low order method.
          'RK6': Runge-Kutta 6th order.

        *ti*:
          Initial VAR file index for tracer time sequences. Overrides
          'varfile'.

        *tf*:
          Final VAR file index for tracer time sequences. Overrides 'varfile'.

        *datadir*:
          Directory where the data is stored.

        *n_proc*:
          Number of cores for multi core computation.
        """

        import numpy as np
        from .. import read
        from .. import math

        # Return the tracers for the specified starting locations.
        def __sub_tracers(queue, var, field, t_idx, i_proc, n_proc):

            xx = np.zeros([(self.x0.shape[0]+n_proc-1-i_proc)/n_proc,
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
                    stream = Stream(field, self.params, interpolation=interpolation,
                                    h_min=h_min, h_max=h_max, len_max=len_max, tol=tol,
                                    iter_max=iter_max, xx=xx[int(ix/n_proc), iy, :])
                    sub_x1[int(ix/n_proc), iy] = stream.tracers[stream.stream_len, 0]
                    sub_y1[int(ix/n_proc), iy] = stream.tracers[stream.stream_len, 1]
                    sub_z1[int(ix/n_proc), iy] = stream.tracers[stream.stream_len, 2]
                    sub_l[int(ix/n_proc), iy] = stream.len
                    if any(np.array(self.params.int_q) == 'curly_A'):
                        for l in range(stream.stream_len):
                            aaInt = vec_int((stream.tracers[l+1] + stream.tracers[l])/2,
                                             var, aa, interpolation=self.params.interpolation)
                            sub_curly_A[int(ix/n_proc), iy] += \
                                np.dot(aaInt, (stream.tracers[l+1] - stream.tracers[l]))
                    if any(np.array(self.params.int_q) == 'ee'):
                        for l in range(stream.stream_len):
                            eeInt = vec_int((stream.tracers[l+1] + stream.tracers[l])/2,
                                             var, ee, interpolation=self.params.interpolation)
                            sub_ee[int(ix/n_proc), iy] += \
                                np.dot(eeInt, (stream.tracers[l+1] - stream.tracers[l]))

                    # Create the color mapping.
                    if (sub_z1[int(ix/n_proc), iy] > self.params.Oz+self.params.Lz-self.params.dz*4):
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
                    else:
                        sub_mapping[int(ix/n_proc), iy, :] = [1, 1, 1]

            queue.put((i_proc, sub_x1, sub_y1, sub_z1, sub_l, sub_mapping,
                       sub_curly_A, sub_ee))


        # Write the tracing parameters.
        self.params.trace_field = trace_field
        self.params.h_min = h_min
        self.params.h_max = h_max
        self.params.len_max = len_max
        self.params.tol = tol
        self.params.interpolation = interpolation
        self.params.trace_sub = trace_sub
        self.params.int_q = int_q
        self.params.varfile = varfile
        self.params.ti = ti
        self.params.tf = tf
        self.params.integration = integration
        self.params.datadir = datadir
        self.params.n_proc = n_proc

        # Multi core setup.
        if not(np.isscalar(n_proc)) or (n_proc%1 != 0):
            print("error: invalid processor number")
            return -1
        queue = mp.Queue()

        # Convert int_q string into list.
        if not isinstance(int_q, list):
            int_q = [int_q]

        # Read the data.
        magic = []
        if trace_field == 'bb':
            magic.append('bb')
        if trace_field == 'jj':
            magic.append('jj')
        if trace_field == 'vort':
            magic.append('vort')
        if any(np.array(int_q) == 'ee'):
            magic.append('bb')
            magic.append('jj')
        dim = read.dim(datadir=datadir)

        # Check if user wants a tracer time series.
        if (ti%1 == 0) and (tf%1 == 0) and (ti >= 0) and (tf >= ti):
            series = True
            nTimes = tf-ti+1
        else:
            series = False
            nTimes = 1

        # Initialize the arrays.
        self.x0 = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny), nTimes])
        self.y0 = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny), nTimes])
        self.x1 = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny), nTimes])
        self.y1 = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny), nTimes])
        self.z1 = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny), nTimes])
        self.l = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny), nTimes])
        if any(np.array(int_q) == 'curly_A'):
            self.curly_A = np.zeros([int(trace_sub*dim.nx),
                                     int(trace_sub*dim.ny), nTimes])
        if any(np.array(int_q) == 'ee'):
            self.ee = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny),
                                nTimes])
        self.mapping = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny),
                                 nTimes, 3])
        self.t = np.zeros(nTimes)

        for t_idx in range(ti, tf+1):
            if series:
                varfile = 'VAR' + str(t_idx)

            # Read the data.
            var = read.var(var_file=varfile, datadir=datadir, magic=magic,
                           quiet=True, trimall=True)
            grid = read.grid(datadir=datadir, quiet=True, trim=True)
            param2 = read.param(datadir=datadir, param2=True, quiet=True)
            self.t[t_idx] = var.t

            # Extract the requested vector trace_field.
            field = getattr(var, trace_field)
            if any(np.array(int_q) == 'curly_A'):
                aa = var.aa
            if any(np.array(int_q) == 'ee'):
                ee = var.jj*param2.eta - math.cross(var.uu, var.bb)

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
            for ix in range(int(trace_sub*dim.nx)):
                for iy in range(int(trace_sub*dim.ny)):
                    self.x0[ix, iy, t_idx] = grid.x[0] + grid.dx/trace_sub*ix
                    self.x1[ix, iy, t_idx] = self.x0[ix, iy, t_idx].copy()
                    self.y0[ix, iy, t_idx] = grid.y[0] + grid.dy/trace_sub*iy
                    self.y1[ix, iy, t_idx] = self.y0[ix, iy, t_idx].copy()
                    self.z1[ix, iy, t_idx] = grid.z[0]

            proc = []
            sub_data = []
            for i_proc in range(n_proc):
                proc.append(mp.Process(target=__sub_tracers, args=(queue, var, field, t_idx, i_proc, n_proc)))
            for i_proc in range(n_proc):
                proc[i_proc].start()
            for i_proc in range(n_proc):
                sub_data.append(queue.get())
            for i_proc in range(n_proc):
                proc[i_proc].join()
            for i_proc in range(n_proc):
                # Extract the data from the single cores. Mind the order.
                sub_proc = sub_data[i_proc][0]
                self.x1[sub_proc::n_proc, :, t_idx] = sub_data[i_proc][1]
                self.y1[sub_proc::n_proc, :, t_idx] = sub_data[i_proc][2]
                self.z1[sub_proc::n_proc, :, t_idx] = sub_data[i_proc][3]
                self.l[sub_proc::n_proc, :, t_idx] = sub_data[i_proc][4]
                self.mapping[sub_proc::n_proc, :, t_idx, :] = sub_data[i_proc][5]
                if any(np.array(int_q) == 'curly_A'):
                    self.curly_A[sub_proc::n_proc, :, t_idx] = sub_data[i_proc][6]
                if any(np.array(int_q) == 'ee'):
                    self.ee[sub_proc::n_proc, :, t_idx] = sub_data[i_proc][7]
            for i_proc in range(n_proc):
                proc[i_proc].terminate()


    def write(self, datadir='./data', destination='tracers.hdf5'):
        """
        Write the tracers into a file.

        call signature::

        write(self, datadir='./data', destination='tracers.hdf5')

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *destination*:
          Destination file.
        """

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
            set_q = []
            for q in self.params.int_q:
                if not q == '':
                    set_q.append(f.create_dataset(q, getattr(self, q).shape,
                                                  dtype=getattr(self, q).dtype))
                    set_q[-1][...] = getattr(self, q)[...]
            set_t = f.create_dataset("t", self.t.shape, dtype=self.l.dtype)
            set_m = f.create_dataset("mapping", self.mapping.shape,
                                     dtype=self.mapping.dtype)
            set_t[...] = self.t[...]
            set_m[...] = self.mapping[...]
            # Write the parameters into their own group.
            group_params = f.create_group('params')
            for key in dir(self.params):
                if not key.startswith('_'):
                    value = getattr(self.params, key)
                    group_params.attrs[key] = value
            f.close()
        else:
            print("error: empty destination file")


    def read(self, datadir='./data', file_name='tracers.hdf5'):
        """
        Read the tracers from a file.

        call signature::

        read(self, datadir='./data', file_name='tracers.hdf5')

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *file_name*:
          File with the tracer data.
        """

        import numpy as np

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
        if any(np.array(f.keys()) == 'curly_A'):
            self.curly_A = f['curly_A'].value
        if any(np.array(f.keys()) == 'ee'):
            self.ee = f['ee'].value

        # Extract parameters.
        params = f['params']
        self.params = TracersParameterClass()
        for param in params.attrs.keys():
            setattr(self.params, param, params.attrs[param])

        f.close()


## Need some updating of this routine.
#def tracer_movie(datadir='data/', tracer_file='tracers.hf5',
#                 fixedFile='fixed_points.dat', zlim=False,
#                 head_size=3, hm=1,
#                 imageDir='./', movieFile='fixed_points.mpg',
#                 fps=5.0, bitrate=1800):
#    """
#    Plots the color mapping together with the fixed points.
#    Creates a movie file.
#
#    call signature::
#
#      tracer_movie(datadir='data/', tracerFile='tracers.dat',
#                   tracer_file='fixed_points.dat', zlim = [],
#                   head_size=3, hm=1,
#                   imageDir='./', movieFile='fixed_points.mpg',
#                   fps=5.0, bitrate=1800)
#
#    Plots the field line mapping together with the fixed points and creates
#    a movie file.
#
#      *datadir*:
#        Data directory.
#
#      *tracerFile*:
#        Name of the tracer file.
#
#      *fixedFile*:
#        Name of the fixed points file.
#
#      *zlim*:
#        The upper limit for the field line mapping at which a field line is
#        considered to have reached the upper boundary.
#
#      *head_size*:
#        Size of the fortran header in binary data. Most of the time this is 3.
#        For the St Andrews cluster it is 5.
#
#      *hm*:
#        Header multiplication factor in case Fortran's binary data writes
#        extra large header. For most cases hm = 1 is sufficient.
#        For the cluster in St Andrews use hm = 2.
#
#      *imageDir*:
#        Directory with the temporary png files.
#
#      *movieFile*:
#        Output file for the movie. Ending should be 'mpg', since the compression
#        format is mpg.
#
#      *fps*:
#        Frames per second of the animation.
#
#      *bitrate*:
#        Bitrate of the movie file. Set to higher value for higher quality.
#    """
#
#
#    # Read the mapping and the fixed point positions.
#    tracers, mapping, t = pc.read_tracers(datadir=datadir, fileName=tracer_file,
#                                          zlim=zlim, head_size=head_size)
#    fixed = pc.read_fixed_points(datadir=datadir, fileName=fixedFile, hm=hm)
#
#    # Read the parameters for the domain boundaries.
#    params = pc.read_param(quiet=True)
#    domain = [params.xyz0[0], params.xyz1[0], params.xyz0[1], params.xyz1[1]]
#
#    # Determine how much faster the fixed pints have been written out than
#    # the color mapping.
#    advance = np.ceil(float(len(fixed.t))/len(mapping[0, 0, :, 0]))
#
#    # Determine the colors for the fixed points.
#    colors = np.zeros(np.shape(fixed.q) + (3,))
#    colors[:, :, :] = 0.
#    print np.shape(colors)
#    for j in range(len(colors[:, 0, 0])):
#        for k in range(len(colors[0, :, 0])):
#            if fixed.q[j, k] >= 0:
#                colors[j, k, 1] = colors[j, k, 2] = (1-fixed.q[j, k]/\
#                                  np.max(np.abs(fixed.q[:, k])))
#                colors[j, k, 0] = fixed.q[j, k]/np.max(np.abs(fixed.q[:, k]))
#            else:
#                colors[j, k, 0] = colors[j, k, 1] = (1+fixed.q[j, k]/\
#                                  np.max(np.abs(fixed.q[:, k])))
#                colors[j, k, 2] = -fixed.q[j, k]/np.max(np.abs(fixed.q[:, k]))
#
#    # Prepare the plot.
#    width = 6
#    height = 6
#    plt.rc("figure.subplot", left=(60/72.27)/width)
#    plt.rc("figure.subplot", right=(width-20/72.27)/width)
#    plt.rc("figure.subplot", bottom=(50/72.27)/height)
#    plt.rc("figure.subplot", top=(height-20/72.27)/height)
#    figure = plt.figure(figsize=(width, height))
#
#    for k in range(len(fixed.x[0, :])):
#        dots = plt.plot(fixed.x[0, k], fixed.y[0, k], 'o', c=colors[0, k, :])
#    image = plt.imshow(zip(*mapping[:, ::-1, 0, :]), interpolation='nearest',
#                       extent=domain)
#    j = 0
#    frameName = imageDir + 'images%06d.png'%j
#    imageFiles = []
#    imageFiles.append(frameName)
#    figure.savefig(frameName)
#
#    for j in range(1, len(fixed.t)):
#        #time.sleep(0.5)
#        figure.clear()
#        for k in range(len(fixed.x[j, :])):
#            dots = plt.plot(fixed.x[j, k], fixed.y[j, k], 'o',
#                            c=colors[j, k, :])
#        image = plt.imshow(zip(*mapping[:, ::-1, np.floor(j/advance), :]),
#                           interpolation='nearest', extent=domain)
#        frameName = imageDir + 'images%06d.png'%j
#        imageFiles.append(frameName)
#        figure.savefig(frameName)
#
#    # Convert the images into a mpg file.
#    mencodeCommand = "mencoder 'mf://"+imageDir+"images*.png'\
#        -mf type=png:fps="+np.str(fps)+" -ovc lavc \
#        -lavcopts vcodec=mpeg4:vhq:vbitrate="+np.str(bitrate)+" -ffourcc MP4S \
#        -oac copy -o "+movieFile
#    os.system(mencodeCommand)
#
#    # Remove the image files.
#    for fname in imageFiles:
#        os.remove(fname)


# Class containing simulation and tracing parameters.
class TracersParameterClass(object):
    """
    __TracersParameterClass -- Holds the simulation and tracing parameters.
    """
    def __init__(self):
        """
        Initialize the parameters.
        """
        self.dx = 0; self.dy = 0; self.dz = 0
        self.Ox = 0; self.Oy = 0; self.Oz = 0
        self.Lx = 0; self.Ly = 0; self.Lz = 0
        self.nx = 0; self.ny = 0; self.nz = 0
        self.trace_field = ''
        self.h_min = 2e-3
        self.h_max = 2e4
        self.len_max = 500
        self.tol = 1e-2
        self.iter_max = 1e3
        self.interpolation = 'trilinear'
        self.trace_sub = 1
        self.int_q = ''
        self.varfile = 'VAR0'
        self.ti = -1
        self.tf = -1
        self.integration = 'simple'
        self.datadir = 'data/'
        self.destination = 'tracers.hdf5'
        self.n_proc = 1
