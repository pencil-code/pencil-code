# fixed_points.py
# Written by Simon Candelaresi (iomsn1@gmail.com)
"""
Creates the fixed point values.
"""

#import numpy as np
import math as m
import multiprocessing as mp
import os
try:
    import h5py
except:
    print('!! ERR in diag/fixed_points.py: Dependency of h5py not fullfilled.')
from ..diag.tracers import TracersParameterClass
from ..diag.tracers import Tracers
from ..calc.streamlines import Stream
from ..math.interpolation import vec_int


class FixedPoint(object):
    """
    FixedPoint -- Holds the fixed points and additional integrated quantities.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.params = TracersParameterClass()
        self.t = []
        self.fidx = []
        self.fixed_points = []
        self.fixed_sign = []
        self.poincare = []
        self.tracers = []


    def find_fixed(self, datadir='data/', destination='fixed_points.hf5',
                   varfile='VAR0', ti=-1, tf=-1, trace_field='bb', h_min=2e-3,
                   h_max=2e4, len_max=500, tol=1e-2, interpolation='trilinear',
                   trace_sub=1, integration='simple', int_q=[''], n_proc=1,
                   tracer_file_name=''):
        """
        Find the fixed points.

        call signature::

        find_fixed(datadir='data/', destination='fixed_points.hf5',
                   varfile='VAR0', ti=-1, tf=-1, trace_field='bb', h_min=2e-3,
                   h_max=2e4, len_max=500, tol=1e-2, interpolation='trilinear',
                   trace_sub=1, integration='simple', int_q=[''], n_proc=1):

        Finds the fixed points. Returns the fixed points positions.

        Keyword arguments:

          *datadir*:
            Data directory.

          *destination*:
            Name of the fixed points file.

         *varfile*:
           Varfile to be read.

          *ti*:
            Initial VAR file index for tracer time sequences.

          *tf*:
            Final VAR file index for tracer time sequences.

         *trace_field*:
           Vector field used for the streamline tracing.

         *h_min*:
           Minimum step length for and underflow to occur.

         *h_max*:
           Parameter for the initial step length.

         *len_max*:
           Maximum length of the streamline. Integration will stop if
           l >= len_max.

         *tol*:
           Tolerance for each integration step.
           Reduces the step length if error >= tol.

         *interpolation*:
           Interpolation of the vector field.
           'mean': takes the mean of the adjacent grid point.
           'trilinear': weights the adjacent grid points according to
                        their distance.

         *trace_sub*:
           Number of sub-grid cells for the seeds for the initial mapping.

          *integration*:
            Integration method.
            'simple': low order method.
            'RK6': Runge-Kutta 6th order.

         *int_q*:
           Quantities to be integrated along the streamlines.

         *n_proc*:
           Number of cores for multi core computation.

         *tracer_file_name*
           Name of the tracer file to be read.
           If equal to '' it will compute the tracers.
        """


        import numpy as np

        # Return the fixed points for a subset of the domain.
        def __sub_fixed(queue, ix0, iy0, field, tracers, tidx, var, i_proc):
            diff = np.zeros((4, 2))
            fixed = []
            fixed_sign = []
            fidx = 0
            poincare_array = np.zeros((tracers.x0[i_proc::self.params.n_proc].shape[0],
                                       tracers.x0.shape[1]))

            for ix in ix0[i_proc::self.params.n_proc]:
                for iy in iy0:
                    # Compute Poincare index around this cell (!= 0 for potential fixed point).
                    diff[0, :] = np.array([tracers.x1[ix, iy, tidx] - tracers.x0[ix, iy, tidx],
                                           tracers.y1[ix, iy, tidx] - tracers.y0[ix, iy, tidx]])
                    diff[1, :] = np.array([tracers.x1[ix+1, iy, tidx] - tracers.x0[ix+1, iy, tidx],
                                           tracers.y1[ix+1, iy, tidx] - tracers.y0[ix+1, iy, tidx]])
                    diff[2, :] = np.array([tracers.x1[ix+1, iy+1, tidx] - tracers.x0[ix+1, iy+1, tidx],
                                           tracers.y1[ix+1, iy+1, tidx] - tracers.y0[ix+1, iy+1, tidx]])
                    diff[3, :] = np.array([tracers.x1[ix, iy+1, tidx] - tracers.x0[ix, iy+1, tidx],
                                           tracers.y1[ix, iy+1, tidx] - tracers.y0[ix, iy+1, tidx]])
                    if sum(np.sum(diff**2, axis=1) != 0):
                        diff = np.swapaxes(np.swapaxes(diff, 0, 1) /
                               np.sqrt(np.sum(diff**2, axis=1)), 0, 1)
                    poincare = __poincare_index(field, tracers.x0[ix:ix+2, iy, tidx],
                                                tracers.y0[ix, iy:iy+2, tidx], diff)
                    poincare_array[ix/n_proc, iy] = poincare

                    if abs(poincare) > 5: # Use 5 instead of 2*pi to account for rounding errors.
                        # Subsample to get starting point for iteration.
                        nt = 4
                        xmin = tracers.x0[ix, iy, tidx]
                        ymin = tracers.y0[ix, iy, tidx]
                        xmax = tracers.x0[ix+1, iy, tidx]
                        ymax = tracers.y0[ix, iy+1, tidx]
                        xx = np.zeros((nt**2, 3))
                        tracers_part = np.zeros((nt**2, 5))
                        i1 = 0
                        for j1 in range(nt):
                            for k1 in range(nt):
                                xx[i1, 0] = xmin + j1/(nt-1.)*(xmax-xmin)
                                xx[i1, 1] = ymin + k1/(nt-1.)*(ymax-ymin)
                                xx[i1, 2] = self.params.Oz
                                i1 += 1
                        for it1 in range(nt**2):
                            stream = Stream(field, self.params,
                                            h_min=self.params.h_min,
                                            h_max=self.params.h_max,
                                            len_max=self.params.len_max,
                                            tol=self.params.tol,
                                            interpolation=self.params.interpolation,
                                            integration=self.params.integration,
                                            xx=xx[it1, :])
                            tracers_part[it1, 0:2] = xx[it1, 0:2]
                            tracers_part[it1, 2:] = stream.tracers[stream.stream_len-1, :]
                        min2 = 1e6
                        minx = xmin
                        miny = ymin
                        i1 = 0
                        for j1 in range(nt):
                            for k1 in range(nt):
                                diff2 = (tracers_part[i1+k1*nt, 2] - \
                                         tracers_part[i1+k1*nt, 0])**2 + \
                                        (tracers_part[i1+k1*nt, 3] - \
                                         tracers_part[i1+k1*nt, 1])**2
                                if diff2 < min2:
                                    min2 = diff2
                                    minx = xmin + j1/(nt-1.)*(xmax - xmin)
                                    miny = ymin + k1/(nt-1.)*(ymax - ymin)
                                it1 += 1

                        # Get fixed point from this starting position using Newton's method.
                        point = np.array([minx, miny])
                        fixed_point = __null_point(point, var)

                        # Check if fixed point lies inside the cell.
                        if ((fixed_point[0] < tracers.x0[ix, iy, tidx]) or
                            (fixed_point[0] > tracers.x0[ix+1, iy, tidx]) or
                            (fixed_point[1] < tracers.y0[ix, iy, tidx]) or
                            (fixed_point[1] > tracers.y0[ix, iy+1, tidx])):
                                pass
                        else:
                            fixed.append(fixed_point)
                            fixed_sign.append(np.sign(poincare))
                            fidx += np.sign(poincare)

            queue.put((i_proc, fixed, fixed_sign, fidx, poincare_array))


        # Find the Poincare index of this grid cell.
        def __poincare_index(field, sx, sy, diff):
            poincare = 0
            poincare += __edge(field, [sx[0], sx[1]], [sy[0], sy[0]],
                               diff[0, :], diff[1, :], 0)
            poincare += __edge(field, [sx[1], sx[1]], [sy[0], sy[1]],
                               diff[1, :], diff[2, :], 0)
            poincare += __edge(field, [sx[1], sx[0]], [sy[1], sy[1]],
                               diff[2, :], diff[3, :], 0)
            poincare += __edge(field, [sx[0], sx[0]], [sy[1], sy[0]],
                               diff[3, :], diff[0, :], 0)
            return poincare


        # Compute rotation along one edge.
        def __edge(field, sx, sy, diff1, diff2, rec):
            phiMin = np.pi/8.
            dtot = m.atan2(diff1[0]*diff2[1] - diff2[0]*diff1[1],
                           diff1[0]*diff2[0] + diff1[1]*diff2[1])
            if (abs(dtot) > phiMin) and (rec < 4):
                xm = 0.5*(sx[0]+sx[1])
                ym = 0.5*(sy[0]+sy[1])

                # Trace the intermediate field line.
                stream = Stream(field, self.params, h_min=self.params.h_min,
                                h_max=self.params.h_max,
                                len_max=self.params.len_max,
                                tol=self.params.tol,
                                interpolation=self.params.interpolation,
                                integration=self.params.integration,
                                xx=np.array([xm, ym, self.params.Oz]))
                stream_x0 = stream.tracers[0, 0]
                stream_y0 = stream.tracers[0, 1]
                stream_x1 = stream.tracers[stream.stream_len-1, 0]
                stream_y1 = stream.tracers[stream.stream_len-1, 1]
                stream_z1 = stream.tracers[stream.stream_len-1, 2]

                # Discard any streamline which does not converge or hits the boundary.
#                if ((stream.len >= len_max) or
#                (stream_z1 < self.params.Oz+self.params.Lz-10*self.params.dz)):
#                    dtot = 0.
                if False:
                    pass
                else:
                    diffm = np.array([stream_x1 - stream_x0, stream_y1 - stream_y0])
                    if sum(diffm**2) != 0:
                        diffm = diffm/np.sqrt(sum(diffm**2))
                    dtot = __edge(field, [sx[0], xm], [sy[0], ym], diff1, diffm, rec+1) + \
                           __edge(field, [xm, sx[1]], [ym, sy[1]], diffm, diff2, rec+1)
            return dtot


        # Finds the null point of the mapping, i.e. fixed point, using Newton's method.
        def __null_point(point, var):
            dl = np.min(var.dx, var.dy)/100.
            it = 0
            # Tracers used to find the fixed point.
            tracers_null = np.zeros((5, 4))
            while True:
                # Trace field lines at original point and for Jacobian.
                # (second order seems to be enough)
                xx = np.zeros((5, 3))
                xx[0, :] = np.array([point[0], point[1], self.params.Oz])
                xx[1, :] = np.array([point[0]-dl, point[1], self.params.Oz])
                xx[2, :] = np.array([point[0]+dl, point[1], self.params.Oz])
                xx[3, :] = np.array([point[0], point[1]-dl, self.params.Oz])
                xx[4, :] = np.array([point[0], point[1]+dl, self.params.Oz])
                for it1 in range(5):
                    stream = Stream(field, self.params, h_min=self.params.h_min,
                                    h_max=self.params.h_max, len_max=self.params.len_max,
                                    tol=self.params.tol, interpolation=self.params.interpolation,
                                    integration=self.params.integration, xx=xx[it1, :])
                    tracers_null[it1, :2] = xx[it1, :2]
                    tracers_null[it1, 2:] = stream.tracers[stream.stream_len-1, 0:2]

                # Check function convergence.
                ff = np.zeros(2)
                ff[0] = tracers_null[0, 2] - tracers_null[0, 0]
                ff[1] = tracers_null[0, 3] - tracers_null[0, 1]
                if sum(abs(ff)) <= 1e-3*np.min(self.params.dx, self.params.dy):
                    fixed_point = np.array([point[0], point[1]])
                    break

                # Compute the Jacobian.
                fjac = np.zeros((2, 2))
                fjac[0, 0] = ((tracers_null[2, 2] - tracers_null[2, 0]) -
                              (tracers_null[1, 2] - tracers_null[1, 0]))/2./dl
                fjac[0, 1] = ((tracers_null[4, 2] - tracers_null[4, 0]) -
                              (tracers_null[3, 2] - tracers_null[3, 0]))/2./dl
                fjac[1, 0] = ((tracers_null[2, 3] - tracers_null[2, 1]) -
                              (tracers_null[1, 3] - tracers_null[1, 1]))/2./dl
                fjac[1, 1] = ((tracers_null[4, 3] - tracers_null[4, 1]) -
                              (tracers_null[3, 3] - tracers_null[3, 1]))/2./dl

                # Invert the Jacobian.
                fjin = np.zeros((2, 2))
                det = fjac[0, 0]*fjac[1, 1] - fjac[0, 1]*fjac[1, 0]
                if abs(det) < dl:
                    fixed_point = point
                    break
                fjin[0, 0] = fjac[1, 1]
                fjin[1, 1] = fjac[0, 0]
                fjin[0, 1] = -fjac[0, 1]
                fjin[1, 0] = -fjac[1, 0]
                fjin = fjin/det
                dpoint = np.zeros(2)
                dpoint[0] = -fjin[0, 0]*ff[0] - fjin[0, 1]*ff[1]
                dpoint[1] = -fjin[1, 0]*ff[0] - fjin[1, 1]*ff[1]
                point += dpoint

                # Check root convergence.
                if sum(abs(dpoint)) < 1e-3*np.min(self.params.dx, self.params.dy):
                    fixed_point = point
                    break

                if it > 20:
                    fixed_point = point
                    break

                it += 1

            return fixed_point


        # Find the fixed point using Newton's method, starting at previous fixed point.
        def __sub_fixed_series(queue, t_idx, field, var, i_proc):
            fixed = []
            fixed_sign = []
            for i, point in enumerate(self.fixed_points[t_idx-1][i_proc::self.params.n_proc]):
                fixed_tentative = __null_point(point, var)
                # Check if the fixed point lies outside the domain.
                if fixed_tentative[0] >= self.params.Ox and \
                fixed_tentative[1] >= self.params.Oy and \
                fixed_tentative[0] <= self.params.Ox+self.params.Lx and \
                fixed_tentative[1] <= self.params.Oy+self.params.Ly:
                    fixed.append(fixed_tentative)
                    fixed_sign.append(self.fixed_sign[t_idx-1][i_proc+i*n_proc])
            queue.put((i_proc, fixed, fixed_sign))


        # Discard fixed points which are too close to each other.
        def __discard_close_fixed_points(fixed, fixed_sign, var):

            fixed_new = []
            fixed_sign_new = []
            if len(fixed) > 0:
                fixed_new.append(fixed[0])
                fixed_sign_new.append(fixed_sign[0])

                dx = fixed[:, 0] - np.reshape(fixed[:, 0], (fixed.shape[0], 1))
                dy = fixed[:, 1] - np.reshape(fixed[:, 1], (fixed.shape[0], 1))
                mask = (abs(dx) > var.dx/2) + (abs(dy) > var.dy/2)

                for idx in range(1, fixed.shape[0]):
                    if all(mask[idx, :idx]):
                        fixed_new.append(fixed[idx])
                        fixed_sign_new.append(fixed_sign[idx])

            return np.array(fixed_new), np.array(fixed_sign_new)

        from .. import read
        from .. import math

        # Convert int_q string into list.
        if not isinstance(int_q, list):
            int_q = [int_q]
        self.params.int_q = int_q
        if any(np.array(self.params.int_q) == 'curly_A'):
            self.curly_A = []
        if any(np.array(self.params.int_q) == 'ee'):
            self.ee = []

        # Multi core setup.
        if not(np.isscalar(n_proc)) or (n_proc%1 != 0):
            print("error: invalid processor number")
            return -1
        queue = mp.Queue()

        # Write the tracing parameters.
        self.params = TracersParameterClass()
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
        self.params.destination = destination
        self.params.n_proc = n_proc

        # Make sure to read the var files with the correct magic.
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
            varfile = 'VAR' + str(ti)
            n_times = tf-ti+1
        else:
            series = False
            n_times = 1
        self.t = np.zeros(n_times)

        # Read the initial field.
        var = read.var(var_file=varfile, datadir=datadir, magic=magic,
                       quiet=True, trimall=True)
        self.t[0] = var.t
        grid = read.grid(datadir=datadir, quiet=True, trim=True)
        field = getattr(var, trace_field)
        param2 = read.param(datadir=datadir, param2=True, quiet=True)
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

        tracers = Tracers()
        # Create the mapping for all times.
        if not tracer_file_name:
            tracers.find_tracers(trace_field=trace_field, h_min=h_min, h_max=h_max,
                                 len_max=len_max, tol=tol,
                                 interpolation=interpolation,
                                 trace_sub=trace_sub, varfile=varfile, ti=ti, tf=tf,
                                 integration=integration, datadir=datadir,
                                 int_q=int_q, n_proc=n_proc)
        else:
            tracers.read(datadir=datadir, file_name=tracer_file_name)
        self.tracers = tracers

        # Set some default values.
        self.t = np.zeros((tf-ti+1)*series + (1-series))
        self.fidx = np.zeros((tf-ti+1)*series + (1-series))
        self.poincare = np.zeros([int(trace_sub*dim.nx),
                                  int(trace_sub*dim.ny), n_times])
        ix0 = range(0, int(self.params.nx*trace_sub)-1)
        iy0 = range(0, int(self.params.ny*trace_sub)-1)

        # Start the parallelized fixed point finding.
        for tidx in range(n_times):
            if tidx > 0:
                var = read.var(var_file='VAR'+str(tidx+ti), datadir=datadir,
                               magic=magic, quiet=True, trimall=True)
                field = getattr(var, trace_field)
                self.t[tidx] = var.t

            proc = []
            sub_data = []
            fixed = []
            fixed_sign = []
            for i_proc in range(n_proc):
                proc.append(mp.Process(target=__sub_fixed,
                                       args=(queue, ix0, iy0, field, self.tracers,
                                             tidx, var, i_proc)))
            for i_proc in range(n_proc):
                proc[i_proc].start()
            for i_proc in range(n_proc):
                sub_data.append(queue.get())
            for i_proc in range(n_proc):
                proc[i_proc].join()
            for i_proc in range(n_proc):
                # Extract the data from the single cores. Mind the order.
                sub_proc = sub_data[i_proc][0]
                fixed.extend(sub_data[i_proc][1])
                fixed_sign.extend(sub_data[i_proc][2])
                self.fidx[tidx] += sub_data[i_proc][3]
                self.poincare[sub_proc::n_proc, :, tidx] = sub_data[i_proc][4]
            for i_proc in range(n_proc):
                proc[i_proc].terminate()

            # Discard fixed points which lie too close to each other.
            fixed, fixed_sign = __discard_close_fixed_points(np.array(fixed),
                                                             np.array(fixed_sign),
                                                             var)
            self.fixed_points.append(np.array(fixed))
            self.fixed_sign.append(np.array(fixed_sign))

        # Compute the traced quantities along the fixed point streamlines.
        if any(np.array(self.params.int_q) == 'curly_A') or \
        any(np.array(self.params.int_q) == 'ee'):
            for t_idx in range(0, n_times):
                if any(np.array(self.params.int_q) == 'curly_A'):
                    self.curly_A.append([])
                if any(np.array(self.params.int_q) == 'ee'):
                    self.ee.append([])
                for fixed in self.fixed_points[t_idx]:
                    # Trace the stream line.
                    xx = np.array([fixed[0], fixed[1], self.params.Oz])
                    stream = Stream(field, self.params,
                                    h_min=self.params.h_min,
                                    h_max=self.params.h_max,
                                    len_max=self.params.len_max,
                                    tol=self.params.tol,
                                    interpolation=self.params.interpolation,
                                    integration=self.params.integration,
                                    xx=xx)
                    # Do the field line integration.
                    if any(np.array(self.params.int_q) == 'curly_A'):
                        curly_A = 0
                        for l in range(stream.stream_len-1):
                            aaInt = vec_int((stream.tracers[l+1] + stream.tracers[l])/2,
                                             var, var.aa, interpolation=self.params.interpolation)
                            curly_A += np.dot(aaInt, (stream.tracers[l+1] - stream.tracers[l]))
                        self.curly_A[-1].append(curly_A)
                    if any(np.array(self.params.int_q) == 'ee'):
                        ee_p = 0
                        for l in range(stream.stream_len-1):
                            eeInt = vec_int((stream.tracers[l+1] + stream.tracers[l])/2,
                                             var, ee, interpolation=self.params.interpolation)
                            ee_p += np.dot(eeInt, (stream.tracers[l+1] - stream.tracers[l]))
                        self.ee[-1].append(ee_p)
                if any(np.array(self.params.int_q) == 'curly_A'):
                    self.curly_A[-1] = np.array(self.curly_A[-1])
                if any(np.array(self.params.int_q) == 'ee'):
                    self.ee[-1] = np.array(self.ee[-1])


    def write(self, datadir='./data', destination='fixed_points.hdf5'):
        """
        Write the fixed points into a file.

        call signature::

        write(self, datadir='./data', destination='fixed_points.hdf5')

        Keyword arguments:

        *datadir*:
          Directory where the data is stored.

        *destination*:
          Destination file.
        """

        import numpy as np

        self.params.destination = destination

        # Write the results into hdf5 file.
        if destination != '':
            f = h5py.File(os.path.join(datadir, destination), 'w')
            # Write main data arrays.
            set_t = f.create_dataset("t", self.t.shape, dtype=self.t.dtype)
            set_t[...] = self.t[...]
            set_fidx = f.create_dataset("fidx", self.fidx.shape,
                                        dtype=self.fidx.dtype)
            set_fidx[...] = self.fidx[...]
            set_poincare = f.create_dataset("poincare", self.poincare.shape,
                                            dtype=self.poincare.dtype)
            set_poincare[...] = self.poincare[...]

            # Write the parameters into their own group.
            group_params = f.create_group('params')
            for key in dir(self.params):
                if not key.startswith('_'):
                    value = getattr(self.params, key)
                    group_params.attrs[key] = value

            # Create a new group for each time step.
            fixed_groups = []
            for t_idx in range(len(self.t)):
                fixed_groups.append(f.create_group('{0}'.format(t_idx)))
                set_fixed_points = fixed_groups[-1].create_dataset("fixed_points", self.fixed_points[t_idx].shape, dtype=self.fixed_points[t_idx].dtype)
                set_fixed_points[...] = self.fixed_points[t_idx]
                set_fixed_sign = fixed_groups[-1].create_dataset("fixed_sign", self.fixed_sign[t_idx].shape, dtype=self.fixed_sign[t_idx].dtype)
                set_fixed_sign[...] = self.fixed_sign[t_idx]
                if any(np.array(self.params.int_q) == 'curly_A'):
                    set_curly_A = fixed_groups[-1].create_dataset("curly_A", self.curly_A[t_idx].shape, dtype=self.curly_A[t_idx].dtype)
                    set_curly_A[...] = self.curly_A[t_idx]
                if any(np.array(self.params.int_q) == 'ee'):
                    set_ee = fixed_groups[-1].create_dataset("ee", self.ee[t_idx].shape, dtype=self.ee[t_idx].dtype)
                    set_ee[...] = self.ee[t_idx]
            f.close()

            # Write the tracers into their own file.
            tracer_file_name = self.params.destination[:-5] + '_tracers' + \
                               self.params.destination[-5:]
            self.tracers.write(datadir=self.params.datadir,
                               destination=tracer_file_name)
        else:
            print("error: empty destination file")


    def read(self, datadir='./data', file_name='fixed_points.hdf5'):
        """
        Read the fixed points from a file.

        call signature::

        read(self, datadir='./data', file_name='fixed_points.hdf5')

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
        self.fidx = f['fidx'].value
        self.poincare = f['poincare'].value

        # Extract parameters.
        params = f['params']
        self.params = TracersParameterClass()
        for param in params.attrs.keys():
            setattr(self.params, param, params.attrs[param])

        # Read the time series.
        self.fixed_points = []
        self.fixed_sign = []
        if any(np.array(self.params.int_q) == 'curly_A'):
            self.curly_A = []
        if any(np.array(self.params.int_q) == 'ee'):
            self.ee = []
        for t_idx in range(len(self.t)):
            group = f['{0}'.format(t_idx)]
            self.fixed_points.append(group['fixed_points'].value)
            self.fixed_sign.append(group['fixed_sign'].value)
            if any(np.array(self.params.int_q) == 'curly_A'):
                self.curly_A.append(group['curly_A'].value)
            if any(np.array(self.params.int_q) == 'ee'):
                self.ee.append(group['ee'].value)
        f.close()

        # Read the tracers from their own file.
        tracer_file_name = self.params.destination[:-5] + '_tracers' + \
                           self.params.destination[-5:]
        self.tracers = Tracers()
        self.tracers.read(datadir=datadir, file_name=tracer_file_name)
