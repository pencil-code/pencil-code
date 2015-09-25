# fixed_points.py
# Written by Simon Candelaresi (iomsn1@gmail.com)
"""
Creates the fixed point values.
"""

import numpy as np
import pencil as pc
import math as m
import threading
from pencilnew.diag.tracers import TracersParameterClass
from pencilnew.diag.tracers import Tracers
from pencilnew.calc.streamlines import Stream


class FixedPoint(object):
    """
    fixed_points -- Holds the fixed points and additional integrated quantities.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.params = TracersParameterClass()
        self.t = []
        self.fidx = []
        self.fixed_points = []
        self.poincare = []


    def find_fixed(self, data_dir='data/', destination='fixed_points.hf5',
                   varfile='VAR0', ti=-1, tf=-1, trace_field='bb', h_min=2e-3,
                   h_max=2e4, len_max=500, tol=1e-2, interpolation='weighted',
                   trace_sub=1, integration='simple', int_q=[''], n_proc=1):
        """
        Find the fixed points.

        call signature::

        find_fixed(data_dir='data/', destination='fixed_points.hf5',
                   varfile='VAR0', ti=-1, tf=-1, trace_field='bb', h_min=2e-3,
                   h_max=2e4, len_max=500, tol=1e-2, interpolation='weighted',
                   trace_sub=1, integration='simple', int_q=[''], n_proc=1):

        Finds the fixed points. Returns the fixed points positions.

        Keyword arguments:

          *data_dir*:
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
           'weighted': weights the adjacent grid points according to
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
        """


        # Return the fixed points for a subset of the domain for the initial fixed points.
        def __sub_fixed_init(ix0, iy0, field, tracers, var, fixed, i_proc):
            diff = np.zeros((4, 2))

            for ix in ix0[i_proc::self.params.n_proc]:
                for iy in iy0:
                    # Compute Poincare index around this cell (!= 0 for potential fixed point).
                    diff[0, :] = np.array([tracers.x1[ix, iy, 0] - tracers.x0[ix, iy, 0],
                                           tracers.y1[ix, iy, 0] - tracers.y0[ix, iy, 0]])
                    diff[1, :] = np.array([tracers.x1[ix+1, iy, 0] - tracers.x0[ix+1, iy, 0],
                                           tracers.y1[ix+1, iy, 0] - tracers.y0[ix+1, iy, 0]])
                    diff[2, :] = np.array([tracers.x1[ix+1, iy+1, 0] - tracers.x0[ix+1, iy+1, 0],
                                           tracers.y1[ix+1, iy+1, 0] - tracers.y0[ix+1, iy+1, 0]])
                    diff[3, :] = np.array([tracers.x1[ix, iy+1, 0] - tracers.x0[ix, iy+1, 0],
                                           tracers.y1[ix, iy+1, 0] - tracers.y0[ix, iy+1, 0]])
                    if sum(np.sum(diff**2, axis=1) != 0):
                        diff = np.swapaxes(np.swapaxes(diff, 0, 1) /
                               np.sqrt(np.sum(diff**2, axis=1)), 0, 1)
                    poincare = __poincare_index(field, tracers.x0[ix:ix+2, iy, 0],
                                                tracers.y0[ix, iy:iy+2, 0], diff)
                    self.poincare[ix, iy] = poincare

                    if abs(poincare) > 5: # Use 5 instead of 2*pi to account for rounding errors.
                        # Subsample to get starting point for iteration.
                        nt = 4
                        xmin = tracers.x0[ix, iy, 0]
                        ymin = tracers.y0[ix, iy, 0]
                        xmax = tracers.x0[ix+1, iy, 0]
                        ymax = tracers.y0[ix, iy+1, 0]
                        xx = np.zeros((nt**2, 3))
                        tracers_sub = np.zeros((nt**2, 5))
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
                            tracers_sub[it1, 0:2] = xx[it1, 0:2]
                            tracers_sub[it1, 2:] = stream.tracers[stream.stream_len-1, :]
                        min2 = 1e6
                        minx = xmin
                        miny = ymin
                        i1 = 0
                        for j1 in range(nt):
                            for k1 in range(nt):
                                diff2 = (tracers_sub[i1+k1*nt, 2] - tracers_sub[i1+k1*nt, 0])**2 + \
                                        (tracers_sub[i1+k1*nt, 3] - tracers_sub[i1+k1*nt, 1])**2
                                if diff2 < min2:
                                    min2 = diff2
                                    minx = xmin + j1/(nt-1.)*(xmax - xmin)
                                    miny = ymin + k1/(nt-1.)*(ymax - ymin)
                                it1 += 1

                        # Get fixed point from this starting position using Newton's method.
                        point = np.array([minx, miny])
                        fixed_point = __null_point(point, var)

                        # Check if fixed point lies inside the cell.
                        if ((fixed_point[0] < tracers.x0[ix, iy, 0]) or
                            (fixed_point[0] > tracers.x0[ix+1, iy, 0]) or
                            (fixed_point[1] < tracers.y0[ix, iy, 0]) or
                            (fixed_point[1] > tracers.y0[ix, iy+1, 0])):
                            print "warning: fixed point lies outside the cell"
                        else:
                            fixed.append(fixed_point)
                            self.fidx += 1


        # Finds the Poincare index of this grid cell.
        def __poincare_index(field, sx, sy, diff):
            poincare = 0
            poincare += __edge(field, [sx[0], sx[1]], [sy[0], sy[0]], diff[0, :], diff[1, :], 0)
            poincare += __edge(field, [sx[1], sx[1]], [sy[0], sy[1]], diff[1, :], diff[2, :], 0)
            poincare += __edge(field, [sx[1], sx[0]], [sy[1], sy[1]], diff[2, :], diff[3, :], 0)
            poincare += __edge(field, [sx[0], sx[0]], [sy[1], sy[0]], diff[3, :], diff[0, :], 0)
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
                stream = Stream(field, self.params, h_min=self.params.h_min, h_max=self.params.h_max,
                                len_max=self.params.len_max, tol=self.params.tol,
                                interpolation=self.params.interpolation,
                                integration=self.params.integration, xx=np.array([xm, ym, self.params.Oz]))
                stream_x0 = stream.tracers[0, 0]
                stream_y0 = stream.tracers[0, 1]
                stream_x1 = stream.tracers[stream.stream_len-1, 0]
                stream_y1 = stream.tracers[stream.stream_len-1, 1]
                stream_z1 = stream.tracers[stream.stream_len-1, 2]

                # Discard any streamline which does not converge or hits the boundary.
                if ((stream.len >= len_max) or
                (stream_z1 < self.params.Oz+self.params.Lz-self.params.dz)):
                    dtot = 0.
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
                    print "warning: Newton did not converged"
                    break

                it += 1
                
            return fixed_point


        # Find the fixed point using Newton's method, starting at previous fixed point.
        def __sub_fixed_series(t_idx, field, var, fixed, i_proc):
            for point in self.fixed_points[t_idx-1][i_proc::self.params.n_proc]:
                fixed_tentative = __null_point(point, var)
                # Check if the fixed point lies outside the domain.
                if fixed_tentative[0] >= self.params.Ox and \
                fixed_tentative[1] >= self.params.Oy and \
                fixed_tentative[0] <= self.params.Ox+self.params.Lx and \
                fixed_tentative[1] <= self.params.Oy+self.params.Ly:
                    fixed.append(fixed_tentative)
                    
                    
        # Convert int_q string into list.
        if not isinstance(int_q, list):
            int_q = [int_q]

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
        self.params.data_dir = data_dir
        self.params.destination = destination
        self.params.n_proc = n_proc

        # Multi core setup.
        if not(np.isscalar(n_proc)) or (n_proc%1 != 0):
            print "error: invalid processor number"
            return -1

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
        dim = pc.read_dim(datadir=data_dir)

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
        var = pc.read_var(varfile=varfile, datadir=data_dir, magic=magic,
                          quiet=True, trimall=True)
        self.t[0] = var.t
        grid = pc.read_grid(datadir=data_dir, quiet=True, trim=True)
        field = getattr(var, trace_field)

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

        # Create the initial mapping.
        tracers = Tracers()
        tracers.find_tracers(trace_field='bb', h_min=h_min, h_max=h_max,
                             len_max=len_max, tol=tol,
                             interpolation=interpolation,
                             trace_sub=trace_sub, varfile=varfile,
                             integration=integration, data_dir=data_dir,
                             int_q=int_q, n_proc=n_proc)
        self.tracers = tracers

        # Set some default values.
        self.t = np.zeros((tf-ti+1)*series + (1-series))
        self.fidx = np.zeros((tf-ti+1)*series + (1-series))
        if any(np.array(int_q) == 'curlyA'):
            self.curly_A = np.zeros([int(trace_sub*dim.nx),
                                     int(trace_sub*dim.ny), n_times])
        if any(np.array(int_q) == 'ee'):
            self.ee = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny),
                                n_times])
        self.poincare = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny)])
        ix0 = range(0, int(self.params.nx*trace_sub)-1)
        iy0 = range(0, int(self.params.ny*trace_sub)-1)

        # Start the parallelized fixed point finding for the initial time.
        threads = []
        fixed = []
        for i_proc in range(n_proc):
            threads.append(threading.Thread(target=__sub_fixed_init,
                                            args=(ix0, iy0, field, tracers,
                                                  var, fixed, i_proc)))
            threads[-1].start()
        for i_proc in range(n_proc):
            threads[i_proc].join()
        self.fixed_points.append(np.array(fixed))

        # Find the fixed points for the remaining times.
        for t_idx in range(1, n_times):
            # Read the data.
            varfile = 'VAR' + str(t_idx+ti)
            var = pc.read_var(varfile=varfile, datadir=data_dir, magic=magic,
                              quiet=True, trimall=True)
            field = getattr(var, trace_field)
            self.t[t_idx] = var.t
            
            # Find the new fixed points.
            threads = []
            fixed = []
            for i_proc in range(int(np.min((n_proc, self.fidx[t_idx-1])))):
                threads.append(threading.Thread(target=__sub_fixed_series,
                                                args=(t_idx, field, var,
                                                      fixed, i_proc)))
                threads[-1].start()
            for i_proc in range(n_proc):
                threads[i_proc].join()
            self.fixed_points.append(np.array(fixed))
            self.fidx[t_idx] = len(fixed)


def read_fixed_points(data_dir='data/', fileName='fixed_points.dat', hm=1):
    """
    Reads the fixed points files.

    call signature::

      fixed = read_fixed_points(data_dir='data/', fileName='fixed_points.dat', hm=1)

    Reads from the fixed points files. Returns the fixed points positions.

    Keyword arguments:

      *data_dir*:
        Data directory.

      *fileName*:
        Name of the fixed points file.

      *hm*:
        Header multiplication factor in case Fortran's binary data writes extra large
        header. For most cases hm = 1 is sufficient. For the cluster in St Andrews use hm = 2.
    """


    # read the cpu structure
    dim = pc.read_dim(data_dir=data_dir)
    if dim.n_procz > 1:
        print "error: number of cores in z-direction > 1"

    data = []

    # read the data
    fixed_file = open(data_dir+fileName, 'rb')
    tmp = fixed_file.read(4*hm)

    data = pc.fixed_struct()
    eof = 0
    # length of longest array of fixed points
    fixedMax = 0
    if tmp == '':
        eof = 1
    while eof == 0:
        data.t.append(struct.unpack("<"+str(hm+1)+"f", fixed_file.read(4*(hm+1)))[0])
        n_fixed = int(struct.unpack("<"+str(2*hm+1)+"f", fixed_file.read(4*(2*hm+1)))[1+hm/2])

        x = list(np.zeros(n_fixed))
        y = list(np.zeros(n_fixed))
        q = list(np.zeros(n_fixed))
        for j in range(n_fixed):
            x[j] = struct.unpack("<"+str(hm+1)+"f", fixed_file.read(4*(hm+1)))[-1]
            y[j] = struct.unpack("<f", fixed_file.read(4))[0]
            q[j] = struct.unpack("<"+str(hm+1)+"f", fixed_file.read(4*(hm+1)))[0]
        data.x.append(x)
        data.y.append(y)
        data.q.append(q)
        data.fidx.append(n_fixed)

        tmp = fixed_file.read(4*hm)
        if tmp == '':
            eof = 1

        if fixedMax < len(x):
            fixedMax = len(x)

    fixed_file.close()

    # add NaN to fill up the times with smaller number of fixed points
    fixed = pc.fixed_struct()
    for i in range(len(data.t)):
        annex = list(np.zeros(fixedMax - len(data.x[i])) + np.nan)
        fixed.t.append(data.t[i])
        fixed.x.append(data.x[i] + annex)
        fixed.y.append(data.y[i] + annex)
        fixed.q.append(data.q[i] + annex)
        fixed.fidx.append(data.fidx[i])

    fixed.t = np.array(fixed.t)
    fixed.x = np.array(fixed.x)
    fixed.y = np.array(fixed.y)
    fixed.q = np.array(fixed.q)
    fixed.fidx = np.array(fixed.fidx)

    return fixed


class fixed_struct:
    def __init__(self):
        self.t = []
        self.fidx = [] # number of fixed points at this time
        self.x = []
        self.y = []
        self.q = []

