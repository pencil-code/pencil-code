# streamlines.py
#
# Trace field lines for a given vector field.
#
# Written by Simon Candelaresi (iomsn1@gmail.com)
"""
Traces streamlines of a vector field from z0 to z1, similar to
'streamlines.f90'.
"""

'''
Test case:
nxyz = np.array([100, 100, 100])
oxyz = np.array([-1, -1, -1])
x = np.linspace(-1, 1, 100)
y = np.linspace(-1, 1, 100)
z = np.linspace(-1, 1, 100)
dxyz = np.array([x[1]-x[0], y[1]-y[0], z[1]-z[0]])
xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
bb = np.array([yy, -xx, np.zeros_like(zz)+.1])
#bb = np.array([np.zeros_like(zz), np.zeros_like(zz), np.zeros_like(zz)+1])
bb = bb.swapaxes(1, 3)
params = pc.diag.TracersParameterClass()
params.dx = dxyz[0]
params.dy = dxyz[1]
params.dz = dxyz[2]
params.Ox = x[0]
params.Oy = y[0]
params.Oz = z[0]
params.nx = len(x)
params.ny = len(y)
params.nz = len(z)
params.Lx = x[-1] - x[0]
params.Ly = y[-1] - y[0]
params.Lz = z[-1] - z[0]
params.interpolation = 'trilinear'
time = np.linspace(0, 20, 100)
stream = pc.calc.Stream(bb, params, xx=[0.5, 0.5, -1], time=time)
'''
class Stream(object):
    """
    Contains the methods and results for the streamline tracing for a field on a grid.
    """

    def __init__(self, field, params, xx=(0, 0, 0), time=(0, 1), metric=None, splines=None):
        """
        Trace a field starting from xx in any rectilinear coordinate system
        with constant dx, dy and dz and with a given metric.

        call signature:

          tracer(field, params, xx=[0, 0, 0], time=[0, 1], metric=None, splines=None):

        Keyword arguments:

        *field*:
          Vector field which is integrated over with shape [n_vars, nz, ny, nx].
          Its elements are the components of the field using unnormed
          unit-coordinate vectors.

        *params*:
          Simulation and tracer parameters.

        *xx*:
          Starting point of the field line integration with starting time.

        *time*:
            Time array for which the tracer is computed.

        *metric*:
            Metric function that takes a point [x, y, z] and an array
            of shape [3, 3] that has the comkponents g_ij.
            Use 'None' for Cartesian metric.

        *splines*:
            Spline interpolation functions for the tricubic interpolation.
            This can speed up the calculations greatly for repeated streamline tracing on the same data.
            Accepts a list of the spline functions for the three vector components.
        """

        import numpy as np
        from scipy.integrate import solve_ivp
        from ..math.interpolation import vec_int

        if params.interpolation == 'tricubic':
            try:
                import warnings

                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=Warning)
                    from eqtools.trispline import Spline
            except:
                print('Warning: Could not import eqtools.trispline.Spline for tricubic interpolation.\n')
                print('Warning: Fall back to trilinear.')
                params.interpolation = 'trilinear'

        if not metric:
            metric = lambda xx: np.eye(3)

        dxyz = np.array([params.dx, params.dy, params.dz])
        oxyz = np.array([params.Ox, params.Oy, params.Oz])
        nxyz = np.array([params.nx, params.ny, params.nz])

        # Redefine the derivative y for the scipy ode integrator using the given parameters.
        if (params.interpolation == 'mean') or (params.interpolation == 'trilinear'):
            odeint_func = lambda t, xx: vec_int(xx, field, dxyz, oxyz, nxyz, params.interpolation)
        if params.interpolation == 'tricubic':
            x = np.linspace(params.Ox, params.Ox+params.Lx, params.nx)
            y = np.linspace(params.Oy, params.Oy+params.Ly, params.ny)
            z = np.linspace(params.Oz, params.Oz+params.Lz, params.nz)
            if splines is None:
                field_x = Spline(z, y, x, field[0, ...])
                field_y = Spline(z, y, x, field[1, ...])
                field_z = Spline(z, y, x, field[2, ...])
            else:
                field_x = splines[0]
                field_y = splines[1]
                field_z = splines[2]
            odeint_func = lambda t, xx: self.trilinear_func(xx, field_x, field_y, field_z, params)
            del(x)
            del(y)
            del(z)
        # Set up the ode solver.
        self.tracers = solve_ivp(odeint_func, (time[0], time[-1]), xx,
                                 t_eval=time, rtol=params.rtol, atol=params.atol,
                                 jac=metric, method=params.method).y.T

        # Remove points that lie outside the domain and interpolation on the boundary.
        cut_mask = ((self.tracers[:, 0] > params.Ox+params.Lx) + \
                    (self.tracers[:, 0] < params.Ox))*(not params.periodic_x) + \
                   ((self.tracers[:, 1] > params.Oy+params.Ly) + \
                    (self.tracers[:, 1] < params.Oy))*(not params.periodic_y) + \
                   ((self.tracers[:, 2] > params.Oz+params.Lz) + \
                    (self.tracers[:, 2] < params.Oz))*(not params.periodic_z)
        if np.sum(cut_mask) > 0:
            # Find the first point that lies outside.
            idx_outside = np.min(np.where(cut_mask))
            # Interpolate.
            p0 = self.tracers[idx_outside-1, :]
            p1 = self.tracers[idx_outside, :]
            lam = np.zeros([6])
            if p0[0] == p1[0]:
                lam[0] = np.inf
                lam[1] = np.inf
            else:
                lam[0] = (params.Ox + params.Lx - p0[0])/(p1[0] - p0[0])
                lam[1] = (params.Ox - p0[0])/(p1[0] - p0[0])
            if p0[1] == p1[1]:
                lam[2] = np.inf
                lam[3] = np.inf
            else:
                lam[2] = (params.Oy + params.Ly - p0[1])/(p1[1] - p0[1])
                lam[3] = (params.Oy - p0[1])/(p1[1] - p0[1])
            if p0[2] == p1[2]:
                lam[4] = np.inf
                lam[5] = np.inf
            else:
                lam[4] = (params.Oz + params.Lz - p0[2])/(p1[2] - p0[2])
                lam[5] = (params.Oz - p0[2])/(p1[2] - p0[2])
            lam_min = np.min(lam[lam >= 0])
            if abs(lam_min) == np.inf:
                lam_min = 0
            self.tracers[idx_outside, :] = p0 + lam_min*(p1-p0)
            # We don't want to cut the interpolated point (was first point outside).
            cut_mask[idx_outside] = False
            cut_mask[idx_outside+1:] = True
            # Remove outside points.
            self.tracers = self.tracers[~cut_mask, :].copy()

        self.params = params
        self.xx = xx
        self.time = time

        # Compute the length of the line segments.
        middle_point_list = (self.tracers[1:, :]+self.tracers[:-1, :])/2
        diff_vectors = (self.tracers[1:, :]-self.tracers[:-1, :])
        self.section_l = np.zeros(diff_vectors.shape[0])
        for idx, middle_point in enumerate(middle_point_list):
            self.section_l[idx] = np.sqrt(np.sum(diff_vectors[idx]*np.sum(metric(middle_point)*diff_vectors[idx], axis=1)))
        self.total_l = np.sum(self.section_l)

        self.iterations = len(time)
        self.section_dh = time[1:] - time[:-1]
        self.total_h = time[-1] - time[0]


    def trilinear_func(self, xx, field_x, field_y, field_z, params):
        """
        Trilinear spline interpolation like eqtools.trispline.Spline
        but return 0 if the point lies outside the box.

        call signature:

        trilinear_func(xx, field_x, field_y, field_z, params)

        Keyword arguments:

        *xx*:
          The zyx coordinates of the point to interpolate the data.

        *field_xyz*:
          The Spline objects for the velocity fields.

        *params*:
          Parameter object.
        """

        import numpy as np

        if (xx[0] < params.Ox) + (xx[0] > params.Ox + params.Lx) + \
           (xx[1] < params.Oy) + (xx[1] > params.Oy + params.Ly) + \
           (xx[2] < params.Oz) + (xx[2] > params.Oz + params.Lz):
            return np.zeros(3)
        return np.array([field_x.ev(xx[2], xx[1], xx[0]),
                         field_y.ev(xx[2], xx[1], xx[0]),
                         field_z.ev(xx[2], xx[1], xx[0])])[:, 0]



#    """
#    Stream -- Holds the traced streamline.
#    """
#
#    def __init__(self, xx, field, params):
#        """
#        Creates the traced streamline for a specified vector field.
#
#        call signature:
#
#          Stream(field, params, xx):
#
#        Keyword arguments:
#
#         *xx*:
#            Seed position.
#
#         *field*:
#            Vector field which is integrated over.
#
#         *params*:
#           Simulation and tracer parameters.
#        """
#
#        import numpy as np
#        from ..math.interpolation import vec_int
#
#        # Tentative streamline length.
#        self.tracers = np.zeros([int(params.iter_max), 3], dtype='float32')
#
#        tol2 = params.tol**2
#        dh = np.sqrt(params.h_max*params.h_min) # Initial step size.
#
#        # Declare some vectors.
#        xMid = np.zeros(3)
#        xSingle = np.zeros(3)
#        xHalf = np.zeros(3)
#        xDouble = np.zeros(3)
#        dxyz = np.array([params.dx, params.dy, params.dz])
#        oxyz = np.array([params.Ox, params.Oy, params.Oz])
#        nxyz = np.array([params.nx, params.ny, params.nz])
#
#        # Initialize the coefficient for the 6th order adaptive time step RK.
#        a = np.zeros(6); b = np.zeros((6, 5)); c = np.zeros(6); cs = np.zeros(6)
#        k = np.zeros((6, 3))
#        a[1] = 0.2; a[2] = 0.3; a[3] = 0.6; a[4] = 1; a[5] = 0.875
#        b[1, 0] = 0.2
#        b[2, 0] = 3/40.; b[2, 1] = 9/40.
#        b[3, 0] = 0.3; b[3, 1] = -0.9; b[3, 2] = 1.2
#        b[4, 0] = -11/54.; b[4, 1] = 2.5; b[4, 2] = -70/27.; b[4, 3] = 35/27.
#        b[5, 0] = 1631/55296.; b[5, 1] = 175/512.; b[5, 2] = 575/13824.
#        b[5, 3] = 44275/110592.; b[5, 4] = 253/4096.
#        c[0] = 37/378.; c[2] = 250/621.; c[3] = 125/594.; c[5] = 512/1771.
#        cs[0] = 2825/27648.; cs[2] = 18575/48384.; cs[3] = 13525/55296.
#        cs[4] = 277/14336.; cs[5] = 0.25
#
#        # Do the streamline tracing.
#        self.tracers[0, :] = xx
#        outside = False
#        stream_len = 0
#        length = 0
#
#        if params.integration == 'simple':
#            while ((length < params.len_max) and (stream_len < params.iter_max-1) and
#                   (not np.isnan(xx[0])) and (outside == False)):
#                # (a) single step (midpoint method)
#                xMid = xx + 0.5*dh*vec_int(xx, field, dxyz, oxyz, nxyz, params.interpolation)
#                xSingle = xx + dh*vec_int(xMid, field, dxyz, oxyz, nxyz, params.interpolation)
#
#                # (b) two steps with half stepsize
#                xMid = xx + 0.25*dh*vec_int(xx, field, dxyz, oxyz, nxyz, params.interpolation)
#                xHalf = xx + 0.5*dh*vec_int(xMid, field, dxyz, oxyz, nxyz, params.interpolation)
#                xMid = xHalf + 0.25*dh*vec_int(xHalf, field, dxyz, oxyz, nxyz, params.interpolation)
#                xDouble = xHalf + 0.5*dh*vec_int(xMid, field, dxyz, oxyz, nxyz, params.interpolation)
#
#                # (c) Check error (difference between methods).
#                dist2 = np.sum((xSingle-xDouble)**2)
#                if dist2 > tol2:
#                    dh = 0.5*dh
#                    if abs(dh) < params.h_min:
#                        print("Error: stepsize underflow")
#                        break
#                else:
#                    length += np.sqrt(np.sum((xx-xDouble)**2))
#                    xx = xDouble.copy()
#                    if abs(dh) < params.h_min:
#                        dh = 2*dh
#                    stream_len += 1
#                    self.tracers[stream_len, :] = xx.copy()
#                    if (dh > params.h_max) or (np.isnan(dh)):
#                        dh = params.h_max
#                    # Check if this point lies outside the domain.
#                    if ((xx[0] < params.Ox-params.dx) or
#                            (xx[0] > params.Ox+params.Lx+params.dx) or
#                            (xx[1] < params.Oy-params.dy) or
#                            (xx[1] > params.Oy+params.Ly+params.dy) or
#                            (xx[2] < params.Oz) or (xx[2] > params.Oz+params.Lz)):
#                        outside = True
#
#        if params.integration == 'RK6':
#            while ((length < params.len_max) and (stream_len < params.iter_max-1) and
#                   (not np.isnan(xx[0])) and (outside == False)):
#                k[0, :] = dh*vec_int(xx, field, dxyz, oxyz, nxyz, params.interpolation)
#                k[1, :] = dh*vec_int(xx + b[1, 0]*k[0, :], field, dxyz, oxyz, nxyz,
#                                     params.interpolation)
#                k[2, :] = dh*vec_int(xx + b[2, 0]*k[0, :] + b[2, 1]*k[1, :],
#                                     field, dxyz, oxyz, nxyz, params.interpolation)
#                k[3, :] = dh*vec_int(xx + b[3, 0]*k[0, :] + b[3, 1]*k[1, :] +
#                                     b[3, 2]*k[2, :], field, dxyz, oxyz, nxyz, params.interpolation)
#                k[4, :] = dh*vec_int(xx + b[4, 0]*k[0, :] + b[4, 1]*k[1, :] +
#                                     b[4, 2]*k[2, :] + b[4, 3]*k[3, :],
#                                     field, dxyz, oxyz, nxyz, params.interpolation)
#                k[5, :] = dh*vec_int(xx + b[5, 0]*k[0, :] + b[5, 1]*k[1, :] +
#                                     b[5, 2]*k[2, :] + b[5, 3]*k[3, :] +
#                                     b[5, 4]*k[4, :], field, dxyz, oxyz, nxyz, params.interpolation)
#
#                xNew = xx + c[0]*k[0, :] + c[1]*k[1, :]  + c[2]*k[2, :]  + \
#                       c[3]*k[3, :] + c[4]*k[4, :] + c[5]*k[5, :]
#                xNewS = xx + cs[0]*k[0, :] + cs[1]*k[1, :] + cs[2]*k[2, :] + \
#                        cs[3]*k[3, :] + cs[4]*k[4, :] + cs[5]*k[5, :]
#
#                delta2 = np.dot((xNew-xNewS), (xNew-xNewS))
#                delta = np.sqrt(delta2)
#
#                if delta2 > tol2:
#                    dh = dh*(0.9*abs(params.tol/delta))**0.2
#                    if abs(dh) < params.h_min:
#                        print("Error: step size underflow")
#                        break
#                else:
#                    length += np.sqrt(np.sum((xx-xNew)**2))
#                    xx = xNew
#                    if abs(dh) < params.h_min:
#                        dh = 2*dh
#                    stream_len += 1
#                    self.tracers[stream_len, :] = xx
#                    if (dh > params.h_max) or (np.isnan(dh)):
#                        dh = params.h_max
#                    # Check if this point lies outside the domain.
#                    if ((xx[0] < params.Ox-params.dx) or
#                            (xx[0] > params.Ox+params.Lx+params.dx) or
#                            (xx[1] < params.Oy-params.dy) or
#                            (xx[1] > params.Oy+params.Ly+params.dy) or
#                            (xx[2] < params.Oz) or (xx[2] > params.Oz+params.Lz)):
#                        outside = True
#                if (dh > params.h_max) or (delta == 0) or (np.isnan(dh)):
#                    dh = params.h_max
#
#                print(xNew, xNewS)
#        # Linearly interpolate if the last point lies above.
#        if outside and (xx[2] > params.Oz+params.Lz):
#            weight = (params.Oz+params.Lz - self.tracers[stream_len-1, 2])/(self.tracers[stream_len, 2] - self.tracers[stream_len-1, 2])
#            self.tracers[stream_len, :] = weight*(self.tracers[stream_len, :] - self.tracers[stream_len-1, :]) + self.tracers[stream_len-1, :]
#
#        self.tracers = np.resize(self.tracers, (stream_len+1, 3))
#        self.len = length
#        self.stream_len = stream_len
#        self.params = params
