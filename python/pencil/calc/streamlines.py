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
nxyz = np.array([10, 10, 10])
oxyz = np.array([-1, -1, -1])
x = np.linspace(-1, 1, 10)
y = np.linspace(-1, 1, 10)
z = np.linspace(-1, 1, 10)
dxyz = np.array([x[1]-x[0], y[1]-y[0], z[1]-z[0]])
xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
bb = np.array([yy, -xx, np.zeros_like(zz)])
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
params.h_max = 1
stream = pc.calc.Stream([0.5, 0.5, -0.9], bb, params)
'''
class Stream(object):
    """
    Stream -- Holds the traced streamline.
    """

    def __init__(self, xx, field, params):
        """
        Creates the traced streamline for a specified vector field.

        call signature:

          Stream(field, params, xx):

        Keyword arguments:

         *xx*:
            Seed position.

         *field*:
            Vector field which is integrated over.

         *params*:
           Simulation and tracer parameters.
        """

        import numpy as np
        from ..math.interpolation import vec_int

        # Tentative streamline length.
        self.tracers = np.zeros([int(params.iter_max), 3], dtype='float32')

        tol2 = params.tol**2
        dh = np.sqrt(params.h_max*params.h_min) # Initial step size.

        # Declare some vectors.
        xMid = np.zeros(3)
        xSingle = np.zeros(3)
        xHalf = np.zeros(3)
        xDouble = np.zeros(3)
        dxyz = np.array([params.dx, params.dy, params.dz])
        oxyz = np.array([params.Ox, params.Oy, params.Oz])
        nxyz = np.array([params.nx, params.ny, params.nz])

        # Initialize the coefficient for the 6th order adaptive time step RK.
        a = np.zeros(6); b = np.zeros((6, 5)); c = np.zeros(6); cs = np.zeros(6)
        k = np.zeros((6, 3))
        a[1] = 0.2; a[2] = 0.3; a[3] = 0.6; a[4] = 1; a[5] = 0.875
        b[1, 0] = 0.2
        b[2, 0] = 3/40.; b[2, 1] = 9/40.
        b[3, 0] = 0.3; b[3, 1] = -0.9; b[3, 2] = 1.2
        b[4, 0] = -11/54.; b[4, 1] = 2.5; b[4, 2] = -70/27.; b[4, 3] = 35/27.
        b[5, 0] = 1631/55296.; b[5, 1] = 175/512.; b[5, 2] = 575/13824.
        b[5, 3] = 44275/110592.; b[5, 4] = 253/4096.
        c[0] = 37/378.; c[2] = 250/621.; c[3] = 125/594.; c[5] = 512/1771.
        cs[0] = 2825/27648.; cs[2] = 18575/48384.; cs[3] = 13525/55296.
        cs[4] = 277/14336.; cs[5] = 0.25

        # Do the streamline tracing.
        self.tracers[0, :] = xx
        outside = False
        stream_len = 0
        length = 0

        if params.integration == 'simple':
            while ((length < params.len_max) and (stream_len < params.iter_max-1) and
                   (not np.isnan(xx[0])) and (outside == False)):
                # (a) single step (midpoint method)
                xMid = xx + 0.5*dh*vec_int(xx, field, dxyz, oxyz, nxyz, params.interpolation)
                xSingle = xx + dh*vec_int(xMid, field, dxyz, oxyz, nxyz, params.interpolation)

                # (b) two steps with half stepsize
                xMid = xx + 0.25*dh*vec_int(xx, field, dxyz, oxyz, nxyz, params.interpolation)
                xHalf = xx + 0.5*dh*vec_int(xMid, field, dxyz, oxyz, nxyz, params.interpolation)
                xMid = xHalf + 0.25*dh*vec_int(xHalf, field, dxyz, oxyz, nxyz, params.interpolation)
                xDouble = xHalf + 0.5*dh*vec_int(xMid, field, dxyz, oxyz, nxyz, params.interpolation)

                # (c) Check error (difference between methods).
                dist2 = np.sum((xSingle-xDouble)**2)
                if dist2 > tol2:
                    dh = 0.5*dh
                    if abs(dh) < params.h_min:
                        print("Error: stepsize underflow")
                        break
                else:
                    length += np.sqrt(np.sum((xx-xDouble)**2))
                    xx = xDouble.copy()
                    if abs(dh) < params.h_min:
                        dh = 2*dh
                    stream_len += 1
                    self.tracers[stream_len, :] = xx.copy()
                    if (dh > params.h_max) or (np.isnan(dh)):
                        dh = params.h_max
                    # Check if this point lies outside the domain.
                    if ((xx[0] < params.Ox-params.dx) or
                            (xx[0] > params.Ox+params.Lx+params.dx) or
                            (xx[1] < params.Oy-params.dy) or
                            (xx[1] > params.Oy+params.Ly+params.dy) or
                            (xx[2] < params.Oz) or (xx[2] > params.Oz+params.Lz)):
                        outside = True

        if params.integration == 'RK6':
            while ((length < params.len_max) and (stream_len < params.iter_max-1) and
                   (not np.isnan(xx[0])) and (outside == False)):
                k[0, :] = dh*vec_int(xx, field, dxyz, oxyz, nxyz, params.interpolation)
                k[1, :] = dh*vec_int(xx + b[1, 0]*k[0, :], field, dxyz, oxyz, nxyz,
                                     params.interpolation)
                k[2, :] = dh*vec_int(xx + b[2, 0]*k[0, :] + b[2, 1]*k[1, :],
                                     field, dxyz, oxyz, nxyz, params.interpolation)
                k[3, :] = dh*vec_int(xx + b[3, 0]*k[0, :] + b[3, 1]*k[1, :] +
                                     b[3, 2]*k[2, :], field, dxyz, oxyz, nxyz, params.interpolation)
                k[4, :] = dh*vec_int(xx + b[4, 0]*k[0, :] + b[4, 1]*k[1, :] +
                                     b[4, 2]*k[2, :] + b[4, 3]*k[3, :],
                                     field, dxyz, oxyz, nxyz, params.interpolation)
                k[5, :] = dh*vec_int(xx + b[5, 0]*k[0, :] + b[5, 1]*k[1, :] +
                                     b[5, 2]*k[2, :] + b[5, 3]*k[3, :] +
                                     b[5, 4]*k[4, :], field, dxyz, oxyz, nxyz, params.interpolation)

                xNew = xx + c[0]*k[0, :]  + c[1]*k[1, :]  + c[2]*k[2, :]  + \
                       c[3]*k[3, :]  + c[4]*k[4, :]  + c[5]*k[5, :]
                xNewS = xx + cs[0]*k[0, :] + cs[1]*k[1, :] + cs[2]*k[2, :] + \
                        cs[3]*k[3, :] + cs[4]*k[4, :] + cs[5]*k[5, :]

                delta2 = np.dot((xNew-xNewS), (xNew-xNewS))
                delta = np.sqrt(delta2)

                if delta2 > tol2:
                    dh = dh*(0.9*abs(params.tol/delta))**0.2
                    if abs(dh) < params.h_min:
                        print("Error: step size underflow")
                        break
                else:
                    length += np.sqrt(np.sum((xx-xNew)**2))
                    xx = xNew
                    if abs(dh) < params.h_min:
                        dh = 2*dh
                    stream_len += 1
                    self.tracers[stream_len, :] = xx
                    if (dh > params.h_max) or (np.isnan(dh)):
                        dh = params.h_max
                    # Check if this point lies outside the domain.
                    if ((xx[0] < params.Ox-params.dx) or
                            (xx[0] > params.Ox+params.Lx+params.dx) or
                            (xx[1] < params.Oy-params.dy) or
                            (xx[1] > params.Oy+params.Ly+params.dy) or
                            (xx[2] < params.Oz) or (xx[2] > params.Oz+params.Lz)):
                        outside = True
                if (dh > params.h_max) or (delta == 0) or (np.isnan(dh)):
                    dh = params.h_max

        # Linearly interpolate if the last point lies above.
        if outside and (xx[2] > params.Oz+params.Lz):
            weight = (params.Oz+params.Lz - self.tracers[stream_len-1, 2])/(self.tracers[stream_len, 2] - self.tracers[stream_len-1, 2])
            self.tracers[stream_len, :] = weight*(self.tracers[stream_len, :] - self.tracers[stream_len-1, :]) + self.tracers[stream_len-1, :]

        self.tracers = np.resize(self.tracers, (stream_len+1, 3))
        self.len = length
        self.stream_len = stream_len
        self.params = params
