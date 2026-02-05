# streamlines.py
#
# Trace field lines for a given vector field.
#
# Written by Simon Candelaresi (iomsn1@gmail.com)
"""
Traces streamlines of a vector field from z0 to z1, similar to
'streamlines.f90'.
"""

"""
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
"""


class Stream(object):
    """
    Contains the methods and results for the streamline tracing for a field on a grid.
    """

    def __init__(
        self, field, params, xx=(0, 0, 0), time=(0, 1), metric=None, splines=None
    ):
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
            This can speed up the calculations greatly for repeated streamline
            tracing on the same data.
            Accepts a list of the spline functions for the three vector components.
        """

        import numpy as np
        from scipy.integrate import ode
        from scipy.integrate import solve_ivp
        from pencil.math.interpolation import vec_int

        if params.interpolation == "tricubic":
            try:
                from scipy.ndimage import map_coordinates
            except (ImportError, ModuleNotFoundError) as e:
                print(
                    f"Warning: Could not import scipy.ndimage.map_coordinates "
                    f"for tricubic interpolation: {e}\n"
                )
                print("Warning: Fall back to trilinear.")
                params.interpolation = "trilinear"

        if not metric:
            metric = lambda xx: np.eye(3)

        dxyz = np.array([params.dx, params.dy, params.dz])
        oxyz = np.array([params.Ox, params.Oy, params.Oz])
        nxyz = np.array([params.nx, params.ny, params.nz])

        # Redefine the derivative y for the scipy ode integrator using the given parameters.
        if (params.interpolation == "mean") or (params.interpolation == "trilinear"):
            odeint_func = lambda t, xx: vec_int(
                xx, field, dxyz, oxyz, nxyz, params.interpolation
            )
        if params.interpolation == "tricubic":
            # For map_coordinates, we pass field data directly (not interpolator objects)
            if splines is None:
                field_data = field
            else:
                # splines here would be the field data array
                field_data = splines
            odeint_func = lambda t, xx: self.tricubic_func(
                xx, field_data, params, map_coordinates
            )
        # Set up the ode solver.
        methods_ode = ["vode", "zvode", "lsoda", "dopri5", "dop853"]
        methods_ivp = ["RK45", "RK23", "Radau", "BDF", "LSODA"]
        if params.method in methods_ode:
            solver = ode(odeint_func, jac=metric)
            solver.set_initial_value(xx, time[0])
            solver.set_integrator(params.method, rtol=params.rtol, atol=params.atol)
            self.tracers = np.zeros([len(time), 3])
            self.tracers[0, :] = xx
            for i, t in enumerate(time[1:]):
                self.tracers[i + 1, :] = solver.integrate(t)
        if params.method in methods_ivp:
            self.tracers = solve_ivp(
                odeint_func,
                (time[0], time[-1]),
                xx,
                t_eval=time,
                rtol=params.rtol,
                atol=params.atol,
                method=params.method,
            ).y.T

        # Remove points that lie outside the domain and interpolation on the boundary.
        cut_mask = (
            (
                (self.tracers[:, 0] > params.Ox + params.Lx)
                + (self.tracers[:, 0] < params.Ox)
            )
            * (not params.periodic_x)
            + (
                (self.tracers[:, 1] > params.Oy + params.Ly)
                + (self.tracers[:, 1] < params.Oy)
            )
            * (not params.periodic_y)
            + (
                (self.tracers[:, 2] > params.Oz + params.Lz)
                + (self.tracers[:, 2] < params.Oz)
            )
            * (not params.periodic_z)
        )
        if np.sum(cut_mask) > 0:
            # Find the first point that lies outside.
            idx_outside = np.min(np.where(cut_mask))
            # Interpolate.
            p0 = self.tracers[idx_outside - 1, :]
            p1 = self.tracers[idx_outside, :]
            lam = np.zeros([6])
            if p0[0] == p1[0]:
                lam[0] = np.inf
                lam[1] = np.inf
            else:
                lam[0] = (params.Ox + params.Lx - p0[0]) / (p1[0] - p0[0])
                lam[1] = (params.Ox - p0[0]) / (p1[0] - p0[0])
            if p0[1] == p1[1]:
                lam[2] = np.inf
                lam[3] = np.inf
            else:
                lam[2] = (params.Oy + params.Ly - p0[1]) / (p1[1] - p0[1])
                lam[3] = (params.Oy - p0[1]) / (p1[1] - p0[1])
            if p0[2] == p1[2]:
                lam[4] = np.inf
                lam[5] = np.inf
            else:
                lam[4] = (params.Oz + params.Lz - p0[2]) / (p1[2] - p0[2])
                lam[5] = (params.Oz - p0[2]) / (p1[2] - p0[2])
            lam_min = np.min(lam[lam >= 0])
            if abs(lam_min) == np.inf:
                lam_min = 0
            self.tracers[idx_outside, :] = p0 + lam_min * (p1 - p0)
            # We don't want to cut the interpolated point (was first point outside).
            cut_mask[idx_outside] = False
            cut_mask[idx_outside + 1 :] = True
            # Remove outside points.
            self.tracers = self.tracers[~cut_mask, :].copy()

        self.params = params
        self.xx = xx
        self.time = time

        # Compute the length of the line segments.
        middle_point_list = (self.tracers[1:, :] + self.tracers[:-1, :]) / 2
        diff_vectors = self.tracers[1:, :] - self.tracers[:-1, :]
        self.section_l = np.zeros(diff_vectors.shape[0])
        for idx, middle_point in enumerate(middle_point_list):
            self.section_l[idx] = np.sqrt(
                np.sum(
                    diff_vectors[idx]
                    * np.sum(metric(middle_point) * diff_vectors[idx], axis=1)
                )
            )
        self.total_l = np.sum(self.section_l)

        self.iterations = len(time)
        self.section_dh = time[1:] - time[:-1]
        self.total_h = time[-1] - time[0]

    def tricubic_func(self, xx, field_data, params, map_coordinates):
        """
        Tricubic spline interpolation using scipy.ndimage.map_coordinates.
        Returns 0 if the point lies outside the box.

        call signature:

        tricubic_func(xx, field_data, params, map_coordinates)

        Keyword arguments:

        *xx*:
          The xyz coordinates of the point to interpolate the data.

        *field_data*:
          The vector field data array with shape [3, nz, ny, nx].

        *params*:
          Parameter object.

        *map_coordinates*:
          The scipy.ndimage.map_coordinates function.
        """

        import numpy as np

        if (
            (xx[0] < params.Ox)
            + (xx[0] > params.Ox + params.Lx)
            + (xx[1] < params.Oy)
            + (xx[1] > params.Oy + params.Ly)
            + (xx[2] < params.Oz)
            + (xx[2] > params.Oz + params.Lz)
        ):
            return np.zeros(3)

        # Convert physical coordinates to array indices
        # Data is stored as [nz, ny, nx], so indices are [z_idx, y_idx, x_idx]
        x_idx = (xx[0] - params.Ox) / params.dx
        y_idx = (xx[1] - params.Oy) / params.dy
        z_idx = (xx[2] - params.Oz) / params.dz

        coords = np.array([[z_idx], [y_idx], [x_idx]])

        return np.array([
            map_coordinates(field_data[0], coords, order=3, mode='constant', cval=0.0)[0],
            map_coordinates(field_data[1], coords, order=3, mode='constant', cval=0.0)[0],
            map_coordinates(field_data[2], coords, order=3, mode='constant', cval=0.0)[0],
        ])
