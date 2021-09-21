# poisson_solver.py
"""
This code contains various functions to solve the vector form of the Poisson
equation in various coordinate systems (Cartesian, cylindrical and spherical)
using finite differences.
"""


def poisson_vector_cartesian(bx, by, bz, x, y, z, hx, hy, hz, niter=1000):
    """
    Solve the vector form of the Poisson equation in 3D Cartesian coordinates,
    $\nabla^2 u = h$, using finite differences.

    Signature:

    poisson_vector_cartesian(bx, by, bz, x, y, z, hx, hy, hz, niter=1000)

    Parameters
    ----------
    *bx, by, bz*: ndarray of shape [nz, ny, nx]
        Boundary conditions on exterior points for each component.
        Keep the inner points 0.

    *x, y, z*: ndarrays of shape [nx], [ny] and [nz]
        Coordinate arrays to calculate grid spacing.

    *hx, hy, hz*: ndarray of shape [nz, ny, nx]
        Representing the components of the known function h.

    *niter*:
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [3, nz, ny, nx], representing solution to the Poisson equation.
    """

    import numpy as np

    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    m = 1 / (2 / dx ** 2 + 2 / dy ** 2 + 2 / dz ** 2)

    ux = bx
    uy = by
    uz = bz
    iteration = 0
    while iteration < niter:
        iteration += 1
        Aux = ux.copy()
        Aux = m * (
            1 / dx ** 2 * (np.roll(ux, -1, 2) + np.roll(ux, 1, 2))
            + 1 / dy ** 2 * (np.roll(ux, -1, 1) + np.roll(ux, 1, 1))
            + 1 / dz ** 2 * (np.roll(ux, -1, 0) + np.roll(ux, 1, 0))
            - hx
        )
        ux[1:-1, 1:-1, 1:-1] = Aux[1:-1, 1:-1, 1:-1]

        Auy = uy.copy()
        Auy = m * (
            1 / dx ** 2 * (np.roll(uy, -1, 2) + np.roll(uy, 1, 2))
            + 1 / dy ** 2 * (np.roll(uy, -1, 1) + np.roll(uy, 1, 1))
            + 1 / dz ** 2 * (np.roll(uy, -1, 0) + np.roll(uy, 1, 0))
            - hy
        )
        uy[1:-1, 1:-1, 1:-1] = Auy[1:-1, 1:-1, 1:-1]

        Auz = uz.copy()
        Auz = m * (
            1 / dx ** 2 * (np.roll(uz, -1, 2) + np.roll(uz, 1, 2))
            + 1 / dy ** 2 * (np.roll(uz, -1, 1) + np.roll(uz, 1, 1))
            + 1 / dz ** 2 * (np.roll(uz, -1, 0) + np.roll(uz, 1, 0))
            - hz
        )
        uz[1:-1, 1:-1, 1:-1] = Auz[1:-1, 1:-1, 1:-1]

    return np.array([ux, uy, uz])


def poisson_scalar_cylindrical(bc, r, theta, z, h, niter=200):
    """
    Solve the scalar form of the Poisson equation in cylindical coordinates using finite differences.

    Signature:

    poisson_scalar_cylindrical(bc, r, theta, z, h, niter=200)

    Parameters
    ----------
    *bc*: ndarray of shape [nz, ny, nx]
        Boundary conditions on exterior points.
        Keep the inner points 0.

    *r, theta, z*: ndarrays of shape [nx], [ny] and [nz]
        Coordinate arrays to calculate grid spacing.

    *h*: ndarray of shape [nz, ny, nx]
        Representing the known function h.

    *niter*:
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [nz, ny, nx], representing solution to the Poisson equation.
    """

    import numpy as np

    radius_matrix = np.swapaxes(np.meshgrid(theta, r, z, indexing="ij")[1], 0, 2)
    u = bc
    dx = r[1] - r[0]
    dy = theta[1] - theta[0]
    dz = z[1] - z[0]
    m = 1 / (2 * (1 / dx ** 2 + 1 / (dy ** 2 * radius_matrix ** 2) + 1 / dz ** 2))

    iteration = 0
    while iteration < niter:
        iteration += 1
        Au = u.copy()
        Au = m * (
            (1 / dx ** 2 + 1 / (2 * dx * radius_matrix)) * np.roll(u, -1, 2)
            + (1 / dx ** 2 - 1 / (2 * dx * radius_matrix)) * np.roll(u, 1, 2)
            + 1
            / (dy ** 2 * radius_matrix ** 2)
            * (np.roll(u, 1, 1) + np.roll(u, -1, 1))
            + 1 / (dz ** 2) * (np.roll(u, 1, 0) + np.roll(u, -1, 0))
            - h
        )
        u[:, :, 1:-1] = Au[:, :, 1:-1]

    return u


def poisson_vector_cylindrical(br, btheta, bz, r, theta, z, hr, htheta, hz, niter=1000):
    """
    Solve the vector form of the Poisson equation, $\nabla^2 u = h$, in cylindrical
    coordinates using finite diffferences.

    Signature:

    poisson_vector_cylindrical(br, btheta, bz, r, theta, z, hr, htheta, hz, niter=1000)

    Parameters*
    ----------
    *br, btheta, bz*: ndarray of shape [nz, ny, nx]
        Boundary conditions on exterior points for each component.
        Keep the inner points 0.

    *r, theta, z*: ndarrays of shape [nx], [ny] and [nz]
        Coordinate arrays to calculate grid spacing.

    *hr, htheta, hz*: ndarray of shape [nz, ny, nx]
        Representing the components of the known function h.

    *niter*:
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [3, nz, ny, nx], representing solution to the Poisson equation.
    """

    import numpy as np

    theta_matrix, radius_matrix, z_matrix = np.meshgrid(theta, r, z, indexing="ij")
    theta_matrix = np.swapaxes(theta_matrix, 0, 2)
    radius_matrix = np.swapaxes(radius_matrix, 0, 2)
    z_matrix = np.swapaxes(z_matrix, 0, 2)
    dx = r[1] - r[0]
    dy = theta[1] - theta[0]
    dz = z[1] - z[0]
    m = 1 / (
        2 / dx ** 2
        + 2 / (radius_matrix ** 2 * dy ** 2)
        + 2 / dz ** 2
        + 1 / radius_matrix ** 2
    )

    iteration = 0
    while iteration < niter:
        iteration += 1
        R = br.copy()
        Theta = btheta.copy()
        R = m * (
            (1 / dx ** 2 + 1 / (2 * dx * radius_matrix))
            * (np.roll(R, -1, 2) + np.roll(R, 1, 2))
            + 1
            / (radius_matrix ** 2 * dy ** 2)
            * (np.roll(R, -1, 1) + np.roll(R, 1, 1))
            + 1 / dz ** 2 * (np.roll(R, -1, 0) + np.roll(R, 1, 0))
            - 1
            / (dy * radius_matrix ** 2)
            * (np.roll(Theta, -1, 1) + np.roll(Theta, 1, 1))
            - hr
        )
        Theta = m * (
            (1 / dx ** 2 + 1 / (2 * dx * radius_matrix))
            * (np.roll(Theta, -1, 2) + np.roll(Theta, 1, 2))
            + 1
            / (radius_matrix ** 2 * dy ** 2)
            * (np.roll(Theta, -1, 1) + np.roll(Theta, 1, 1))
            + 1 / dz ** 2 * (np.roll(Theta, -1, 0) + np.roll(Theta, 1, 0))
            + 1 / (dy * radius_matrix ** 2) * (np.roll(R, -1, 1) + np.roll(R, 1, 1))
            - htheta
        )

        br[:, :, 1:-1] = R[:, :, 1:-1]
        btheta[:, 1:-1, :] = Theta[:, 1:-1, :]

    return np.array([R, Theta, poisson_scalar_cylindrical(bz, r, theta, z, hz)])


def poisson_vector_spherical(
    br, btheta, bphi, r, theta, phi, hr, htheta, hphi, niter=200
):
    """
    Solve the vector form of the Poisson equation, $\nabla^2 u = h$,
    in spherical coordinates using finite differences.

    Signature:

    poisson_vector_spherical(br, btheta, bphi, r, theta, phi, hr, htheta, hphi, niter=200)

    Parameters*
    ----------
    *br, btheta, bphi*:  ndarray of shape [nz, ny, nx]
        Boundary conditions on exterior points for each component.
        Keep the inner points 0.

    *r, theta, phi*:  ndarrays of shape [nx], [ny] and [nz]
        Coordinate arrays, where the convention is taken that theta ranges from 0 to pi,
        phi ranges from 0 to 2pi.
        Note that singularities arise in the equation when theta = 0.
        This may produce unintended results.

    *hr, htheta, hphi*: ndarray of shape [nz, ny, nx]
        Representing the components of the known function h.

    *niter*:
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [3, nz, ny, nx], representing solution to the Poisson equation.
    """

    import numpy as np

    theta_matrix, radius_matrix, phi_matrix = np.meshgrid(theta, r, phi, indexing="ij")
    theta_matrix = np.swapaxes(theta_matrix, 0, 2)
    radius_matrix = np.swapaxes(radius_matrix, 0, 2)
    phi_matrix = np.swapaxes(phi_matrix, 0, 2)
    dx = r[1] - r[0]
    dy = theta[1] - theta[0]
    dz = phi[1] - phi[0]

    fraction_r = 1 / (
        2
        * (
            1 / dx ** 2
            + 1 / (dy * radius_matrix) ** 2
            + 1 / (dz * radius_matrix * np.sin(theta_matrix)) ** 2
            + 1 / radius_matrix ** 2
        )
    )
    fraction_theta = 1 / (
        2
        * (
            1 / dx ** 2
            + 1 / (dy * radius_matrix) ** 2
            + 1 / (dz * radius_matrix * np.sin(theta_matrix)) ** 2
            + 1 / (radius_matrix * np.sin(theta_matrix)) ** 2
        )
    )
    fraction_phi = fraction_theta

    def pre(coord):
        """
        This function simply returns the first five terms in each of the differential equations.
        The parameter coord is simply the most recent iteration of the unknown.
        """

        return np.nan_to_num(
            (
                1
                / (dx * radius_matrix)
                * (np.roll(coord, -1, 2) - np.roll(coord, 1, 2))
                + 1 / dx ** 2 * (np.roll(coord, -1, 2) + np.roll(coord, 1, 2))
                + 1
                / (dy * radius_matrix) ** 2
                * (np.roll(coord, -1, 1) + np.roll(coord, 1, 1))
                + 1
                / (dz * radius_matrix * np.sin(theta_matrix)) ** 2
                * (np.roll(coord, -1, 0) + np.roll(coord, 1, 0))
                + 1
                / (2 * dy * radius_matrix ** 2 * np.tan(theta_matrix))
                * (np.roll(coord, -1, 1) + np.roll(coord, 1, 1))
            )
        )

    iteration = 0
    while iteration <= niter:
        iteration += 1
        R = br.copy()
        Theta = btheta.copy()
        Phi = bphi.copy()

        R = fraction_r * np.nan_to_num(
            (
                pre(R)
                - 1
                / (dy * radius_matrix ** 2)
                * (np.roll(Theta, -1, 1) - np.roll(Theta, 1, 1))
                - 1
                / (dz * radius_matrix ** 2 * np.sin(theta_matrix))
                * (np.roll(Phi, -1, 1) - np.roll(Phi, 1, 1))
                - 2 / (radius_matrix ** 2 * np.nan_to_num(np.tan(theta_matrix))) * Theta
                - hr
            )
        )
        Theta = fraction_theta * np.nan_to_num(
            (
                pre(Theta)
                - 1
                / (
                    np.nan_to_num(np.tan(theta_matrix))
                    * dz
                    * radius_matrix ** 2
                    * np.sin(theta_matrix)
                )
                * (np.roll(Phi, -1, 0) - np.roll(Phi, 1, 0))
                + 1 / (radius_matrix ** 2 * dy) * (np.roll(R, -1, 1) - np.roll(R, 1, 1))
                - htheta
            )
        )
        Phi = fraction_phi * np.nan_to_num(
            (
                pre(Phi)
                - 1
                / (dz * radius_matrix ** 2 * np.sin(theta_matrix))
                * (np.roll(R, -1, 0) - np.roll(R, 1, 0))
                + 1
                / (
                    np.nan_to_num(np.tan(theta_matrix))
                    * radius_matrix ** 2
                    * dz
                    * np.sin(theta_matrix)
                )
                * (np.roll(Theta, -1, 1) - np.roll(Theta, 1, 1))
                - hphi
            )
        )

        br[:, :, 1:-1, :, :] = R[:, :, 1:-1]
        btheta[:, 1:-1, :] = Theta[:, 1:-1, :]
        bphi[1:-1, :, :] = Phi[1:-1, :, :]

    return np.array([R, Theta, Phi])
