# laplace_solver.py
"""
This code contains various functions to solve the vector form of the Laplace
equation in various coordinate systems (Cartesian, cylindrical and spherical)
using finite differences.
Additionally, one can find functions to solve the scalar Laplace equation
in Cartesian, cylidrindical and spherical coordinate systems, since these are used
(at least in the case of Cartesian and cylindrical) in the vector solvers.
"""


def laplace_scalar_cartesian(bc, dx, dy, dz, niter=200):
    """
    Solve the scalar Laplace equation in Cartesian coordinates in 3 dimensions
    using finite differences.

    Signature:

    laplace_scalar_cartesian(bc, dx, dy, dz, niter=100)

    Parameters
    ----------
     *bc*: ndarray of shape [nz, ny, nx]
         Boundary conditions on exterior points.
         Keep the inner points 0.

    *dx, dy, dz*: float
        Grid spacing in eaach direction.

    *niter*: int
        Number of iterations.

    Returns
    ----------
    ndarray with the same shape as bc, representing solution to the Laplace equation.
    """

    import numpy as np

    m = 1 / (2 / dx ** 2 + 2 / dy ** 2 + 2 / dz ** 2)

    uu = bc
    iteration = 0
    while iteration < niter:
        iteration += 1
        Au = uu.copy()
        Au = m * (
            1 / dx ** 2 * (np.roll(uu, -1, 2) + np.roll(uu, 1, 2))
            + 1 / dy ** 2 * (np.roll(uu, -1, 1) + np.roll(uu, 1, 1))
            + 1 / dz ** 2 * (np.roll(uu, -1, 0) + np.roll(uu, 1, 0))
        )
        uu[1:-1, 1:-1, 1:-1] = Au[1:-1, 1:-1, 1:-1]

    return uu


def laplace_vector_cartesian(bx, by, bz, dx, dy, dz, niter=200):
    """
    Solve the vector Laplace equation in Cartesian coordinates in 3 dimensions
    using finite differences. This function simply applies the scalar function
    to the three components.

    Signature:

    laplace_vector_cartesian(bx, by, bz, dx, dy, dz, niter=1000)

    Parameters
    ----------
     *bx, by, bz*: ndarrays of shape [nz, ny, nx]
         Boundary conditions on exterior points.
         Keep the inner points 0.

    *dx, dy, dz*: float
        Grid spacing in eaach direction.

    *niter*: int
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [3, nz, ny, nx], representing solution to the Laplace equation.
    """

    import numpy as np

    return np.array(
        [
            laplace_scalar_cartesian(bx, dx, dy, dz, niter=niter),
            laplace_scalar_cartesian(by, dx, dy, dz, niter=niter),
            laplace_scalar_cartesian(bz, dx, dy, dz, niter=niter),
        ]
    )


def laplace_scalar_cylindrical(bc, r, theta, z, niter=200):
    """
    Solve the scalar Laplace equation in cylindical coordinates in 3 dimensions
    using finite differences.

    Signature:

    laplace_scalar_cylindrical(bc, r, theta, z, niter=1000)

    Parameters
    ----------
     *bc*: ndarray of shape [nz, ny, nx]
         Boundary conditions on exterior points.
         Keep the inner points 0.

    r, theta, z: ndarrays of shape [nx], [ny] and [nz]
        1D coordinate arrays of r, theta and z.

    *niter*: int
        Number of iterations.

    Returns
    ----------
    ndarray with the same shape as bc, representing solution to the Laplace equation.
    """

    import numpy as np

    radius_matrix = np.swapaxes(np.meshgrid(theta, r, z, indexing="ij")[1], 0, 2)
    uu = bc
    dr = r[1] - r[0]
    dtheta = theta[1] - theta[0]
    dz = z[1] - z[0]
    m = 1 / (2 * (1 / dr ** 2 + 1 / (dtheta ** 2 * radius_matrix ** 2) + 1 / dz ** 2))

    iteration = 0
    while iteration < niter:
        iteration += 1
        Au = uu.copy()
        Au = m * (
            (1 / dr ** 2 + 1 / (2 * dr * radius_matrix)) * np.roll(uu, -1, 2)
            + (1 / dr ** 2 - 1 / (2 * dr * radius_matrix)) * np.roll(uu, 1, 2)
            + 1
            / (dtheta ** 2 * radius_matrix ** 2)
            * (np.roll(uu, 1, 1) + np.roll(uu, -1, 1))
            + 1 / dz ** 2 * (np.roll(uu, 1, 0) + np.roll(uu, -1, 0))
        )
        uu[1:-1, 1:-1, 1:-1] = Au[1:-1, 1:-1, 1:-1]

    return uu


def laplace_vector_cylindrical(br, btheta, bz, r, theta, z, niter=200):
    """
    Solve the vector Laplace equation in cylindrical coordinates in 3 dimensions
    using finite differences.

    Signature:

    laplace_vector_cylindrical(br, btheta, bz, r, theta, z, niter=200)

    Parameters
    ----------
     *br, btheta, bz*: ndarrays of shape [nz, ny, nx]
         Boundary conditions on exterior points.
         Keep the inner points 0.

    r, theta, z: ndarrays of shape [nx], [ny] and [nz]
        1D coordinate arrays of r, theta and z.

    *niter*: int
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [3, nz, ny, nx], representing solution to the Laplace equation.
    """

    import numpy as np

    radius_matrix = np.swapaxes(np.meshgrid(theta, r, z, indexing="ij")[1], 0, 2)
    dr = r[1] - r[0]
    dtheta = theta[1] - theta[0]
    dz = z[1] - z[0]
    m = 1 / (
        2 / dr ** 2
        + 2 / (radius_matrix ** 2 * dtheta ** 2)
        + 2 / dz ** 2
        + 1 / radius_matrix ** 2
    )

    iteration = 0
    while iteration < niter:
        iteration += 1
        R = br.copy()
        Theta = btheta.copy()
        R = m * (
            (1 / dr ** 2 + 1 / (2 * dr * radius_matrix))
            * (np.roll(R, -1, 2) + np.roll(R, 1, 2))
            + 1
            / (radius_matrix ** 2 * dtheta ** 2)
            * (np.roll(R, -1, 1) + np.roll(R, 1, 1))
            + 1 / dz ** 2 * (np.roll(R, -1, 0) + np.roll(R, 1, 0))
            - 1
            / (dtheta * radius_matrix ** 2)
            * (np.roll(Theta, -1, 1) + np.roll(Theta, 1, 1))
        )
        Theta = m * (
            (1 / dr ** 2 + 1 / (2 * dr * radius_matrix))
            * (np.roll(Theta, -1, 2) + np.roll(Theta, 1, 2))
            + 1
            / (radius_matrix ** 2 * dtheta ** 2)
            * (np.roll(Theta, -1, 1) + np.roll(Theta, 1, 1))
            + 1 / dz ** 2 * (np.roll(Theta, -1, 0) + np.roll(Theta, 1, 0))
            + 1 / (dtheta * radius_matrix ** 2) * (np.roll(R, -1, 1) + np.roll(R, 1, 1))
        )

        br[:, :, 1:-1] = R[:, :, 1:-1]
        btheta[:, 1:-1, :] = Theta[:, 1:-1, :]
        bz = laplace_scalar_cylindrical(bz, r, theta, z, niter=niter)

    return np.array([R, Theta, bz])


def laplace_scalar_spherical(bc, r, theta, phi, niter=200):
    """
    Solve the scalar Laplace equation in spherical coordinates in 3 dimensions
    using finite differences.

    Signature:

    laplace_scalar_spherical(bc, r, theta, phi, niter=200)

    Parameters
    ----------
     *bc*: ndarray of shape [nz, ny, nx]
         Boundary conditions on exterior points.
         Keep the inner points 0.

    r, theta, phi: ndarrays of shape [nx], [ny] and [nz]
        1D coordinate arrays of r, theta and phi.
        The convention is taken that theta ranges from 0 to pi and
        phi ranges from 0 to 2pi.
        Note that singularities arise in the equation when theta = 0 and r = 0.
        This may produce unintended results.

    *niter*: int
        Number of iterations.

    Returns
    ----------
    ndarray with the same shape as bc, representing solution to the Laplace equation.
    """

    import numpy as np

    radius_matrix, theta_matrix, phi_matrix = np.meshgrid(r, theta, phi, indexing="ij")
    radius_matrix = np.swapaxes(radius_matrix, 0, 2)
    theta_matrix = np.swapaxes(theta_matrix, 0, 2)

    uu = bc
    dr = r[1] - r[0]
    dtheta = theta[1] - theta[0]
    dphi = phi[1] - phi[0]
    m = np.nan_to_num(
        1
        / (
            2
            * (
                1 / dr ** 2
                + 1 / (dtheta ** 2 * radius_matrix ** 2)
                + 1 / (dphi ** 2 * np.sin(theta_matrix) ** 2 * radius_matrix ** 2)
            )
        )
    )

    iteration = 0
    while iteration < niter:
        iteration += 1
        Au = uu.copy()
        Au = m * (
            1 / (dr * radius_matrix) * (np.roll(uu, -1, 2) - np.roll(uu, 1, 2))
            + 1 / dr ** 2 * (np.roll(uu, -1, 2) + np.roll(uu, 1, 2))
            + np.nan_to_num(
                1 / (2 * dtheta * r ** 2 * np.nan_to_num(np.tan(theta_matrix)))
            )
            * (np.roll(uu, -1, 1) - np.roll(uu, 1, 1))
            + 1 / (r ** 2 * dtheta ** 2) * (np.roll(uu, -1, 1) + np.roll(uu, 1, 1))
            + np.nan_to_num(1 / (r ** 2 * (np.sin(theta_matrix))) ** 2 * dphi ** 2)
            * (np.roll(uu, 1, 0) + np.roll(uu, -1, 0))
        )
        uu[1:-1, 1:-1, 1:-1] = Au[1:-1, 1:-1, 1:-1]

    return uu


def laplace_vector_spherical(br, btheta, bphi, r, theta, phi, niter=200):
    """
    Solve the scalar Laplace equation in spherical coordinates in 3 dimensions
    using finite differences.

    Signature:

    laplace_vector_spherical(br, btheta, bphi, r, theta, phi, niter=200)

    Parameters
    ----------
     *bc, btheta, bphi*: ndarray of shape [nz, ny, nx]
         Boundary conditions on exterior points.
         Keep the inner points 0.

    r, theta, phi: ndarrays of shape [nx], [ny] and [nz]
        1D coordinate arrays of r, theta and phi.
        The convention is taken that theta ranges from 0 to pi and
        phi ranges from 0 to 2pi.
        Note that singularities arise in the equation when theta = 0 and r = 0.
        This may produce unintended results.

    *niter*: int
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [3, nz, ny, nx], representing solution to the Laplace equation.
    """

    import numpy as np

    radius_matrix, theta_matrix, phi_matrix = np.meshgrid(r, theta, phi, indexing="ij")
    radius_matrix = np.swapaxes(radius_matrix, 0, 2)
    theta_matrix = np.swapaxes(theta_matrix, 0, 2)
    phi_matrix = np.swapaxes(phi_matrix, 0, 2)
    dr = r[1] - r[0]
    dtheta = theta[1] - theta[0]
    dphi = phi[1] - phi[0]

    fraction_r = 1 / (
        2
        * (
            1 / dr ** 2
            + 1 / (dtheta * radius_matrix) ** 2
            + 1 / (dphi * radius_matrix * np.sin(theta_matrix)) ** 2
            + 1 / radius_matrix ** 2
        )
    )
    fraction_theta = 1 / (
        2
        * (
            1 / dr ** 2
            + 1 / (dtheta * radius_matrix) ** 2
            + 1 / (dphi * radius_matrix * np.sin(theta_matrix)) ** 2
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
                / (dr * radius_matrix)
                * (np.roll(coord, -1, 0) - np.roll(coord, 1, 0))
                + 1 / dr ** 2 * (np.roll(coord, -1, 0) + np.roll(coord, 1, 0))
                + 1
                / (dtheta * radius_matrix) ** 2
                * (np.roll(coord, -1, 1) + np.roll(coord, 1, 1))
                + 1
                / (dphi * radius_matrix * np.sin(theta_matrix)) ** 2
                * (np.roll(coord, -1, 2) + np.roll(coord, 1, 2))
                + 1
                / (2 * dtheta * radius_matrix ** 2 * np.tan(theta_matrix))
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
                / (dtheta * radius_matrix ** 2)
                * (np.roll(Theta, -1, 1) - np.roll(Theta, 1, 1))
                - 1
                / (dphi * radius_matrix ** 2 * np.sin(theta_matrix))
                * (np.roll(Phi, -1, 1) - np.roll(Phi, 1, 1))
                - 2 / (radius_matrix ** 2 * np.nan_to_num(np.tan(theta_matrix))) * Theta
            )
        )
        Theta = fraction_theta * np.nan_to_num(
            (
                pre(Theta)
                - 1
                / (
                    np.nan_to_num(np.tan(theta_matrix))
                    * dphi
                    * radius_matrix ** 2
                    * np.sin(theta_matrix)
                )
                * (np.roll(Phi, -1, 2) - np.roll(Phi, 1, 2))
                + 1 / (radius_matrix ** 2 * dr) * (np.roll(R, -1, 1) - np.roll(R, 1, 1))
            )
        )
        Phi = fraction_phi * np.nan_to_num(
            (
                pre(Phi)
                - 1
                / (dphi * radius_matrix ** 2 * np.sin(theta_matrix))
                * (np.roll(R, -1, 2) - np.roll(R, 1, 2))
                + 1
                / (
                    np.nan_to_num(np.tan(theta_matrix))
                    * radius_matrix ** 2
                    * dphi
                    * np.sin(theta_matrix)
                )
                * (np.roll(Theta, -1, 1) - np.roll(Theta, 1, 1))
            )
        )

        br[:, 1:-1, :] = R[:, 1:-1, :]
        btheta[:, 1:-1, :] = Theta[:, 1:-1, :]
        bphi[1:-1, :, :] = Phi[1:-1, :, :]

    return np.array([R, Theta, Phi])
