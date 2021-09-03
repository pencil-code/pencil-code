# laplace_solver.py
'''
This code contains various functions to solve the vector form of the Laplace
equation in various coordinate systems (Cartesian, cylindrical and spherical)
using finite differences. The adaptions to the Poisson equation have also been
added. Additionally, one can find functions to solve the scalar Laplace equation
in Cartesian, cylidrindical and spherical coordinate systems, since these are used
(at least in the case of Cartesian and cylindrical) in the vector solvers.
'''


def laplace_scalar_cartesian(bc, dx, dy, dz, niter=200):
    '''
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
    '''

    import numpy as np

    m = 1/(2/dx**2 + 2/dy**2 + 2/dz**2)

    uu = bc
    iteration = 0
    while iteration < niter:
        iteration += 1
        Au = uu.copy()
        Au = m*(1/dx**2 * (np.roll(uu, -1, 2) + np.roll(uu, 1, 2)) +
                1/dy**2 * (np.roll(uu, -1, 1) + np.roll(uu, 1, 1)) +
                1/dz**2 * (np.roll(uu, -1, 0) + np.roll(uu, 1, 0)))
        uu[1:-1, 1:-1, 1:-1] = Au[1:-1, 1:-1, 1:-1]

    return uu


def laplace_vector_cartesian(bx, by, bz, dx, dy, dz, niter=200):
    '''
    Solve the vector Laplace equation in Cartesian coordinates in 3 dimensions
    using finite differences. This function simply applies the scalar function
    to the three components.

    Signature:

    laplace_vector_cartesian(bx, by, bz, dx, dy, dz, niter=1000)

    Parameters
    ----------
     *bx, by, bz*: ndarray of shape [nz, ny, nx]
         Boundary conditions on exterior points.
         Keep the inner points 0.

    *dx, dy, dz*: float
        Grid spacing in eaach direction.

    *niter*: int
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [3, nz, ny, nx], representing solution to the Laplace equation.
    '''

    import numpy as np

    return np.array([laplace_scalar_cartesian(bx, dx, dy, dz, niter=niter),
                     laplace_scalar_cartesian(by, dx, dy, dz, niter=niter),
                     laplace_scalar_cartesian(bz, dx, dy, dz, niter=niter)])


def laplace_scalar_cylindrical(bc, r, theta, z, niter=200):
    '''
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
    '''

    import numpy as np

    radius_matrix = np.swapaxes(np.meshgrid(theta, r, z, indexing='ij')[1], 0, 2)
    uu = bc
    dr = r[1] - r[0]
    dtheta = theta[1] - theta[0]
    dz = z[1] - z[0]
    m = 1/(2 * (1/dr**2 + 1/(dtheta**2 * radius_matrix**2) + 1/dz**2))

    iteration = 0
    while iteration < niter:
        iteration += 1
        Au = uu.copy()
        Au = m*((1/dr**2 + 1/(2*dr*radius_matrix))*np.roll(uu, -1, 2) +
                (1/dr**2 - 1/(2*dr*radius_matrix))*np.roll(uu, 1, 2) +
                1/(dtheta**2 * radius_matrix**2)*(np.roll(uu, 1, 1) + np.roll(uu, -1, 1)) +
                1/dz**2*(np.roll(uu, 1, 0) + np.roll(uu, -1, 0)))
        uu[1:-1, 1:-1, 1:-1] = Au[1:-1, 1:-1, 1:-1]

    return uu


def laplace_vector_cylindrical(br, btheta, bz, r, theta, z, niter=200):
    '''
    Solve the vector Laplace equation in cylindrical coordinates in 3 dimensions
    using finite differences.

    Signature:

    laplace_vector_cylindrical(br, btheta, bz, r, theta, z, niter=200)

    Parameters
    ----------
     *br, btheta, bz*: ndarray of shape [nz, ny, nx]
         Boundary conditions on exterior points.
         Keep the inner points 0.

    r, theta, z: ndarrays of shape [nx], [ny] and [nz]
        1D coordinate arrays of r, theta and z.

    *niter*: int
        Number of iterations.

    Returns
    ----------
    ndarray with the shape [3, nz, ny, nx], representing solution to the Laplace equation.
    '''

    import numpy as np

    radius_matrix = np.swapaxes(np.meshgrid(theta, r, z, indexing='ij')[1], 0, 2)
    dr = r[1] - r[0]
    dtheta = theta[1] - theta[0]
    dz = z[1] - z[0]
    m = 1/(2/dr**2 + 2/(radius_matrix**2*dtheta**2) + 2/dz**2 + 1/radius_matrix**2)

    iteration = 0
    while iteration < niter:
        iteration += 1
        R = br.copy()
        Theta = btheta.copy()
        R = m*((1/dr**2 + 1/(2*dr*radius_matrix))*(np.roll(R, -1, 2) + np.roll(R, 1, 2)) +
               1/(radius_matrix**2*dtheta**2)*(np.roll(R, -1, 1) + np.roll(R, 1, 1)) +
               1/dz**2*(np.roll(R, -1, 0) + np.roll(R, 1, 0)) -
               1/(dtheta*radius_matrix**2)*(np.roll(Theta, -1, 1) + np.roll(Theta, 1, 1)))
        Theta = m*((1/dr**2 + 1/(2*dr*radius_matrix))*(np.roll(Theta, -1, 2) + np.roll(Theta, 1, 2)) +
                   1/(radius_matrix**2*dtheta**2)*(np.roll(Theta, -1, 1) + np.roll(Theta, 1, 1)) +
                   1/dz**2 * (np.roll(Theta, -1, 0) + np.roll(Theta, 1, 0)) +
                   1/(dtheta*radius_matrix**2)*(np.roll(R, -1, 1) + np.roll(R, 1, 1)))

        br[:, :, 1:-1] = R[:, :, 1:-1]
        btheta[:, 1:-1, :] = Theta[:, 1:-1, :]
        bz = laplace_scalar_cylindrical(bz, r, theta, z, niter=niter)

    return np.array([R, Theta, bz])



def laplace_scalar_spherical(bc, r, theta, phi, niter=200):
    '''
    Solve the vector Laplace equation in spherical coordinates in 3 dimensions
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
    '''

    '''
    x = np.linspace(1, 10, 101)
    y = np.linspace(0.1, 1, 101)
    z = np.linspace(0, 2*np.pi, 101)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    bx = np.zeros_like(xx)
    by = np.zeros_like(xx)
    bz = np.exp(-(xx-2)**2-81*(yy-0.5)**2) - np.exp(-(xx-8)**2-81*(yy-0.5)**2)
    bz = np.swapaxes(bz, 0, 2)
    bz[1:, :, :] = 0
    solution = pc.math.laplace_scalar_spherical(bz, x, y, z, niter=100)
    '''

    import numpy as np

    radius_matrix, theta_matrix, phi_matrix = np.meshgrid(r, theta, phi, indexing='ij')
    radius_matrix = np.swapaxes(radius_matrix, 0, 2)
    theta_matrix = np.swapaxes(theta_matrix, 0, 2)

    uu = bc
    dr = r[1] - r[0]
    dtheta = theta[1] - theta[0]
    dphi = phi[1] - phi[0]
    m = np.nan_to_num(1/(2*(1/dr**2 + 1/(dtheta**2*radius_matrix**2) +
                            1/(dphi**2*np.sin(theta_matrix)**2*radius_matrix**2))))

    iteration = 0
    while iteration < niter:
        iteration += 1
        Au = uu.copy()
        Au = m*(1/(dr*radius_matrix)*(np.roll(uu, -1, 2) - np.roll(uu, 1, 2)) +
                1/dr**2*(np.roll(uu, -1, 2) + np.roll(uu, 1, 2)) +
                np.nan_to_num(1/(2*dtheta*r**2*np.nan_to_num(np.tan(theta_matrix))))*(np.roll(uu, -1, 1) - np.roll(uu, 1, 1)) +
                1/(r**2*dtheta**2)*(np.roll(uu, -1, 1) + np.roll(uu, 1, 1)) +
                np.nan_to_num(1/(r**2*(np.sin(theta_matrix)))**2*dphi**2)*(np.roll(uu, 1, 0) + np.roll(uu, -1, 0)))
        uu[1:-1, 1:-1, 1:-1] = Au[1:-1, 1:-1, 1:-1]

    return uu


#def laplace_spherical_vector(br, btheta, bphi, r, theta, phi):
#    '''
#    Solve the vector form of the Laplace equation in spherical coordinates using finite differences.
#    
#    Call signature:
#    
#    br, btheta, bphi: 3D ndarrays containing boundary conditions on exterior points. If the interior points are non-zero, the
#              solution may not be correct. 0th index corresponds to smallest radis; theta increases as index increases.
#              
#    r, theta, phi: 1D coordinate arrays of r, theta and phi where the convention is taken that theta ranges from 0 to pi, phi
#                   ranges from 0 to 2pi. Note that singularities arise in the equation when theta = 0. This may produce 
#                   unintended results. 
#                   
#    Returns 3 3D ndarrays representing the components of the solution to the vector Laplace equation in spherical coordinates.
#    '''
#    
#    import numpy as np
#    
#    theta_matrix, radius_matrix, phi_matrix = np.meshgrid(theta, r, phi)
#    h = r[1] - r[0]
#    h_prime = theta[1] - theta[0]
#    h_doubleprime = phi[1] - phi[0]
#    
#    fraction_r = 1/( 
#                    2 * (1/h**2 + 1/(h_prime * radius_matrix)**2
#                          + 1/(h_doubleprime * radius_matrix * np.sin(theta_matrix))**2 + 1/radius_matrix**2)
#                   )
#    
#    fraction_theta = 1/(
#                        2 * (1/h**2 + 1/(h_prime * radius_matrix)**2 + 1/(h_doubleprime * radius_matrix * np.sin(theta_matrix))**2
#                        + 1/(radius_matrix * np.sin(theta_matrix))**2)
#                        )
#    
#    fraction_phi = fraction_theta
#        
#    def pre(C):
#        '''
#        This function simply returns the first five terms in each of the differential equations to increase readability.
#        The parameter C is simply the most recent iteration of the unknown.
#        '''
#        return np.nan_to_num((
#                1/(h * radius_matrix) * (np.roll(C, -1, 0) - np.roll(C, 1, 0)) 
#                + 1/h**2 * (np.roll(C, -1, 0) + np.roll(C, 1, 0))
#                + 1/(h_prime * radius_matrix)**2 * (np.roll(C, -1, 1) + np.roll(C, 1, 1))
#                + 1/(h_doubleprime * radius_matrix * np.sin(theta_matrix))**2 * (np.roll(C, -1, 2) + np.roll(C, 1, 2))
#                + 1/(2 * h_prime * radius_matrix**2 * np.tan(theta_matrix)) * (np.roll(C, -1, 1) + np.roll(C, 1, 1))
#                ))
#    
#    
#    
#    iteration = 0
#    while iteration <= 200:
#        iteration += 1
#        R = br.copy()
#        Theta = btheta.copy()
#        Phi = bphi.copy()
#    
#        R = fraction_r * np.nan_to_num((
#                                        pre(R)
#                                        - 1/(h_prime * radius_matrix**2) * (np.roll(Theta, -1, 1) - np.roll(Theta, 1, 1))
#                                        - 1/(h_doubleprime * radius_matrix**2 * np.sin(theta_matrix)) 
#                                        * (np.roll(Phi, -1, 1) - np.roll(Phi, 1, 1))
#                                        - 2/(radius_matrix**2 * np.nan_to_num(np.tan(theta_matrix))) * Theta
#                                        ))
#                        
#        Theta = fraction_theta * np.nan_to_num ((
#                                                  pre(Theta)
#                                                  - 1/(np.nan_to_num(np.tan(theta_matrix)) * h_doubleprime * radius_matrix**2 * np.sin(theta_matrix)) * (np.roll(Phi, -1, 2) - np.roll(Phi, 1, 2))
#                                                  + 1/(radius_matrix**2 * h_prime) * (np.roll(R, -1, 1) - np.roll(R, 1, 1))
#                                 ))
#        
#        Phi = fraction_phi * np.nan_to_num((
#                                              pre(Phi)
#                                              - 1/(h_doubleprime * radius_matrix**2 * np.sin(theta_matrix)) * (np.roll(R,-1,2)-np.roll(R,1,2))
#                                              + 1/(np.nan_to_num(np.tan(theta_matrix)) * radius_matrix**2 * h_doubleprime * np.sin(theta_matrix)) * (np.roll(Theta, -1, 1) - np.roll(Theta, 1, 1))
#                             ))
#    
#        br[1:-1, :, :] = R[1:-1, :, :]
#        btheta[:, 1:-1, :] = Theta[:, 1:-1, :]
#        bphi[:, 1:-1, :] = Phi[:, 1:-1, :]
#        
#    return R, Theta, Phi
#
#
#
#
#def poisson_cartesian_vector(bx, by, bz, x, y, z, hx, hy, hz):
#    '''
#    Solves the (vector) form of the Poisson equation in 3D Cartesian coordinates, $\nabla^2 u = h$, using finite differences.
#    
#    Call signature:
#    
#    poisson_cartesian_vector(bx, by, bz, x, y, z, hx, hy, hz)
#    
#    Arguments:
#    bx, by, bz: 3D ndarrays containing boundary conditions on exterior points for each component.
#                If the interior points are non-zero, the solution may not be correct.
#    x, y, z:    1D coordinate arrays. These are only used to calculate grid spacing.
#    hx, hy, hz: 3D ndarrays representing the components of the known function h. 
#    
#    Returns the 3 components of the solution to the vector Poisson equation in Cartesian coordinates.
#    '''
#     
#    import numpy as np
#   
#    h = x[1] - x[0]
#    h_prime = y[1] - y[0]
#    h_doubleprime = z[1] - z[0]
#    m = 1/(2 * 1/(h**2) + 2 * 1/(h_prime**2) + 2 * 1/(h_doubleprime**2))
#    
#    ux = bx
#    uy = by
#    uz = bz
#    iteration = 0
#    while iteration < 1000:
#        iteration += 1
#        Aux = ux.copy()
#        Aux = m * (
#                    1/h**2 * (np.roll(ux, -1, 0) + np.roll(ux, 1, 0)) 
#                    + 1/h_prime**2 * (np.roll(ux, -1, 1) + np.roll(ux, 1, 1))
#                    + 1/h_doubleprime**2 * (np.roll(ux, -1, 2) + np.roll(ux, 1, 2))
#                    - hx
#                    )
#        ux[1:-1, 1:-1, 1:-1] = Aux[1:-1,1:-1, 1:-1]
#        
#        Auy = uy.copy()
#        Auy = m * (
#                    1/h**2 * (np.roll(uy, -1, 0) + np.roll(uy, 1, 0)) 
#                    + 1/h_prime**2 * (np.roll(uy, -1, 1) + np.roll(uy, 1, 1))
#                    + 1/h_doubleprime**2 * (np.roll(uy, -1, 2) + np.roll(uy, 1, 2))
#                    -hy
#                    )
#        uy[1:-1, 1:-1, 1:-1] = Auy[1:-1,1:-1, 1:-1]
#        
#        Auz = uz.copy()
#        Auz = m * (
#                    1/h**2 * (np.roll(uz, -1, 0) + np.roll(uz, 1, 0)) 
#                    + 1/h_prime**2 * (np.roll(uz, -1, 1) + np.roll(uz, 1, 1))
#                    + 1/h_doubleprime**2 * (np.roll(uz, -1, 2) + np.roll(uz, 1, 2))
#                    -hz
#                    )
#        uz[1:-1, 1:-1, 1:-1] = Auz[1:-1,1:-1, 1:-1]
#        
#        
#
#    return ux, uy, uz
#
#
#
#def poisson_cylindrical_scalar(bc, r, theta, z, h):
#    '''
#    Solve the scalar form of the Poisson equation in cylindical coordinates using finite differences. 
#    
#    Call signature:
#    
#    poisson_cylindrical_scalar(bc, r, theta, z, h)
#    
#    Arguments:
#    bc: 3D ndarray containing boundary conditions on exterior points. If the interior points are non-zero, the
#        solution may not be correct. 0th index corresponds to minimum radius; theta increases from left to right.
#        
#    r, theta, z: 1D coordinate arrays of r, theta and z.
#    
#    h: ndarray representing known function h.
#    
#    Returns the 3 components of the solution to the scalar Poisson equation in cylindrical coordinates.
#    '''
#    
#    import numpy as np
#
#    theta_matrix, radius_matrix, z_matrix = np.meshgrid(theta, r, z)
#    u = bc
#    h = r[1] - r[0]
#    h_prime = theta[1] - theta[0]
#    h_doubleprime = z[1] - z[0]
#    m = 1/(2 * (1/h**2 + 1/(h_prime**2 * radius_matrix**2) + 1/h_doubleprime**2))
#    
#    iteration = 0
#    while iteration < 200:
#        iteration += 1
#        Au = u.copy()
#        Au =  m * (
#                    (1/h**2 + 1/(2 * h * radius_matrix)) * np.roll(u, -1, 0) 
#                    + (1/h**2 - 1/(2 * h * radius_matrix)) * np.roll(u, 1, 0) 
#                    + 1/(h_prime**2 * radius_matrix**2) * (np.roll(u, 1, 1) + np.roll(u, -1, 1))
#                    + 1/(h_doubleprime**2) * (np.roll(u, 1, 2) + np.roll(u, -1, 2))
#                    -h
#                    )
#        
#        u[1:L-1, :, :] = Au[1:L-1, :, :]
#        
#    return u
#
#def poisson_cylindrical_vector(br, btheta, bz, r, theta, z, hr, htheta, hz):
#    '''
#    Solve the vector form of the Poisson equation, $\nabla^2 u = h$, in cylindrical coordinates using finite diffferences.
#    
#    Call signature:
#    
#    poisson_cylindrical_vector(br, btheta, bz, r, theta, z, hr, htheta, hz)
#    
#    Arguments:
#    br, btheta, bz:  3D ndarrays containing boundary conditions on exterior points for each component. 
#                     If the interior points are non-zero, the solution may not be correct. See PDF for
#                     derivation.
#                       
#    r, theta, z: 1D coordinate arrays. 
#    
#    hr, htheta, hz: 3d ndarrays representing the components of the known function h.
#    
#    Returns the 3 components of the solution to the vector Poisson equation in cylindrical coordinates.
#    
#    '''
#    
#    import numpy as np
#
#    theta_matrix, radius_matrix, z_matrix = np.meshgrid(theta, r,z)
#    h = r[1] - r[0]
#    h_prime = theta[1] - theta[0]
#    h_doubleprime = z[1] - z[0]
#    m = 1/(2/h**2 + 2/(radius_matrix**2 * h_prime**2) + 2/h_doubleprime**2 + 1/radius_matrix**2)
#    
#    iteration = 0
#    while iteration < 1000:
#        iteration += 1
#        R = br.copy()
#        Theta = btheta.copy()
#        R = m * (
#                  (1/h**2 + 1/(2 * h * radius_matrix)) * (np.roll(R, -1, 0) + np.roll(R, 1, 0))
#                  + 1/(radius_matrix**2 * h_prime**2) * (np.roll(R, -1, 1) + np.roll(R, 1, 1))
#                  + 1/h_doubleprime**2 * (np.roll(R, -1, 2) + np.roll(R, 1, 2))
#                  - 1/(h_prime * radius_matrix**2) * (np.roll(Theta, -1, 1) + np.roll(Theta, 1, 1))
#                  -hr
#                )
#        Theta = m * (
#                      (1/h**2 + 1/(2 * h * radius_matrix)) * (np.roll(Theta, -1, 0) + np.roll(Theta, 1, 0))
#                      + 1/(radius_matrix**2 * h_prime**2) * (np.roll(Theta, -1, 1) + np.roll(Theta, 1, 1))
#                      + 1/h_doubleprime**2 * (np.roll(Theta, -1, 2) + np.roll(Theta, 1, 2))
#                      + 1/(h_prime * radius_matrix**2) * (np.roll(R, -1, 1) + np.roll(R, 1, 1))
#                      -htheta
#                    )
#        
#        br[1:-1,:,:] = R[1:-1,:,:]
#        btheta[:,1:-1,:] = Theta[:,1:-1,:]
#        
#    return R, Theta, poisson_cylindrical_scalar(bz, r, theta, z, hz)
#
#
#
#def poisson_spherical_vector(br, btheta, bphi, r, theta, phi, hr, htheta, hphi):
#    '''
#    Solve the vector form of the Poisson equation, $\nabla^2 u = h$, in spherical coordinates using finite differences.
#    
#    Call signature:
#    
#    poisson_spherical_vector(br, btheta, bphi, r, theta, phi, hr, htheta, hphi)
#    
#    Arguments:
#    
#    br, btheta, bphi: 3D ndarrays containing boundary conditions on exterior points. If the interior points are non-zero, the
#              solution may not be correct. 0th index corresponds to smallest radis; theta increases as index increases.
#              
#    r, theta, phi: 1D coordinate arrays of r, theta and phi where the convention is taken that theta ranges from 0 to pi, phi
#                   ranges from 0 to 2pi. Note that singularities arise in the equation when theta = 0. This may produce 
#                   unintended results. 
#    hr, htheta, hphi: ndarrays in the same shape as br, btheta, bphi which represent the components of the known function h. 
#    
#    Returns the 3 components of the solution to the vector Poisson equation in spherical coordinates.
#    '''
#        
#    import numpy as np
#
#    theta_matrix, radius_matrix, phi_matrix = np.meshgrid(theta, r, phi)
#    h = r[1] - r[0]
#    h_prime = theta[1] - theta[0]
#    h_doubleprime = phi[1] - phi[0]
#    
#    fraction_r = 1/( 
#                    2 * (1/h**2 + 1/(h_prime * radius_matrix)**2
#                          + 1/(h_doubleprime * radius_matrix * np.sin(theta_matrix))**2 + 1/radius_matrix**2)
#                   )
#    
#    fraction_theta = 1/(
#                        2 * (1/h**2 + 1/(h_prime * radius_matrix)**2 + 1/(h_doubleprime * radius_matrix * np.sin(theta_matrix))**2
#                        + 1/(radius_matrix * np.sin(theta_matrix))**2)
#                        )
#    
#    fraction_phi = fraction_theta
#        
#    def pre(C):
#        '''
#        This function simply returns the first five terms in each of the differential equations to increase readability.
#        The parameter C is simply the most recent iteration of the unknown.
#        '''
#        return np.nan_to_num((
#                1/(h * radius_matrix) * (np.roll(C, -1, 0) - np.roll(C, 1, 0)) 
#                + 1/h**2 * (np.roll(C, -1, 0) + np.roll(C, 1, 0))
#                + 1/(h_prime * radius_matrix)**2 * (np.roll(C, -1, 1) + np.roll(C, 1, 1))
#                + 1/(h_doubleprime * radius_matrix * np.sin(theta_matrix))**2 * (np.roll(C, -1, 2) + np.roll(C, 1, 2))
#                + 1/(2 * h_prime * radius_matrix**2 * np.tan(theta_matrix)) * (np.roll(C, -1, 1) + np.roll(C, 1, 1))
#                ))
#    
#    
#    
#    iteration = 0
#    while iteration <= 200:
#        iteration += 1
#        R = br.copy()
#        Theta = btheta.copy()
#        Phi = bphi.copy()
#    
#        R = fraction_r * np.nan_to_num((
#                                        pre(R)
#                                        - 1/(h_prime * radius_matrix**2) * (np.roll(Theta, -1, 1) - np.roll(Theta, 1, 1))
#                                        - 1/(h_doubleprime * radius_matrix**2 * np.sin(theta_matrix)) 
#                                        * (np.roll(Phi, -1, 1) - np.roll(Phi, 1, 1))
#                                        - 2/(radius_matrix**2 * np.nan_to_num(np.tan(theta_matrix))) * Theta
#                                        - hr
#                                        ))
#                        
#        Theta = fraction_theta * np.nan_to_num ((
#                                                  pre(Theta)
#                                                  - 1/(np.nan_to_num(np.tan(theta_matrix)) * h_doubleprime * radius_matrix**2 * np.sin(theta_matrix)) * (np.roll(Phi, -1, 2) - np.roll(Phi, 1, 2))
#                                                  + 1/(radius_matrix**2 * h_prime) * (np.roll(R, -1, 1) - np.roll(R, 1, 1))
#                                                  - htheta
#                                 ))
#        
#        Phi = fraction_phi * np.nan_to_num((
#                                              pre(Phi)
#                                              - 1/(h_doubleprime * radius_matrix**2 * np.sin(theta_matrix)) * (np.roll(R,-1,2)-np.roll(R,1,2))
#                                              + 1/(np.nan_to_num(np.tan(theta_matrix)) * radius_matrix**2 * h_doubleprime * np.sin(theta_matrix)) * (np.roll(Theta, -1, 1) - np.roll(Theta, 1, 1))
#                                              - hphi
#                             ))
#    
#        br[1:-1, :, :] = R[1:-1, :, :]
#        btheta[:, 1:-1, :] = Theta[:, 1:-1, :]
#        bphi[:, 1:-1, :] = Phi[:, 1:-1, :]
#        
#    return R, Theta, Phi
#
