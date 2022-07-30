#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test the functions that calculate derivatives.
"""

import numpy as np
import pencil as pc
from test_utils import (
    test,
    _assert_close_arr,
)
from numpy import exp, sin, cos

pcd = pc.math.derivatives


@test
def partial_derivatives() -> None:
    """Partial derivatives of scalar fields"""
    z, y, x = generate_xyz_grid()
    dx = x[1, 1, 1] - x[0, 0, 0]
    dy = y[1, 1, 1] - y[0, 0, 0]
    dz = z[1, 1, 1] - z[0, 0, 0]

    f = sin(x)
    df_ana = cos(x)
    df_num = pcd.xder(f, dx)
    check_arr_close(df_ana, df_num)

    f = sin(2 * x) * cos(y) * exp(z)
    df_ana = 2 * cos(2 * x) * cos(y) * exp(z)
    df_num = pcd.xder(f, dx)
    check_arr_close(df_ana, df_num)

    f = sin(2 * x) * cos(y) * exp(z)
    df_ana = -sin(2 * x) * sin(y) * exp(z)
    df_num = pcd.yder(f, dy)
    check_arr_close(df_ana, df_num)

    f = sin(2 * x) * cos(y) * exp(z)
    df_ana = sin(2 * x) * cos(y) * exp(z)
    df_num = pcd.zder(f, dz)
    check_arr_close(df_ana, df_num)

    f = sin(2 * x) * cos(y) * exp(2 * z)
    df_ana = 2 * sin(2 * x) * cos(y) * exp(2 * z)
    df_num = pcd.zder(f, dz)
    check_arr_close(df_ana, df_num)


@test
def derivatives_lap_grad() -> None:
    """Laplacian and gradient of scalar fields"""
    z, y, x = generate_xyz_grid()
    dx = x[1, 1, 1] - x[0, 0, 0]
    dy = y[1, 1, 1] - y[0, 0, 0]
    dz = z[1, 1, 1] - z[0, 0, 0]

    f = sin(x) * exp(y) * cos(z)
    grad_f_x = exp(y) * cos(x) * cos(z)
    grad_f_y = exp(y) * sin(x) * cos(z)
    grad_f_z = -exp(y) * sin(x) * sin(z)
    grad_f = np.stack([grad_f_x, grad_f_y, grad_f_z], axis=0)
    Lap_f = -exp(y) * sin(x) * cos(z)

    check_arr_close(grad_f, pcd.grad(f, dx, dy, dz))
    check_arr_close(Lap_f, pcd.del2(f, dx, dy, dz))


@test
def derivatives_vector() -> None:
    """Derivatives of vector fields"""
    z, y, x = generate_xyz_grid()
    dx = x[1, 1, 1] - x[0, 0, 0]
    dy = y[1, 1, 1] - y[0, 0, 0]
    dz = z[1, 1, 1] - z[0, 0, 0]

    # Build the vector field to be tested
    vx = sin(x) * exp(y) * cos(z)
    vy = sin(x + y + z)
    vz = exp(x) * cos(y + z)
    v = np.stack([vx, vy, vz], axis=0)

    div_v = exp(y) * cos(x) * cos(z) + cos(x + y + z) - exp(x) * sin(y + z)
    curl_v_x = -exp(x) * sin(y + z) - cos(x + y + z)
    curl_v_y = -exp(x) * cos(y + z) - exp(y) * sin(x) * sin(z)
    curl_v_z = -exp(y) * sin(x) * cos(z) + cos(x + y + z)
    curl_v = np.stack([curl_v_x, curl_v_y, curl_v_z], axis=0)

    check_arr_close(div_v, pcd.div(v, dx, dy, dz))
    check_arr_close(curl_v, pcd.curl(v, dx, dy, dz))


@test
def derivatives_cylindrical() -> None:
    """Derivatives in cylindrical coordinates"""
    z, theta, r = generate_cylindrical_grid()
    dr = r[1, 1, 1] - r[0, 0, 0]
    dtheta = theta[1, 1, 1] - theta[0, 0, 0]
    dz = z[1, 1, 1] - z[0, 0, 0]
    r_1d = r[0, 0, :]

    # Define a scalar field
    f = exp(r) * sin(theta) * cos(z)

    # Define a vector field
    vx = sin(r) * exp(theta) * cos(z)
    vy = sin(r + theta + z)
    vz = exp(r) * cos(theta + z)
    v = np.stack([vx, vy, vz], axis=0)

    # Laplacian
    df_ana = (r - 1) * exp(r) * sin(theta) * cos(z) / r**2
    df_num = pcd.del2(f, dr, dtheta, dz, x=r_1d, coordinate_system="cylindrical")
    check_arr_close(df_ana, df_num)

    # Gradient
    df_ana_x = exp(r) * sin(theta) * cos(z)
    df_ana_y = exp(r) * cos(theta) * cos(z) / r
    df_ana_z = -exp(r) * sin(theta) * sin(z)
    df_ana = np.stack([df_ana_x, df_ana_y, df_ana_z], axis=0)
    df_num = pcd.grad(f, dr, dtheta, dz, x=r_1d, coordinate_system="cylindrical")
    check_arr_close(df_ana, df_num)

    # Curl
    curl_r = -cos(r + theta + z) - exp(r) * sin(theta + z) / r
    curl_theta = -exp(r) * cos(theta + z) - exp(theta) * sin(r) * sin(z)
    curl_z = (
        r * cos(r + theta + z) - exp(theta) * sin(r) * cos(z) + sin(r + theta + z)
    ) / r
    df_ana = np.stack([curl_r, curl_theta, curl_z], axis=0)
    df_num = pcd.curl(v, dr, dtheta, dz, x=r_1d, coordinate_system="cylindrical")
    check_arr_close(df_ana, df_num)

    # Divergence
    df_ana = (
        -exp(r) * sin(theta + z)
        + exp(theta) * cos(r) * cos(z)
        + exp(theta) * sin(r) * cos(z) / r
        + cos(r + theta + z) / r
    )
    df_num = pcd.div(v, dr, dtheta, dz, x=r_1d, coordinate_system="cylindrical")
    check_arr_close(df_ana, df_num)


@test
def derivatives_spherical() -> None:
    """Derivatives in spherical coordinates"""
    phi, theta, r = generate_spherical_grid()
    dr = r[1, 1, 1] - r[0, 0, 0]
    dtheta = theta[1, 1, 1] - theta[0, 0, 0]
    dphi = phi[1, 1, 1] - phi[0, 0, 0]
    r_1d = r[0, 0, :]
    theta_1d = theta[0, :, 0]

    # Define a scalar field
    f = exp(r) * sin(theta) * cos(phi)

    # Define a vector field
    vx = sin(r) * exp(theta) * cos(phi)
    vy = sin(r + theta + phi)
    vz = exp(r) * cos(theta + phi)
    v = np.stack([vx, vy, vz], axis=0)

    # Laplacian
    df_ana = (r**2 + 2 * r - 2) * exp(r) * sin(theta) * cos(phi) / r**2
    df_num = pcd.del2(
        f, dr, dtheta, dphi, x=r_1d, y=theta_1d, coordinate_system="spherical"
    )
    check_arr_close(df_ana, df_num)

    # Gradient
    df_ana_x = exp(r) * sin(theta) * cos(phi)
    df_ana_y = exp(r) * cos(phi) * cos(theta) / r
    df_ana_z = -exp(r) * sin(phi) / r
    df_ana = np.stack([df_ana_x, df_ana_y, df_ana_z], axis=0)
    df_num = pcd.grad(
        f, dr, dtheta, dphi, x=r_1d, y=theta_1d, coordinate_system="spherical"
    )
    check_arr_close(df_ana, df_num)

    # Curl
    curl_r = (
        -exp(r) * sin(theta) * sin(phi + theta)
        + exp(r) * cos(theta) * cos(phi + theta)
        - cos(phi + r + theta)
    ) / (r * sin(theta))
    curl_theta = (
        -(r + 1) * exp(r) * sin(theta) * cos(phi + theta)
        - exp(theta) * sin(phi) * sin(r)
    ) / (r * sin(theta))
    curl_phi = (
        r * cos(phi + r + theta) - exp(theta) * sin(r) * cos(phi) + sin(phi + r + theta)
    ) / r
    df_ana = np.stack([curl_r, curl_theta, curl_phi], axis=0)
    df_num = pcd.curl(
        v, dr, dtheta, dphi, x=r_1d, y=theta_1d, coordinate_system="spherical"
    )
    check_arr_close(df_ana, df_num)

    # Divergence
    df_ana = (
        (r * cos(r) + 2 * sin(r)) * exp(theta) * sin(theta) * cos(phi)
        - exp(r) * sin(phi + theta)
        + sin(theta) * cos(phi + r + theta)
        + sin(phi + r + theta) * cos(theta)
    ) / (r * sin(theta))
    df_num = pcd.div(
        v, dr, dtheta, dphi, x=r_1d, y=theta_1d, coordinate_system="spherical"
    )
    check_arr_close(df_ana, df_num)


def check_arr_close(a, b):
    _assert_close_arr(trim(a), trim(b), "max abs difference")


def trim(arr):
    """
    Given a pencil code array, return a view without the ghost zones
    """
    nghost = 3
    ind = [slice(None) for i in range(arr.ndim)]
    for i in [-1, -2, -3]:
        ind[i] = slice(3, -3)
    return arr[tuple(ind)]


def generate_xyz_grid():
    """
    Generate the grid on which we can construct arrays to test various derivative operators.
    """
    x = np.linspace(0, 1, 15)
    y = np.linspace(0, 1, 15)
    z = np.linspace(0, 1, 15)
    z_grid, y_grid, x_grid = np.meshgrid(z, y, x, indexing="ij")
    return z_grid, y_grid, x_grid


def generate_cylindrical_grid():
    """
    Generate the grid on which we can construct arrays to test various derivative operators.
    """
    r = np.linspace(1, 2, 15)
    theta = np.linspace(0, 1, 15)
    z = np.linspace(0, 1, 15)
    z_grid, theta_grid, r_grid = np.meshgrid(z, theta, r, indexing="ij")
    return z_grid, theta_grid, r_grid


def generate_spherical_grid():
    """
    Generate the grid on which we can construct arrays to test various derivative operators.
    """
    r = np.linspace(1, 2, 15)
    theta = np.linspace(1, 2, 15)
    phi = np.linspace(1, 2, 15)
    phi_grid, theta_grid, r_grid = np.meshgrid(phi, theta, r, indexing="ij")
    return phi_grid, theta_grid, r_grid
