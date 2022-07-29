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

pcd = pc.math.derivatives


@test
def partial_derivatives() -> None:
    """Partial derivatives of scalar fields"""
    z, y, x = generate_xyz_grid()
    dx = x[1, 1, 1] - x[0, 0, 0]
    dy = y[1, 1, 1] - y[0, 0, 0]
    dz = z[1, 1, 1] - z[0, 0, 0]

    f = np.sin(x)
    df_ana = np.cos(x)
    df_num = pcd.xder(f, dx)
    check_arr_close(df_ana, df_num)

    f = np.sin(2 * x) * np.cos(y) * np.exp(z)
    df_ana = 2 * np.cos(2 * x) * np.cos(y) * np.exp(z)
    df_num = pcd.xder(f, dx)
    check_arr_close(df_ana, df_num)

    f = np.sin(2 * x) * np.cos(y) * np.exp(z)
    df_ana = -np.sin(2 * x) * np.sin(y) * np.exp(z)
    df_num = pcd.yder(f, dy)
    check_arr_close(df_ana, df_num)

    f = np.sin(2 * x) * np.cos(y) * np.exp(z)
    df_ana = np.sin(2 * x) * np.cos(y) * np.exp(z)
    df_num = pcd.zder(f, dz)
    check_arr_close(df_ana, df_num)

    f = np.sin(2 * x) * np.cos(y) * np.exp(2 * z)
    df_ana = 2 * np.sin(2 * x) * np.cos(y) * np.exp(2 * z)
    df_num = pcd.zder(f, dz)
    check_arr_close(df_ana, df_num)


@test
def derivatives_lap_grad() -> None:
    """Laplacian and gradient of scalar fields"""
    z, y, x = generate_xyz_grid()
    dx = x[1, 1, 1] - x[0, 0, 0]
    dy = y[1, 1, 1] - y[0, 0, 0]
    dz = z[1, 1, 1] - z[0, 0, 0]

    f = np.sin(x) * np.exp(y) * np.cos(z)
    grad_f_x = np.exp(y) * np.cos(x) * np.cos(z)
    grad_f_y = np.exp(y) * np.sin(x) * np.cos(z)
    grad_f_z = -np.exp(y) * np.sin(x) * np.sin(z)
    grad_f = np.stack([grad_f_x, grad_f_y, grad_f_z], axis=0)
    Lap_f = -np.exp(y) * np.sin(x) * np.cos(z)

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
    vx = np.sin(x) * np.exp(y) * np.cos(z)
    vy = np.sin(x + y + z)
    vz = np.exp(x) * np.cos(y + z)
    v = np.stack([vx, vy, vz], axis=0)

    div_v = (
        np.exp(y) * np.cos(x) * np.cos(z)
        + np.cos(x + y + z)
        - np.exp(x) * np.sin(y + z)
    )
    curl_v_x = -np.exp(x) * np.sin(y + z) - np.cos(x + y + z)
    curl_v_y = -np.exp(x) * np.cos(y + z) - np.exp(y) * np.sin(x) * np.sin(z)
    curl_v_z = -np.exp(y) * np.sin(x) * np.cos(z) + np.cos(x + y + z)
    curl_v = np.stack([curl_v_x, curl_v_y, curl_v_z], axis=0)

    check_arr_close(div_v, pcd.div(v, dx, dy, dz))
    check_arr_close(curl_v, pcd.curl(v, dx, dy, dz))


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
