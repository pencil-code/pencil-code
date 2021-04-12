#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Test reading data files from Python"""


import numpy as np
import os
from typing import Any, Tuple

try:
    from proboscis import test
    from proboscis.asserts import assert_equal, assert_true
except ImportError:
    from proboscis_dummy import test, assert_equal, assert_true


from pencil.read.timeseries import ts
from pencil.read.dims import dim
from pencil.read.varfile import var
from pencil.read.params import param


DATA_DIR = os.path.realpath(
    os.path.join(
        __file__, *[os.path.pardir] * 3, "tests", "input", "serial-1"
    )  # ../../tests/input/serial-1
)


def data_file(file_name: str) -> str:
    path = os.path.join(DATA_DIR, file_name)
    if os.path.exists(path):
        return path
    else:
        raise Exception("File {} not found.".format(path))


def _assert_close(
    expected: float, actual: float, eps: float = 1.0e-6
) -> None:
    """Assert that actual and expected differ by at most eps."""
    assert_true(
        abs(actual - expected) <= eps,
        "|{} - {}| > {}".format(expected, actual, eps),
    )


def _assert_equal_tuple(
    expected: Tuple[int, ...], actual: Tuple[int, ...]
) -> None:
    """Assert that actual and expected differ by at most eps."""
    assert_true(expected == actual, "{} ≠ {}".format(expected, actual))


@test()
def test_read_time_series() -> None:
    time_series = ts(data_file("time-series-1.dat"))
    expected = {
        "it": np.array([0, 50, 100, 150]),
        "t": np.array([0.000, 0.441, 0.939, 1.480]),
        "dt": np.array([8.63e-3, 9.18e-3, 1.10e-2, 1.02e-2]),
        "urms": np.array([0.7071, 0.5917, 0.3064, 0.26]),
        "rhom": np.array([1.0, 1.0, 1.0, 1.0]),
        "ecrm": np.array([1.0, 1.020, 1.059, 1.058]),
        "ecrmax": np.array([1.000, 1.930, 2.492, 1.835]),
    }
    for key, val in expected.items():
        assert_true(np.allclose(val, getattr(time_series, key)))


@test()
def test_read_dim() -> None:
    global_dim = dim(DATA_DIR)
    assert_equal(global_dim.mx, 10)
    assert_equal(global_dim.my, 12)
    assert_equal(global_dim.mz, 11)
    assert_equal(global_dim.mvar, 5)
    assert_equal(global_dim.precision, "S")
    assert_equal(global_dim.nghostx, 3)
    assert_equal(global_dim.nghosty, 3)
    assert_equal(global_dim.nghostz, 3)
    assert_equal(global_dim.nprocx, 1)
    assert_equal(global_dim.nprocy, 1)
    assert_equal(global_dim.nprocz, 1)

    proc_dim = dim(DATA_DIR, 0)
    # As we don't have a Dim.__eq__() method:
    attributes = [
        "mx",
        "my",
        "mz",
        "mvar",
        "precision",
        "nghostx",
        "nghosty",
        "nghostz",
    ]
    for attr in attributes:
        assert_equal(
            getattr(global_dim, attr),
            getattr(proc_dim, attr),
            "global_dim.{0} = {1} ≠ proc_dim.{0} = {2}".format(
                attr, getattr(global_dim, attr), getattr(proc_dim, attr)
            ),
        )


@test()
def test_read_param() -> None:
    params = param(DATA_DIR)
    assert_equal(params.coord_system, "cartesian")
    assert_equal(params.lcollective_io, False)
    assert_equal(params.gamma, 1.666_666_6)
    assert_equal(params.kx_uu, 1.0)
    assert_equal(params.cs2top, 1.0)
    assert_equal(params.lhydro, True)
    assert_equal(params.ldensity, True)
    assert_equal(params.lentropy, True)
    assert_equal(params.ltemperature, False)


@test()
def test_read_var() -> None:
    data = var("var.dat", DATA_DIR, proc=0, quiet=True)
    import sys

    print("t({}) = {:.20g}".format(type(data.t), data.t), file=sys.stderr)
    _assert_close(data.t, 3.865971)
    _assert_close(data.dx, 1.333333)
    _assert_close(np.mean(data.x), 0.0)
    _assert_close(data.dx, 1.3333334)
    _assert_close(np.mean(data.z), 1.6332519)
    _assert_close(np.std(data.z), 5.918408)
    _assert_equal_tuple(data.f.shape, (5, 11, 12, 10))
    _assert_close(np.mean(data.f[0, :, :, :]), 0.0)
    _assert_close(np.mean(data.f[1, :, :, :]), -1.668_489e-16, 1.0e-22)
    _assert_close(np.mean(data.f[2, :, :, :]), -7.817_168e-11, 1.0e-17)
    _assert_close(np.mean(data.f[3, :, :, :]), 1.763_629e-9, 1.0e-15)
    _assert_close(np.mean(data.f[4, :, :, :]), 2.544_411e-19, 1.0e-25)
    _assert_close(np.std(data.f[0, :, :, :]), 0.0)
    _assert_close(np.std(data.f[1, :, :, :]), 1.705_128e-9, 1.0e-15)
    _assert_close(np.std(data.f[2, :, :, :]), 1.171_468e-9, 1.0e-15)
    _assert_close(np.std(data.f[3, :, :, :]), 2.497_441e-9, 1.0e-15)
    _assert_close(np.std(data.f[4, :, :, :]), 2.047_645e-19, 1.0e-25)
