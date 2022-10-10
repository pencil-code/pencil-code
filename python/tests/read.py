#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Test reading data files from Python"""


import numpy as np
import os
from typing import Any, Tuple

from test_utils import (
    test,
    assert_equal,
    assert_true,
    _assert_close,
    _assert_equal_tuple,
    test_extracted,
)

from pencil.read.timeseries import ts
from pencil.read.dims import dim
from pencil.read.varfile import var
from pencil.read.params import param
from pencil.read.powers import power


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


@test
def test_read_time_series() -> None:
    """Read time series."""
    time_series = ts(data_file("time-series-1.dat"), quiet=True)
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
        expect = val
        actual = getattr(time_series, key)
        assert_true(
            np.allclose(expect, actual),
            "time_series.{}: expected {}, got {}".format(
                key, expect, actual
            ),
        )
    _assert_close(time_series.rhom[2], 1.0, "rhom[2]")
    _assert_close(time_series.urms[3], 0.26, "urms[3]")
    _assert_close(time_series.ecrmax[3], 1.835, "ecrmax[3]")


@test
def test_read_dim() -> None:
    """Read dim.dat file."""
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


@test
def test_read_param() -> None:
    """Read param.nml file."""
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


@test
def test_read_var() -> None:
    """Read var.dat (data cube) file."""
    data = var("var.dat", DATA_DIR, proc=0, quiet=True)
    _assert_equal_tuple(data.f.shape, (5, 11, 12, 10))

    def ident(x: Any) -> Any:
        return x

    expected = [
        # (key, extractor, expected, eps)
        ("t", ident, 3.865971, 1.0e-6),
        ("dx", ident, 1.333333, 1.0e-6),
        ("x", np.mean, 0.0, 1.0e-6),
        ("dx", ident, 1.3333334, 1.0e-6),
        ("z", np.mean, 1.6332519, 1.0e-6),
        ("z", lambda z: np.std(z), 5.918408, 1.0e-6),
        ("f", lambda f: np.mean(f[0, :, :, :]), 0.0, 1.0e-6),
        ("f", lambda f: np.mean(f[1, :, :, :]), -1.668_489e-16, 1.0e-22),
        ("f", lambda f: np.mean(f[2, :, :, :]), -7.817_168e-11, 1.0e-17),
        ("f", lambda f: np.mean(f[3, :, :, :]), 1.763_629e-9, 1.0e-15),
        ("f", lambda f: np.mean(f[4, :, :, :]), 2.544_411e-19, 1.0e-25),
        ("f", lambda f: np.std(f[0, :, :, :]), 0.0, 1.0e-6),
        ("f", lambda f: np.std(f[1, :, :, :]), 1.705_128e-9, 1.0e-15),
        ("f", lambda f: np.std(f[2, :, :, :]), 1.171_468e-9, 1.0e-15),
        ("f", lambda f: np.std(f[3, :, :, :]), 2.497_441e-9, 1.0e-15),
        ("f", lambda f: np.std(f[4, :, :, :]), 2.047_645e-19, 1.0e-25),
    ]
    for (key, extract, expect, eps) in expected:
        test_extracted(getattr(data, key), extract, expect, key, eps)


@test
def test_read_power() -> None:
    """Read power spectra"""
    ps = power(datadir=DATA_DIR, quiet=True)

    expected = {
        "t": np.array([1.0477389, 2.0494874]),
        "krms": np.array([0.0, 1.29]),
        "kin": np.array([[1.88e-10, 1.41e-07], [8.16e-10, 1.40e-06]]),
        "hel_kin": np.array([[-2.08e-15, 1.14e-08], [2.26e-15, 2.66e-07]]),
        # TODO: test reading complex 'spectra' as well.
    }
    for key, val in expected.items():
        expect = val
        actual = getattr(ps, key)
        assert_true(
            np.allclose(expect, actual),
            "power.{}: expected {}, got {}".format(key, expect, actual),
        )
