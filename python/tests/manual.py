#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Examples from the manual

Currently, this test requires that there is a run directory with very
specific data in samples/conv-slab-noequi/ .

"""


import numpy as np
import os
import re
from typing import Any, List, Tuple

import pencil as pc

from test_utils import (
    test,
    assert_true,
    fail,
    _assert_close,
    _assert_equal_tuple,
    standalone_test,
    test_extracted,
)


@test
def read_var() -> None:
    """Read var.dat (data cube) file."""
    # Fixme: we shouldn't change directories here, as this can influence
    # other tests
    pwd = os.getcwd()
    os.chdir(get_run_dir())
    var = pc.read.var(trimall=True)
    _assert_equal_tuple(var.f.shape, (5, 32, 32, 32))

    def ident(x: Any) -> Any:
        return x

    expected = [
        # (key, extractor, expected, eps)
        ("t", ident, 0.354_875, 1.0e-6),
        ("dx", ident, 0.03125, 1.0e-6),
        ("x", np.mean, 0.0, 1.0e-6),
        ("z", np.mean, 0.550_700, 1.0e-6),
        ("z", lambda z: np.std(z), 0.574_175, 1.0e-6),
        ("f", lambda f: np.mean(f[0, :, :, :]), 0.0, 1.0e-6),
        ("f", lambda f: np.std(f[0, :, :, :]), 0.001_606_503, 1.0e-9),
        ("f", lambda f: np.mean(f[1, :, :, :]), -5.983_958e-11, 1.0e-17),
        ("f", lambda f: np.std(f[1, :, :, :]), 0.002_251_669, 1.0e-9),
    ]
    for (key, extract, expect, eps) in expected:
        test_extracted(getattr(var, key), extract, expect, key, eps)
    os.chdir(pwd)


@test
def get_help() -> None:
    """Get doc strings of imported functions."""
    math_dot_help = pc.math.dot.__doc__
    assert_true(
        re.search(
            r"Take dot product of two pencil-code vectors", math_dot_help
        ),
        "Unexpected docstring for pc.math.dot: {}".format(math_dot_help),
    )


@test
def get_help_standalone() -> None:
    """Get doc string in a separate Python process."""
    math_dot_help = standalone_test(["print(pc.math.dot.__doc__)"])
    marker = r"Take dot product of two pencil-code vectors"
    for line in math_dot_help:
        if re.search(marker, line):
            return
    fail(
        "Line '{}' not found in pc.math.dot.__doc__:\n  {}".format(
            marker, "\n  ".join(math_dot_help)
        )
    )


@test
def read_ts_standalone() -> None:
    """Read time series in a separate Python process."""
    read_and_check_type(
        [
            "time_series = pc.read.ts(",
            "    file_name='time_series.dat',",
            "    datadir='{}'".format(get_data_dir()),
            ")",
        ],
        variable="time_series",
        expected_type="pc.read.timeseries.TimeSeries",
        method="pc.read_ts()",
    )


@test
def read_var_standalone() -> None:
    """Read data cube in a separate Python process."""
    read_and_check_type(
        [
            "data_cube = pc.read.var(",
            "    datadir='{}'".format(get_data_dir()),
            ")",
        ],
        variable="data_cube",
        expected_type="pc.read.varfile.DataCube",
        method="pc.read_var()",
    )


@test
def read_slices_standalone() -> None:
    """Read slices in a separate Python process."""
    read_and_check_type(
        [
            "slices = pc.read.slices(",
            "    field='bb1',",
            "    extension='xy',",
            "    datadir='{}'".format(get_data_dir2()),
            ")",
        ],
        variable="slices",
        expected_type="pc.read.allslices.SliceSeries",
        method="pc.read.slices()",
    )


@test
def remesh() -> None:
    """Remeshing: [Not yet implemented]."""
    pass


def read_and_check_type(
    python_code: List[str], variable: str, expected_type: str, method: str
) -> None:
    """Read quantity in a spearate Python process and check type."""
    standalone_test(
        python_code
        + [
            "assert isinstance({}, {}), \\".format(
                variable, expected_type
            ),
            "    '{} result: expected a {}, got a {})'"
            + ".format('{}', '{}', type({}))".format(
                method, expected_type, variable
            ),
        ]
    )


def get_data_dir() -> str:
    return os.path.join(get_run_dir(), "data")


def get_data_dir2() -> str:
    return os.path.join(get_run_dir2(), "data")


def get_run_dir() -> str:
    pencil_home = os.getenv("PENCIL_HOME")
    assert pencil_home is not None
    run_dir = os.path.join(pencil_home, "samples", "conv-slab-noequi")
    if not os.path.isdir(run_dir):
        raise Exception("Run directory {} does not exist".format(run_dir))
    return run_dir


def get_run_dir2() -> str:
    pencil_home = os.getenv("PENCIL_HOME")
    assert pencil_home is not None
    run_dir = os.path.join(pencil_home, "samples", "conv-slab")
    if not os.path.isdir(run_dir):
        raise Exception("Run directory {} does not exist".format(run_dir))
    return run_dir
