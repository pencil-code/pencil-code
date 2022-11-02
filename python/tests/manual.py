#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Examples from the manual

Currently, this test requires that there are run directories with very
specific data in samples/conv-slab-noequi/ and samples/conv-slab/ .

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
    get_docstring_standalone,
    read_and_check_type,
    standalone_test,
    test_extracted,
)


@test
def read_var() -> None:
    """Read var.dat (data cube) file."""
    var = pc.read.var(trimall=True, datadir=get_data_dir())
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


@test
def get_help() -> None:
    """Get doc strings of imported functions."""
    math_dot_help = pc.math.dot.__doc__
    assert_true(
        re.search(r"Take dot product of two pencil-code vectors", math_dot_help),
        "Unexpected docstring for pc.math.dot: {}".format(math_dot_help),
    )


@test
def get_help_standalone() -> None:
    """Get doc string in a separate Python process."""
    get_docstring_standalone(
        "pc.math.dot", r"Take dot product of two pencil-code vectors"
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


def get_data_dir() -> str:
    sim = pc.sim.get(get_run_dir(), quiet=True)
    datadir = sim.datadir
    if not os.path.exists(os.path.join(datadir, "time_series.dat")):
        #TODO: Is it a good idea to print stuff like this during a test? Is this information needed at all?
        print("Compiling and running {}. This may take some time.".format(sim.path))
        sim.compile(bashrc=False, cleanall=False)
        sim.run(bashrc=False)
    return datadir


def get_data_dir2() -> str:
    sim = pc.sim.get(get_run_dir2(), quiet=True)
    datadir = sim.datadir
    if not os.path.exists(os.path.join(datadir, "time_series.dat")):
        #TODO: Is it a good idea to print stuff like this during a test? Is this information needed at all?
        print("Compiling and running {}. This may take some time.".format(sim.path))
        sim.compile(bashrc=False, cleanall=False)
        sim.run(bashrc=False)
    if not os.path.exists(os.path.join(datadir, "slice_position.dat")):
        sim.bash("pc_build -t read_all_videofiles", bashrc=False, verbose=False)
        sim.bash("src/read_all_videofiles.x", bashrc=False, verbose=False)
    return datadir


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
