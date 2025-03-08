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
import pytest

import pencil as pc

from test_utils import (
    assert_true,
    fail,
    _assert_close,
    _assert_equal_tuple,
    get_docstring_standalone,
    read_and_check_type,
    standalone_test,
    cmp_extracted,
)


@pytest.mark.integration
def test_read_var(datadir_conv_slab_noequi) -> None:
    """Read var.dat (data cube) file."""
    var = pc.read.var(trimall=True, datadir=datadir_conv_slab_noequi)
    _assert_equal_tuple(var.f.shape, (5, 32, 32, 32))

    def ident(x: Any) -> Any:
        return x

    expected = [
        # (key, extractor, expected, eps)
        ("t", ident, 0.291_004, 1.0e-6),
        ("dx", ident, 0.03125, 1.0e-6),
        ("x", np.mean, 0.0, 1.0e-6),
        ("z", np.mean, 0.550_700, 1.0e-6),
        ("z", lambda z: np.std(z), 0.574_175, 1.0e-6),
        ("f", lambda f: np.mean(f[0, :, :, :]), 0.0, 1.0e-6),
        ("f", lambda f: np.std(f[0, :, :, :]), 0.001_821_54, 1.0e-8),
        ("f", lambda f: np.mean(f[1, :, :, :]), 0, 1.0e-8),
        ("f", lambda f: np.std(f[1, :, :, :]), 0.002_363_680, 1.0e-8),
    ]
    for (key, extract, expect, eps) in expected:
        cmp_extracted(getattr(var, key), extract, expect, key, eps)


def test_get_help() -> None:
    """Get doc strings of imported functions."""
    math_dot_help = pc.math.dot.__doc__
    assert_true(
        re.search(r"Take dot product of two pencil-code vectors", math_dot_help),
        "Unexpected docstring for pc.math.dot: {}".format(math_dot_help),
    )


def test_get_help_standalone() -> None:
    """Get doc string in a separate Python process."""
    get_docstring_standalone(
        "pc.math.dot", r"Take dot product of two pencil-code vectors"
    )


@pytest.mark.integration
def test_read_ts_standalone(datadir_conv_slab_noequi) -> None:
    """Read time series in a separate Python process."""
    read_and_check_type(
        [
            "time_series = pc.read.ts(",
            "    file_name='time_series.dat',",
            "    datadir='{}'".format(datadir_conv_slab_noequi),
            ")",
        ],
        variable="time_series",
        expected_type="pc.read.timeseries.TimeSeries",
        method="pc.read_ts()",
    )


@pytest.mark.integration
def test_read_var_standalone(datadir_conv_slab_noequi) -> None:
    """Read data cube in a separate Python process."""
    read_and_check_type(
        [
            "data_cube = pc.read.var(",
            "    datadir='{}'".format(datadir_conv_slab_noequi),
            ")",
        ],
        variable="data_cube",
        expected_type="pc.read.varfile.DataCube",
        method="pc.read_var()",
    )


@pytest.mark.integration
def test_read_slices_standalone(datadir_conv_slab) -> None:
    """Read slices in a separate Python process."""
    read_and_check_type(
        [
            "slices = pc.read.slices(",
            "    field='bb1',",
            "    extension='xy',",
            "    datadir='{}'".format(datadir_conv_slab),
            ")",
        ],
        variable="slices",
        expected_type="pc.read.allslices.SliceSeries",
        method="pc.read.slices()",
    )


def test_remesh() -> None:
    """Remeshing: [Not yet implemented]."""
    pass
