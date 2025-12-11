#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Utility functions fro Python tests."""


import numpy as np
import re
import subprocess
from typing import Any, Callable, List, Tuple
import os
import pytest
from decorator import decorate

import pencil as pc

from proboscis_dummy import (
    TestProgram,
    make_test,
    assert_true,
    assert_equal,
    fail,
    )


def _assert_close(
    expected: float, actual: float, variable: str, eps: float = 1.0e-6
) -> None:
    """Assert that actual and expected differ by at most eps."""
    assert_true(
        abs(actual - expected) <= eps,
        "{}: |{} - {}| > {}".format(variable, expected, actual, eps),
    )

def _assert_close_arr(
    expected: np.ndarray, actual: np.ndarray, variable: str, eps: float = 1.0e-6
    ) -> None:
    """Assert that the maximum difference between arrays actual and expected is at most eps."""
    err = np.max(np.abs(actual - expected))
    assert_true(
        err <= eps,
        "{}: error > {}".format(err, eps),
    )

def _assert_equal_tuple(
    expected: Tuple[int, ...], actual: Tuple[int, ...]
) -> None:
    """Assert that actual and expected differ by at most eps."""
    assert_true(expected == actual, "{} ≠ {}".format(expected, actual))


def _pretty_print(value: Any) -> str:
    if isinstance(value, str):
        if '"' in value:
            return "'{}'".format(value)
        else:
            return '"{}"'.format(value)
    else:
        return str(value)


def standalone_test(python_code: List[str]) -> List[str]:
    """Run the given lines of Python code in a separate Python process.

    Prepend the import line
        import pencil as pc

    Arguments:
        python_code:
            The Python code to run, as list of lines.
    Returns:
        The output (stdout) from the Python code, as list of lines.
    """
    p = subprocess.Popen(
        ["python3"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    python_code = ["import pencil as pc"] + python_code
    stdout, stderr = p.communicate("\n".join(python_code))
    status = p.wait()
    if status != 0:
        fail("Failed to run Python code:\n  " + "\n  ".join(python_code))
    return stdout.splitlines()


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


def cmp_extracted(
    value: np.ndarray,
    extract: Callable[[np.ndarray], float],
    expected: float,
    tag: str,
    eps: float,
) -> None:
    """Compare numeric value extracted from a variable to reference value.

    Arguments:
        value:
            The value to investigate. Can be a float or an ndarray.
        extract:
            A function mapping 'value' to a float value.
        expected:
            the expected value for that float
        tag:
            A string describing the value.
        eps:
            Absolute accuracy for the comparison expected vs. actual.
    """
    actual = extract(value)
    _assert_close(
        expected, actual, "{}(var.{})".format(extract.__name__, tag), eps
    )


def get_docstring_standalone(symbol: str, marker: str) -> None:
    """In separate Python process, get doc string and compare to marker"""
    doc_lines = standalone_test(["print({}.__doc__)".format(symbol)])
    for line in doc_lines:
        if re.search(marker, line):
            return
    fail(
        "Line '{}' not found in pc.math.dot.__doc__:\n  {}".format(
            marker, "\n  ".join(doc_lines)
        )
    )

def get_rundir(path):
    pencil_home = os.getenv("PENCIL_HOME")
    assert pencil_home is not None
    run_dir = os.path.join(pencil_home, path)
    if not os.path.isdir(run_dir):
        raise Exception("Run directory {} does not exist".format(run_dir))
    return run_dir

def compile_and_run_sample(path):
    """
    Compiles and runs the sample at `path`. Returns the corresponding simulation object.
    """
    rundir = get_rundir(path)
    sim = pc.sim.get(rundir, quiet=True)

    if sim is False:
        raise RuntimeError(f"Could not get simulation in {rundir}")

    sim.compile(
        bashrc=False,
        cleanall=False,
        autoclean=True,
        raise_errors=True,
        )
    sim.run(bashrc=False, cleardata=True, raise_errors=True)

    return sim

def _require_sample_markers(sample):
    return [
        pytest.mark.xdist_group(f"requires_datadir_{sample}"),
        pytest.mark.integration,
        ]

def require_sample(sample):
    def decorator(func):
        def wrapper(func, *args, **kwargs):
            return func(*args, **kwargs)
        func = decorate(func, wrapper)

        for d in _require_sample_markers(sample):
            func = d(func)

        return func
    return decorator
