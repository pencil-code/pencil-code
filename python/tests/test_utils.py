#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Utility functions fro Python tests."""


import numpy as np
import subprocess
from typing import Any, Callable, List, Tuple

try:
    from proboscis import TestProgram
    from proboscis import test
    from proboscis.asserts import assert_true, assert_equal, fail
except ImportError:
    print("These tests work best with Proboscis installed:")
    print("  pip3 install proboscis")
    print("Continuing with dummy implementation")
    from proboscis_dummy import (
        TestProgram,
        test,
        assert_true,
        assert_equal,
        fail,
    )  # noqa


def _assert_close(
    expected: float, actual: float, variable: str, eps: float = 1.0e-6
) -> None:
    """Assert that actual and expected differ by at most eps."""
    assert_true(
        abs(actual - expected) <= eps,
        "{}: |{} - {}| > {}".format(variable, expected, actual, eps),
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


def test_extracted(
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
