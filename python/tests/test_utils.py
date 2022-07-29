#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Utility functions fro Python tests."""


import numpy as np
import re
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
