#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""Utility functions fro Python tests."""


import subprocess
from typing import List

try:
    from proboscis import TestProgram
    from proboscis import test
    from proboscis.asserts import assert_true, assert_equal, fail
except ImportError:
    print("These tests work best with Proboscis installed:")
    print("  pip3 install proboscis")
    print("Continuing with dummy implementation")
    from proboscis_dummy import TestProgram, test, assert_true, assert_equal, fail  # noqa


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
