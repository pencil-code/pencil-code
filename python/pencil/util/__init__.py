#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""General utility functions.

This is the bottom layer of the pencil package hierarchy, so modules
defined here are not allowed to import any other pencil packages.

"""

import os
import re
import warnings
import pathlib

MARKER_FILES = ["run.in", "start.in", "src/cparam.local", "src/Makefile.local"]


def is_sim_dir(path="."):
    """Decide if a path is pointing at a pencil code simulation directory.

    The heuristics used is to check for the existence of start.in, run.in,
    src/ cparam.local and src/Makefile.local .

    """
    return all([os.path.exists(os.path.join(path, f)) for f in MARKER_FILES])


def ffloat(x):
    """
    Numbers are read from fortran code, which has a specific lenght, in this case 8 char
    If we have scientific notation, it cuts the e and the number doesn't make sense.
    Example:
    Instead of 3.76e-291 it will write 3.76-291

    This function checks and converts all numbers to scientific notation in this case

    KG (2024-Apr-20): it is unclear why this function is needed at all.
    float() seems to correctly handle both "3.76e-291" and "3.76E-291",
    which are what the Fortran code outputs. This function does the
    conversion '3.76-291' -> 3.76e-291, but are there any scenarios where
    the Fortran code produces incorrectly formatted numbers like that?
    """

    try:
        return float(x)

    except:
        warnings.warn("This usage of pc.util.ffloat will be removed soon. If you believe your use-case is legitimate, please email <pencil-code-python@googlegroups.com> describing it.")
        val = re.sub(r"(-?\d+\.?\d*)([+-]\d+)", r"\1E\2", x)
        return float(val)

class PathWrapper(pathlib.WindowsPath if os.name == 'nt' else pathlib.PosixPath):
    """
    See documentation of pathlib.Path.

    This wrapper tries to avoid immediately breaking user code which assumes
    paths are always strings

    KG (2024-Oct-10): added
    KG (2024-Nov-09): fixed usage with Python<3.12 (see https://stackoverflow.com/a/78471242 )
    """
    def _add_warning(self):
        warnings.warn("Adding paths to strings will not work in the future; please change your code before it breaks. If you believe your use-case is legitimate, please email <pencil-code-python@googlegroups.com> describing it.")

    def __add__(self, other):
        self._add_warning()
        return str(self) + str(other)

    def __radd__(self, other):
        self._add_warning()
        return str(other) + str(self)
