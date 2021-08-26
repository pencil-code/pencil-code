#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""General utility functions.

This is the bottom layer of the pencil package hierarchy, so modules
defined here are not allowed to import any other pencil packages.

"""

import os
import re

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
    """

    try:
        return float(x)

    except:
        val = re.sub(r"(-?\d+\.?\d*)([+-]\d+)", r"\1E\2", x)
        return float(val)
