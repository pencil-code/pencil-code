#!/usr/bin/python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""General utility functions.

This is the bottom layer of the pencil package hierarchy, so modules
defined here are not allowed to import any other pencil packages.

"""

import os

MARKER_FILES = [
    "run.in",
    "start.in",
    "src/cparam.local",
    "src/Makefile.local",
]


def is_sim_dir(path="."):
    """Decide if a path is pointing at a pencil code simulation directory.

    The heuristics used is to check for the existence of start.in, run.in,
    src/ cparam.local and src/Makefile.local .

    """
    return all(
        [os.path.exists(os.path.join(path, f)) for f in MARKER_FILES]
    )
