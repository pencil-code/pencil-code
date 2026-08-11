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
import functools
import numpy as np

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

    except Exception:
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

class SinglePrinter:
    """
    Instances of this class will print their output only on the MPI root process.

    KG (2025-Feb-07): added
    KG (2025-Dec-03): optimized to import mpi4py only when really needed
    """

    @property
    @functools.lru_cache()
    def print(self):
        """
        This logic really belongs in __init__, but we do things this way because
        importing mpi4py.MPI is very slow. If mpi4py.MPI were imported in
        __init__, that alone would account for half the initialization time of
        the Pencil module
        """
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
        except ImportError:
            rank = 0

        if (rank == 0):
            return True
        else:
            return False

    def __call__(self, message):
        """
        message: str
        """
        if self.print:
            print(message)

pc_print = SinglePrinter()

def copy_docstring(original):
    """
    Decorator, to be used for wrapper functions, that makes the docstring of a
    particular function the same as another function.
    
    Rationale: if you consider pc.read.aver and pc.read.averages.Averages.read,
    both the docstrings are currently independently defined, but the former is
    simply a wrapper for the other. When new functionality is added to the
    latter, the docstring for the former is almost never updated, leading to
    the users being shown outdated help text.
    
    Copied from https://softwareengineering.stackexchange.com/a/386758
    
    NOTE: since sphinx-autoapi (used to generated the readthedocs pages) does
    not know about this decorator, it is best to manually add a minimal
    docstring of the form
    \"\"\"
    Wrapper for :py:meth:`Averages.read`
    \"\"\"
    to the wrapper function; this will allow Sphinx to create a link to the
    original function.
    """
    def wrapper(target):
        target.__doc__ = original.__doc__
        return target
    return wrapper

class DotDict(dict):
    """A dict subclass that also supports attribute-style access.

    This allows sim.param to be used both as a dict (sim.param['key'])
    and with attribute access (sim.param.key), so that it is compatible
    with the Param objects returned by pc.read.param() and accepted by
    all reading routines.
    """

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError("{} does not exist".format(key))

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError(key)


class PencilArray(np.ndarray):
    """
    A numpy.ndarray subclass carrying an `axis_order` attribute: a tuple of
    strings naming what each array axis represents, e.g. ('t', 'z', 'y', 'x').

    This is used by the pencil.read routines to make the (historically
    inconsistent) axis ordering of the returned arrays self-documenting, and
    to let users reorder axes by name instead of remembering numeric axis
    positions.

    axis_order is tracked automatically through basic indexing/slicing
    (integer indices drop the corresponding label, slices/Ellipsis keep it,
    np.newaxis inserts an unnamed axis). It is NOT tracked through fancy
    indexing, boolean masks, .transpose(), .reshape(), or ufuncs -- in those
    cases the result's axis_order falls back to None. Use `.reorder()` to
    permute axes by label; it updates axis_order correctly.

    KG (2026-08-08): added
    """

    def __new__(cls, input_array, axis_order=None):
        obj = np.asarray(input_array).view(cls)
        if axis_order is None:
            axis_order = getattr(input_array, "axis_order", None)
        obj.axis_order = None if axis_order is None else tuple(axis_order)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.axis_order = getattr(obj, "axis_order", None)

    def __getitem__(self, key):
        result = super().__getitem__(key)
        if isinstance(result, PencilArray) and self.axis_order is not None:
            result.axis_order = self._reindex_axis_order(key)
        return result

    def _reindex_axis_order(self, key):
        """
        Best-effort computation of the axis_order after applying `key` to an
        array whose axis_order is self.axis_order. Falls back to None
        whenever the effect of `key` on axis identity is not straightforward
        to determine (fancy/boolean indexing).
        """

        if not isinstance(key, tuple):
            key = (key,)

        # Bail out on anything where axis-dropping is not simply
        # "int drops, slice/ellipsis keeps, newaxis inserts".
        for entry in key:
            if entry is Ellipsis or entry is None or entry is np.newaxis:
                continue
            if isinstance(entry, (int, np.integer, slice)):
                continue
            return None

        if sum(1 for entry in key if entry is Ellipsis) > 1:
            return None

        old_order = list(self.axis_order)

        # Expand Ellipsis into explicit full slices so the rest of the logic
        # only has to deal with int / slice / None.
        n_explicit = sum(
            1 for entry in key if isinstance(entry, (int, np.integer, slice))
        )
        n_fill = len(old_order) - n_explicit
        if n_fill < 0:
            return None
        expanded_key = []
        for entry in key:
            if entry is Ellipsis:
                expanded_key.extend([slice(None)] * n_fill)
            else:
                expanded_key.append(entry)

        new_order = []
        old_pos = 0
        for entry in expanded_key:
            if entry is None or entry is np.newaxis:
                new_order.append(None)
            elif isinstance(entry, (int, np.integer)):
                old_pos += 1
            else:  # slice
                new_order.append(old_order[old_pos])
                old_pos += 1
        new_order.extend(old_order[old_pos:])
        return tuple(new_order)

    def reorder(self, *order):
        """
        Return a view of this array with axes permuted to match `order`,
        a sequence of axis labels (e.g. `arr.reorder('x', 'y', 'z')` or
        `arr.reorder(['t', 'x', 'y', 'z'])`). `order` must contain exactly
        the labels in self.axis_order, in the desired new order.
        """

        if len(order) == 1 and not isinstance(order[0], str):
            order = tuple(order[0])

        if self.axis_order is None:
            raise ValueError(
                "This array has no axis_order set; cannot reorder by label."
            )
        if len(order) != len(self.axis_order) or set(order) != set(self.axis_order):
            raise ValueError(
                "Requested order {} does not contain the same labels as "
                "the current axis_order {}.".format(order, self.axis_order)
            )

        permutation = [self.axis_order.index(label) for label in order]
        result = np.transpose(self, axes=permutation).view(PencilArray)
        result.axis_order = tuple(order)
        return result
