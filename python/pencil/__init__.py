"""
The __init__ file is used not only to import the sub-modules, but also to
set everything up properly.
"""

import lazy_loader as lazy

from .util import (
    pc_print,
    is_sim_dir,
    )

try:
    import h5py
except ImportError:
    pc_print(
        "Error: You need to install h5py library doing 'pip3 install h5py' (Python 3) \
            or 'pip install h5py' (Python 2)."
    )

try:
    import os
    os.environ['HDF5_USE_FILE_LOCKING']='FALSE'
except ImportError:
    #Kishore/2025-Jun-05: what kind of situation is this expected to fail in? `os` is part of the Python standard library.
    pc_print("os could not be imported -> HDF5 file locking still in effect.")

#The following will only be imported when they are first accessed
submodules = [
    "io",
    "diag",
    "visu",
    "calc",
    "math",
    "sim",
    "read",
    "tool_kit",
    "export",
    "backpack",
    "ism_dyn",
    "pipelines",
    ]

__getattr__, __dir__, _ = lazy.attach(__name__, submodules)

# Internal routines.
def get_sim(path=".", quiet=True):
    """
    Return simulation object from 'path, if already existing, or creates new
    simulation object from path, if its as simulation.

    Args:
        path:   Base directory where to look for simulation from.
        quiet:  Switches out the output of the function. Default: False.
    """
    from pencil import sim #needed because of the lazy loading above
    return sim.get(path, quiet=quiet)


def get_sims(path_root=".", depth=1, unhide_all=False, quiet=True):
    """
    Returns all found simulations as object list from all subdirs, not
    following symbolic links.

    Args:
        path_root:  Base directory where to look for simulation from.
        depth:      Depth of searching for simulations, default is 1,
                    i.e. only one level deeper directories will be scanned.
        unhide_all: Unhides all simulation found if True, if False (default)
                    hidden sim will stay hidden.
        quiet:      Switches out the output of the function. Default: True.
    """
    from pencil import sim #needed because of the lazy loading above
    return sim.get_sims(
        path_root=path_root, depth=depth, unhide_all=unhide_all, quiet=quiet
    )


def check_dependencies():
    """
    Check if optional dependencies are fullfilled for pencil.
    """

    import importlib
    from itertools import compress

    dependencies = ["vtk", "tqdm"]

    not_found = [importlib.util.find_spec(dep) is None for dep in dependencies]
    missing_dependencies = list(compress(dependencies, not_found))

    pc_print(
        "WARNING: The following python modules have not been found. \
            Full functionallity may not be granted!"
    )

    if "vtk" in missing_dependencies:
        pc_print("Warning: vtk missing. Try to install the python-vtk or pyevtk module.")
    if "tqdm" in missing_dependencies:
        pc_print("Warning: tqdm missing. Check out https://github.com/tqdm/tqdm.")
