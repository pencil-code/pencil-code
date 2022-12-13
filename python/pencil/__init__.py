"""
The __init__ file is used not only to import the sub-modules, but also to
set everything up properly.
"""

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    l_mpi = True
    l_mpi = l_mpi and (size != 1)
except ImportError:
    rank = 0
    size = 1
    comm = None
    l_mpi = False

if rank == 0:
    print("Warning: pencilnew has moved to pencil.")
    print("         pencil has moved to pencil_old.")
    print("To change your scripts accordingly:")
    print("import pencilnew as pc -> import pencil as pc")
    print("import pencil as pc -> import pencil_old as pc")

try:
    import h5py
# except ImportError:
except:
    if rank == 0:
        print(
            "Error: You need to install h5py library doing 'pip3 install h5py' (Python 3) \
               or 'pip install h5py' (Python 2)."
        )

try:
    import os
    os.environ['HDF5_USE_FILE_LOCKING']='FALSE'
except:
    if rank == 0:
        print("os could not be imported -> HDF5 file locking still in effect.")
    
# Load sub-modules.
from . import io
from . import diag
from . import visu
from . import calc
from . import math
from . import sim
from . import read
from . import tool_kit
from . import export
from . import backpack
from . import ism_dyn
from . import pipelines
from pencil.util import is_sim_dir

# Internal routines.
def get_sim(path=".", quiet=True):
    """
    Return simulation object from 'path, if already existing, or creates new
    simulation object from path, if its as simulation.

    Args:
        path:   Base directory where to look for simulation from.
        quiet:  Switches out the output of the function. Default: False.
    """

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

    return sim.get_sims(
        path_root=path_root, depth=depth, unhide_all=unhide_all, quiet=quiet
    )


def check_dependencies():
    """
    Check if dependencies are fullfilled for pencil.
    """

    import importlib
    from itertools import compress

    dependencies = ["vtk", "tqdm"]

    not_found = [importlib.util.find_spec(dep) is None for dep in dependencies]
    missing_dependencies = list(compress(dependencies, not_found))

    if rank == 0:
        print(
            "WARNING: The following python modules have not been found. \
              Full functionallity may not be granted!"
        )

        if "vtk" in missing_dependencies:
            print("Warning: vtk missing. Try to install the python-vtk or pyevtk module.")
        if "tqdm" in missing_dependencies:
            print("Warning: tqdm missing. Check out https://github.com/tqdm/tqdm.")
