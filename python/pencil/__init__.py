'''
The __init__ file is used not only to import the sub-modules, but also to
set everything up properly.
'''

print("Warning: pencilnew has moved to pencil.")
print("         pencil has moved to pencil_old.")
print("To change your scripts accordingly:")
print("import pencilnew as pc -> import pencil as pc")
print("import pencil as pc -> import pencil_old as pc")

try:
    import h5py
#except ImportError:
except:
    print("Error: You need to install h5py library doing 'pip3 install h5py' (Python 3) \
           or 'pip3 install h5py' (Python 2).")
    import h5py

# Load sub-modules.
from . import io           # input und output functions, like save data or call IDL scripts
from . import diag         # diagnostic scripts and functions
from . import visu         # visualization routines
from . import calc         # math functions and further calculations
from . import math         # basic math functions, like products and derivatives
from . import sim          # handling simulations as python objects
from . import read         # read data and parameters from pencil code directory
from . import tool_kit     # all nice workarounds get stored here (e.g., resubmit script)
from . import export       # exporter (e.g., vtk, xml)
from . import backpack     # third party modules, tribute to the author!
from . import ism_dyn      # diagnostics for ism dynamo simulations


# Internal routines.
def __is_sim_dir__(path='.'):
    """
    Check if a path is pointing at a pencil code simulation.
    """

    return sim.is_sim_dir(path)


def get_sim(path='.', quiet=True):
    """
    Return simulation object from 'path, if already existing, or creates new
    simulation object from path, if its as simulation.

    Args:
        path:   Base directory where to look for simulation from.
        quiet:  Switches out the output of the function. Default: False.
    """

    return sim.get(path, quiet=quiet)


def get_sims(path_root='.', depth=1, unhide_all=False, quiet=True):
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

    return sim.get_sims(path_root=path_root, depth=depth, unhide_all=unhide_all, quiet=quiet)


def check_dependencies():
    """
    Check if dependencies are fullfilled for pencil.
    """

    import importlib
    from itertools import compress

    dependencies = ['vtk', 'tqdm']

    not_found = [importlib.util.find_spec(dep) is None for dep in dependencies]
    missing_dependencies = list(compress(dependencies, not_found))

    print('WARNING: The following python modules have not been found. \
          Full functionallity may not be granted!')

    if 'vtk' in missing_dependencies:
        print('Warning: vtk missing. Try to install the python-vtk or pyevtk module.')
    if 'tqdm' in missing_dependencies:
        print('Warning: tqdm missing. Check out https://github.com/tqdm/tqdm.')
