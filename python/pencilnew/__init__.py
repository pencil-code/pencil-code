#
# The __init__ file is used not only to import the sub-modules, but also to
# set everything up properly.
#

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
from . import mycode

# internal routines
def __is_sim_dir__(path='.'):
    """Checks if a path is pointing at a pencil code simulation."""
    return sim.is_sim_dir(path)

def get_sim(path='.', quiet=True):
    """Returns simulation object from 'path, if already existing, or creates new simulation object from path, if its as simulation."""
    return sim.get(path, quiet=quiet)

def get_sims(path_root='.', depth=1, unhide_all=False, quiet=True):
    """Returns all found simulations as object list from all subdirs, not
       following symbolic links.

    Args:
        depth   depth of searching for simulations, default is 1,
                i.e. only one level deeper directories will be scanned.
        unhide  unhides all simulation found if True, if False (default)
                hidden sim will stay hidden.
    """
    return sim.get_sims(path_root=path_root, depth=1, unhide_all=unhide_all, quiet=quiet)

def check_dependencies():
    """ Check if dependencies are fullfilled for pencilnew. """
    import importlib
    from itertools import compress

    dependencies = ['vtk' ,'tqdm']

    not_found = [importlib.util.find_spec(dep) is None for dep in dependencies]
    missing_dependencies = list(compress(dependencies, not_found))

    print('? WARNING: The following python modules have not been found. Full functionallity may not be granted!')

    if 'vtk' in missing_dependencies: print('? Try to install the python-vtk module. But this is outdated anyway, check out pyevtk, its much better anyway!')
    if 'tqdm' in missing_dependencies: print('? Check out https://github.com/tqdm/tqdm')
