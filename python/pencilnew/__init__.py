#
# The __init__ file is used not only to import the sub-modules, but also to set everything up properly.

# externals
from os.path import isdir as __isdir__
from os.path import join as __join__
from os.path import exists as __exists__

# sub-modules
import io		   # input und output functions, like save data or call IDL scripts
import diag		   # diagnostic scripts and functions
import visu		   # visualisation routines
import calc		   # math functions and further calculations
import vtk		   #
import math
import read		   # read data and parameters from pencil code directory
import tool_kit	   # all nice workarounds get stored here (e.g., resubmit script)
import export	   # exporter (e.g., vtk, xml)

# shortcuts
from pencilnew.sim.class_sim import Simulation as __Simulation__

# internal routines
def __is_sim__(path='.'):
    """ Checks if a path is pointing at a pencil code simulation."""
    from pencilnew.header import *
    if __isdir__(__join__(path,'data')): return True
    if __exists__(__join__(path,'run.in')) and __exists__(__join__(path,'start.in')) and __exists__(__join__(path,'src/cparam.local')) and __exists__(__join__(path,'src/Makefile.local')):
        return True
    return False

def get_sim(path=''):
    """Returns simulation object from 'path/.pc/' if already existing."""
    if __exists__(__join__(path,'.pc/sim.pkl')):
        return io.load('sim', folder='.pc')
    else:
        return False

# Startup and init. processes
if __is_sim__('.'):
    print '~ Pencil Code Simulation found here! Creating Simulation object, accessible via pc.get_sim().'
    __sim__ = __Simulation__('.')
