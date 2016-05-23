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
from pencilnew.sim.simulation import Simulation as __Simulation__

# internal routines
def __is_sim_dir__(path='.'):
    """ Checks if a path is pointing at a pencil code simulation."""
    if __isdir__(__join__(path,'data')) and __exists__(__join__(path,'run.in')) and __exists__(__join__(path,'start.in')) and __exists__(__join__(path,'src/cparam.local')) and __exists__(__join__(path,'src/Makefile.local')):
        return True
    return False

def get_sim(path='.'):
    """Returns simulation object from 'path/.pc/' if already existing."""
    if __exists__(__join__(path,'.pc/sim.pkl')):
        return io.load('sim', folder='.pc')
    else:
        from pencilnew import __is_sim_dir__
        if __is_sim_dir__(path):
            return sim.Simulation(path)
        else:
            print('?? WARNING: No simulation found in '+path+'>. Try get_sims maybe?')
            return False

def get_sims(path='.', depth=1, nuhide=False):
    """Returns all found simulations as object list from all subdirs, not following symbolic links.

    Args:
        depth   depth of searching for simulations, default is only one level higher directories will be searched
        unhide  unhides all simulation found if True, if False hidden sim will stay hidden
    """
    import os
    import numpy as np

    from pencilnew.io import load
    from pencilnew.io import save
    from pencilnew.sim import Simulation
    #from pen.intern.class_simdict import Simdict
    #from intern import get_simdict
    from pencilnew.io import walklevel
    #import intern.debug_breakpoint as debug_breakpoint

    new_sim_list = []
    old_sim_list = False
    print '~ A list of pencil code simulations is generated from this dir downwards, this may take some time..'
    print '~ Symbolic links will not be followed, since this can lead to infinit recursion'

    # get overview of simulations in all lower dirs
    sim_paths = []
    for path,dirs,files in walklevel('.', depth):
      if 'start.in' in files and 'run.in' in files:
        print '# Found Simulation in '+path
        sim_paths.append(path)

    # take care for each simulation found
    # generate new simulation object for all found simulations and append on new_sim_list
    for path in sim_paths:
        sim = Simulation(path)

        # check if sim.name is already existing as name for a different simulation
        for s in new_sim_list:			# check for double names
            if sim.name == s.name:
                sim.name = sim.name+'#'		# add # to dublicate
                print "?? Warning: Found two simulatoins with the same name: "+sim.path+' and '+s.path
                print "?? Changed name of "+sim.path+' to '+sim.name

        if old_sim_list:
            if (not unhide):    # take over hidden status from previouse simDICT if existing and unhide = False (which is default!)
                for old_sim in old_sim_list:
                    if (sim.path == old_sim.path) and (old_sim.hidden==True): sim.hide()

        ## ADD: auto hide simulations which are not runned yet, i.e. without param.nml
        new_sim_list.append(sim)       # add new sim object to new_sim_list
        sim.export()

    # is new_sim_list empty?
    if new_sim_list == []: print '?? WARNING: no simulations found!'
    save(new_sim_list, 'sim_list', folder='.pen')
    return new_sim_list


# Startup and init. processes
if __is_sim_dir__('.'):
    print '~ Pencil Code Simulation found here! Creating Simulation object, accessible via pc.get_sim().'
    __sim__ = __Simulation__('.')
