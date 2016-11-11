def get(path='.'):
    """Returns simulation object from 'path, if already existing, or creates new simulation object from path, if its as simulation."""

    from os.path import isdir
    from os.path import join
    from os.path import exists
    from pencilnew.io import load
    from . simulation import simulation

    if exists(join(path, '.pc/sim.dill')):
        sim = load('sim', folder=join(path, '.pc'))
        sim.update()
        return sim
    else:
        from pencilnew import __is_sim_dir__
        if __is_sim_dir__(path):
            return simulation(path)
        else:
            print('? WARNING: No simulation found in '+path+' -> try get_sims maybe?')
            return False

def get_sims(path='.', depth=1, unhide_all=False):
    """Returns all found simulations as object list from all subdirs, not
       following symbolic links.

    Args:
        depth   depth of searching for simulations, default is 1,
                i.e. only one level deeper directories will be scanned.
        unhide  unhides all simulation found if True, if False (default)
                hidden sim will stay hidden.
    """
    from os.path import join
    import numpy as np

    from pencilnew.io import load
    from pencilnew.io import save
    from pencilnew.sim import simulation
    from pencilnew.io import walklevel
    from is_sim_dir import is_sim_dir

    #from pen.intern.class_simdict import Simdict
    #from intern import get_simdict
    #import intern.debug_breakpoint as debug_breakpoint

    print('~ A list of pencil code simulations is generated from this dir downwards, this may take some time..')
    print('~ (Symbolic links will not be followed, since this can lead to infinit recursion.)')

    # get overview of simulations in all lower dirs
    sim_paths = []
    for path, dirs in walklevel('.', depth):
        #if 'start.in' in files and 'run.in' in files:
        # print('path: '+str(path))
        for dir in dirs:
            # print('dirs: '+str(dir))
            sd = join(path, dir)
            if is_sim_dir(sd): print '# Found Simulation in '+sd; sim_paths.append(sd)

    # take care of each simulation found, i.e.
    # generate new simulation object for each and append the sim.-object on sim_list
    sim_list = []
    for path in sim_paths:
        sim = get(path)

        # check if sim.name is already existing as a name for a different simulation (name conflict?)
        for s in sim_list:			# check for double names
            if sim.name == s.name:
                sim.name = sim.name+'#'		# add # to dublicate
                print "? Warning: Found two simulatoins with the same name: "+sim.path+' and '+s.path
                print "? Changed name of "+sim.path+' to '+sim.name+' -> rename simulation and re-export manually'

        if unhide_all: sim.unhide()
        sim.export()
        sim_list.append(sim)

    # is sim_list empty?
    if sim_list == []: print '? WARNING: no simulations found!'
    return sim_list
