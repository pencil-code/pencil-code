def get(path=".", quiet=False):
    """
    Return simulation object from 'path, if already existing, or creates new
    simulation object from path, if its as simulation.

    Args:
        path:   base directory where to look for simulation from.
        quiet:  Switches out the output of the function. Default: False.
    """

    from os.path import isdir, join, exists, basename

    from pencil.io import load
    from pencil.sim.simulation import simulation
    from pencil import is_sim_dir

    if exists(join(path, "pc/sim.dill")):
        try:
            sim = load("sim", folder=join(path, "pc"))
            sim.update(quiet=quiet)
            return sim
        except:
            import os

            print(
                "? Warning: sim.dill in "
                + path
                + " is not up to date, recreating simulation object.."
            )
            os.system("rm " + join(path, "pc/sim.dill"))

    if is_sim_dir(path):
        if quiet == False:
            print(
                "~ Found simulation in "
                + path
                + " and simulation object is created for the first time. May take some time.. "
            )
        return simulation(path, quiet=quiet)
    else:
        print("? WARNING: No simulation found in " + path + " -> try get_sims maybe?")
        return False


def get_sims(path_root=".", depth=0, unhide_all=True, quiet=False):
    """
    Returns all found simulations as object list from all subdirs, not
    following symbolic links.

    Args:
        path_root:  base directory where to look for simulation from.
        depth       depth of searching for simulations, default is 1,
                    i.e. only one level deeper directories will be scanned.
        unhide_all  unhides all simulation found if True, if False (default)
                    hidden sim will stay hidden.
        quiet:      Switches out the output of the function. Default: False.
    """
    from os.path import join, basename
    import numpy as np

    from pencil.io import load
    from pencil.io import save
    from pencil.sim import simulation
    from pencil.io import walklevel
    from pencil import is_sim_dir

    # from pen.intern.class_simdict import Simdict
    # from intern import get_simdict
    # import intern.debug_breakpoint as debug_breakpoint

    if not quiet:
        print(
            "~ A list of pencil code simulations is generated from this dir downwards, this may take some time.."
        )
        print(
            "~ (Symbolic links will not be followed, since this can lead to infinit recursion.)"
        )

    # get overview of simulations in all lower dirs
    sim_paths = []
    for path, dirs in walklevel(path_root, depth):

        for sdir in dirs:
            if sdir.startswith("."):
                continue
            sd = join(path, sdir)
            if is_sim_dir(sd) and not basename(sd).startswith("."):
                if not quiet:
                    print("# Found Simulation in " + sd)
                sim_paths.append(sd)
    if is_sim_dir("."):
        sim_paths.append(".")

    # take care of each simulation found, i.e.
    # generate new simulation object for each and append the sim.-object on sim_list
    sim_list = []
    for path in sim_paths:
        sim = get(path, quiet=quiet)

        # check if sim.name is already existing as a name for a different simulation (name conflict?)
        for s in sim_list:  # check for double names
            if sim.name == s.name:
                sim.name = sim.name + "#"  # add # to dublicate
                if not quiet:
                    print(
                        "? Warning: Found two simulations with the same name: "
                        + sim.path
                        + " and "
                        + s.path
                    )
                    print(
                        "? Changed name of "
                        + sim.path
                        + " to "
                        + sim.name
                        + " -> rename simulation and re-export manually"
                    )

        if unhide_all:
            sim.unhide()
        sim.export()
        sim_list.append(sim)

    # is sim_list empty?
    if sim_list == [] and not quiet:
        print("? WARNING: no simulations found!")
    return sim_list
