from pathlib import Path
from os.path import (
    isdir,
    join,
    exists,
    basename,
    abspath,
    getmtime,
    )
import numpy as np

from pencil.io import (
    load,
    save,
    walklevel,
    )
from pencil.util import is_sim_dir
from pencil.sim import simulation

def get(path=".", quiet=False):
    """
    get(path=".", quiet=False)

    Return simulation object from 'path, if already existing, or creates new
    simulation object from path, if its as simulation.

    Parameters
    ----------
    path : string
        Base directory where to look for simulation from.

    quiet : bool
        Switches out the output of the function. Default: False.
    """
    if exists(join(path, "pc/sim.dill")):
        try:
            sim = load("sim", folder=join(path, "pc"))
            sim.update(quiet=quiet)

            if sim.path != Path(path).absolute():
                # The name of the directory containing the simulation has somehow changed (maybe the user renamed it). Better to just try to reload the sim from scratch.
                raise RuntimeError

            if getmtime(join(path, "pc/sim.dill")) < getmtime(sim.datadir):
                #Without this, sim.param sometimes does not pick up changes made by the user.
                raise RuntimeError

            return sim
        except:
            import os

            print(
                f"? Warning: sim.dill in {path} is not up to date, recreating simulation object.."
            )
            os.system(f"rm {join(path, 'pc/sim.dill')}")

    if is_sim_dir(path):
        if not quiet:
            print(
                f"~ Found simulation in {path} and simulation object is created for the first time. May take some time.. "
            )
        return simulation(path, quiet=quiet)
    else:
        print(f"? WARNING: No simulation found in {path} -> try get_sims maybe?")
        return False


def get_sims(path_root=".", depth=0, unhide_all=True, quiet=False):
    """
    get_sims(path_root=".", depth=0, unhide_all=True, quiet=False)

    Returns all found simulations as object list from all subdirs, not
    following symbolic links.

    Parameters
    ----------
    path_root : string
        Base directory where to look for simulation from.

    depth : int
        depth of searching for simulations, default is 1,
        i.e. only one level deeper directories will be scanned.

    unhide_all : bool
        Unhides all simulation found if True, if False (default)
        hidden sim will stay hidden.

    quiet : bool
        Switches out the output of the function. Default: False.
    """

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
                    print(f"# Found Simulation in {sd}")
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

def get_cparam(filepath):
    """
    Read contents of src/cparam.local into a dictionary

    filepath : string
        path to cparam.local file relative or absolute
    """
    cpars = dict()

    lines = open(filepath, "r").readlines()

    for aline in lines:
        if not aline[0]=="!":
            for subline in aline.strip().split():
                if "=" in subline:
                    for ssubline in subline.split(","):
                        cpars[ssubline.split("=")[0]]=int(ssubline.split("=")[1])

    return cpars
