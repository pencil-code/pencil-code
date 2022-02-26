from pencil.io import mkdir, change_value_in_file
from pencil.sim import get as get_sim
import subprocess
import os
from os.path import join
import f90nml

def clone_sims(simset, simsdir=None):
    """
    Create a set of simulation directories from a dictionary object, simset,
    which contains the list of parameters to combine across set of clones.

    Signature:

    clone_sims(simset, simsdir=None)

    Parameters
    ----------
    *simset*:  Dictionary of parameters to combine across set of clones.

    *simsdir*:  Root directory for collection of clone simulations.

    Returns
    -------
    Set of uncompiled and unexecuted simulation run directories with
    parameters updated.

    Notes
    -----
    It is assumed that the user is working in the compiled source directory.

    Examples
    --------
    >>> simsdir = '/path/to/set_of_clones'
    >>> params = pencil.pipelines.parameter_table('example_filename.txt')
    >>> params = user_own_trim_table(params)#See pc.pipelines.trim_table
    >>> simset = pencil.pipelines.make_sims_dict(params)
    >>> clone_sims(simset,simsdir=simsdir)
    """
    # If user provides no clone path
    if not simsdir:
        simsdir = os.getcwd().strip(os.getcwd().split("/")[-1])
    mkdir(simsdir)
    # For each set of simulation parameters create new simulation subdirectory
    sourced = False
    for sim in simset:
        newdir = join(simsdir, sim)
        cmd = ["pc_newrun", "-s", newdir]
        # Only compile if makefile.local or cparam.local change
        if "compile" in simset[sim].keys():
            if not sourced:
                moduleinfo = "src/.moduleinfo"
                cmd = ["source " + moduleinfo]
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                try:
                    outs, errs = process.communicate()
                except TimeoutExpired:
                    process.kill()
                    outs, errs = process.communicate()
            if simset[sim]["compile"]:
                cmd = ["pc_newrun", newdir]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        process.communicate()
        for filename in simset[sim]:
            #'compile' flag only used above
            if not filename == "compile":
                # Files which are f90nml-compatible
                if "run.in" in filename or "start.in" in filename:
                    pars = f90nml.read(filename)
                    newpars = pars.copy()
                    for group in simset[sim][filename]:
                        for item in simset[sim][filename][group]:
                            newpars[group][item] = simset[sim][filename][group][item]
                    newpars.write(join(newdir, filename), force=True)
                else:
                    file_path = newdir
                    if "local" in filename:
                        file_path = join(file_path, "src")
                    for item in simset[sim][filename]:
                        change_value_in_file(
                            filename,
                            item,
                            simset[sim][filename][item],
                            filepath=file_path,
                        )

def clone_sims_from_obj(simset, simsdir="..", template_sim=None, specify_nml=True):
    """
    Create a set of simulation directories from a dictionary object, simset,
    which contains the list of parameters to combine across set of clones.
    This differs from the clone_sims function in that it is more configurable, 
    and uses the native Python way of copying the sim.

    Signature:

    clone_sims_from_obj(simset, simsdir="..", template_sim=None, specify_nml=True):

    Parameters
    ----------
    *simset*:
        Dictionary of parameters to combine across set of clones.

    *simsdir*:
        Root directory for collection of clone simulations.
    
    *template_sim*:
        A pencil simulation object (returned by pc.get_sim). 
        If not specified, the simulation in the current directory will be used.
    
    *specify_nml*:
        Whether, for files like run.in, you also specify the namelist in simset. 
        If so, f90nml will be directly used to write the value. If not, the 
        change_value_in_file function will be used even for these files.

    Returns
    -------
    Set of unexecuted simulation run directories with
    parameters updated.

    Example
    --------
    >>> template_sim = pc.get_sim("simulation_to_clone")
    >>> template_sim.optionals.append("job.pbs") #A job submission script which
    ... #one wants to copy along with the other simulation files.
    >>> params = pc.pipelines.parameter_table("example_filename.txt")
    >>> simset = pc.pipelines.make_sims_dict(params)
    >>> pc.pipelines.clone_sims_from_obj(simset, template_sim=template_sim, simsdir=".")
    """
    if not template_sim:
        template_sim = get_sim()
    if not os.path.isdir(simsdir):
        if os.path.exists(simsdir):
            raise RuntimeError("simsdir ({}) exists but is not a directory!".format(simsdir))
        else:
            mkdir(simsdir)
    # For each set of simulation parameters create new simulation subdirectory
    for sim in simset:
        newdir = join(simsdir, sim)
        out = template_sim.copy(path_root=simsdir, name=sim)
        if out is False:
            raise RuntimeError("Copying sim failed")
        for filename in simset[sim]:
            if filename == "compile":
                if simset[sim]["compile"]:
                    print("Warning: clone_sims_from_obj: compilation not implemented yet, so not compiling.")
                    ## KG: I am not sure if the following works as intended; it just hangs for me.
                    #new_sim = get_sim(join(simsdir, newdir))
                    #new_sim.compile()
            elif specify_nml and ("run.in" in filename or "start.in" in filename):
                #Use f90nml for these files
                pars = f90nml.read(filename)
                newpars = pars.copy()
                for group in simset[sim][filename]:
                    for item in simset[sim][filename][group]:
                        newpars[group][item] = simset[sim][filename][group][item]
                newpars.write(join(newdir, filename), force=True)
            else:
                file_path = newdir
                if "local" in filename:
                    file_path = join(file_path, "src")
                for item in simset[sim][filename]:
                    change_value_in_file(
                        filename,
                        item,
                        simset[sim][filename][item],
                        filepath=file_path,
                    )
