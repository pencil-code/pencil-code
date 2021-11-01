from pencil.io import mkdir,change_value_in_file
import subprocess
from os.path import join
import f90nml

def clone_sims(simset,simsdir=None):
    """
    Create a set of simulation directories from a dictionary object, simset,
    which contains the list of parameters to combine across set of clones

    Signature:

    clone_sims(simset, simsdir=None)

    Parameters
    ----------
    *simset*:  Dictionary of parameters to combine across set of clones.

    *simsdir*:  Root directory for collection of clone simulations.

    Returns
    -------
    Set of uncompiled and unexecuted simulation run directories with 
    parameters updated

    Notes
    -----

    Examples
    --------
    >>> simsdir = '/path/to/set_of_clones'
    >>> params = pencil.pipelines.parameter_table('example_filename.txt')
    >>> params = user_own_trim_table(params)#See pc.pipelines.trim_table
    >>> simset = pencil.pipelines.make_sims_dict(params)
    >>> clone_sims(simset,simsdir=simsdir)
    """
    #If user provides no clone path 
    if not simsdir:
        simsdir = os.getcwd().strip(os.getcwd().split('/')[-1])
    mkdir(simsdir)
    #For each set of simulation parameters create new simulation subdirectory
    for sim in simset:
        newdir = join(simsdir,sim)
        cmd = ['pc_newrun',"-s",newdir]
        #Only compile if makefile.local or cparam.local change
        if "compile" in simset[sim].keys():
            if simset[sim]['compile']:
                cmd = ['pc_newrun',newdir]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        process.communicate()
        for filename in simset[sim]:
            #'compile' flag only used above
            if not filename == 'compile':
                #Files which are f90nml-compatible
                if "run.in" in filename or "start.in" in filename:
                    pars = f90nml.read(filename)
                    newpars = pars.copy()
                    for group in simset[sim][filename]:
                        for item in simset[sim][filename][group]:
                            newpars[group][item] =  simset[sim][filename][group][item]
                    newpars.write(join(newdir,filename),force=True)
                else:
                    file_path = newdir
                    if "local" in filename:
                        file_path = join(file_path,'src')
                    for item in simset[sim][filename]:
                        change_value_in_file(filename, item,
                                simset[sim][filename][item],filepath=file_path)
