"""
Contains the simulation class which can be used to directly create, access and
manipulate simulations.

"""

import os
from os.path import join, exists, split, islink, realpath, abspath, basename
import numpy as np

from pencil.util import PathWrapper, pc_print

class CommandFailedError(RuntimeError):
    pass

def simulation(*args, **kwargs):
    """
    simulation(*args, **kwargs)

    Generate simulation object from parameters.
    Simulation objects are containers for simulations. pencil can work with
    several of them at once if stored in a simulations object.

    Parameters
    ----------
    path : string
        Path to simulation, default = '.'.

    hidden : bool
        Set True to set hidden flag, default is False.

    quiet : bool
        Suppress irrelevant output, default is False.

    hard : bool
        Force update, default is False

    Properties
    ----------
    self.name:             name of
    self.path:             path to simulation
    self.datadir:          path to simulation data-dir (./data/)
    self.pc_dir:           path to simulation pc-dir (./pc/)
    self.pc_datadir:       path to simulation pendir in datadir (data/pc/)
    self.components:       list of files which are necessary components of
                           the simulation
    self.optionals:        list of files which are optional components of
                           the simulation
    self.start_components: list of files provided at the simulation start
    self.start_optionals:  list of optional files after the simulation start
    self.hidden:           Default is False, if True this simulation will
                               be ignored by pencil
    self.param:            list of param file
    self.grid:             grid object
    self.index             index object
    self.dim:              dim object
    self.tmp_dict:         temporal dictionary of stuff, will not be saved
    """

    return __Simulation__(*args, **kwargs)


class __Simulation__(object):
    """
    Simulation object.
    """

    def __init__(self, path=".", hidden=False, hard=False, quiet=False):
        # from pen.intern.hash_sim import hash_sim

        self.path = PathWrapper(path).absolute()
        self.name = self.path.name  # find out name and store it

        if not quiet:
            print(f"# Creating Simulation object for {self.path}")
        self.datadir = self.path/"data"
        self.pc_dir = self.path/"pc"
        self.pc_datadir = self.path/"data"/"pc"

        # core files of a simulation run
        self.components = [
            "src/cparam.local",
            "src/Makefile.local",
            "start.in",
            "run.in",
            "print.in",
        ]
        # files in which quanitities can be searched
        self.quantity_searchables = ["src/cparam.local", "start.in", "run.in"]
        # optional files that should stick with the simulation when copied
        self.optionals = ["*.in", "*.py", "submit*"]
        # files required from a simulation start
        self.start_components = [
            "index.pro",
            "param.nml",
            #                                 'jobid.dat',
            "pencils.list",
        ]
        # optional files that should stick with the simulation when copied
        self.start_optionals = ["time_series.dat","*.pro", "*.h5"]

        self.hidden = hidden  # hidden is default False
        self.param = False
        self.grid = False
        self.index = False
        self.dim = False
        self.ghost_grid = False
        self.tmp_dict = {}
        self = self.update(hard=hard, quiet=quiet)  # auto-update (read param.nml)
        # Done

    def copy(
        self,
        path_root=".",
        name=False,
        start_optionals=True,
        optionals=True,
        quiet=True,
        rename_submit_script=False,
        OVERWRITE=False,
        lfs=False,
        MB=1,
        count=1,
    ):
        """
        This method does a copy of the simulation object by creating a new
        directory 'name' in 'path_root' and copy all simulation components and
        optiona)
                ls to its directory.
        This method neither links/compiles the simulation.
        If start_optionals it creates data dir.
        It does not overwrite anything, unless OVERWRITE is True.

        Submit Script Rename:
            Name in submit scripts will be renamed if possible!
            Submit scripts will be identified by submit* plus appearenace of old
            simulation name inside, latter will be renamed!

        Parameters
        ----------
        path_root : string
            Path to new sim.-folder(sim.-name). This folder will
            be created if not existing! Relative paths are
            thought to be relative to the python current workdir

        name : string
            Name of new simulation, will be used as folder name.
            Rename will also happen in submit script if found.
            Simulation folders is not allowed to preexist!

        optionals : bool
            Add list of further files to be copied. Wildcasts
            allowed according to glob module!
            Set True to use self.optionals.

        start optionals :
            Add list of further files to be copied.
            Wildcasts allowed according to glob module!
            Set True to use self.optionals.

        quiet : bool
            Set True to suppress output.

        rename_submit_script : bool
            Set False if no renames shall be performed in submit* files

        OVERWRITE : bool
            Set True to overwrite no matter what happens!
        """
        from glob import glob
        from numpy import size
        from os import listdir, symlink
        from shutil import copyfile

        from pencil import get_sim
        from pencil.io import (
            mkdir,
            get_systemid,
            rename_in_submit_script,
            debug_breakpoint,
        )
        from pencil import is_sim_dir

        # set up paths
        if path_root == False or type(path_root) != type("string"):
            print("! ERROR: No path_root specified to copy the simulation to.")
            return False
        path_root = abspath(path_root)  # simulation root dir

        # name and folder of new simulation but keep name of old if sim with old
        # name is NOT existing in NEW directory
        if name == False:
            name = self.name
        if exists(join(path_root, name)) and OVERWRITE == False:
            name = name + "_copy"
            if exists(join(path_root, name)):
                name = name + str(
                    size([f for f in listdir(path_root) if f.startswith(name)])
                )
            print(
                "? Warning: No name specified and simulation with that name "
                + "already found! New simulation name now "
                + name
            )
        path_newsim = join(path_root, name)  # simulation abspath
        path_newsim_src = join(path_newsim, "src")
        if islink(join(path_root, self.name, "data")):
            link_data = True
            oldtmp = os.path.realpath(join(path_root, self.name, "data"))
            newtmp = oldtmp.replace(self.name, name)
            if exists(newtmp) and OVERWRITE == False:
                raise ValueError("Data directory {0} already exists".format(newtmp))
            else:
                path_newsim_data = newtmp
                path_newsim_data_link = join(path_newsim, "data")
        else:
            link_data = False
            path_newsim_data = join(path_newsim, "data")

        path_initial_condition = join(self.path, "initial_condition")
        if exists(path_initial_condition):
            has_initial_condition_dir = True
            path_newsim_initcond = join(path_newsim, "initial_condition")
        else:
            has_initial_condition_dir = False

        if isinstance(optionals, list):
            optionals = self.optionals + optionals
        elif isinstance(optionals, str):
            #Kishore: why is it that only in this case, self.optionals is not appended?
            optionals = [optionals]
        elif optionals is True:
            optionals = self.optionals
        elif optionals is False:
            optionals = []
        else:
            raise TypeError("optionals must be bool, string, or list")

        tmp = []
        for opt in optionals:
            files = glob(join(self.path, opt))
            for f in files:
                tmp.append(basename(f))
        optionals = tmp

        # optional files to be copied
        if isinstance(start_optionals, list):
            start_optionals = self.start_optionals + start_optionals
        elif start_optionals is False:
            start_optionals = []
        elif start_optionals is True:
            start_optionals = self.start_optionals
        elif isinstance(start_optionals, str):
            start_optionals = [start_optionals]
        else:
            raise TypeError("start_optionals must be of type list, str, or bool!")

        tmp = []
        for opt in start_optionals:
            files = glob(join(self.datadir, opt))
            for f in files:
                tmp.append(basename(f))
        start_optionals = tmp
        ## check if the copy was already created
        if is_sim_dir(path_newsim) and OVERWRITE == False:
            if not quiet:
                print(
                    "? WARNING: Simulation already exists."
                    + " Returning with existing simulation."
                )
            return get_sim(path_newsim, quiet=quiet)

        ## expand list of optionals wildcasts

        # check existence of path_root+name, a reason to stop and not overwrite
        if OVERWRITE == False and exists(path_newsim):
            print(
                f"! ERROR: Folder to copy simulation to already exists!\n! -> {path_newsim}"
            )
            return False

        # check existance of self.components
        for comp in self.components:
            if not exists(join(self.path, comp)):
                print(
                    "! ERROR: Couldnt find component "
                    + comp
                    + " from simulation "
                    + self.name
                    + " at location "
                    + join(self.path, comp)
                )
                return False

        # check existance of optionals
        for opt in optionals:
            if not exists(join(self.path, opt)):
                print(
                    "! ERROR: Couldnt find optional component "
                    + opt
                    + " from simulation "
                    + self.name
                    + " at location "
                    + join(self.path, opt)
                )
                return False

        if self.started():
            start_components = self.start_components
        else:
            start_components = []

        # check existance of self.start_components
        for comp in start_components:
            if not exists(join(self.datadir, comp)):
                print(
                    "! ERROR: Couldnt find start_component "
                    + comp
                    + " from simulation "
                    + self.name
                    + " at location "
                    + join(self.path, comp)
                )
                return False

        # check existance of start_optionals
        for opt in start_optionals:
            if not exists(join(self.datadir, opt)):
                print(
                    "! ERROR: Couldnt find optional component "
                    + opt
                    + " from simulation "
                    + self.name
                    + " at location "
                    + join(self.datadir, opt)
                )
                return False

        # create folders
        if mkdir(path_newsim) == False and OVERWRITE == False:
            print(
                f"! ERROR: Couldnt create new simulation directory {path_newsim} !!"
            )
            return False

        if mkdir(path_newsim_src) == False and OVERWRITE == False:
            print(
                f"! ERROR: Couldnt create new simulation src directory {path_newsim_src} !!"
            )
            return False

        if mkdir(path_newsim_data) == False and OVERWRITE == False:
            print(
                f"! ERROR: Couldnt create new simulation data directory {path_newsim_data} !!"
            )
            return False
        if link_data:
            symlink(path_newsim_data, path_newsim_data_link)

        # copy files
        files_to_be_copied = []
        for f in self.components + optionals:
            f_path = abspath(join(self.path, f))
            copy_to = abspath(join(path_newsim, f))
            if f_path == copy_to:
                print(
                    "!! ERROR: file path f_path equal to destination "
                    + "copy_to. Debug this line manually!"
                )
                debug_breakpoint()
            copyfile(f_path, copy_to)

        for f in start_components + start_optionals:
            f_path = abspath(join(self.datadir, f))
            copy_to = abspath(join(path_newsim_data, f))
            if f_path == copy_to:
                print(
                    "!! ERROR: file path f_path equal to destination "
                    + "copy_to. Debug this line manually!"
                )
                debug_breakpoint()
            copyfile(f_path, copy_to)

        # Organizes any personalized initial conditions
        if has_initial_condition_dir:
            if mkdir(path_newsim_initcond) == False and OVERWRITE == False:
                print(
                    f"! ERROR: Couldnt create new simulation initial_condition directory {path_newsim_initcond} !!"
                )
                return False

            for f in listdir(path_initial_condition):
                f_path = abspath(join(path_initial_condition, f))
                copy_to = abspath(join(path_newsim_initcond, f))

                if f_path == copy_to:
                    print(
                        "!! ERROR: file path f_path equal to destination "
                        + "copy_to. Debug this line manually!"
                    )
                    debug_breakpoint()
                copyfile(f_path, copy_to)

        # modify name in submit script files
        if rename_submit_script != False:
            if type(rename_submit_script) == type("STRING"):
                rename_in_submit_script(
                    new_name=rename_submit_script, sim=get_sim(path_newsim)
                )
            else:
                print(
                    "!! ERROR: Could not understand rename_submit_script="
                    + str(rename_submit_script)
                )

        # done
        return get_sim(path_newsim)

    def resume_from_var(self, sim_source, varno, DEBUG=False):
        """
        Copies everything to resume a run from an older state.

        It uses VAR-file number >varno< as new VAR0 and var.dat.
        Does copy PVAR as well if available.

        Parameters
        ----------
        sim_source : string
            Simulation from where to copy all the files.

        varno : int
            var-file number # from which to copy (VAR#).
        """

        from os import listdir
        from os.path import exists, join, isdir
        from pencil.math import is_int
        from pencil.io import mkdir

        def copyfile(src, dst, DEBUG=False):
            from shutil import copy2
            from os.path import exists

            if not exists(src):
                return False
            if DEBUG:
                print("< " + src)
            if DEBUG:
                print("> " + dst)
            copy2(src, dst)

        src = sim_source.datadir
        dst = self.datadir
        if is_int(varno):
            varno = "VAR" + str(int(varno))

        if not exists(src):
            print("! ERROR: Source data directory does not exist: " + str(src))
            return False
        if not exists(dst):
            print("! ERROR: Destination data directory does not exist: " + str(dst))
            return False
        if not varno in sim_source.get_varlist():
            print(
                "! ERROR: Could not find "
                + varno
                + " in procX folder of sim_source: "
                + sim_source.name
            )
            return False

        data_folders = [p for p in listdir(src) if isdir(join(src, p))]
        procX_folder = [p for p in data_folders if p.startswith("proc")]
        for p in data_folders:
            mkdir(join(dst, p))

        # data/
        files = [
            "def_var.pro",
            "dim.dat",
            "index.pro",
            "move-me.list",
            "particles_stalker_header.dat",
            "params.log",
            "pc_constants.pro",
            "pdim.dat",
            "pencils.list",
            "pvarname.dat",
            "svnid.dat",
            "var.general",
            "variables.pro",
            "varname.dat",
        ]
        for f in files:
            copyfile(join(src, f), dst, DEBUG=DEBUG)

        # data/allprocs/
        files = ["grid.dat"]
        for f in files:
            copyfile(join(src, "allprocs", f), join(dst, "allprocs/"), DEBUG=DEBUG)

        # data/procX
        files = ["dim.dat", "grid.dat", "proc_bounds.dat"]
        for X in procX_folder:
            for f in files:
                copyfile(join(src, X, f), join(dst, X), DEBUG=DEBUG)
            copyfile(join(src, X, varno), join(dst, X, "VAR0"), DEBUG=DEBUG)
            copyfile(join(src, X, "P" + varno), join(dst, X, "PVAR0"), DEBUG=DEBUG)
            copyfile(join(src, X, varno), join(dst, X, "var.dat"), DEBUG=DEBUG)
            copyfile(join(src, X, "P" + varno), join(dst, X, "pvar.dat"), DEBUG=DEBUG)

        print("? WARNING: KNOWN ERRORS:")
        print(
            "? RUN MIGHT NOT START BECAUSE data/param.nml can get damaged in"
            + " a run that crashes. This is not fixed by this routine."
        )
        print(
            "? TRY AND START A SINGLE CORE RUN WITH THIS SETUP AND USE THE"
            + " CREATED param.nml FOR YOUR PURPOSE INSTEAD."
        )
        print("? SAME FOR: - tstalk.dat")

        return True

    def update(self, hard=False, quiet=True):
        """Update simulation object:
        if not read in:
            - read param.nml
            - read grid and ghost grid

        Set hard=True to force update.
        """
        from os.path import exists
        from os.path import join
        from pencil.read import param, grid, dim, index

        REEXPORT = False

        if hard == True:
            self.param = False
            self.grid = False
            self.ghost_grid = False
            self.dim = False
            REEXPORT = True

        if self.param == False:
            try:
                if exists(join(self.datadir, "param.nml")):
                    if not quiet: print("~ Reading param.nml.. ")
                    param = param(quiet=quiet, datadir=self.datadir)
                    self.param = {}
                    # read params into Simulation object
                    for key in dir(param):
                        if key.startswith("_") or key == "read":
                            continue
                        if type(getattr(param, key)) in [bool, list, float, int, str]:
                            self.param[key] = getattr(param, key)
                        else:
                            try:
                                # allow for nested param objects
                                self.param[key] = {}
                                for subkey in dir(getattr(param, key)):
                                    if subkey.startswith("_") or subkey == "read":
                                        continue
                                    if type(getattr(getattr(param, key), subkey)) in [
                                        bool,
                                        list,
                                        float,
                                        int,
                                        str,
                                    ]:
                                        self.param[key][subkey] = getattr(
                                            getattr(param, key), subkey
                                        )
                            except:
                                # not nested param objects
                                continue
                    REEXPORT = True
                else:
                    if not quiet:
                        print(
                            f"? WARNING: for {self.path}",
                            + "? Simulation has not run yet! Meaning: No param.nml found!",
                            sep='\n',
                        )
                    REEXPORT = True
            except:
                print(f"! ERROR: while reading param.nml for {self.path}")
                self.param = False
                REEXPORT = True

        if self.param != False and (self.grid == False or self.ghost_grid == False):
            # read grid only if param is not False
            try:
                if not quiet: print("~ Reading grid.. ")
                self.grid = grid(datadir=self.datadir, trim=True, quiet=True)
                if not quiet: print("~ Reading ghost_grid.. ")
                self.ghost_grid = grid(datadir=self.datadir, trim=False, quiet=True)
                if not quiet: print("~ Reading index.. ")
                self.index = index(datadir=self.datadir)
                if not quiet: print("~ Reading dim.. ")
                self.dim = dim(datadir=self.datadir)
                if not quiet:
                    print("# Updating grid and ghost_grid succesfull")
                REEXPORT = True
                # adding lx, dx etc to params
                self.param["Lx"] = self.grid.Lx
                self.param["Ly"] = self.grid.Ly
                self.param["Lz"] = self.grid.Lz
                self.param["lx"] = self.grid.Lx
                self.param["ly"] = self.grid.Ly
                self.param["lz"] = self.grid.Lz
                self.param["dx"] = self.grid.dx
                self.param["dy"] = self.grid.dy
                self.param["dz"] = self.grid.dz
            except:
                if not quiet:
                    print(
                        "? WARNING: Updating grid and ghost_grid "
                        + "was not successfull, since run has not yet started."
                    )
                if self.started() or (not quiet):
                    print(f"? WARNING: Couldnt load grid for {self.path}")
                self.grid = False
                self.ghost_grid = False
                self.index = False
                self.dim = False
                REEXPORT = True
        elif self.param == False:
            if not quiet:
                print(
                    "? WARNING: Updating grid and ghost_grid "
                    + "was not successfull, since run has not yet started."
                )
            self.grid = False
            self.ghost_grid = False
            self.index = False
            self.dim = False
            REEXPORT = True

        if REEXPORT == True:
            self.export()
        return self

    def hide(self):
        """Set hide flag True for this simulation."""
        self.hidden = True
        self.export()

    def unhide(self):
        """Set hide flag False for this simulation."""
        self.hidden = False
        self.export()

    def export(self):
        """Export simulation object to its root/.pc-dir"""
        from pencil.io import save

        if self == False:
            #Kishore (2024-Nov-28): is it ever possible for this code path to be reached? I think not.
            raise RuntimeError("Simulation object is bool object and False!")

        # clean self.tmp_dict
        tmp_dict = self.tmp_dict
        self.tmp_dict = {}

        try:
            save(self, name="sim", folder=self.pc_dir)
            # restore self.tmp_dict
            self.tmp_dict = tmp_dict
        except Exception as e:
            pc_print("Warning: pencil.io.save failed. If dill is not installed, try:")
            pc_print("'pip3 install dill' (Python 3) or 'pip install dill' (Python 2).")
            pc_print(f"The raised error was: {e}")

    def started(self):
        """Returns whether simulation has already started.
        This is indicated by existing time_series.dat in data directory."""
        from os.path import exists, realpath, join

        return exists(join(self.path, "data", "time_series.dat"))

    def compile(
        self,
        cleanall=False,
        fast=False,
        verbose=False,
        hostfile=None,
        autoclean=True,
        **kwargs,
        ):
        """Compiles the simulation. Per default the linking is done before the
        compiling process is called. This method will use your settings as
        defined in your .bashrc-file.

        Parameters
        ----------
        cleanall : bool
            Before calling pc_build, pc_build --cleanall is called.

        verbose : bool
            Activate for verbosity.

        fast : bool
            Set True for fast compilation.

        autoclean : bool
            If compilation fails, automatically set cleanall=True and retry.

        Accepts all other keywords accepted by self.bash
        """

        from pencil import io
        from os.path import join

        timestamp = io.timestamp()

        command = []
        command.append("pc_build")

        if cleanall:
            self.cleanall(verbose=verbose, hostfile=hostfile, **kwargs)
        if fast == True:
            command.append(" --fast")
        if hostfile:
            command.append(" -f " + hostfile)
        if verbose != False:
            print(f"! Compiling {self.path}")

        try:
            ret = self.bash(
                command=" ".join(command),
                verbose=verbose,
                logfile=join(self.pc_dir, "compilelog_" + timestamp),
                **kwargs,
                )
        except CommandFailedError:
            if autoclean:
                ret = None
            else:
                raise
        finally:
            if (ret is not True) and autoclean and (not cleanall):
                #If cleanall was already passed, no point in cleaning again and retrying.
                return self.compile(
                    cleanall=True,
                    autoclean=False, #prevent infinite recursion
                    fast=fast,
                    verbose=verbose,
                    hostfile=hostfile,
                    **kwargs,
                    )
            else:
                return ret

    def build(self, **kwargs):
        """Same as compile()"""
        return self.compile(**kwargs)

    def bash(
        self,
        command,
        verbose="last100",
        logfile=False,
        bashrc=True,
        raise_errors=False,
        ):
        """Executes command in simulation directory.
        This method will use your settings as defined in your .bashrc-file.
        A log file will be produced within 'self.path/pc'-folder

        Parameters
        ----------
        command : string
            Command to be executed, can be a list of commands.

        verbose : bool
            lastN = show last N lines of output afterwards
             False = no output
             True = all output

        bashrc : bool
            True: source bashrc in the subprocess
            False: don't source bashrc in the subprocess. Instead, only pass along the environment variables from the current session.
        
        raise_errors: bool
            If True, a nonzero return code for command will be raised as a Python error.
        """

        import subprocess
        from pencil import io
        from os.path import join
        import os

        timestamp = io.timestamp()
        io.mkdir(self.pc_dir)
        if not type(logfile) == type("string"):
            logfile = join(self.pc_dir, "bash_log_" + timestamp)

        commands = ["cd " + realpath(self.path)]
        # commands.append('source ~/.bashrc')
        # commands.append('shopt -s expand_aliases')

        if type(command) == type(["list"]):
            for c in command:
                commands.append(c)
        elif type(command) == type("string"):
            commands.append(command)
        else:
            print("! ERROR: Couldnt understand the command parameter: " + str(command))

        if bashrc:
            shellcmd = ["/bin/bash", "-i", "-c"]
            env = None
        else:
            shellcmd = ["/bin/bash", "-c"]
            env = os.environ

        shellcmd.append(";".join(commands))

        with open(logfile, "w") as f:
            rc = subprocess.call(shellcmd, stdout=f, stderr=f, env=env)

        if type(verbose) == type("string"):
            outputlength = -int(verbose.split("last")[-1])
            with open(logfile, "r") as f:
                strList = f.read().split("\n")[outputlength:]
                print("\n".join([s for s in strList if not s == ""]))
        elif verbose == True:
            with open(logfile, "r") as f:
                print(f.read())

        if rc == 0:
            return True
        else:
            message = (
                f"! ERROR: Execution ended with error code {rc}"
                + "!\n! Please check log file in"
                + f"\n! {logfile}"
                )
            if raise_errors:
                raise CommandFailedError(message)
            else:
                print(message)
                return rc

    def cleanall(self, verbose=False, hostfile=None, **kwargs):
        """Runs `pc_build --cleanall` in the simulation directory

        Parameters
        ----------
        verbose : bool
            Activate for verbosity.

        Accepts all other keywords accepted by self.bash
        """

        from pencil import io
        from os.path import join

        timestamp = io.timestamp()

        command = []
        command.append("pc_build")
        command.append(" --cleanall")

        if hostfile:
            command.append(" -f " + hostfile)
        if verbose != False:
            print(f"! Cleaning {self.path}")

        return self.bash(
            command=" ".join(command),
            verbose=verbose,
            logfile=join(self.pc_dir, "cleanlog_" + timestamp),
            **kwargs,
            )

    def clear_src(self, do_it=False, do_it_really=False):
        """This method clears the src directory of the simulation!
        All files in src get deleted, except of whats in components and optionals!
        By default, everything except Makefile.local and cparam.local gets erased!

        Args:
            to activate pass        True, True
        """
        from os import listdir
        from os.path import join, exists
        from pencil.io import remove_files as remove

        folder = join(self.path, "src")
        keeps = [f.split("/")[-1] for f in self.components + self.optionals]

        if not exists(folder):
            print("? Warning: No src directory found!")
            return True

        filelist = listdir(folder)  # remove everything INSIDE
        for f in filelist:
            if f in keeps:
                continue
            remove(join(folder, f), do_it=do_it, do_it_really=do_it_really)
        return True

    def clear_data(self, do_it=False, do_it_really=False):
        """This method clears the data directory of the simulation!
        All files in data get deleted!

        Args:
            to activate pass        True, True
        """
        from os import listdir
        from os.path import join, exists
        from pencil.io import remove_files as remove

        folder = join(self.path, "data")
        keeps = []

        if not exists(folder):
            print("? Warning: No data directory found!")
            return True

        filelist = listdir(folder)  # remove everything INSIDE
        for f in filelist:
            if f in keeps:
                continue
            remove(join(folder, f), do_it=do_it, do_it_really=do_it_really)
        return True

    def remove(self, do_it=False, do_it_really=False, remove_data=False):
        """This method removes the WHOLE simulation,
        but NOT the DATA directory per default.
        Do remove_data=True to delete data dir as well.

        Args:
            to activate pass        True, True
            remove_data:            also clear data directory
        """
        from os import listdir
        from os.path import join
        from pencil.io import remove_files as remove

        self.clear_src(do_it=do_it, do_it_really=do_it_really)
        if remove_data:
            self.clear_data(do_it=do_it, do_it_really=do_it_really)

        filelist = listdir(self.path)
        for f in filelist:
            remove(join(self.path, f), do_it=do_it, do_it_really=do_it_really)
        return True

    def get_T_last(self):
        """Returns ts.t[-1] WITHOUTH reading the whole time series!"""

        if self.started() != True:
            return 0

        from os.path import join

        with open(join(self.datadir, "time_series.dat"), "rb") as fh:
            first = next(fh).decode()
            fh.seek(-1024, 2)
            last = fh.readlines()[-1].decode()

        header = [
            i
            for i in first.split("-")
            if not i == "" and not i == "#" and not i == "\n"
        ]
        values = [i for i in last.split(" ") if not i == ""]

        if len(header) != len(values):
            return self.get_ts().t[-1]

        return float(dict(zip(header, values))["t"])

    def get_extent(self, dimensions="xy"):
        """Returns extent as [xmin, xmax, ymin, ymax], as needed by e.g. imshow.

        Arguments:
            dimensions: specify here if you want x, y or z dimensions.
        """

        a = getattr(self.grid, dimensions[0])
        b = getattr(self.grid, dimensions[1])

        da = getattr(self.grid, "d" + dimensions[0])
        db = getattr(self.grid, "d" + dimensions[1])

        return [a[0] - da / 2, a[-1] + da / 2, b[0] - db / 2, b[-1] + db / 2]

    def get_varlist(self, pos=False, particle=False, down=False):
        """
        Get a list of all existing VAR# file names.

        pos = False:                 give full list
        pos = 'last'/'first':        give latest/first var file
        pos = 'lastXXX' / 'firstXXX' give last/first XXX varfiles
        pos = list of numbers:       give varfiles at this positions
        particle = True:             return PVAR- instead of VAR-list
        down = True:                 return list of downsampled snapshots
        """

        import glob
        from os.path import join as join
        from os.path import basename
        from pencil.math import natural_sort

        if particle:
            key = "PVAR"
        elif down:
            key = "VARd"
        else:
            key = "VAR"

        if isinstance(self.param, dict) and self.param['io_strategy'] == 'HDF5':
            proc = "allprocs"
        else:
            proc = "proc0"

        varlist = natural_sort(
            [
                basename(i)
                for i in glob.glob(join(self.datadir, proc) + "/" + key + "*")
            ]
        )

        if pos == False:
            return varlist
        if pos == "first":
            return [varlist[0]]
        if pos == "last":
            return [varlist[-1]]
        if pos.startswith("last"):
            return varlist[-int(pos[4:]) :]
        if pos.startswith("first"):
            return varlist[: int(pos[4:])]
        if type(pos) == type([]):
            if pos[0].startswith("VAR"):
                pos = [i[3:] for i in pos]
            if pos[0].startswith("PVAR"):
                pos = [i[3:] for i in pos]
            return [varlist[int(i)] for i in pos]
        return varlist

    def get_var_time(self, var_file):
        """
        Read varN.list to find the time corresponding to a varfile

        Arguments:
            var_file: string or list of strings

        Returns:
            float or list of floats depending on the type of var_file
        """

        if isinstance(self.param, dict) and self.param['io_strategy'] == 'HDF5':
            proc = "allprocs"

            def fixname(name):
                if name[-3:] == ".h5":
                    name = name[:-3]
                return name
        else:
            proc = "proc0"
            fixname = lambda name: name

        fname = os.path.join(self.datadir, proc, "varN.list")
        if not os.path.isfile(fname):
            raise FileNotFoundError(fname)
        var_dict = {name: float(time) for name, time in np.loadtxt(fname, dtype=str, ndmin=2)}

        if isinstance(var_file, list):
            return [var_dict[fixname(k)] for k in var_file]
        else:
            return var_dict[fixname(var_file)]

    def get_pvarlist(self, pos=False):
        """Same as get_varfiles(pos, particles=True)."""
        return self.get_varlist(pos=pos, particle=True)

    def get_lastvarfilename(self, particle=False, id=False):
        """Returns las varfile name as string."""
        if id == False:
            return self.get_varlist(pos="last", particle=particle)
        return int(self.get_varlist(pos="last", particle=particle)[0].split("VAR")[-1])

    def get_value(self, quantity, DEBUG=False):
        """
        Optimized version of get_value_from_file. Just state quantity for
        simulation and param-list together with searchable components will be
        searched."""

        if DEBUG:
            print("~ DEBUG: Updating simulation.")
        self.update()

        if DEBUG:
            print("~ DEBUG: Searching through simulation.params and dim ...")
        if type(self.param) == type({"dictionary": "with_values"}):
            if quantity in self.param.keys():
                if DEBUG:
                    print("~ DEBUG: " + quantity + " found in simulation.params ...")
                q = self.param[quantity]
                return q
            elif quantity in [d for d in dir(self.dim) if not d.startswith("__")]:
                if DEBUG:
                    print("~ DEBUG: " + quantity + " found in simulation.params ...")
                q = getattr(self.dim, quantity)
                return q

        if DEBUG:
            print("~ DEBUG: Searching through simulation.quantity_searchables ...")
        from pencil.io import get_value_from_file

        for filename in self.quantity_searchables:
            q = get_value_from_file(
                filename, quantity, sim=self, DEBUG=DEBUG, silent=True
            )
            if q is not None:
                if DEBUG:
                    print("~ DEBUG: " + quantity + " found in " + filename + " ...")
                return q
            else:
                if DEBUG:
                    print("~ DEBUG: Couldnt find quantity here.. continue searching")

        print("! ERROR: Couldnt find " + quantity + "!")
        return None

    def get_ts(self, unique_clean=True):
        """Returns time series object.
        Args:
            unique_clean:  set True, np.unique is used to clean up the ts,
                           e.g. remove errors at the end of crashed runs"""
        from pencil.read import ts

        # check if already loaded
        if (
            "ts" in self.tmp_dict.keys()
            and self.tmp_dict["ts"].t[-1] == self.get_T_last()
        ):
            return self.tmp_dict["ts"]

        if self.started():
            ts = ts(sim=self, quiet=True, unique_clean=unique_clean)
            self.tmp_dict["ts"] = ts
            return ts
        else:
            print(
                "? WARNING: Simulation "
                + self.name
                + " has not yet been started. No timeseries available!"
            )
            return False

    def change_value_in_file(
        self, filename, quantity, newValue, filepath=False, DEBUG=False
    ):
        """Same as pencil.io.change_value_in_file."""
        from pencil.io import change_value_in_file

        return change_value_in_file(
            filename, quantity, newValue, sim=self, filepath=filepath, DEBUG=DEBUG
        )

    def run(self, verbose=False, hostfile=None, cleardata=False, **kwargs):
        """Runs the simulation.

        Parameters
        ----------
        verbose : bool
            Whether to print progress and detailed output

        hostfile : string
            Pencil config file to use

        cleardata : bool
            Whether to clear existing data

        Accepts all other keywords accepted by self.bash
        """

        from pencil import io
        from os.path import join
        import os

        timestamp = io.timestamp()

        command = []
        command.append("pc_run")

        if hostfile:
            command.append(" -f " + hostfile)
        if verbose:
            print(f"! Running {self.path}")

        if cleardata:
            if verbose:
                print("Clearing existing data")
            self.clear_data(True, True)

        if not os.path.isdir(self.datadir):
            if os.path.exists(self.datadir):
                if verbose:
                    print("datadir exists but is not a directory; removing it.")
                os.remove(self.datadir)
            if verbose:
                print("Creating datadir")
            os.mkdir(self.datadir)

        return self.bash(
            command=" ".join(command),
            verbose=verbose,
            logfile=join(self.pc_dir, "runlog_" + timestamp),
            **kwargs,
            )
