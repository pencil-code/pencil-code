def dill_load(name, folder=False, sim=False, quiet=True):
    """This scripts loads an dill-file. It automatically checks known folders if no folder is specified.
    Args:
        name:        Name of dill-file  (<name>.dill)
        folder:        Folder containing dill file
        sim:        Simulation for checking automatically folders for

    Example:
       to read ".pc/sim.dill" use: dill_load('sim', '.pc')
       or simply dill_load('sim'), since dill_load checks for following folders automatically: '.pc', 'data/.pc'
    """

    import dill
    from os.path import join, exists

    if folder == "pc" and name.startswith("pc/"):
        name = name[3:]
    if not name.endswith(".dill"):
        name = name + ".dill"  # add .dill to name if not already ending with

    # if folder is not defined try to find the dill-file at typical places
    sim_path = "."
    if sim:
        sim_path = sim.path
    if not folder:
        if exists(join(sim_path, "pc", name)):
            folder = join(sim_path, "pc")
            if not quiet:
                print("~ Found " + name + " in " + folder)
        elif exists(join(sim_path, "data/pc", name)):
            folder = join(sim_path, "data/pc")
            if not quiet:
                print("~ Found " + name + " in " + folder)
        elif exists(join(sim_path, ".", name)):
            folder = join(sim_path, ".")
            if not quiet:
                print("~ Found " + name + " in " + folder)
        else:
            print("!! ERROR: Couldnt find file " + name)
            return False

    # open file
    filepath = join(folder, name)
    if not quiet:
        print(filepath)
    # from pc.io import debug_breakpoint; debug_breakpoint()
    try:  # check on existance
        if not exists(filepath) or not exists(join(sim_path, filepath)):
            print("!! ERROR: dill_load couldnt load " + filepath)
            return False
        # try:                                               # open file and return it
        with open(filepath, "rb") as f:
            obj = dill.load(f)
        return obj
        # except:
        # with open(join(sim_path, filepath), 'rb') as f:
        # obj = dill.load(f)
        # return obj

    except:  # if anything goes wrong, try dry importing, i.e. if python2 and python3 usage was mixed
        print("? Something went wrong with the dill importer, trying backup solution..")
        try:
            import pickle

            with open(filepath, "rb") as f:
                u = pickle._Unpickler(f)
                u.encoding = "latin1"
                data = u.load()
                print("? Success!")
                return data
        except:
            print(
                "!! ERROR: Something went wrong while importing dill-file: " + filepath
            )
            return False
