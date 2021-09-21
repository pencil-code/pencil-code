def dill_exists(name, folder=False, sim=False):
    """This scripts checks if a certain dill-file already exists.

    Args:
      name:		Name of dill file  (<name>.dill)
      folder:		Folder containing dill file
      sim:        specific simulation where the dill file should be
    """

    from pencil import sim
    from os.path import join, exists

    if not name.endswith(".dill"):
        name = name + ".dill"

    if folder == False:
        if type(sim) == sim.__Simulation__:
            folder = sim.pc_datadir
        else:
            # if folder is not defined try to find file at typical places
            if exists(join("pc", name)):
                folder = "pc"
            elif exists(join("data/pc", name)):
                folder = "data/pc"
            else:
                return False

    file = join(folder, name)
    try:  # check on existance
        if not exists(file):
            return False
        return True

    except:  # if anything goes wrong
        print("!! ERROR: Something went wrong when checking for the dill file!")
        return False
