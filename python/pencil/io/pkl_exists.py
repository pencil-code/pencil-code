def pkl_exists(name, folder=False):
    """This scripts checks if a certain pkl-file already exists.

    Args:
      name:        Name of pkl file  (<name>.pkl)
      folder:        Folder containing pkl file
    """

    if not name.endswith(".pkl"):
        name = name + ".pkl"

    from os.path import exists, join

    # if folder is not defined try to find file at typical places
    if not folder:
        if exists(join(".pc", name)):
            folder = ".pc"
        elif exists(join("data/.pc", name)):
            folder = "data/.pc"
        else:
            return False

    file = join(folder, name)
    try:  # check on existance
        if not exists(file):
            return False
        return True

    except:  # if anything goes wrong
        print("!! ERROR: Something went wrong when checking for the pkl file!")
        return False
