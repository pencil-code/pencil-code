def dill_save(obj, name, folder="pc"):
    """This scripts saves any kind of object as a dill file in a folder.

    Args:
        obj:		object you want to save in an pkl file
        name:		name of pkl file, '.pkl' will be added automatically if missing
    """
    from pencil.io.mkdir import mkdir
    from os import remove
    from os.path import join, exists
    import dill

    mkdir(folder)  ## prepare folder

    if not name.endswith(".dill"):
        name = name + ".dill"
    if folder == "pc" and name.startswith("pc/"):
        name = name[3:]

    full_path = join(folder, name)

    if exists(full_path):
        remove(full_path)

    with open(join(folder, name), "wb") as f:
        dill.dump(obj, f)

    return True
