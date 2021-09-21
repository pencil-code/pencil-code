def exists_file(file, folder=False):
    """Better version of exists, compared with os.path.exists!

    Args:
        file:       filename to look for in folder, or complete filepath string
        folder:     optinal, folder in which file should be

    Returns:
        True if file in folder of path specified in file exists.
    """
    from os.path import join, exists, split

    if not folder:
        folder = split(file)[0]
        if folder == "":
            folder = "."
    if exists(join(folder, file)):
        return True
    else:
        return file in [
            i for sublist in [j.split(".") for j in listdir(folder)] for i in sublist
        ]
