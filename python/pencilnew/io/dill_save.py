
def dill_save(obj, name, folder='.pc'):
    """This scripts saves any kind of object as a dill file in a folder.

    Args:
        obj:		object you want to save in an pkl file
        name:		name of pkl file, '.pkl' will be added automatically if missing
    """
    from pencilnew.io.mkdir import mkdir as __mkdir__
    from os.path import join as __join__
    import pencilnew.backpack.dill as dill

    __mkdir__(folder)        ## prepare folder

    if (not name.endswith('.dill')): name = name+'.dill'

    with open(__join__(folder, name), 'wb') as f:
        dill.dump(obj, f)
        return True

    return False
