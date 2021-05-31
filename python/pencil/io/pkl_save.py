
def pkl_save(obj, name, folder='.pc', ascii=True):
    """This scripts saves any kind of object as a pkl file in folder.

    Args:
    obj:		object you want to save in an pkl file
    name:		name of pkl file, '.pkl' will be added automatically if missing
    ascii:		if you dont want to save pkl file in ascii set this to false
    """
    from pencil.io.mkdir import mkdir
    from os.path import join
    import pickle

    mkdir(folder)        ## prepare folder

    if (not name.endswith('.pkl')): name = name+'.pkl'

    with open(join(folder, name), 'wb') as f:
        if ascii:
            pickle.dump(obj, f, 0)
            return True
        else:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
            return True

    return False
