
def pkl_save(obj, name, folder='.pc', ascii=True):
    """This scripts saves any kind of object as a pkl file in folder.

    Args:
    obj:		object you want to save in an pkl file
    name:		name of pkl file, '.pkl' will be added automatically if missing
    ascii:		if you dont want to save pkl file in ascii set this to false
    """
    from pencilnew.io.mkdir import mkdir as __mkdir__
    from os.path import join as __join__
    import pickle

    __mkdir__(folder)        ## prepare folder

    if (not name.endswith('.pkl')): name = name+'.pkl'

    with open(__join__(folder, name), 'wb') as f:
        if ascii:
            pickle.dump(obj, f, 0)
            return True
        else:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
            return True

    return False
