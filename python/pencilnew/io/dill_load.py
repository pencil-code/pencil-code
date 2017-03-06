def dill_load(name, folder=False, sim=False):
    """This scripts loads an dill-file. It automatically checks known folders if no folder is specified.
    Args:
        name:        Name of dill-file  (<name>.dill)
        folder:        Folder containing dill file
        sim:        Simulation for checking automatically folders for

    Example:
       to read ".pc/sim.dill" use: dill_load('sim', '.pc')
       or simply dill_load('sim'), since dill_load checks for following folders automatically: '.pc', 'data/.pc'
    """

    import pencilnew.backpack.dill as dill
    from os.path import join as __join__
    from os.path import exists as  __exists__

    if (not name.endswith('.dill')): name = name+'.dill' # add .dill to name if not already ending with

    # if folder is not defined try to find the dill-file at typical places
    sim_path = '.'
    if sim: sim_path = sim.path
    if not folder:
        if __exists__(__join__(sim_path, '.pc', name)):
            folder = __join__(sim_path, '.pc');       print('~ Found '+name+' in '+folder)
        elif __exists__(__join__(sim_path, 'data/.pc', name)):
            folder = __join__(sim_path, 'data/.pc');  print('~ Found '+name+' in '+folder)
        elif __exists__(__join__(sim_path, '.', name)):
            folder = __join__(sim_path, '.');         print('~ Found '+name+' in '+folder)
        else:
            print('~ Couldnt find file '+name);        return False

    # open file
    file = __join__(folder, name)
    try:                                                   # check on existance
        if not __exists__(file) or not __exists__(__join__(sim_path, file)):
            print('!! ERROR: dill_load couldnt load '+file); return False
        try:                                               # open file and return it
            with open(file, 'r') as f: return dill.load(f)
        except:
            with open(__join__(sim_path, file), 'r') as f: return dill.load(f)

    except: # if anything goes wrong
       print('!! ERROR: Something went wrong while importing dill-file!'); return False
