
def pkl_load(name, folder=False, sim=False, ascii=True):
  """This scripts loads an pkl-file. It automatically checks known folders if no folder is specified.

  Args:
    name:        Name of pkl-file  (<name>.pkl)
    folder:        Folder containing pkl file
    sim:        Simulation for checking automatically folders for
    ascii:        Switch to False if pkl file was NOT written in ascii

  Example:
    to read ".pc/sim.pkl" use: pkl_load('sim', '.pc')
    or simply pkl_load('sim'), since pkl_load checks for following folders automatically: '.pc', 'data/.pc'
  """

  import pickle
  from os.path import join as __join__
  from os.path import exists as  __exists__

  if (not name.endswith('.pkl')): name = name+'.pkl' # add .pkl to name if not already ending with

  # if folder is not defined try to find the pkl-file at typical places
  sim_path = '.'
  if sim: sim_path = sim.path
  if not folder:
    if __exists__(__join__(sim_path, '.pc', name)):
      folder = __join__(sim_path, '.pc'); print('~ Found '+name+' in '+folder)
    elif __exists__(__join__(sim_path, 'data/.pc', name)):
      folder = __join__(sim_path, 'data/.pc'); print('~ Found '+name+' in '+folder)
    elif __exists__(__join__(sim_path, '.', name)):
      folder = __join__(sim_path, '.'); print('~ Found '+name+' in '+folder)
    else:
      print('~ Couldnt find file '+name); return False

  # open file
  file = __join__(folder, name)
  try:                                                   # check on existance
    if not __exists__(file) or not __exists__(__join__(sim_path, file)):
      print('!! ERROR: pkl_load couldnt load '+file); return False
    try:                                               # open file and return it
      with open(file, 'r') as f: return pickle.load(f)
    except:
      with open(__join__(sim_path, file), 'r') as f: return pickle.load(f)

  except:                        # if anything goes wrong
    print('!! ERROR: Something went wrong while importing pkl-file!'); return False
