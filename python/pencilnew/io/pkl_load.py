
def pkl_load(name, folder=False, sim=False):
  """This scripts loads an pkl-file. It automatically checks known folders if no folder is specified.

  Args:
    name:        Name of pkl-file  (<name>.pkl)
    folder:        Folder containing pkl file
    sim:        Simulation for checking automatically folders for

  Example:
    to read ".pc/sim.pkl" use: pkl_load('sim', '.pc')
    or simply pkl_load('sim'), since pkl_load checks for following folders automatically: '.pc', 'data/.pc'
  """

  import pickle
  from os.path import join, exists

  if (not name.endswith('.pkl')): name = name+'.pkl' # add .pkl to name if not already ending with

  # if folder is not defined try to find the pkl-file at typical places
  sim_path = '.'
  if sim: sim_path = sim.path
  if not folder:
    if exists(join(sim_path, '.pc', name)):
      folder = join(sim_path, '.pc'); print('~ Found '+name+' in '+folder)
    elif exists(join(sim_path, 'data/.pc', name)):
      folder = join(sim_path, 'data/.pc'); print('~ Found '+name+' in '+folder)
    elif exists(join(sim_path, '.', name)):
      folder = join(sim_path, '.'); print('~ Found '+name+' in '+folder)
    else:
      print('~ Couldnt find file '+name); return False

  # open file
  file = join(folder, name)
  try:                                                   # check on existance
    if not exists(file) or not exists(join(sim_path, file)):
      print('!! ERROR: pkl_load couldnt load '+file); return False
    try:                                               # open file and return it
      with open(file, 'rb') as f: return pickle.load(f)
    except:
      with open(join(sim_path, file), 'rb') as f: return pickle.load(f)

  except:                        # if anything goes wrong
    print('!! ERROR: Something went wrong while importing pkl-file!'); return False
