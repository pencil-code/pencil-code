
def group(simulations, groupby, sort=True):
  """Group simulation by a quantity. Each Simulation object can only be part of one group.

  Args:
    simulations:       put here a list of Simulation objects or list of simulation [sim1, sim2, ...]
    groupby:        put here the heyword after which the grouping shall happen
    sort:           set True to sort returned dictionary naturally

  Return:
    a dictionary with keywords are the group entries and values are lists of simulations in that group
  """

  from collections import OrderedDict
  from pencilnew.math import natural_sort

  sim_dict_grouped = {}

  # case the groupby-keyword can be found in param.keys()
  if groupby in sim_list[0].param.keys():
    for sim in sim_list:
      q = str(sim.param[groupby])
      if (not q in sim_dict_grouped.keys()):
        sim_dict_grouped[q] = [sim]
      else:
        sim_dict_grouped[q].append(sim)

  # special cases:
  elif groupby in ['Lx', 'Ly', 'Lz']:
    for sim in sim_list:
      q = str(sim.param['lxyz'][0])
      if (not q in sim_dict_grouped.keys()):
        sim_dict_grouped[q] = [sim]
      else:
        sim_dict_grouped[q].append(sim)

  else:
    print '!! ERROR: Coudnt group simulations, no fitting groupby-keyword has been found to match "'+groupby+'"!'
    return False

  if sort:
    sim_dict_grouped_n_sorted = OrderedDict()
    for key in natural_sort(sim_dict_grouped.keys()):
      sim_dict_grouped_n_sorted[key] = sim_dict_grouped[key]
    sim_dict_grouped = sim_dict_grouped_n_sorted

  return sim_dict_grouped
