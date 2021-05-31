
def sort(simulations, sortby, only_started=False, reverse=False):
  """Sort simulations by a quantity.
  Based on group.py

  Args:
    simulations:    put here a Simulations object or a list of simulations [sim1, sim2, ...]
    sortby:         put here the heyword after which the sorting shall happen
    only_started:   only sort simulations that already has started
    reverse:        reverse order

  Return:
    a list with simulation objects in their sorted order
  """

  from pencil.sim import group

  def flatten(l):
      al = []
      for el in l:
          if type(el) == type(['list']):
              al.extend(flatten(el))
          else:
              al.append(el)

      return al

  sim_dict = group(simulations=simulations,
                   groupby=sortby,
                   sort=True,
                   only_started=only_started,
                   reverse=reverse)

  if reverse:
      return flatten(sim_dict.values())[::-1]
  else:
      return flatten(sim_dict.values())
