def group(simulations, groupby, sort=True, only_started=False, reverse=False):
    """
    group(simulations, groupby, sort=True, only_started=False, reverse=False)

    Group simulation by a quantity. Each Simulation object can only be part of one group.

    Parameters
    ----------
    simulations : list of sim objects
        Put here a Simulations object or a list of simulations [sim1, sim2, ...].

    groupby : string
        Put here the heyword after which the grouping shall happen.

    sort : bool
        Set True to sort returned dictionary naturally.

    only_started : bool
        Only group simulations that already has started.

    reverse : bool
        Flag for reverse grouping.

    Returns
    -------
    Dictionary with keywords are the group entries and values are lists of simulations in that group.
    """

    from collections import OrderedDict
    from pencil.math import natural_sort

    sim_dict_grouped = {}

    if type(simulations) == type(["list"]):
        sim_list = simulations
    # elif type(simulations) == Simulations:
    #      sim_list = simulations.sims
    else:
        print("!! ERROR: Dont know how to interprated simulations argument..")
        return False

    # sort out simulations that has not started
    if only_started == True:
        sim_list = [s for s in sim_list if s.started()]

    # special cases:
    if groupby in ["Lx", "Ly", "Lz"]:
        if groupby[-1] == "x":
            ii = 0
        elif groupby[-1] == "y":
            ii = 1
        elif groupby[-1] == "z":
            ii = 2
        for sim in sim_list:
            q = str(sim.param["lxyz"][ii])
            if not q in sim_dict_grouped.keys():
                sim_dict_grouped[q] = [sim]
            else:
                sim_dict_grouped[q].append(sim)

    elif groupby in ["nx", "ny", "nz"]:
        for sim in sim_list:
            q = str(getattr(sim.dim, groupby))
            if not q in sim_dict_grouped.keys():
                sim_dict_grouped[q] = [sim]
            else:
                sim_dict_grouped[q].append(sim)

    # case the groupby-keyword can be found via __simulation__.get_value
    elif sim_list[0].get_value(groupby) != None:
        for sim in sim_list:
            q = str(sim.get_value(groupby))
            if not q in sim_dict_grouped.keys():
                sim_dict_grouped[q] = [sim]
            else:
                sim_dict_grouped[q].append(sim)

    else:
        print(
            '!! ERROR: Coudnt group simulations, no fitting groupby-keyword has been found to match "'
            + groupby
            + '"!'
        )
        return False

    if sort:
        sim_dict_grouped_n_sorted = OrderedDict()
        keys = sim_dict_grouped.keys()
        for key in natural_sort(keys, reverse=reverse):
            sim_dict_grouped_n_sorted[key] = sim_dict_grouped[key]

        return sim_dict_grouped_n_sorted

    return sim_dict_grouped
