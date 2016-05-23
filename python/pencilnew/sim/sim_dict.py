# $Id$
#
# read a dicionary storing you simulations from all subfolders
#
# Author: A. Schreiber (aschreiber@mpia.de)

def sim_dict(*args, **kwargs):
    """
    Read simulation setup
    """

  return sim_dict(*args, **kwargs)

class sim_dict(object):
    from pencilnew.io import load
    from pencilnew.io import save

    def __init__(self, datadir='./data'):
    """
    Read simulation data and store it in a dictionary object.
    """
    print '## Loading simulation parameters..'
    sim_dict = load('sim_dict')
    if not sim_dict:
        sim_dict = {}
        save(sim_dict, sim_dict_name)
