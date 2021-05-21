#
# simulations.py
#
# Create simulations object which stores simulation objects.
#
# Authors:
# A. Schreiber (aschreiber@mpia.de)
#
"""
Contains the simulations class which can be used to store and perform actions on
multiple simulations at once.
"""

def simulations(*args, **kwargs):
    """
    Generate simulations object, which is a container for simulation objects.

    Ways to use the Constructor:
        - no args:          create empty object, simulation objects can be added
                            with add()
        - list, tuple:      provide a list of simulation objects or paths or names
                            if one of the latter is provide, pc.get_sim()
                            is used to generate simulation object from path or
                            name

    Properties:
        self.sims:          direct access on simulations list

    Methods:
        self.add:           add a simulation object, provide simulation object
                            or name or path (pc.get_sim() is then used)
        self.sort           sort self.sims list by, default by name, but also
                            different sorting algorithm are provided
    """

    return __Simulations__(*args, **kwargs)

class __Simulations__(object):
    """
    Simulations object.
    """

    def __init__(self, *args, **kwargs):

        self.sims = []          # list of all simulation stored

        for arg in args:
            #print('\n __init__ : '+ str(arg))
            self.add(args)

        for kw, arg in kwargs:
            print('!! ERROR: Not prepared for kwargs yet!!')
            print(kw+': '+str(arg))

        if type(self.sims) == type(False) and self.sims == False: return False

        # sort self.sims list by name
        import re
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key.name) ]
        self.sims = sorted(self.sims, key = alphanum_key)

        # Done

    def add(self, *args):
        """Add simulation(s) to simulations object.

        Args:
            - simulation object
            - list of simulation objects
            - iterable of simulation objects
            - simulations object
        """

        from pencil.math import is_iterable
        from .simulation import __Simulation__
        from .simulation import simulation
        import numpy as np
        from ..pio import debug_breakpoint
        from .. import get_sim

        for arg in args:

            if isinstance(arg, __Simulation__):
                #print('\n self.add: __Simulation__')
                #print(arg.path)
                self.sims.append(arg)
                return True

            elif isinstance(arg, str):
                #print('\n self.add: string')
                self.sims.append(get_sim(arg))
                return True

            elif is_iterable(arg):
                #print('\n self.add: iterable: '+str(args))
                for ar in arg:
                    #print('n self.add: iterable add: '+str(a))
                    self.add(ar)
                return True

            else:
                print('!! ERROR: Couldnt add to simulations object: '+str(args))

        for kw, arg in kwargs:
            print('!! ERROR: Not prepared for kwargs yet!!')
            print(kw+': '+str(arg))

        return False
