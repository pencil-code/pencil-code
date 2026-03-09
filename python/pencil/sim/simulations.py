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

from pencil.util import copy_docstring

class Simulations:
    """
    A container for multiple `Simulation` objects.

    Ways to use the constructor:
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

    def __init__(self, *args, **kwargs):

        self.sims = []  # list of all simulation stored

        self.add(*args, **kwargs)

        if type(self.sims) == type(False) and self.sims == False:
            return False

        # sort self.sims list by name
        import re

        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key.name)]
        self.sims = sorted(self.sims, key=alphanum_key)

    def add(self, *args):
        """Add simulation(s) to simulations object.

        Args:
            - simulation object
            - list of simulation objects
            - iterable of simulation objects
            - simulations object
        """

        from pencil.math import is_iterable
        from pencil.sim.simulation import Simulation
        import numpy as np
        from pencil import get_sim

        for arg in args:

            if isinstance(arg, Simulation):
                self.sims.append(arg)
                return True

            elif isinstance(arg, str):
                self.sims.append(get_sim(arg))
                return True

            elif is_iterable(arg):
                for ar in arg:
                    self.add(ar)
                return True

            else:
                print("!! ERROR: Couldnt add to simulations object: " + str(args))

        for kw, arg in kwargs:
            print("!! ERROR: Not prepared for kwargs yet!!")
            print(kw + ": " + str(arg))

        return False

@copy_docstring(Simulations)
def simulations(*args, **kwargs):
    """
    Wrapper for :py:class:`Simulations`
    """
    #2026-03-09/Kishore: TODO: this seems a redundant wrapper; why not just ask the user to call Simulations(*args, **kwargs)?
    return Simulations(*args, **kwargs)
