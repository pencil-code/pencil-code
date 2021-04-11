# draglift.py
#
# Routines for computing drag and lift coefficients, as well as non-dimensional
# shedding frequency (Strouhal number) of unsteady flow past a circular cylinder.
#
# Author:
# J. Aarnes (jorgenaarnes@gmail.com)
#
"""
Contains classes and methods to compute drag, lift and Strouhal number
for a cylinder in a cross flow
"""

import numpy as np
from .. import sim

def draglift(simulations, sortby='dx', **kwargs):
    """
    Compute mean drag coefficient, rms lift coefficient, and Strouhal number
    for flow past a circular cylinder.
    If the flow is steady, only drag and lift coefficients are computed.

    Call signature:
      draglift(simulations, d_cylinder=0.1, u_0=1.0, flow_dir='y',
               t_start=-1, sortby='dx'):

    Keyword arguments:

    *simulations*
      array of simulation names to be included in the computations

    *d_cylinder*:
      diameter of the cylinder

    *u_0*
      velocity at the inlet

    *flow_dir*:
      direction of the flow

    *t_start*
      time to start the drag computations from
      should be where the a steady vortex shedding has developed

    *sortby*
      property to sort the arrays by
      typical choices for parametric studies are grid size, length of domain, etc.

    Returns
      object with information: sim-name, drag (mean), lift (rms), and st (non-dim shedding frequency)
    """
    sims = []
    for simulation in simulations:
        sims.append(sim.get(simulation,quiet=True))

    # Sort simulations
    sims = sim.sort(sims,sortby,reverse=True)

    draglift_tmp = []
    for i,thissim in enumerate(sims):
        draglift_tmp.append(Draglift())
        draglift_tmp[i].set_name(thissim)
        draglift_tmp[i].compute(thissim, **kwargs)

    return draglift_tmp


class Draglift(object):
    """
    Grid -- holds pencil code time grid data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.name = ''
        self.drag = 0
        self.lift = 0
        self.st = 0

    def set_name(self,simulation):
        self.name=simulation.name

    #def compute(self,simulation, **kwargs):
        #self.drag,self.lift,self.st = compute_drag(thissim.get_ts(), **kwargs)

    def compute(self, simulation, d_cylinder=0.1, u_0=1.0, flow_dir='y', t_start=-1):
        """
        Compute the drag coefficienc, rms drag, lift coefficient, rms lift and Strouhal number of a given time series
        Requires time series T as input, on the form given in the pencil code where T.t, T.c_dragx and T.c_dragy are avaliable quantities
        Cylinder diameter, fluid mean velocity , flow direction and start time of drag computations velocity should also be given as input.
        """

        time_series = simulation.get_ts()
        time = time_series.t
        if(flow_dir=='y'):
            lift = time_series.c_dragx
            drag = time_series.c_dragy
        elif(flow_dir=='x'):
            lift = time_series.c_dragy
            drag = time_series.c_dragx
        else:
            print('ERROR: Incorrect flow direction specified. Must be x- or y-direction')
            return False

        if(t_start<0):
        # Midpoint as standart start time for drag computations if nothing else is given
            t_0_ind = int(np.ceil(drag.size*0.5))
        else:
            t_0_ind = np.argmax(time_series.t>t_start)

        maximum_drag = find_extrema(drag[t_0_ind:],1)
        minimum_drag = find_extrema(drag[t_0_ind:],-1)

        maximum_lift = find_extrema(lift[t_0_ind:],1)
        minimum_lift = find_extrema(lift[t_0_ind:],-1)

        if(maximum_drag==[]):
            print('ERROR: No vortex shedding detected in flow.')
            return False

        C_d = np.mean(drag[t_0_ind+maximum_drag[0]:t_0_ind+minimum_drag[-1]])
        L_d = np.mean(lift[t_0_ind+maximum_lift[0]:t_0_ind+minimum_lift[-1]])

        rms_drag = np.sqrt(np.mean(np.square(drag[t_0_ind+maximum_drag[0]:t_0_ind+minimum_drag[-1]])))
        rms_lift = np.sqrt(np.mean(np.square(lift[t_0_ind+maximum_lift[0]:t_0_ind+minimum_lift[-1]])))

        freq = (maximum_lift.size-1)/(time[t_0_ind+maximum_lift[-1]]-time[t_0_ind+maximum_lift[0]])
        strouhal= freq*d_cylinder/u_0

        # mean lift and rms drag are computed, but not returned at the moment
        # can perhaps be useful someday, so not removed from compute_drag()-routine
        self.drag = C_d
        self.lift = rms_lift
        self.st = strouhal

def find_extrema(series, maxmin):
    """
    Input: Data series, extrama of intereset (max or min)
    Use scipy to find this
    scipy.signal.argrelextrema(array type of np,comparison operator eg np.greater or np.less)
    """
    from scipy.signal import argrelextrema # use to find local minima or maxima in numpy array
    if (maxmin == 1):
        return argrelextrema(series,np.greater)[0]
    elif(maxmin == -1):
        return argrelextrema(series,np.less)[0] # Tuple, return only nparray
