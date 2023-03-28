# gravity.py
#
# 05-may-20
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
""" Derive gravity profile from parameters
"""

import numpy as np
from pencil import read
import os

def grav_profile(
    grav, x, y, z, par=None
):
    if not par:
        par=pc.read.par()
    if grav == "gravz_profile":
        if par.__getattribute__(grav) == "zero":
            return np.zeros_like(z)
        elif par.__getattribute__(grav) == "Ferriere":
            g_A_cgs=4.4e-9; g_C_cgs=1.7e-9; g_B_cgs=6.172e20; g_D_cgs=3.086e21
            g_A =par.g_a_factor*g_A_cgs/par.unit_velocity*par.unit_time
            g_C =par.g_c_factor*g_C_cgs/par.unit_velocity*par.unit_time
            g_B =par.g_b_factor*g_B_cgs/par.unit_length
            g_D =par.g_d_factor*g_D_cgs/par.unit_length
            print(g_A,g_B,g_C,g_D)
            return -(g_A*z/np.sqrt(z**2+g_B**2) + g_C*z/g_D)
        else:
            if rank == 0:
                print("gravz_profile not defined: setting to zero")
            return np.zeros_like(z)
    if grav == "gravy_profile":
        if rank == 0:
            print("gravy_profile not defined: setting to zero")
        return np.zeros_like(y)
    if grav == "gravx_profile":
        if rank == 0:
            print("gravx_profile not defined: setting to zero")
        return np.zeros_like(x)
