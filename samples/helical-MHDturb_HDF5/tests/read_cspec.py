#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading cross-spectra
"""

# Set up Python load path and configure a matplotlib backend that does not
# need X11. This needs to happen before importing the pencil module.
import sys
sys.path.append('../../../python')
import matplotlib
matplotlib.use('agg')
import numpy as np
import pencil as pc

sim = pc.sim.get(path="..", quiet=True)

p_uru = pc.read.power(
    datadir=sim.datadir,
    file_name = "power_u_ru.dat",
    )
p_uu = pc.read.power(
    datadir=sim.datadir,
    file_name = "power_u_u.dat",
    )

def write(f, qty):
    f.write(f"{qty} :")
    for xp in eval(qty):
        f.write(f' {xp:g}')
    f.write('\n')

with open(f'{__file__[:-3]}.out', 'w') as f:
    # This should be consistent with what is calculated by p.kin
    # I have manually checked that
    # np.all(np.isclose(p.u_u/p.kin, 2, rtol=1e-2))
    # is True.
    write(f, "p_uu.t")
    write(f, "p_uu.u_u[2,1:6]")

    write(f, "p_uru.t")
    write(f, "p_uru.u_ru[2,1:6]")
