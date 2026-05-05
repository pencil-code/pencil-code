#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading power spectra
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

p = pc.read.power(
    datadir=sim.datadir,
    quiet=True,
    )

def run_and_write(s, f):
    f.write(f"{s} :")
    for a in eval(s):
        #NOTE: the numbers in power*.dat have only 3 significant figures.
        f.write(f" {a:.2e}")
    f.write('\n')

with open(f'{__file__[:-3]}.out', 'w') as f:
    run_and_write("p.t[:]", f)
    run_and_write("np.real(p.kin[-1,2:6])", f)
