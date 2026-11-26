#!/usr/bin/env python3
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading phiaver
"""

# Set up Python load path and configure a matplotlib backend that does not
# need X11. This needs to happen before importing the pencil module.
import sys
sys.path.append('../../../../python')
import matplotlib
matplotlib.use('agg')
import numpy as np
import pencil as pc

sim = pc.sim.get(path="..", quiet=True)

av = pc.read.aver(
    plane_list=['phi'],
    datadir=sim.datadir,
    simdir=sim.path,
    )

def run_and_write(s, f):
    f.write(f"{s} :")
    for a in eval(s):
        f.write(f" {a}")
    f.write('\n')

with open(f'{__file__[:-3]}.out', 'w') as f:
    run_and_write("av.t", f)
    run_and_write("av.phi.u2mphi[2,15:20,1]", f)
