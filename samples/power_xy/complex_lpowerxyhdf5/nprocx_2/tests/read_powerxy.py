#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading power_xy
"""

# Set up Python load path and configure a matplotlib backend that does not
# need X11. This needs to happen before importing the pencil module.
import sys
sys.path.append('../../../../../python')
import matplotlib
matplotlib.use('agg')
import numpy as np
import pencil as pc

sim = pc.sim.get(path="..", quiet=True)

p = pc.read.power(
    datadir=sim.datadir,
    file_name="poweruz_xy.h5",
    quiet=True,
    )

def run_and_write(s, f):
    f.write(f"{s} :")
    for a in eval(s):
        f.write(f" {a:g}")
    f.write('\n')

with open(f'{__file__[:-3]}.out', 'w') as f:
    run_and_write("p.t[:5]", f)
    run_and_write("np.real(p.uz_xy[0,16,1,:3])", f)
    run_and_write("np.real(p.uz_xy[3,7:15,1,1])", f)
    run_and_write("np.real(p.uz_xy[3,16:20,1,1])", f)
    run_and_write("np.abs(p.uz_xy[3,9,26:,1])", f)
    run_and_write("np.imag(p.uz_xy[3,9,29,-5:])", f)
