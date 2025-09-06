#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading xyaver
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

av = pc.read.aver(
    datadir=sim.datadir,
    simdir=sim.path,
    plane_list=['xy'],
    quiet=True,
    )

def write(f, qty):
    f.write(f"{qty} :")
    for xp in eval(qty):
        f.write(f' {xp:g}')
    f.write('\n')

with open(f'{__file__[:-3]}.out', 'w') as f:
    write(f, "av.t[:5]")
    write(f, "av.xy.u2mz[3,:4]")
    write(f, "av.xy.rhoupmz[3,1:5]")
    write(f, "av.xy.ss2downmz[3,1:5]")
    write(f, "av.xy.dtvmaxz[2,-4:]")
