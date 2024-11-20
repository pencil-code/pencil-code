#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading yaver
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
    plane_list=['y'],
    quiet=True,
    )

with open(f'{__file__[:-3]}.out', 'w') as f:
    f.write('av.t :')
    for t in av.t[:2]:
        f.write(f' {t:g}')
    f.write('\n')

    f.write('av.y.ux2mxz[0,10,4:9] :')
    for a in av.y.ux2mxz[0,10,4:9]:
        f.write(f' {a:g}')
    f.write('\n')

    f.write('av.y.ux2mxz[1,10,4:9] :')
    for a in av.y.ux2mxz[1,10,4:9]:
        f.write(f' {a:g}')
    f.write('\n')
