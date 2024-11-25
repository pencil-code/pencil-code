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

with open(f'{__file__[:-3]}.out', 'w') as f:
    f.write('av.t :')
    for t in av.t[:5]:
        f.write(f' {t:g}')
    f.write('\n')

    f.write('av.xy.u2mz[3,:4] :')
    for a in av.xy.u2mz[3,:4]:
        f.write(f' {a:g}')
    f.write('\n')

