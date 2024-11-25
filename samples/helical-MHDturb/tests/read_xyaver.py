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

datadir = '../data'
simdir = '..'

av = pc.read.aver(
    datadir=datadir,
    simdir=simdir,
    plane_list=['xy'],
    quiet=True,
    )

# Now write to file
file = open(f'{__file__[:-3]}.out', 'w')

file.write('av.t :')
for t in av.t[:5]:
    file.write(f' {t:g}')
file.write('\n')

file.write('av.xy.bymz[3,:4] :')
for a in av.xy.bymz[3,:4]:
    file.write(f' {a:g}')
file.write('\n')

file.close()
