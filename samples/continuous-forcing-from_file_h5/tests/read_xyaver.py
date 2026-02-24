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

def write(f, qty):
    f.write(f"{qty} :")
    for xp in eval(qty):
        f.write(f' {xp:g}')
    f.write('\n')

with open(f'{__file__[:-3]}.out', 'w') as f:
    write(f, "av.t")
    write(f, "av.xy.rufmz[-1,:]")
