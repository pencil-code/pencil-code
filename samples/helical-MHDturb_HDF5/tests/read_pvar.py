#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading pvar.dat
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

pvar2 = pc.read.pvar(
    datadir=datadir,
    pvarfile="PVAR2",
    precision='d',
    )
pvar0 = pc.read.pvar(
    datadir=datadir,
    pvarfile="PVAR0",
    precision='d',
    )

def write(f, qty):
    f.write(f"{qty}: ")
    for xp in eval(qty):
        f.write(f' {xp:g}')
    f.write('\n')


# Now write to file
with open(f'{__file__[:-3]}.out', 'w') as f:
    write(f, "pvar0.xp[:5]")
    write(f, "np.shape(pvar0.xp)")
    
    write(f, "pvar2.xp[:5]")
    write(f, "np.shape(pvar2.xp)")
