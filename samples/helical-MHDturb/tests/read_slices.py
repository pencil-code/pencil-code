#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading slices
"""

# Set up Python load path and configure a matplotlib backend that does not
# need X11. This needs to happen before importing the pencil module.
import sys
sys.path.append('../../../python')
import matplotlib
matplotlib.use('agg')
import numpy as np
import pencil as pc
import os

sim = pc.sim.get(path="..", quiet=True)

if not os.path.exists(sim.datadir/"slice_position.dat"):
    #Collect the slice files from the individual processors
    sim.bash("pc_build -t read_all_videofiles", bashrc=False, verbose=False)
    sim.bash("src/read_all_videofiles.x", bashrc=False, verbose=False)

sl = pc.read.slices(
    field = ["bb3"],
    extension = ["xy"],
    datadir=sim.datadir,
    quiet=True,
    )

with open(f'{__file__[:-3]}.out', 'w') as f:
    f.write('sl.xy.bb3[1,4:9,13] :')
    for b in sl.xy.bb3[1,4:9,13]:
        f.write(f' {b:g}')
    f.write('\n')
