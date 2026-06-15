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
    bash_opts = {
        'bashrc': False,
        'verbose': False,
        'raise_errors': True,
        }
    sim.compile(
        previous_flags=True,
        additional_options="-t read_all_videofiles",
        **bash_opts,
        )
    sim.bash("src/read_all_videofiles.x", **bash_opts)

sl = pc.read.slices(
    field = ["bb3"],
    datadir=sim.datadir,
    quiet=True,
    )

def write(f, qty):
    f.write(f"{qty} :")
    for xp in eval(qty):
        f.write(f' {xp:g}')
    f.write('\n')

with open(f'{__file__[:-3]}.out', 'w') as f:
    write(f, "sl.xy.bb3[1,4:9,13]")
    write(f, "sl.position.xy")
    write(f, "sl.position.xy2")
    write(f, "sl.position.xz")
    write(f, "sl.position.yz")
