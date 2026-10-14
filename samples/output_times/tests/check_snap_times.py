#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

"""
Test reading power_xy
"""

# Set up Python load path and configure a matplotlib backend that does not
# need X11. This needs to happen before importing the pencil module.
import sys
sys.path.append('../../../python')
import matplotlib
matplotlib.use('agg')

from utils import Sim

sim = Sim(path="..", quiet=True)

def run_and_write(s, f):
    f.write(f"{s} :")
    for a in [float(eval(f"var.{s}")) for var in sim.snaps]:
        f.write(f" {a:g}")
    f.write('\n')

with open(f'{__file__[:-3]}.out', 'w') as f:
    run_and_write("t", f)
    run_and_write("special[0,0,0]", f)
