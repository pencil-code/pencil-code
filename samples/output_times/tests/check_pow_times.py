#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

# Set up Python load path and configure a matplotlib backend that does not
# need X11. This needs to happen before importing the pencil module.
import sys
sys.path.append('../../../python')
import matplotlib
matplotlib.use('agg')
import numpy as np

from utils import Sim

def run_and_write(s, f):
    f.write(f"{s} :")
    for a in eval(s):
        f.write(f" {a:g}")
    f.write('\n')

sim = Sim(path="..", quiet=True)
p = sim.p

with open(f'{__file__[:-3]}.out', 'w') as f:
    run_and_write("p.t", f)
    run_and_write("np.sqrt(sim.p.sp[:,0])", f)
