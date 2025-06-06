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

class Av(pc.read.averages.Averages):
    def _read_h5_aver(self, *args, **kwargs):
        t, plane = super()._read_h5_aver(*args, **kwargs)
        for k in plane.__dict__.keys():
            if k != "t":
                val = getattr(plane, k)
                if val.ndim == 3:
                    """
                    To be consistent with what happens with io_dist.
                    Axis ordering is now [t,x,z] for yaver, or [t,x,y] for zaver.
                    """
                    setattr(plane, k, val.swapaxes(1,2))
        return t, plane

sim = pc.sim.get(path="..", quiet=True)
av = Av()
av.read(
    datadir=sim.datadir,
    simdir=sim.path,
    plane_list=['z'],
    quiet=True,
    )

with open(f'{__file__[:-3]}.out', 'w') as f:
    f.write('av.t :')
    for t in av.t[:2]:
        f.write(f' {t:g}')
    f.write('\n')

    f.write('av.z.bymxy[0,10,4:9] :')
    for a in av.z.bymxy[0,10,4:9]:
        f.write(f' {a:g}')
    f.write('\n')

    f.write('av.z.bymxy[1,10,4:9] :')
    for a in av.z.bymxy[1,10,4:9]:
        f.write(f' {a:g}')
    f.write('\n')
