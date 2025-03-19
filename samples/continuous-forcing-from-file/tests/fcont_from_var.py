#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

import sys
sys.path.append('../../../python') #so that pencil can be imported without sourceme

import pencil as pc
import numpy as np

"""
Check if the file created by generate_forcing_cont.py was correctly read by forcing.f90
"""

var = pc.read.var(trimall=True, datadir="../data")
force_from_var = np.stack([var.fx, var.fy, var.fz])
with open('fcont_from_var.out', 'w') as f:
	for elem, i in zip(np.nditer(force_from_var), range(0, np.size(force_from_var))):
		f.write("{}\n".format(elem))
