#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

# Set up Python load path and configure a matplotlib backend that does not
# need X11. This needs to happen before importing the pencil module.
import sys
sys.path.append('../../../python')
import matplotlib
matplotlib.use('agg')
import pencil as pc
import numpy as np

"""
Write out a continuous forcing file
"""

b = np.arange(0,81).reshape((3,3,3,3))
pc.tool_kit.write_forcing_cont(b, outfile="fcont_write_test.out")
