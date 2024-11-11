#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

import pencil as pc
import numpy as np

"""
Write out a continuous forcing file
"""

b = np.arange(0,81).reshape((3,3,3,3))
pc.tool_kit.write_forcing_cont(b, outfile="fcont_write_test.out")
