#!/usr/bin/python
# -*- coding: utf-8 -*-   vim: set fileencoding=utf-8 :

import pencil as pc


# Description:
#   Read time series and data cube, write a few values

datadir = '../data'

ts = pc.read_ts(datadir=datadir, plot_data=False, quiet=True)
var = pc.read_var(datadir=datadir, trimall=True, quiet=True)


# Now write to file
file = open('read_data.out', 'w')

file.write('ts.times :')
for t in ts.t[0:5]:
    file.write(' %g' % (t, ))
file.write('\n')

file.write('aa(5,5,0:4,1) :')
for a in var.aa[1, 0:5, 5, 5]:
    file.write(' %g' % (a, ))
file.write('\n')

file.close()
