#!/usr/bin/env python

"Python module for everything Pencil Code related."

import re
import numpy

def read_ts(filename=None):

  "Reads in a time series and stores the data in a dictionary."

  if filename==None: filename = 'data/time_series.dat'

  file = open(filename,'rt')
  legend = re.findall('\w+',file.readline())
  data = numpy.loadtxt(file,comments='#')
  file.close()

  list = ["\'"+legend[i]+"\': data[:,"+str(i)+"]" for i in range(len(legend))]
  exec("ts = {"+", ".join(list)+"}")

  return ts
