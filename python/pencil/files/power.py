# $Id$

import numpy as N

def read_power(file):

    infile = open(file, 'r')
    lines = infile.readlines()
    infile.close()

    infile = open(file, 'r')
    t=N.zeros(1, dtype='Float32')
    data=N.zeros(1, dtype='Float32')
    for i in range(0, len(lines), 2):
      st=infile.readline()
      t=N.append(t, float(st))
      st=infile.readline()
      data=N.append(data, N.asarray(st.split()).astype('f'))
    infile.close()

    t=t[1:] ; data=data[1:]
    nt=len(t)
    nk=len(data)/nt
    data=data.reshape(nt, nk)
    return t, data

