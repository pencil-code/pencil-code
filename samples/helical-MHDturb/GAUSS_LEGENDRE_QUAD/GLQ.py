#!/usr/bin/python

import numpy as np
nmax=36
precision=8

file='gauss_legendre_quadrature_n'+str(nmax)+'p'+str(precision)+'.dat'
with open(file,'w') as f:
    f.write("%s\n" % nmax)
t='{:.'+str(precision)+'f}'

for i in range(1,nmax+1):
    glq=np.transpose(np.asarray((np.polynomial.legendre.leggauss(i))))
    with open(file, 'a') as f:
        for j in glq:
            for k in j:
                f.write("%s\n" % t.format(k))
        for j in range(2*nmax-2*i):
            f.write("0.\n")
