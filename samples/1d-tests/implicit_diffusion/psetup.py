#!/usr/bin/env python

#  $Id$
#  20-may-12/dintrans: coded
#  Plot the initial setup (density, temperature, entropy and radia. cond. K)
#  Work for both non-parallel and parallel runs in the z-direction
#

import pylab as P
from pencil import read_dim

dim=read_dim()

z=[] ; rho=[] ; temp=[] ; ss=[] ; hcond=[]
for i in range(dim.nprocz):
    z0,rho0,temp0,ss0,hcond0=P.loadtxt('data/proc%i/setup.dat'%i,skiprows=1,unpack=True)
    z=P.hstack((z,z0))
    rho=P.hstack((rho,rho0))
    temp=P.hstack((temp,temp0))
    ss=P.hstack((ss,ss0))
    hcond=P.hstack((hcond,hcond0))

P.rc("lines", linewidth=2)
P.subplot(221)
P.semilogy(z,rho)
P.title('density')

P.subplot(222)
P.plot(z,temp)
P.title('temperature')

P.subplot(223)
P.plot(z,ss)
P.title('entropy')

P.subplot(224)
P.plot(z,hcond)
P.title('radiative conductivity')

P.show()
