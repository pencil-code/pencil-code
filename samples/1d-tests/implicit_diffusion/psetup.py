#!/usr/bin/env python

#  $Id$
#  20-may-12/dintrans: coded
#  Plot the initial setup (density, temperature, entropy and radia. cond. K)
#

import pylab as P

z,rho,temp,ss,hcond=P.loadtxt('data/proc0/setup.dat',unpack=True,skiprows=1)

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
