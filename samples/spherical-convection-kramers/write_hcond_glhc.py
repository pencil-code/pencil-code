#
#  Primitive python script to generate hcond_glhc.dat for 
#  spherical-convection-kramers sample.
#
#  Author: P. Käpylä (pkaepyl/uni-goettinge[dot]de)
#
import numpy as np
#
# Create x-array
nx=64; x=np.linspace(0.7, 1, nx)
dx=x[1]-x[0]
#
# Compute K such that it would carry total flux at top
Fbot=8.4387005e-07
Ktop=(2./3.)*Fbot # dT/dr = -3, Ftop/Fbot \approx 0.5
#
# tanh profile such that K is non-negligible only near surface
Kprof=Ktop*0.5*(np.tanh((x-.975)/.015)+1.)+1e-12
dKprof=np.gradient(Kprof,dx)
#
hcond_glhc=np.column_stack([Kprof,dKprof])
#
np.savetxt('hcond_glhc.dat', hcond_glhc, delimiter=' ')
