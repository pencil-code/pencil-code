
import numpy as np

nlambda    = 1000
lambda_init  = 0.1
lambda_final = 10.

nstars=1
Teff=4000
radius_star = 2.0
mass_star = 1.0

rsun=6.96e10
msun=1.98892e33

iformat=2
position=np.zeros(3)

f = open('stars.inp','w')

f.write('%d \n'%(iformat))
f.write('%d %d \n'%(nstars,nlambda))
f.write('%13.7e %13.7e %13.7e %13.7e %13.7e \n\n'%(radius_star*rsun, mass_star*msun,position[0],position[1],position[2]))

wavelengths = np.logspace(np.log10(lambda_init),np.log10(lambda_final),nlambda)
for i in range(nlambda):
    f.write('%13.7e\n'%(wavelengths[i]))

f.write('\n%13.6e\n'%(-1*Teff))
    
f.close()
