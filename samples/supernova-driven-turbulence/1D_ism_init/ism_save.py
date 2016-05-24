import pencil as pc
import numpy as np

varfile = 'var.dat' # or specific snaphot as required 'VAR?'
var=pc.read_var(varfile,magic=['tt'],trimall=True,quiet=True)
param=pc.read_param(quiet=True)

filename='init_ism.dat'
f = open(filename, 'w')
#smooth and ensure symmetric about midplane - assumes centred
#convert to cgs - so units can be applied independently
rho = (var.rho[:,0,0]+var.rho[::-1,0,0])/2*param.unit_density
tt  = (var.tt[:,0,0]+var.tt[::-1,0,0])/2*param.unit_temperature
#f.write('#--rho-------------------TT-------------\n')
for i in range(0,var.rho[:,0,0].size):
    f.write(str(rho[i])+'    '+str(tt[i])+'\n')
f.closed




