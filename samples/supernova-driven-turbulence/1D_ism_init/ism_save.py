import pencil as pc
import numpy as np

varfile = 'var.dat' # or specific snaphot as required 'VAR?'
var=pc.read_var(varfile,magic=['tt'],trimall=True,quiet=True)

filename='init_ism.dat'
f = open(filename, 'w')
#f.write('#--rho-------------------TT-------------\n')
for i in range(0,var.rho[:,0,0].size):
    f.write(str(var.rho[i,0,0])+'    '+str(var.tt[i,0,0])+'\n')
f.closed




