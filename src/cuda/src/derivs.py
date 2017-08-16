import numpy as np

#
# Derivative operations 
#

def xder(f,dx):
   # 6th order x dervative 
   dx2 = 1./(60.*dx)
   dfdx = np.zeros_like(f)
   l1 = 3
   l2 = f.shape[0]-3
   dfdx[l1:l2,...] = dx2*( +45.*(f[l1+1:l2+1,...]-f[l1-1:l2-1,...])
                           -9.*(f[l1+2:l2+2,...]-f[l1-2:l2-2,...])
                           +(f[l1+3:l2+3,...]-f[l1-3:l2-3,...]) )
   return dfdx

def yder(f,dy):
   # 6th order y dervative 
   dy2 = 1./(60.*dy)
   dfdy = np.zeros_like(f)
   m1 = 3
   m2 = f.shape[1]-3
   dfdy[:,m1:m2,...] = dy2*( +45.*(f[:,m1+1:m2+1,...]-f[:,m1-1:m2-1,...]) 
                             -9.*(f[:,m1+2:m2+2,...]-f[:,m1-2:m2-2,...]) 
                             +(f[:,m1+3:m2+3,...]-f[:,m1-3:m2-3,...]) ) 
   return dfdy

def zder(f,dz):
   # 6th order z dervative 
   dz2 = 1./(60.*dz)
   dfdz = np.zeros_like(f)
   n1 = 3
   n2 = f.shape[2]-3
   dfdz[:,:,n1:n2,...] = dz2*( +45.*(f[:,:,n1+1:n2+1,...]-f[:,:,n1-1:n2-1,...])
                               -9.*(f[:,:,n1+2:n2+2,...]-f[:,:,n1-2:n2-2,...]) 
                               +(f[:,:,n1+3:n2+3,...]-f[:,:,n1-3:n2-3,...]) )
   return dfdz


def curl(gridx, gridy, gridz, dx, dy, dz):
   # Calculate a curl in the array (adapted from pecin code python stuff)
   print "Curl"
   curl = np.empty_like(gridx)
   curlx = yder(gridz,dy) - zder(gridy,dz)
   curly = zder(gridx,dz) - xder(gridz,dx)
   curlz = xder(gridy,dx) - yder(gridx,dy)
   print "Done. \n"
 
   return curlx, curly, curlz

def dot(gridAx, gridAy, gridAz, gridBx, gridBy, gridBz):
   # Calculate a dot product
   dotp = np.zeros_like(gridAx)
   dotp = gridAx*gridBx + gridAy*gridBy + gridAz*gridBz

   return dotp
   


