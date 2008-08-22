#!/usr/bin/env python

import pencil as pc
P=pc.P

dim=pc.read_dim()
index=pc.read_index()
param=pc.read_param(quiet=True)
grid=pc.read_grid(trim=True,param=param,quiet=True)

P.ion()
P.figure(figsize=(6,6),dpi=64)
frame=grid.x.min(),grid.x.max(),grid.y.min(),grid.y.max()
P.subplots_adjust(bottom=0,top=1,left=0,right=1)
x0=grid.x.mean()
P.axvline(x0+0.5,color='black',linestyle='--')
P.axvline(x0-0.5,color='black',linestyle='--')
P.axhline(0.5,color='black',linestyle='--')
P.axhline(-0.5,color='black',linestyle='--')

for ivar in range(0,8): 
  print "read VAR%d"%ivar
  var=pc.read_var(ivar=ivar,run2D=param.lwrite_2d,param=param,dim=dim,index=index,quiet=True,trimall=True)
  f=var.lnrho[dim.nz/2,...]

# acceleration using an handle
  if (ivar==0):
    im=P.imshow(f,extent=frame,origin='lower',aspect='auto')
  else:
    im.set_data(f)
    im.set_clim(f.min(),f.max())
    P.draw()

# filename='movie/img%04d.png'%ivar
# print 'write', filename
# P.savefig(filename,dpi=64)

P.show()
