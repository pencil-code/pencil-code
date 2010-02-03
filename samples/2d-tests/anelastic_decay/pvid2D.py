#!/usr/bin/python

# $Id: pvid2D.py,v 1.1.1.1 2009-12-16 17:37:16 dintrans Exp $

import numpy as N
import pylab as P
import pencil as pc
from os import system

system('cat video.in')
field=raw_input('which field? ')
f, t = pc.read_slices(field=field, proc=0, extension='xz')
ux, t = pc.read_slices(field='uu1', proc=0, extension='xz')
uz, t = pc.read_slices(field='uu3', proc=0, extension='xz')

dim=pc.read_dim()
grid=pc.read_grid(trim=True)
param=pc.read_param(quiet=True)
nt=len(t)
f=f.reshape(nt, dim.nz, dim.nx)
ux=ux.reshape(nt, dim.nz, dim.nx)
uz=uz.reshape(nt, dim.nz, dim.nx)

P.ion()
frame=param.xyz0[0],param.xyz1[0],param.xyz0[2],param.xyz1[2]
qs1=N.random.random_integers(0,dim.nx-1, 1000)
qs2=N.random.random_integers(0,dim.nz-1, 1000)
xx,zz=P.meshgrid(grid.x, grid.z)

im=P.imshow(f[0,...], extent=frame, origin='lower', aspect='auto')
a=ux[0, qs2, qs1]**2+uz[0, qs2, qs1]**2 ; norm=N.sqrt(a.max())
ux[0, qs2, qs1]=ux[0, qs2, qs1]/norm
uz[0, qs2, qs1]=uz[0, qs2, qs1]/norm
vel=P.quiver(xx[qs2, qs1], zz[qs2, qs1], ux[0, qs2, qs1],
uz[0, qs2, qs1], scale=7)
st=P.figtext(0.8,0.2,'t=%.1f'%t[0], color='w')
P.xlabel('x')
P.ylabel('z')

for i in range(1, nt):
    im.set_data(f[i, ...])
    a=ux[i, qs2, qs1]**2+uz[i, qs2, qs1]**2 ; norm=N.sqrt(a.max())
    ux[i, qs2, qs1]=ux[i, qs2, qs1]/norm
    uz[i, qs2, qs1]=uz[i, qs2, qs1]/norm
    vel.set_UVC(ux[i, qs2, qs1], uz[i, qs2, qs1])
    st.set_text('t=%.1f'%t[i])
    P.draw()

P.show()

