#!/usr/bin/env python

# $Id$
# 19-mai-2009/dintrans: coded
# Plot the time evolution of the acoustic waves:
#   - Fig 1: animation
#   - Fig 2: summary in the plane (z,t)


import numpy as N
import pylab as P
import pencil as pc

# read the slice of uz in the plane (x,z) 
# [equivalent here to the plane (y,z) for this 1-D problem]
f, t = pc.read_slices(field='uu3', proc=0, extension='xz')

dim=pc.read_dim()
nt=len(t)
ff=f.reshape(nt, dim.nz)

P.ion()
P.subplot(211)
line, = P.plot(ff[0, :])
P.xlim(xmax=dim.nz-1)
P.ylim(ymin=ff.min(), ymax=ff.max())
P.title('velocity')
st=P.figtext(0.2, 0.85, 't=%.1f'%t[0])

for i in range(1, nt):
  line.set_ydata(ff[i, :])
  st.set_text('t=%.1f'%t[i])
  P.draw()

P.plot(ff[i, :],'o')
P.xlim(xmax=dim.nz-1)
P.draw()

P.subplot(212)
frame = 0, dim.nz-1, t[0], t[-1]
P.imshow(ff, origin='lower', extent=frame, aspect='auto')
P.ylabel('time')
P.show()

