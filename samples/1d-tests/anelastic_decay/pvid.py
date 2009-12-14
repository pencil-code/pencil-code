#!/usr/bin/python

# $Id$

import numpy as N
import pylab as P
import pencil as pc
from os import system

system('cat video.in')
field=raw_input('which field? ')
f, t = pc.read_slices(field=field, proc=0, extension='xy')

dim=pc.read_dim()
grid=pc.read_grid(trim=True)
nt=len(t)
ff=f.reshape(nt, dim.nx)

P.ion()
P.subplot(211)
line, = P.plot(grid.x, ff[0, :])
P.xlim(grid.x[0], grid.x[-1])
P.ylim(ymin=ff.min(), ymax=ff.max())
P.title(field)
st=P.figtext(0.2, 0.85, 't=%.1f'%t[0])

for i in range(1, nt):
  line.set_ydata(ff[i, :])
  st.set_text('t=%.1f'%t[i])
  P.draw()

P.plot(grid.x, ff[0, :],'o')
P.xlim(grid.x[0], grid.x[-1])
P.draw()

P.subplot(212)
frame = grid.x[0], grid.x[-1], t[0], t[-1]
P.imshow(ff, origin='lower', extent=frame, aspect='auto')
P.ylabel('time')
P.show()

