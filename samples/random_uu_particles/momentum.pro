old_multi = !p.multi
!p.multi = [ 3, 3, 1 ]
charsize = 2.5

pc_read_const, obj=cst
pc_read_dim, obj=dim
pc_read_pdim, obj=pdim
pc_read_ts, obj=ts

rhopm = cst.rhop_swarm*pdim.npar/double(dim.nx*dim.ny*dim.nz)

window, 12, title='Momentum conservation analysis', retain=2

title = 'ruxm (green), vpxm (red), sum (white)'
plot, ts.t, ts.ruxm, yrange=[-1.0,1.0]*5.0e-5, title=title, charsize=charsize, /nodata
oplot, ts.t, ts.ruxm, col=20020
oplot, ts.t, ts.ruxm, col=20020, psym=2
oplot, ts.t, ts.vpxm*rhopm, col=200
oplot, ts.t, ts.vpxm*rhopm, col=200, psym=2
oplot, ts.t, ts.ruxm+ts.vpxm*rhopm
oplot, ts.t, ts.ruxm+ts.vpxm*rhopm, psym=2
title = 'ruym (green), vpym (red), sum (white)'
plot, ts.t, ts.ruym, yrange=[-1.0,1.0]*2.0e-4, title=title, charsize=charsize, /nodata
oplot, ts.t, ts.ruym, col=20020
oplot, ts.t, ts.ruym, col=20020, psym=2
oplot, ts.t, ts.vpym*rhopm, col=200
oplot, ts.t, ts.vpym*rhopm, col=200, psym=2
oplot, ts.t, ts.ruym+ts.vpym*rhopm
oplot, ts.t, ts.ruym+ts.vpym*rhopm, psym=2
title = 'ruzm (green), vpzm (red), sum (white)'
plot, ts.t, ts.ruzm, yrange=[-1.0,1.0]*1.0e-4, title=title, charsize=charsize, /nodata
oplot, ts.t, ts.ruzm, col=20020
oplot, ts.t, ts.ruzm, col=20020, psym=2
oplot, ts.t, ts.vpzm, col=200
oplot, ts.t, ts.vpzm, col=200, psym=2
oplot, ts.t, ts.ruzm+ts.vpzm*rhopm
oplot, ts.t, ts.ruzm+ts.vpzm*rhopm, psym=2

!p.multi = old_multi

END

