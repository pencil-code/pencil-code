old_multi = !p.multi
!p.multi = [ 1, 1, 1 ]
charsize = 1.25

pc_read_ts, obj=ts

window, 11, title='Energy conservation analysis', retain=2

title = 'ethm (blue), ekin (green), ekinp (red), sum (white)'
plot, ts.t, ts.ethm, yrange=[-1.0,1.0]*2.0e-3, title=title, charsize=charsize, /nodata
oplot, ts.t, ts.ethm-0.9, col=12810020
oplot, ts.t, ts.ethm-0.9, col=12810020, psym=2
oplot, ts.t, ts.ekin, col=20020
oplot, ts.t, ts.ekin, col=20020, psym=2
oplot, ts.t, ts.ekinp, col=200
oplot, ts.t, ts.ekinp, col=200, psym=2
oplot, ts.t, ts.ethm+ts.ekin+ts.ekinp-0.9
oplot, ts.t, ts.ethm+ts.ekin+ts.ekinp-0.9, psym=2

!p.multi = old_multi

END

