old_multi = !p.multi
!p.multi = [ 1, 1, 1 ]
charsize = 1.25
nbins = 40

pc_read_pvar, object=pvars, varfile='pvar.dat'

x_hist = histogram (pvars.xx[*,0], nbins=nbins, locations=xx)
y_hist = histogram (pvars.xx[*,1], nbins=nbins, locations=yy)
z_hist = histogram (pvars.xx[*,2], nbins=nbins, locations=zz)
total_hist = histogram (pvars.xx[*,*], nbins=nbins, locations=xx) / 3.0

window, 13, title='Particles distribution analysis', retain=2

title = 'x-coordinate (red), y-coordinate (green), z-corrdinate (blue), total (white)'
plot, xx, x_hist, yrange=[0.75,1.1]*max(total_hist), title=title, charsize=charsize, /nodata
oplot, zz, z_hist, col=12810020
oplot, zz, z_hist, col=12810020, psym=2
oplot, yy, y_hist, col=20020
oplot, yy, y_hist, col=20020, psym=2
oplot, xx, x_hist, col=200
oplot, xx, x_hist, col=200, psym=2
oplot, xx, total_hist
oplot, xx, total_hist, psym=2

!p.multi = old_multi

END

