;$Id$
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4
end
;
;  reads in time series file
;
!p.charsize=1.6
pc_read_ts,obj=ts
plot_io,ts.t,ts.brms,li=1,yr=[5e-5,1],xtit='t',ytit='urms, brms'
oplot,ts.t,ts.urms,li=0

END
