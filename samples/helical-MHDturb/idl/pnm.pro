;$Id$
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4
end
;
;  reads in time series file
;  mv idl.ps fig/pnm.ps
;
pc_read_ts,obj=ts
plot_io,ts.t,ts.bmx,li=1,yr=[5e-6,1],xtit='t'
oplot,ts.t,ts.bmy,li=2
oplot,ts.t,ts.bmz,li=0
;
xx=250. & dx=100. & siz=1.6
legend,xx,dx,3e-4,1,'x-dir',siz=siz
legend,xx,dx,1e-4,2,'y-dir',siz=siz
legend,xx,dx,3e-5,0,'z-dir',siz=siz
legend
END
