;$Id: add3.tex,v 1.418 2009-07-27 01:19:04 brandenb Exp $
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
;  mv idl.ps
;  mv idl.ps ~/tex/koen/RFP/fig/pfrac.ps
;
!p.charsize=1.6
;
a=rtable('scaling.dat',2,head=1)
N=reform(a(0,*))
pi=reform(a(1,*))
dx=1./(N-1.)
!x.title='!7D!8r!6'
!y.title='!6error of integral'
plot_oo,dx,pi-!pi,xr=[.001,.1],yr=[1e-5,1e-1],ps=6
xx=[.0014,.06] & oplot,xx,10.*xx^2
;
END
