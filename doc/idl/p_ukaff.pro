; $Id$
;
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
;  mv idl.ps ~/boris/Galaxycode/proposals/ukaff/ptimings.ps
;
fact=128.^2*512./1e6
a=rtable('timings.dat',2,head=1)
b=rtable('oldcode.dat',2)
c=rtable('dr.dat',2)
n=reform(a(0,*)) & t=reform(a(1,*))
m=reform(b(0,*)) & s=reform(b(1,*))
k=reform(c(0,*)) & r=reform(c(1,*))
;
!p.multi=[0,1,2]
!p.multi=[0,2,1]
!p.multi=0
!p.charsize=1.6
!x.title='# of procs'
!y.title='!7l!6s/step/point'
!y.title='!6sec/step'
!x.range=[.8,72]
;
plot_oo,n,fact*t,ps=-1,yr=fact*[.1,30]
oplot,m,fact*s,ps=-5,li=1
;oplot,k,fact*r*2,ps=-6,li=2
;
!p.charsize=1.0
!x.range=[.0,72]
!p.position=[.32,.32,.72,.65]
plot,n,fact*t,ps=-1,/noerase,yr=fact*[.0,24],xst=5,yst=5
axis,0,0,xax=0
axis,0,0,yax=0
oplot,m,fact*s,ps=-5,li=1
;
!p.position=0
END
