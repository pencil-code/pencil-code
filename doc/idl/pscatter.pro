; $Id$
;
window,xs=600,ys=480
!p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
;  mv idl.ps ~/boris/Galaxycode/proposals/ukaff/ptimings.ps
;
fact=1.
a=rtable('horseshoe.dat',2,head=1)
b=rtable('horseshoe_scatter.dat',2,head=1)
c=rtable('kabul.dat',2,head=1)
n=reform(a(0,*)) & t=reform(a(1,*))
m=reform(b(0,*)) & s=reform(b(1,*))
k=reform(c(0,*)) & r=reform(c(1,*))
;
!p.multi=[0,1,2]
!p.multi=[0,2,1]
!p.multi=0
!p.charsize=2.0
!p.title='different layouts'
!x.title='# of procs'
!y.title='!7l!6s/step/point'
!x.range=[.8,160]
;
plot_oo,n,fact*t,ps=-1,yr=fact*[.07,15],li=0,back=255,col=1
oplot,m,fact*s,ps=5,li=0,col=122
oplot,k,fact*r,ps=-6,li=2,col=55
;
;xx=5. & dx=20. & siz=2.0
;legend,xx,dx,10^0.9,1,'!6Origin3000',col=1,siz=siz
;legend,xx,dx,10^0.7,0,'!6Horseshoe',col=122,siz=siz
;legend,xx,dx,10^0.5,2,'!6GB switch',col=55,siz=siz
;
; xx=[1,120] & oplot,col=1,xx,4./xx^.7
print,'import ptimings_scatter.jpg'
print,'scp2 ptimings_scatter.jpg $scr/ccp2001'
END
