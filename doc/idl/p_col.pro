; $Id$
;
!p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
  !p.charsize=1.6
  siz=1.6
endif else begin
  window,xs=600,ys=480
  !p.charsize=2.0
  siz=2.0
end
;
;  mv idl.ps ~/boris/Galaxycode/proposals/ukaff/ptimings.ps
;  mv idl.ps ../figs/ptimings.ps
;
fact=1.
a=rtable('timings.dat',2,head=1)
b=rtable('horseshoe.dat',2,head=1)
c=rtable('kabul.dat',2,head=1)
d=rtable('horseshoe_mega.dat',2,head=1)
e=rtable('horseshoe_giga2.dat',2,head=1)
n=reform(a(0,*)) & t=reform(a(1,*))
m=reform(b(0,*)) & s=reform(b(1,*))
k=reform(c(0,*)) & r=reform(c(1,*))
l=reform(d(0,*)) & q=reform(d(1,*))
l2=reform(e(0,*)) & q2=reform(e(1,*))
;
!p.multi=[0,1,2]
!p.multi=[0,2,1]
!p.multi=0
!x.title='# of procs'
!y.title='!7l!6s/step/point'
!x.range=[.8,400]
!y.range=[.01,15]*fact
;
fact2=fact*1.3/1.7
fact2=fact
plot_oo,n,fact*t,ps=-1,li=1,back=255,col=1
oplot,m,fact*s,ps=-5,li=0,col=122
oplot,k,fact*r,ps=-6,li=2,col=55
oplot,l,fact*q,ps=-6,li=3,col=166
oplot,l2,fact2*q2,ps=-6,li=4,col=100
;
xx=[8,200] & oplot,xx,5./xx,col=1,thick=1
;
xx=14. & dx=20. & xx2=7.
legend,xx,dx,10^0.9,1,'!3Origin3000',col=1,siz=siz
legend,xx,dx,10^0.7,0,'!3Horseshoe',col=122,siz=siz
legend,xx,dx,10^0.5,2,'!3Kabul',col=55,siz=siz
legend,xx,dx,10^0.3,3,'!3Giga-queue',col=166,siz=siz
legend,xx,dx,10^0.1,4,'!3Giga2-queue',col=100,siz=siz
;
; xx=[1,120] & oplot,col=1,xx,4./xx^.7
print,'import ptimings.png'
print,'mv ptimings.png ~/MyPictures/PencilCode/'
END
