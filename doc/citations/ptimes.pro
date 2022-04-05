;$Id: add3.tex,v 1.502 2015/02/19 06:35:19 brandenb Exp $
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
;  mv idl.ps fig/ptimes.ps
;  convert fig/ptimes.ps /D/Print/PCUserMeeting2017/ptimes.png
;  convert fig/ptimes.ps /D/Print/PC18/ptimes.png
;
!p.charsize=1.7
!x.margin=[6.1,.5]
!y.margin=[3.2,.5]
!x.title='!6year'
!y.title='!6number of papers'
;
a=rtable('times.txt',4,head=1)
n=a(0,*)
y=a(1,*)
c=a(2,*)
o=a(3,*)
print,n
print
print,total(n)
good=where(y ge 2011.)
nm=mean(n(good))
cm=mean(c(good))
om=mean(o(good))
;
xr=[2002.8,2022.7]
plot,y,n,ps=10,yr=[0,64],xr=xr
oplot,y,c,ps=10,col=122
oplot,y,o,ps=10,col=55
;
ygood=[min(y(good)),xr[1]]
oplot,ygood,ygood*0+nm
oplot,ygood,ygood*0+cm,col=122
oplot,ygood,ygood*0+om,col=55
;
print,'total(n)=',fix(total(n))
print,'total(c)=',fix(total(c))
print,'total(o)=',fix(total(o)),' ',nint(100*total(o)/total(n)),'%'
;
siz=1.7
xyouts,2003.4,58,'w/o Brandenburg',col=55,siz=siz*.9
xyouts,2003.4,52,'comp & ref papers',col=122,siz=siz*.8
;
print,"$mv idl.ps fig/ptimes.ps"
END
