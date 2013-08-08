;$Id: p.pro,v 1.1 2012-06-21 03:09:52 brandenb Exp $
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
;  mv idl.ps pparticipants.ps
;
!p.charsize=1.7
!x.margin=[8.8,.5]
!y.margin=[3.2,.5]
!x.title='!6year'
!x.range=[2004.5,2012.5]
w=1.05
;
a=rtable(head=1,'participants.dat',3)
y=reform(a(0,*))
M=reform(a(1,*))
P=reform(a(2,*))
good=where(y ge 2008.)
;
;  ratios and error estimate
;
PM=P+M
r=M/P
rr=M/PM
rm=mean(r(good))
rrm=mean(rr(good))
print,rm,rrm
;
;  relative errors
;
dr=sqrt(M)/M+sqrt(P)/P
drr=sqrt(M)/M+sqrt(PM)/PM
;
;  plot total number of participants
;
!y.title='!6participants'
plot,y,PM,yr=[0,27]
oplot,y,M,col=122
oplot,y,P,col=55
wait,w
;
;  plot fraction of student
;
s=.5 ;(re-scaled relative error)
!y.title='!6fraction of student'
plot,y,rr,li=1,yr=[0,1.],ps=-6
errplot,y,rr*(1.-s*drr),rr*(1.+s*drr)
xx=[2007.5,2012.5] & oplot,xx,xx*0.+rrm
wait,w
;
;  plot students/non-student ratio
;
!y.title='!6students/non-student ratio'
plot,y,M/P,li=1,yr=[0,3.3],ps=-6
print,y
errplot,y,r*(1.-dr),r*(1.+dr)
xx=[2007.5,2012.5] & oplot,xx,xx*0.+rm
;
END
