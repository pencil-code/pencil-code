$Id: pcomp.pro,v 1.3 2013/10/20 09:35:46 brandenb Exp $
;
;  sample to compare time series form different directories
;
a=rtable('reference.out',head=1,17)
b=rtable('data/time_series.dat',head=1,17)
;
!p.multi=[0,4,4]
;!x.range=[0,300]
;
;plot,ts1.t,ts1.orms
;oplot,ts2.t,ts2.orms,col=122,li=2

for j=2,16 do begin
  plot,a[1,*],a[j,*],ps=-1
  oplot,a[1,*],a[j,*],ps=-1,li=2,col=122,th=4
endfor
;
END
