power,'_kin','_mag',k=k,spec1=spec1,spec2=spec2,i=n,tt=t,/noplot
nt=n-1
;
plot_oo,k[1:*],spec2[1:*,0]
;
for it=0,nt-1 do begin
  oplot,k[1:*],spec2[1:*,it],col=122
  oplot,k[1:*],spec1[1:*,it],col=55
  wait,.4
endfor
;
END
