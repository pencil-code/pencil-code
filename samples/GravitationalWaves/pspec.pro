;
;  Reads spectra and plots GW spectra at k=2 vs time.
;
power,'_kin','_mag',k=k,spec1=kin,spec2=mag,tt=t,/noplo
power,'_GWs','_GWh',k=k,spec1=GWs,spec2=GWh,tt=t,/noplo
power,'hel_GWs','hel_GWh',k=k,spec1=GWshel,spec2=GWhhel,$
  tt=t,/noplot
;
nt=n_elements(t)
plot,k,mag[*,0],yr=[0,20]
for it=1,nt-1 do oplot,k,mag[*,it]
;
;for it=0,nt-1 do oplot,k,kin[*,it]
;for it=0,nt-1 do oplot,k,GWs[*,it],col=188
;for it=0,nt-1 do oplot,k,GWh[*,it],col=155
;for it=0,nt-1 do oplot,k,GWshel[*,it],col=188
;for it=0,nt-1 do oplot,k,GWhhel[*,it],col=155
;
print,'red: S_hdot(k=2,t) = (4pi)^2'
plot,t,GWs[2,*],yr=[0,160]
oplot,t,GWshel[2,*],li=2,col=122
;
fac=1.
print,'blue: S_h(k=2,t) = (4pi)^2'
oplot,t,fac*GWh[2,*]
oplot,t,fac*GWhhel[2,*],li=2,col=55
;
END
