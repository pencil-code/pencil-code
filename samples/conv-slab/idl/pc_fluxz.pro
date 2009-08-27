; $Id$
;
;  for plotting 1-D energy fluxes from convection simulations (conv-slab)
;
ii=i_fmassz-1
ii=i_fkinz-1
ii=i_fradz-1
ii=i_fconvz-1
ii=i_uzmz-1
ii=i_dcoolz-1
;
if i_fmassz ne 0 then fmasst=fmzt(*,*,i_fmassz-1)
if i_fturbz ne 0 then fturbt=fmzt(*,*,i_fturbz-1)
if i_dcoolz ne 0 then dcoolt=fmzt(*,*,i_dcoolz-1)
if i_fconvz ne 0 then fconvt=fmzt(*,*,i_fconvz-1)
if i_fradz ne 0 then fradt=fmzt(*,*,i_fradz-1)
if i_fkinz ne 0 then fkint=fmzt(*,*,i_fkinz-1)
;
;  integrate to get cooling flux
;
fcoolt=make_array(size=size(dcoolt),/nozero)
for it=1,nt-1 do begin
  fcoolt(*,it)=-integr(dcoolt(*,it),x=z)
endfor
;
;  plot result for last time
;
ftott=Fradt+Fconvt+Fkint+Fcoolt+Fturbt
fbot=fradt(0,nt-1)
;
plot,zzz,ftott(*,nt-1),yr=[-.5,1.2]*fbot
oplot,zzz,fkint(*,nt-1),col=55
oplot,zzz,fradt(*,nt-1),col=122
oplot,zzz,fconvt(*,nt-1),col=188
oplot,zzz,fturbt(*,nt-1),col=155
;
;  do the same for all times (if continued)
;
wait,.1
fo='(f8.2)'
default,w,.005
for it=0,nt-1 do begin
  plot,zzz,ftott(*,it),yr=[-.5,1.2]*fbot,tit='t ='+string(tt(it),fo=fo)
  oplot,zzz,fkint(*,it),col=55
  oplot,zzz,fradt(*,it),col=122
  oplot,zzz,fconvt(*,it),col=188
  oplot,zzz,fturbt(*,it),col=155
  wait,w
endfor
;
;help,fmzt
;for it=1,nt-1 do begin
;  plot,zzz,fmzt(*,it,ii),yr=[-1,1]*1e-3
;endfor
;
END
