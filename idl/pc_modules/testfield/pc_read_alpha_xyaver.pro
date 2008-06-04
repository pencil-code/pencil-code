;$Id: pc_read_alpha_xyaver.pro,v 1.5 2008-06-04 11:25:43 brandenb Exp $
;
;  In order to determine the z-dependence of the alpha and eta tensors
;  we have to read the horizontal averages of Epq, i.e. we assume that
;  E111z, ..., E222z has been set in xyaver.in
;
;  30-mar-08/axel
;
pc_read_xyaver,o=xyaver
pc_read_param,o=param,/param2
;
default,use_grid,1
;
if use_grid eq 1 then begin
  pc_read_grid,o=grid
  pc_read_dim,o=dim
  nz=dim.nz & n1=dim.n1 & n2=dim.n2
  zzz=grid.z(n1:n2)
endif else begin
  ;
  ;  map z-array to a (-pi,pi) interval
  ;
  pc_read_dim,o=dim
  default,nz,dim.nz
  n1=3 & n2=nz+2
  z1=4.*!pi & z0=-z1
  dz=(z1-z0)/nz
  zzz=z0+dz*(.5+findgen(nz))
endelse
;
;  prepare relevant array for testfield sine and cosine functions
;
if param.ltestfield_newz then begin
  dz=2.*!pi/nz
  ztestfield=-!pi+dz*(findgen(nz)+.5)
endif else begin
  ztestfield=zzz
endelse
;
;  time array
;
tt=xyaver.t
nt=n_elements(tt)
;
;  prepare sine and cosine functions
;
k=param.ktestfield
cz=cos(k*ztestfield)
sz=sin(k*ztestfield)
;
k1cz=cz/k
k1sz=sz/k
;
alpij=fltarr(nz,nt,2,2)
etaij=fltarr(nz,nt,2,2)
;
;  These formulae are consistent with B05 (AN 326, 787).
;  They are also consistent with BRRK08 and BRS08, except
;  that there p and q are interchanged.
;
for it=0,nt-1 do begin
  alpij(*,it,0,0)=cz*xyaver.E111z(*,it)+sz*xyaver.E121z(*,it)
  alpij(*,it,0,1)=cz*xyaver.E112z(*,it)+sz*xyaver.E122z(*,it)
  alpij(*,it,1,0)=cz*xyaver.E211z(*,it)+sz*xyaver.E221z(*,it)
  alpij(*,it,1,1)=cz*xyaver.E212z(*,it)+sz*xyaver.E222z(*,it)
;
  etaij(*,it,0,0)=-(k1sz*xyaver.E112z(*,it)-k1cz*xyaver.E122z(*,it))
  etaij(*,it,0,1)=+(k1sz*xyaver.E111z(*,it)-k1cz*xyaver.E121z(*,it))
  etaij(*,it,1,0)=-(k1sz*xyaver.E212z(*,it)-k1cz*xyaver.E222z(*,it))
  etaij(*,it,1,1)=+(k1sz*xyaver.E211z(*,it)-k1cz*xyaver.E221z(*,it))
endfor
;
print,'tvscl,alpij(*,*,0,0)'
;
alpijm=total(alpij,2)/nt
etaijm=total(etaij,2)/nt
;
save,file='alpetaij.sav',zzz,alpijm,etaijm
END
