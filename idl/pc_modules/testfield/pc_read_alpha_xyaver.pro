;$Id: pc_read_alpha_xyaver.pro,v 1.1 2008-03-30 16:33:38 brandenb Exp $
;
;  In order to determine the z-dependence of the alpha and eta tensors
;  we have to read the horizontal averages of Epq, i.e. we assume that
;  E111z, ..., E222z has been set in xyaver.in
;
;  30-mar-08/axel
;
pc_read_xyaver,o=xyaver
pc_read_param,o=param,/param2
pc_read_grid,o=grid
pc_read_dim,o=dim
;
tt=xyaver.t
nt=n_elements(tt)
;
;  map z-array to a (-pi,pi) interval
;
nz=dim.nz
dz=2.*!pi/nz
zz=-!pi+dz*(findgen(nz)+.5)
;
;  prepare sine and cosine functions
;
k=param.ktestfield
cz=cos(k*zz)
sz=sin(k*zz)
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
END
