;$Id$
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
;
;  In order to determine the z-dependence of the alpha and eta tensors
;  we have to read the horizontal averages of Epq, i.e. we assume that
;  E111z, ..., E222z has been set in xyaver.in
;
;  Note: if you forgot to say use_grid=0 first time round, you may need
;  to exit and run this again.
;
;  30-mar-08/axel
;
@xder_6th
pc_read_xyaver,o=xyaver
pc_read_param,o=param,/param2
pc_read_param,o=param1
z0=param1.xyz0(2)
z1=param1.xyz1(2)
Lz=z1-z0
;
default,use_grid,1
@data/testfield_info.dat
spawn,'touch parameters.pro'
@parameters
;
;  read 
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
  ;n1=3 & n2=nz+2
  ;z1=!pi & z0=-z1
  ;dz=(z1-z0)/nz
  ;zzz=z0+dz*(.5+findgen(nz))
  dz=Lz/(nz-1)
  print,'assume non-periodic domain with dz=',dz
  zzz=z0+dz*findgen(nz)
endelse
;
;  check whether or not we have sinusoidal test fields
;
if param.itestfield eq 'B11-B22_lin' then sinusoidal=0
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
help,nt
;
;  prepare sine and cosine functions
;  Note that param.ktestfield is only correct if Lz=2*pi.
;  Multiply param.ktestfield by 2.*!pi/Lz, which works in all cases.
;
k=param.ktestfield
kscaled=k*2.*!pi/Lz
cz=cos(k*ztestfield)
sz=sin(k*ztestfield)
zz=ztestfield
;
k1cz=cz/kscaled
k1sz=sz/kscaled
;
alpij=fltarr(nz,nt,2,2)
etaij=fltarr(nz,nt,2,2)
;
;  These formulae are consistent with B05 (AN 326, 787).
;  They are also consistent with BRRK08 and BRS08, except
;  that there p and q are interchanged.
;
default,sinusoidal,1
if sinusoidal eq 1 then begin
  print,'assuming sinusoidal test-fields'
  for it=0L,nt-1L do begin
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
endif else begin
  print,'assuming linear test-fields'
  for it=0L,nt-1L do begin
    alpij(*,it,0,0)=   xyaver.E111z(*,it)
    alpij(*,it,0,1)=   xyaver.E112z(*,it)
    alpij(*,it,1,0)=   xyaver.E211z(*,it)
    alpij(*,it,1,1)=   xyaver.E212z(*,it)
  ;
    etaij(*,it,0,0)=+(-zz*xyaver.E112z(*,it)+xyaver.E122z(*,it))
    etaij(*,it,0,1)=-(-zz*xyaver.E111z(*,it)+xyaver.E121z(*,it))
    etaij(*,it,1,0)=+(-zz*xyaver.E212z(*,it)+xyaver.E222z(*,it))
    etaij(*,it,1,1)=-(-zz*xyaver.E211z(*,it)+xyaver.E221z(*,it))
  endfor
endelse
;
print,'tvscl,alpij(*,*,0,0)'
print,'contour,transpose(alpij(*,*,1,1)),tt(0:*),zzz,nlev=20,/fil'
print,'contour,transpose(alpijc(*,*,1,1)),ttt,zzz,nlev=20,/fil'
;
;tarray=spread(xyaver.t,[0,2,3],[nz,2,2])
;default,good,where(tarray gt 0.)
;spawn,'touch good.pro'
;@good
;ntgood=n_elements(tarray(good))/(nz*2*2)
;
default,it1,0
ntgood=n_elements(xyaver.t(it1:*))
;
alpijmt=total(alpij(*,*,*,*),1)/nz
etaijmt=total(etaij(*,*,*,*),1)/nz
;
bxmz=total(xyaver.bxmz,2)/ntgood
bymz=total(xyaver.bymz,2)/ntgood
;
jxmz=total(xyaver.jxmz,2)/ntgood
jymz=total(xyaver.jymz,2)/ntgood
;
alpijm=total(alpij(*,it1:*,*,*),2)/ntgood
etaijm=total(etaij(*,it1:*,*,*),2)/ntgood
;
alpij_end=reform(alpij(*,nt-1,*,*))
etaij_end=reform(etaij(*,nt-1,*,*))
print,tt(nt-1),nt
tt=xyaver.t
;
;  coarsegrain data by factor nevery
;
default,nevery,40
ntout=nt/nevery
print,'coarsegrain data: nevery,nt,ntout=',nevery,nt,ntout
print,'nz=',nz
it2=nt-1
it1=nt-nevery*ntout
;
;  coarsegrain using averaging
;
pc_coarsegrain,tt,alpij,nevery,ttt,alpijc,/aver
pc_coarsegrain,tt,etaij,nevery,ttt,etaijc,/aver
save,file='alpetaijc.sav',zzz,alpijc,etaijc,alpij_end,etaij_end,tt,ttt
save,file='alpetaijm.sav',zzz,alpijm,etaijm,bxmz,bymz,jxmz,jymz,tt,ttt
;
END
