; $Id$
;
;  read global sizes
;
;  calculate x and z arrays
;
dx=2*!pi/nx & x=-!pi+dx*(.5+dindgen(nx))
dz=2*!pi/nz & z=-!pi+dz*(.5+dindgen(nz))
;
;  put data into alpijxz and etaijkxz arrays
;
;bmxz=reform(fmxz(*,*,i_bxmxz-1:i_bzmxz-1),nx,nz,3)
alpijxz=reform(fmxz(*,*,i_alp11xz-1:i_alp33xz-1),nx,nz,3,3)
etaijkxz=reform(fmxz(*,*,i_eta111xz-1:i_eta333xz-1),nx,nz,3,3,2)
;
;  define arrays for final alpha and eta arrays
;
alpijxz1=fltarr(nx,nz,3,3)
etaijkxz1=fltarr(nx,nz,3,3,2)
;
;  prepare matrix for calculating proper alpha and eta tensors
;
sx=sin(x) & cx=cos(x) & ux=replicate(1.,nx)
sz=sin(z) & cz=cos(z) & uz=replicate(1.,nz)
mat=fltarr(nx,nz,3,3)
mat(*,*,0,0)=cx#cz & mat(*,*,0,1)=-sx#cz & mat(*,*,0,2)=-cx#sz
mat(*,*,1,0)=sx#uz & mat(*,*,1,1)=+cx#uz & mat(*,*,1,2)=0.
mat(*,*,2,0)=ux#sz & mat(*,*,2,1)=0.     & mat(*,*,2,2)=+ux#cz
mat1=fltarr(nx,nz,3,3)
;
;  invert
;
for iz=0,nz-1 do begin
for ix=0,nx-1 do begin
  mat1(ix,iz,*,*)=invert(reform(mat(ix,iz,*,*)))
endfor
endfor
;
;  remove alpha part from etaijkxz
;
for j=0,2 do begin
for i=0,2 do begin
  alpijxz1 (*,*,i,j)  =mat1(*,*,0,0)*alpijxz(*,*,i,j)+mat1(*,*,0,1)*etaijkxz(*,*,i,j,0)+mat1(*,*,0,2)*etaijkxz(*,*,i,j,1)
  etaijkxz1(*,*,i,j,0)=mat1(*,*,1,0)*alpijxz(*,*,i,j)+mat1(*,*,1,1)*etaijkxz(*,*,i,j,0)+mat1(*,*,1,2)*etaijkxz(*,*,i,j,1)
  etaijkxz1(*,*,i,j,1)=mat1(*,*,2,0)*alpijxz(*,*,i,j)+mat1(*,*,2,1)*etaijkxz(*,*,i,j,0)+mat1(*,*,2,2)*etaijkxz(*,*,i,j,1)
endfor
endfor
;
;  calculate spatial average
;
fot='(f5.1)'
fo3='(3f10.3)'
alpijxz1m=total(total(alpijxz1,1),1)/(nx*nz)
etaijkxz1m=total(total(etaijkxz1,1),1)/(nx*nz)
print
print,'alpha_ij,  t =',string(t,fo=fot)
print,alpijxz1m,fo=fo3
print
print,'eta_ijx'
print,etaijkxz1m(*,*,0),fo=fo3
print
print,'eta_ijz'
print,etaijkxz1m(*,*,1),fo=fo3
;
END
