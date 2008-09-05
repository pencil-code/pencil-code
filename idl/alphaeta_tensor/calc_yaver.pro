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
alpijexz=reform(fmxz(*,*,i_alp11exz-1:i_alp33exz-1),nx,nz,3,3)
;
;  define arrays for final alpha and eta arrays
;
alpijxz1=fltarr(nx,nz,3,3)
etaijkxz1=fltarr(nx,nz,3,3,2)
alpijexz1=fltarr(nx,nz,3,3)
;
;  prepare matrix for calculating proper alpha and eta tensors
;
sx=sin(x) & cx=cos(x)
sz=sin(z) & cz=cos(z)
mat=fltarr(nz,2,2) & mat1=fltarr(nz,2,2) & mat(*,0,0)=+sz & mat(*,0,1)=+cz & mat(*,1,0)=+cz & mat(*,1,1)=-sz
matx=fltarr(nx,2,2) & matx1=fltarr(nx,2,2) & matx(*,0,0)=+sx & matx(*,0,1)=+cx & matx(*,1,0)=+cx & matx(*,1,1)=-sx
;
;  invert
;
for iz=0,nz-1 do begin
  mat1(iz,*,*)=invert(reform(mat(iz,*,*)))
endfor
;
for ix=0,nx-1 do begin
  matx1(ix,*,*)=invert(reform(matx(ix,*,*)))
endfor
;
;  remove alpha part from etaijkxz
;
for j=0,2 do begin
for i=0,2 do begin
for l=0,nx-1 do begin
  alpijxz1(l,*,i,j)=mat1(*,0,0)*alpijxz(l,*,i,j)+mat1(*,0,1)*etaijkxz(l,*,i,j,1)
  etaijkxz1(l,*,i,j,1)=mat1(*,1,0)*alpijxz(l,*,i,j)+mat1(*,1,1)*etaijkxz(l,*,i,j,1)
  ;etaijkxz1(l,*,i,j,0)=(etaijkxz(l,*,i,j,0)-alpijxz1(l,*,i,j)*cos(x(l)))/(-sin(x(l)))
endfor
endfor
endfor
;
;  remove alpha part from etaijkxz
;
for j=0,2 do begin
for i=0,2 do begin
for n=0,nz-1 do begin
  alpijexz1(*,n,i,j)=matx1(*,0,0)*alpijexz(*,n,i,j)+matx1(*,0,1)*etaijkxz(*,n,i,j,0)
  etaijkxz1(*,n,i,j,0)=matx1(*,1,0)*alpijexz(*,n,i,j)+matx1(*,1,1)*etaijkxz(*,n,i,j,0)
endfor
endfor
endfor
;
;  calculate spatial average
;
fot='(f5.1)'
fo3='(3f10.3)'
alpijxz1m=total(total(alpijxz1,1),1)/(nx*nz)
etaijkxz1m=total(total(etaijkxz1,1),1)/(nx*nz)
alpijexz1m=total(total(alpijexz1,1),1)/(nx*nz)
print
print,'alpha_ij,  t =',string(t,fo=fot)
print,alpijxz1m,fo=fo3
print
print,'eta_ijx'
print,etaijkxz1m(*,*,0),fo=fo3
print
print,'eta_ijz'
print,etaijkxz1m(*,*,1),fo=fo3
print
print,'alphae_ij'
print,alpijexz1m,fo=fo3
;
END
