; $Id: ryaver_alp.pro,v 1.1 2005-06-12 19:21:41 brandenb Exp $
;
;  reads the yaver.dat file, puts the result into fmxz array
;  this routine keeps only the last time
;
t=0.
@data/index
default, datatopdir, 'data'
default, dimfile, 'dim.dat'
default, t1, 0.
icount=0
;
mx=0L & my=0L & mz=0L & nvar=0L & naux=0L
prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
nprocx=0L & nprocy=0L & nprocz=0L
close,1
openr,1,datatopdir+'/'+dimfile
readf,1,mx,my,mz,nvar,naux
readf,1,prec
readf,1,nghostx,nghosty,nghostz
readf,1,nprocx,nprocy,nprocz
close,1
;
;  calculate dimension and define input array
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
fmxz=fltarr(nx,nz,nnamexz)
;
;  calculate x and z arrays
;
dx=2*!pi/nx & x=-!pi+dx*(.5+dindgen(nx))
dz=2*!pi/nz & z=-!pi+dz*(.5+dindgen(nz))
;
;  define arrays for final alpha and eta arrays
;
alpijxz1=fltarr(nx,nz,3,3)
etaijkxz1=fltarr(nx,nz,3,3,2)
alpijexz1=fltarr(nx,nz,3,3)
;
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
;  start reading
;
close,1
openr,1,datatopdir+'/yaverages.dat'
;
fo='(8e12.4)'
fo3='(3f10.2)'
default,w,.01
while not eof(1) do begin
  readf,1,t
  readf,1,fmxz,fo=fo
  ;
  ;  put data into alpijxz and etaijkxz arrays
  ;
  alpijxz=reform(fmxz(*,*,i_alp11xz-1:i_alp33xz-1),nx,nz,3,3)
  etaijkxz=reform(fmxz(*,*,i_eta111xz-1:i_eta333xz-1),nx,nz,3,3,2)
  alpijexz=reform(fmxz(*,*,i_alp11exz-1:i_alp33exz-1),nx,nz,3,3)
  ;
  ;  calculate alpha etc
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
  alpijxz1m=total(total(alpijxz1,1),1)/(nx*nz)
  etaijkxz1m=total(total(etaijkxz1,1),1)/(nx*nz)
  alpijexz1m=total(total(alpijexz1,1),1)/(nx*nz)
  ;
  ;  print
  ;
  print
  print,'-----------------------------------'
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
  if icount eq 0 then begin
    alpijxz1mt=alpijxz1m
    etaijkxz1mt=etaijkxz1m
    alpijexz1mt=alpijexz1m
    tt=t
  endif else begin
    alpijxz1mt=[alpijxz1mt,alpijxz1m]
    etaijkxz1mt=[etaijkxz1mt,etaijkxz1m]
    alpijexz1mt=[alpijexz1mt,alpijexz1m]
    tt=[tt,t]
  endelse
  icount=icount+1
  ;
  wait,w
end
close,1
;
;  reform
;
alpijxz1mt=reform(alpijxz1mt,3,icount,3)
etaijkxz1mt=reform(etaijkxz1mt,3,icount,3,2)
alpijexz1mt=reform(alpijexz1mt,3,icount,3)
;
print
print,'alpijxz1mt,etaijkxz1mt,alpijexz1mt'
help,alpijxz1mt,etaijkxz1mt,alpijexz1mt
END
