; $Id$
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
;  mean vorticity
;
Wx=+cx#sz
Wz=-sx#cz
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
  ; alpha tensor
  ;
  alphayy=alpijxz1(*,*,1,1)
  alpheyy=alpijexz1(*,*,1,1)
  ;
  ;  project gamma vector against vorticity
  ;
  gammax=.5*(alpijxz1(*,*,2,1)-alpijxz1(*,*,1,2))
  gammay=.5*(alpijxz1(*,*,0,2)-alpijxz1(*,*,2,0))
  gammaz=.5*(alpijxz1(*,*,1,0)-alpijxz1(*,*,0,1))
  gammaW=gammax*Wx+gammaz*Wz
  ;
  gammex=.5*(alpijexz1(*,*,2,1)-alpijexz1(*,*,1,2))
  gammey=.5*(alpijexz1(*,*,0,2)-alpijexz1(*,*,2,0))
  gammez=.5*(alpijexz1(*,*,1,0)-alpijexz1(*,*,0,1))
  gammeW=gammex*Wx+gammez*Wz
  ;
  ;  project delta vector against vorticity
  ;
  etaxx=+.5*etaijkxz1(*,*,0,1,1)
  etazz=-.5*etaijkxz1(*,*,2,1,0)
  etayy=.5*(etaijkxz1(*,*,1,2,0)-etaijkxz1(*,*,1,0,1))
  ;
  deltax=.25*(etaijkxz1(*,*,2,0,1)-etaijkxz1(*,*,1,1,0)-etaijkxz1(*,*,2,2,0))
  deltaz=.25*(etaijkxz1(*,*,0,2,0)-etaijkxz1(*,*,0,0,1)-etaijkxz1(*,*,1,1,1))
  deltaW=deltax*Wx+deltaz*Wz
  ;
  ;  total 2-D averages
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
    gammeyt=gammey
    gammeWt=gammeW
    gammayt=gammay
    gammaWt=gammaW
    deltaWt=deltaW
    etaxxt=etaxx
    etayyt=etayy
    etazzt=etazz
    alphayyt=alphayy
    alpheyyt=alpheyy
  endif else begin
    alpijxz1mt=[alpijxz1mt,alpijxz1m]
    etaijkxz1mt=[etaijkxz1mt,etaijkxz1m]
    alpijexz1mt=[alpijexz1mt,alpijexz1m]
    tt=[tt,t]
    gammayt=[gammayt,gammay]
    gammaWt=[gammaWt,gammaW]
    gammeyt=[gammeyt,gammey]
    gammeWt=[gammeWt,gammeW]
    deltaWt=[deltaWt,deltaW]
    etaxxt=[etaxxt,etaxx]
    etayyt=[etayyt,etayy]
    etazzt=[etazzt,etazz]
    alphayyt=[alphayyt,alphayy]
    alpheyyt=[alpheyyt,alpheyy]
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
;
deltax=.25*(etaijkxz1mt(2,*,0,1)-etaijkxz1mt(2,*,2,0))
deltaz=.25*(etaijkxz1mt(0,*,2,0)-etaijkxz1mt(0,*,0,1))
;
gammayt=reform(gammayt,nx,icount,nz)
gammaWt=reform(gammaWt,nx,icount,nz)
gammeyt=reform(gammeyt,nx,icount,nz)
gammeWt=reform(gammeWt,nx,icount,nz)
deltaWt=reform(deltaWt,nx,icount,nz)
etaxxt=reform(etaxxt,nx,icount,nz)
etayyt=reform(etayyt,nx,icount,nz)
etazzt=reform(etazzt,nx,icount,nz)
alphayyt=reform(alphayyt,nx,icount,nz)
alpheyyt=reform(alpheyyt,nx,icount,nz)
;
gammaytm=total(total(gammayt,1),2)/(nx*nz)
gammaWtm=total(total(gammaWt,1),2)/(nx*nz)
gammeytm=total(total(gammeyt,1),2)/(nx*nz)
gammeWtm=total(total(gammeWt,1),2)/(nx*nz)
deltaWtm=total(total(deltaWt,1),2)/(nx*nz)
etaxxtm=total(total(etaxxt,1),2)/(nx*nz)
etayytm=total(total(etayyt,1),2)/(nx*nz)
etazztm=total(total(etazzt,1),2)/(nx*nz)
alphayytm=total(total(alphayyt,1),2)/(nx*nz)
alpheyytm=total(total(alpheyyt,1),2)/(nx*nz)
;
save,file='ryaver_alp.sav',tt,gammaytm,gammaWtm,gammeytm,gammeWtm,deltaWtm,etaxxtm,etayytm,etazztm,alphayytm,alpheyytm,alpijxz1mt,etaijkxz1mt,alpijexz1mt
END
