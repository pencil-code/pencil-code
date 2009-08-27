; $Id$
;
;  read global sizes
;
default, datatopdir, 'data'
default, dimfile, 'dim.dat'
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
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
;
;  calculate z array
;
dz=2*!pi/nz
z=-!pi+dz*(.5+dindgen(nz))
;
;  reads the xyaver.dat file
;
t=0.
bxmz=fltarr(nz)
bymz=fltarr(nz)
alpijz=fltarr(nz,3,3)
etaijkz=fltarr(nz,3,3)
alpijz1=fltarr(nz,3,3)
etaijkz1=fltarr(nz,3,3)
;
;  prepare matrix for calculating proper alpha and eta tensors
;
sz=sin(z)
cz=cos(z)
mat=dblarr(nz,2,2)
mat1=dblarr(nz,2,2)
mat(*,0,0)=+sz
mat(*,0,1)=+cz
mat(*,1,0)=+cz
mat(*,1,1)=-sz
;
;  invert
;
for iz=0,nz-1 do begin
  mat1(iz,*,*)=invert(reform(mat(iz,*,*)))
endfor
;
;  read data
;
close,1
openr,1,datatopdir+'/xyaverages.dat'
print,'read from file: ',datatopdir+'/xyaverages.dat'
;
it=0
fo='(8e12.5)'
default,w,.01
default,tmax,1e33
while not eof(1) do begin
  readf,1,t & print,t
  readf,1,bxmz,bymz,alpijz,etaijkz,fo=fo
  if (t ge tmax) then goto,ending
  ;
  ;  remove alpha part from etaijkz
  ;
  for j=0,2 do begin
  for i=0,2 do begin
    alpijz1(*,i,j)=mat1(*,0,0)*alpijz(*,i,j)+mat1(*,0,1)*etaijkz(*,i,j)
    etaijkz1(*,i,j)=mat1(*,1,0)*alpijz(*,i,j)+mat1(*,1,1)*etaijkz(*,i,j)
  endfor
  endfor
  ;
  ;  print mean tensor on console
  ;
  alpijz1m=total(alpijz1,1)/nz
  etaijkz1m=total(etaijkz1,1)/nz
  ;
  if it eq 0 then begin
    bxmzt=bxmz
    bymzt=bymz
    alpijz1t=alpijz1
    etaijkz1t=etaijkz1
    tt=t
  endif else begin
    bxmzt=[bxmzt,bxmz]
    bymzt=[bymzt,bymz]
    alpijz1t=[alpijz1t,alpijz1]
    etaijkz1t=[etaijkz1t,etaijkz1]
    tt=[tt,t]
  endelse
  it=it+1
end
ending:
close,1
;
nt=n_elements(tt)
bxmzt=reform(bxmzt,nz,nt)
bymzt=reform(bymzt,nz,nt)
alpijz1t=reform(alpijz1t,nz,nt,3,3)
etaijkz1t=reform(etaijkz1t,nz,nt,3,3)
;
END
