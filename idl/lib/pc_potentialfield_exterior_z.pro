;$Id$
PRO pc_potentialfield_exterior_z,aa,aaa=aaa,zzz=zzz,plo=plo
;
;  calculate potential field in the exterior in the z-direction
;
pc_read_dim,obj=dim,/quiet
l1=dim.l1 & l2=dim.l2
m1=dim.m1 & m2=dim.m2
n1=dim.n1 & n2=dim.n2
nx=dim.nx
ny=dim.ny
nz=dim.nz
;
;  read grid
;
pc_read_grid,obj=grid,/quiet
z=grid.z
;
;  expand the new z mesh by a certain factor (currently =3)
;
nznew=3*dim.nz
aaa=fltarr(dim.mx,dim.my,nznew,3)
zzz=fltarr(nznew)
;
;  determine the locations where the interior should be within
;  the new mesh.
;
nn1=nz
nn2=2*nz-1
aaa(*,*,nn1:nn2,*)=aa(*,*,n1:n2,*)
zzz(nn1:nn2)=z(n1:n2)
;
;  calculate wavenumbers
;
dx=grid.dx & dkx=2.*!pi/(nx*dx)
dy=grid.dy & dky=2.*!pi/(ny*dy)
dz=grid.dz
;
kx=shift(rebin(reform(dkx*(findgen(nx)-(nx-1)/2),nx,1,1),nx,ny,nz),-(nx-1)/2,-(ny-1)/2,-(nz-1)/2)
ky=shift(rebin(reform(dky*(findgen(ny)-(ny-1)/2),1,ny,1),nx,ny,nz),-(nx-1)/2,-(ny-1)/2,-(nz-1)/2)
;
kappa=sqrt(kx^2+ky^2)
;
;  calculate exterior in upper part of domain
;
for j=0,2 do begin
  aaak=fft(aa(l1:l2,m1:m2,n1,j),-1)
  for n=0,nn1-1 do begin
    zrel=(n-nn1)*dz
    aaa(l1:l2,m1:m2,n,j)=fft(aaak*exp(+kappa*zrel),1)
    zzz(n)=z(n1)+zrel
  endfor
endfor
;
;  calculate exterior in upper part of domain
;
for j=0,2 do begin
  aaak=fft(aa(l1:l2,m1:m2,n2,j),-1)
  for n=nn2+1,nznew-1 do begin
    zrel=(n-nn2)*dz
    aaa(l1:l2,m1:m2,n,j)=fft(aaak*exp(-kappa*zrel),1)
    zzz(n)=z(n2)+zrel
  endfor
endfor
;
;  set ghost zones in x and y assuming periodicity
;
for n=0,nznew-1 do begin
  aaa(l1-3:l1-1,*,n,*)=aaa(l2-2:l2-0,*,n,*)
  aaa(l2+1:l2+3,*,n,*)=aaa(l1+0:l1+2,*,n,*)
  aaa(*,m1-3:m1-1,n,*)=aaa(*,m2-2:m2-0,n,*)
  aaa(*,m2+1:m2+3,n,*)=aaa(*,m1+0:m1+2,n,*)
endfor
;
if keyword_set(plo) then begin
  x=grid.x
  contour,reform(aaa(*,3,*,1)),x,zzz,nlev=30
  oplot,x,x-x+.5,col=122
  oplot,x,x-x-.5,col=122
endif
;
END
