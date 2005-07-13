;
;  $Id: pc_particles_to_density.pro,v 1.8 2005-07-13 14:15:44 ajohan Exp $
;
;  Convert positions of particles to a number density field.
;
;  Author: Anders Johansen
;
function pc_particles_to_density, xxp, x, y, z, vvp=vvp, sigmap=sigmap

COMMON pc_precision, zero, one

if (n_elements(vvp) ne 0) then lsigma=1

pc_set_precision
if (n_elements(x) gt 1) then dx=x[1]-x[0] else dx=1.0*one
if (n_elements(y) gt 1) then dy=y[1]-y[0] else dy=1.0*one
if (n_elements(z) gt 1) then dz=z[1]-z[0] else dz=1.0*one
dx1=1.0D/dx & dy1=1.0D/dy & dz1=1.0D/dz

npar=0L

npar=n_elements(xxp[*,0])
nx=n_elements(x)
ny=n_elements(y)
nz=n_elements(z)

np=fltarr(nx,ny,nz)*one

if (lsigma) then begin
  vvpm  =fltarr(nx,ny,nz,3)*one
  vvp2m =fltarr(nx,ny,nz,3)*one
  sigmap=fltarr(nx,ny,nz,3)*one
endif

for k=0L,npar-1 do begin
  
  ix = round((xxp[k,0]-x[0])*dx1)
  iy = round((xxp[k,1]-y[0])*dy1)
  iz = round((xxp[k,2]-z[0])*dz1)
  if (ix eq nx) then ix=ix-1
  if (iy eq ny) then iy=iy-1
  if (iz eq nz) then iz=iz-1
  if (ix eq -1) then ix=0
  if (iy eq -1) then iy=0
  if (iz eq -1) then iz=0
  np[ix,iy,iz]=np[ix,iy,iz]+1.0*one
;  Velocity dispersion
  if (lsigma) then begin
    vvpm[ix,iy,iz,*]=vvpm[ix,iy,iz,*]+vvp[k,*]
    vvp2m[ix,iy,iz,*]=vvp2m[ix,iy,iz,*]+vvp[k,*]^2
  endif

endfor

if (lsigma) then begin
  for l=0,nx-1 do begin & for m=0,ny-1 do begin & for n=0,nz-1 do begin
;  Divide by number of particles
    if (np[l,m,n] gt 1.0) then begin
      vvpm[l,m,n,*]=vvpm[l,m,n,*]/np[l,m,n]
      vvp2m[l,m,n,*]=vvp2m[l,m,n,*]/np[l,m,n]
    endif
  endfor & endfor & endfor
  sigmap=vvp2m-vvpm^2
  sigmap[where(sigmap lt 0.0)]=0.0
  sigmap=sqrt(sigmap)
endif

return, reform(np)

end
