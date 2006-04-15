;
;  $Id: pc_particles_to_velocity.pro,v 1.2 2006-04-15 07:53:44 ajohan Exp $
;
;  Convert individual particle velocities to a velocity field.
;
;  Author: Anders Johansen
;
function pc_particles_to_velocity, xxp, vvp, x, y, z, datadir=datadir

COMMON pc_precision, zero, one

pc_set_precision, datadir=datadir
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
vvpsum=fltarr(nx,ny,nz,3)*one
uup=fltarr(nx,ny,nz,3)*one

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

  vvpsum[ix,iy,iz,*]=vvpsum[ix,iy,iz,*]+vvp[k,*]
  np[ix,iy,iz]=np[ix,iy,iz]+1.0*one

endfor

for l=0,nx-1 do begin & for m=0,ny-1 do begin & for n=0,nz-1 do begin
  if (np[l,m,n] ne 0.0) then uup[l,m,n,*]=vvpsum[l,m,n,*]/np[l,m,n]
endfor & endfor & endfor

return, reform(uup)

end
