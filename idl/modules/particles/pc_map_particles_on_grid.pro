;
;  $Id$
;
;  Map the particle positions on the grid.
;
;  Author: Anders Johansen
;
pro pc_map_particles_on_grid, xxp, x, y, z, $
    ineargrid=ineargrid, kneighbour=kneighbour, kshepherd=kshepherd, $
    datadir=datadir

if (n_elements(x) gt 1) then dx=x[1]-x[0] else dx=1.0
if (n_elements(y) gt 1) then dy=y[1]-y[0] else dy=1.0
if (n_elements(z) gt 1) then dz=z[1]-z[0] else dz=1.0
dx1=1.0D/dx
dy1=1.0D/dy
dz1=1.0D/dz

npar=0L

npar=n_elements(xxp[*,0])
nx=n_elements(x)
ny=n_elements(y)
nz=n_elements(z)

ineargrid=lonarr(npar,3)

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

  ineargrid[k,*]=[ix,iy,iz]

endfor

kneighbour=lonarr(npar,/nozero) & kneighbour[*]=-1
kshepherd =lonarr(nx,ny,nz,/nozero) & kshepherd[*,*,*]=-1

for k=0L,npar-1 do begin
  ix=ineargrid[k,0] & iy=ineargrid[k,1] & iz=ineargrid[k,2]
  kneighbour[k]=kshepherd[ix,iy,iz]
  kshepherd[ix,iy,iz]=k
endfor

end
