;
;  $Id: pc_particles_to_density.pro,v 1.1 2005-01-26 10:18:07 ajohan Exp $
;
;  Convert positions of particles to a number density field.
;
;  Author: Anders Johansen
;
function pc_particles_to_density, xxp, x, y, z

npar=n_elements(xxp[*,0])
nx=n_elements(x)
ny=n_elements(y)
nz=n_elements(z)

np=fltarr(nx,ny,nz)
distx=fltarr(nx)
disty=fltarr(ny)
distz=fltarr(nz)

for k=0,npar-1 do begin

  for l=0,nx-1 do begin
    distx[l]=abs(xxp[k,0]-x[l])
  endfor
  ix=where( distx eq min(distx) )

  for m=0,ny-1 do begin
    disty[m]=abs(xxp[k,1]-y[m])
  endfor
  iy=where( disty eq min(disty) )

  for n=0,nz-1 do begin
    distz[n]=abs(xxp[k,2]-x[n])
  endfor
  iz=where( distz eq min(distz) )

  np[ix,iy,iz]=np[ix,iy,iz]+1.0

endfor

return, np

end
