;
;  $Id: pc_closest_particles.pro,v 1.1 2005-02-17 13:24:23 ajohan Exp $
;
;  Find nclost closest particles surrounding the point (x,y,z).
;
;  Author: Anders Johansen
;
function pc_closest_particles, xxp, x, y, z, nclose

Lx=1.32
Ly=1.32
Lz=1.32

npar=0L

npar=n_elements(xxp[*,0])

iip=where( (xxp[*,0]-x) gt  0.5*Lx)
iim=where( (xxp[*,0]-x) lt -0.5*Lx)
if (iip[0] ne -1) then xxp[iip,0]=xxp[iip,0]-Lx
if (iim[0] ne -1) then xxp[iim,0]=xxp[iim,0]+Lx

iip=where( (xxp[*,1]-y) gt  0.5*Ly)
iim=where( (xxp[*,1]-y) lt -0.5*Ly)
if (iip[0] ne -1) then xxp[iip,1]=xxp[iip,1]-Ly
if (iim[0] ne -1) then xxp[iim,1]=xxp[iim,1]+Ly

iip=where( (xxp[*,2]-z) gt  0.5*Lz)
iim=where( (xxp[*,2]-z) lt -0.5*Lz)
if (iip[0] ne -1) then xxp[iip,2]=xxp[iip,2]-Lz
if (iim[0] ne -1) then xxp[iim,2]=xxp[iim,2]+Lz

dist2=fltarr(npar)

for k=0L,npar-1 do begin
  dist2[k]=(xxp[k,0]-x)^2 + (xxp[k,1]-y)^2 + (xxp[k,2]-z)^2
endfor

sort2=sort(dist2)

return, sort2[0:nclose-1]


end
