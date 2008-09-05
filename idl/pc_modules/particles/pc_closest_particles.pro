;
;  $Id$
;
;  Find nclost closest particles surrounding the point (x,y,z).
;
;  Author: Anders Johansen
;
function pc_closest_particles, xxp, xx, nclose, dist2=dist2

Lx=1.32
Ly=1.32
Lz=1.32

npar=0L

npar=n_elements(xxp[*,0])

iip=where( (xxp[*,0]-xx[0]) gt  0.5*Lx)
iim=where( (xxp[*,0]-xx[0]) lt -0.5*Lx)
if (iip[0] ne -1) then xxp[iip,0]=xxp[iip,0]-Lx
if (iim[0] ne -1) then xxp[iim,0]=xxp[iim,0]+Lx

iip=where( (xxp[*,1]-xx[1]) gt  0.5*Ly)
iim=where( (xxp[*,1]-xx[1]) lt -0.5*Ly)
if (iip[0] ne -1) then xxp[iip,1]=xxp[iip,1]-Ly
if (iim[0] ne -1) then xxp[iim,1]=xxp[iim,1]+Ly

iip=where( (xxp[*,2]-xx[2]) gt  0.5*Lz)
iim=where( (xxp[*,2]-xx[2]) lt -0.5*Lz)
if (iip[0] ne -1) then xxp[iip,2]=xxp[iip,2]-Lz
if (iim[0] ne -1) then xxp[iim,2]=xxp[iim,2]+Lz

r2=fltarr(npar)

for k=0L,npar-1 do begin
  r2[k]=(xxp[k,0]-xx[0])^2 + (xxp[k,1]-xx[1])^2 + (xxp[k,2]-xx[2])^2
endfor

sort2=sort(r2)
dist2=r2[sort2[0:nclose-1]]

return, sort2[0:nclose-1]


end
