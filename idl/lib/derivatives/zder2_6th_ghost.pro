;***********************************************************************
function zder2,f
COMPILE_OPT IDL2,HIDDEN
;
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
common cdat_nonequidist,xprim,yprim,zprim,xprim2,yprim2,zprim2,lequidist
;
;  Check if we have a degenerate case (no x-extension)
;
if n_elements(lequidist) ne 3 then lequidist=[1,1,1]
if (nz eq 1) then return,fltarr(nx,ny,nz)
s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
n1=3 & n2=nz-4
;
if lequidist[2] then begin
  dz2=1./(180.*(z[4]-z[3])^2)
endif else begin
  d1=zder(f)
; check if divide by zero:  
  tt = where(zprim eq 0)
  if (tt[0] ne -1) then  zprim[tt] = 1e-10 
  tt = where(zprim2 eq 0)
  if (tt[0] ne -1) then  zprim2[tt] = 1e-10
;
;  nonuniform mesh correction
;  d2f/dz2 = f"/z'^2 - z"*f'/z'^3
;          = f"/z'^2 - z"/z'^2 * df/dz
;
  dz2=spread(spread(1./(180.*zprim^2),0,nx),1,ny)
  dd =d1*spread(spread(zprim2/zprim^2,0,nx),1,ny)
endelse
;
  d[*,*,n1:n2]=dz2*(-490.*f[*,*,n1:n2]$
                   +270.*(f[*,*,n1-1:n2-1]+f[*,*,n1+1:n2+1])$
                    -27.*(f[*,*,n1-2:n2-2]+f[*,*,n1+2:n2+2])$
                     +2.*(f[*,*,n1-3:n2-3]+f[*,*,n1+3:n2+3])$
                   )
;
; apply correction only for nonuniform mesh
;
if not lequidist[2] then d=d-dd
;
return,d
end
