;
;  $Id: zder_6th_ghost.pro,v 1.7 2004-12-25 14:21:23 ajohan Exp $
;
;***********************************************************************
function zder,f
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
  common cdat_nonequidist,xprim,yprim,zprim,xprim2,yprim2,zprim2,lequidist
;
;  Check if we have a degenerate case (no x-extension)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (nz eq 1) then return,fltarr(nx,ny,nz)
  s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
  n1=3 & n2=nz-4
;
  if (lequidist[2]) then begin
    dz2=1./(60.*(z[4]-z[3])) 
  endif else begin
    tt = where(zprim eq 0)
    if (tt[0] ne -1) then zprim[tt] = 1 
    dz2=spread(spread(zprim/60.,0,nx),1,ny)
  endelse
;
  if (s[0] eq 1) then begin
    print, 'error: zder_6th_ghost not implemented for 1-D arrays'
  endif else if (s[0] eq 2) then begin
    print, 'error: zder_6th_ghost not implemented for 2-D arrays'
  endif else if (s[0] eq 3) then begin
;
    if (n2 gt n1) then begin
      d[*,*,n1:n2]=dz2*( +45.*(f[*,*,n1+1:n2+1]-f[*,*,n1-1:n2-1]) $
                          -9.*(f[*,*,n1+2:n2+2]-f[*,*,n1-2:n2-2]) $
                             +(f[*,*,n1+3:n2+3]-f[*,*,n1-3:n2-3]) )
    endif else begin
      d[*,*,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (n2 gt n1) then begin
      d[*,*,n1:n2,*]=dz2*( +45.*(f[*,*,n1+1:n2+1,*]-f[*,*,n1-1:n2-1,*]) $
                            -9.*(f[*,*,n1+2:n2+2,*]-f[*,*,n1-2:n2-2,*]) $
                               +(f[*,*,n1+3:n2+3,*]-f[*,*,n1-3:n2-3,*]) )
    endif else begin
      d[*,*,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: zder_6th_ghost not implemented for ', $
        strtrim(s[0],2), '-D arrays'
  endelse
;
  return, d
;
end
