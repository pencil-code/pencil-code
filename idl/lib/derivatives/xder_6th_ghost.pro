;
;  $Id: xder_6th_ghost.pro,v 1.9 2005-10-24 02:12:21 dobler Exp $
;
;  6th-order first derivative in x direction for date with ghost cells and on
;  potentially non-equidistant grid.
;
;***********************************************************************
function xder,f
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
  common cdat_nonequidist,xprim,yprim,zprim,xprim2,yprim2,zprim2,lequidist
;
;  Check if we have a degenerate case (no x-extension)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (nx eq 1) then return,fltarr(nx,ny,nz)
  s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
  l1=3 & l2=nx-4
;
  if (lequidist[0]) then begin
    dx2=1./(60.*(x[4]-x[3]))
  endif else begin
    tt = where(xprim eq 0)
    if (tt[0] ne -1) then xprim[tt] = 1
    dx2=spread(spread(xprim[n1:n2]/60.,1,ny),2,nz)
  endelse
;
  if (s[0] eq 1) then begin
    print, 'error: xder_6th_ghost not implemented for 1-D arrays'
  endif else if (s[0] eq 2) then begin
    print, 'error: xder_6th_ghost not implemented for 2-D arrays'
  endif else if (s[0] eq 3) then begin
;
    if (l2 gt l1) then begin
      d[l1:l2,*,*]=dx2*( +45.*(f[l1+1:l2+1,*,*]-f[l1-1:l2-1,*,*]) $
                          -9.*(f[l1+2:l2+2,*,*]-f[l1-2:l2-2,*,*]) $
                             +(f[l1+3:l2+3,*,*]-f[l1-3:l2-3,*,*]) )
    endif else begin
      d[l1:l2,*,*]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (l2 gt l1) then begin
      d[l1:l2,*,*,*]=dx2*( +45.*(f[l1+1:l2+1,*,*,*]-f[l1-1:l2-1,*,*,*]) $
                            -9.*(f[l1+2:l2+2,*,*,*]-f[l1-2:l2-2,*,*,*]) $
                               +(f[l1+3:l2+3,*,*,*]-f[l1-3:l2-3,*,*,*]) )
    endif else begin
      d[l1:l2,*,*,*]=0.
    endelse
;
  endif else begin
    print, 'error: xder_6th_ghost not implemented for ', $
        strtrim(s[0],2), '-D arrays'
  endelse
;
  return, d
;
end
