;
;  $Id: yder_6th_ghost.pro,v 1.8 2005-10-21 10:19:13 bingert Exp $
;
;***********************************************************************
function yder,f
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
  common cdat_nonequidist,xprim,yprim,zprim,xprim2,yprim2,zprim2,lequidist
;
;  Check if we have a degenerate case (no y-extension)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (ny eq 1) then return,fltarr(nx,ny,nz)
  s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
  m1=3 & m2=ny-4
;
  if (lequidist[1]) then begin
    dy2=1./(60.*(y[4]-y[3]))
  endif else begin
    tt = where(yprim eq 0)
    if (tt[0] ne -1) then yprim[tt] = 1 
    dy2=spread(spread(yprim[n1:n2]/60.,0,nx),2,nz)
  endelse
;
  if (s[0] eq 1) then begin
    print, 'error: yder_6th_ghost not implemented for 1-D arrays'
  endif else if (s[0] eq 2) then begin
    print, 'error: yder_6th_ghost not implemented for 2-D arrays'
  endif else if (s[0] eq 3) then begin
;
    if (m2 gt m1) then begin
      d[*,m1:m2,*]=dy2*( +45.*(f[*,m1+1:m2+1,*]-f[*,m1-1:m2-1,*]) $
                          -9.*(f[*,m1+2:m2+2,*]-f[*,m1-2:m2-2,*]) $
                             +(f[*,m1+3:m2+3,*]-f[*,m1-3:m2-3,*]) )
    endif else begin
      d[*,m1:m2,*]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (m2 gt m1) then begin
      d[*,m1:m2,*,*]=dy2*( +45.*(f[*,m1+1:m2+1,*,*]-f[*,m1-1:m2-1,*,*]) $
                            -9.*(f[*,m1+2:m2+2,*,*]-f[*,m1-2:m2-2,*,*]) $
                               +(f[*,m1+3:m2+3,*,*]-f[*,m1-3:m2-3,*,*]) )
    endif else begin
      d[*,m1:m2,*,*]=0.
    endelse
;
  endif else begin
    print, 'error: yder_6th_ghost not implemented for ', $
        strtrim(s[0],2), '-D arrays'
  endelse
;
  return, d
;
end
