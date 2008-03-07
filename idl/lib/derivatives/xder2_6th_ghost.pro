;
;  $Id: xder2_6th_ghost.pro,v 1.9 2008-03-07 14:36:14 ajohan Exp $
;
;  Second derivative d^2/dx^2
;  - 6th-order (7-point stencil)
;  - with ghost cells
;  - on potentially non-equidistant grid
;
;***********************************************************************
function xder2,f
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
;
;  calculate nx, ny, and nz, based on the input array size
;
  s=size(f) & d=make_array(size=s)
  nx=s[1] & ny=s[2] & nz=s[3]
;
;  Check for degenerate case (no x-extension)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (nx eq 1) then return, fltarr(nx,ny,nz)
;
;  determine location of ghost zones, assume nghost=3 for now.
;
  l1=3 & l2=nx-4
;
  if (lequidist[0]) then begin
    dx2=1./(180.*(x[4]-x[3])^2)
  endif else begin
    dx2=dx_1[l1:l2]^2/180.
;
;  nonuniform mesh correction
;  d2f/dx2  = f"*xi'^2 + xi"f'   see also the manual
;
    d1=xder(f)
  endelse
;
  if (s[0] eq 3) then begin
    if (l2 gt l1) then begin
      if (lequidist[0] eq 0) then begin
        dx2 =    spread(dx2,     [1,2],[s[2],s[3]])
        dd  = d1*spread(dx_tilde,[1,2],[s[2],s[3]])
      endif
      d[l1:l2,*,*]=dx2*(-490.*f[l1:l2,*,*]$
                        +270.*(f[l1-1:l2-1,*,*]+f[l1+1:l2+1,*,*])$
                         -27.*(f[l1-2:l2-2,*,*]+f[l1+2:l2+2,*,*])$
                          +2.*(f[l1-3:l2-3,*,*]+f[l1+3:l2+3,*,*])$
                       )
    endif else begin
      d[l1:l2,*,*]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (l2 gt l1) then begin
      if (lequidist[0] eq 0) then begin
        dx2 =    spread(dx2,     [1,2,3],[s[2],s[3],s[4]])
        dd  = d1*spread(dx_tilde,[1,2,3],[s[2],s[3],s[4]])
      endif
      d[l1:l2,*,*,*]=dx2*(-490.*f[l1:l2,*,*,*]$
                          +270.*(f[l1-1:l2-1,*,*,*]+f[l1+1:l2+1,*,*,*])$
                           -27.*(f[l1-2:l2-2,*,*,*]+f[l1+2:l2+2,*,*,*])$
                            +2.*(f[l1-3:l2-3,*,*,*]+f[l1+3:l2+3,*,*,*])$
                         )
    endif else begin
      d[l1:l2,*,*,*]=0.
    endelse
;
  endif else begin
    print, 'error: xder2_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
  endelse
;
; apply correction only for nonuniform mesh
;
  if (not lequidist[0]) then d=d+dd
;
  return, d
;
end
