;
;  $Id: xder6_6th_ghost.pro,v 1.1 2007-08-12 12:17:36 ajohan Exp $
;
;  Sixth derivative d^6/dx^6
;  - 6th-order (7-point stencil)
;  - with ghost cells
;
;***********************************************************************
function xder6,f
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
    dx6=1./(x[4]-x[3])^6
  endif else begin
;
;  Nonuniform mesh not implemented.
;
    print, 'Nonuniform mesh not implemented for xder6_6th_ghost.pro'
    stop
  endelse
;
  if (s[0] eq 3) then begin
    if (l2 gt l1) then begin
      d[l1:l2,*,*]=dx6*( -20.*f[l1:l2,*,*]$
                         +15.*(f[l1-1:l2-1,*,*]+f[l1+1:l2+1,*,*])$
                          -6.*(f[l1-2:l2-2,*,*]+f[l1+2:l2+2,*,*])$
                          +1.*(f[l1-3:l2-3,*,*]+f[l1+3:l2+3,*,*])$
                       )
    endif else begin
      d[l1:l2,*,*]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (l2 gt l1) then begin
      d[l1:l2,*,*,*]=dx6*( -20.*f[l1:l2,*,*,*]$
                           +15.*(f[l1-1:l2-1,*,*,*]+f[l1+1:l2+1,*,*,*])$
                            -6.*(f[l1-2:l2-2,*,*,*]+f[l1+2:l2+2,*,*,*])$
                            +1.*(f[l1-3:l2-3,*,*,*]+f[l1+3:l2+3,*,*,*])$
                         )
    endif else begin
      d[l1:l2,*,*,*]=0.
    endelse
;
  endif else begin
    print, 'error: xder6_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
    stop
  endelse
;
  return, d
;
end
