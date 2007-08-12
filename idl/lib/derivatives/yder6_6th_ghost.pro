;
;  $Id: yder6_6th_ghost.pro,v 1.1 2007-08-12 12:17:36 ajohan Exp $
;
;  Sixth derivative d^6/dy^6
;  - 6th-order (7-point stencil)
;  - with ghost cells
;
;***********************************************************************
function yder6,f
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
;
;  Calculate nx, ny, and nz, based on the input array size.
;
  s=size(f) & d=make_array(size=s)
  nx=s[1] & ny=s[2] & nz=s[3]
;
;  Check for degenerate case (no x-extension).
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (ny eq 1) then return,fltarr(nx,ny,nz)
;
;  Determine location of ghost zones, assume nghost=3 for now.
;
  m1=3 & m2=ny-4
;
  if (lequidist[1]) then begin
    dy6=1./(y[4]-y[3])^6
  endif else begin
;
;  Nonuniform mesh not implemented.
;
    print, 'Nonuniform mesh not implemented for yder6_6th_ghost.pro'
    stop
  endelse
;
  if (s[0] eq 3) then begin
    if (m2 gt m1) then begin
      d[*,m1:m2,*]=dy6*( -20.*f[*,m1:m2,*]$
                         +15.*(f[*,m1-1:m2-1,*]+f[*,m1+1:m2+1,*])$
                          -6.*(f[*,m1-2:m2-2,*]+f[*,m1+2:m2+2,*])$
                          +1.*(f[*,m1-3:m2-3,*]+f[*,m1+3:m2+3,*])$
                       )
    endif else begin
      d[*,m1:m2,*]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (m2 gt m1) then begin
      d[*,m1:m2,*,*]=dy6*( -20.*f[*,m1:m2,*,*]$
                           +15.*(f[*,m1-1:m2-1,*,*]+f[*,m1+1:m2+1,*,*])$
                            -6.*(f[*,m1-2:m2-2,*,*]+f[*,m1+2:m2+2,*,*])$
                            +1.*(f[*,m1-3:m2-3,*,*]+f[*,m1+3:m2+3,*,*])$
                         )
    endif else begin
      d[*,m1:m2,*,*]=0.
    endelse
;
  endif else begin
    print, 'error: yder6_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
    stop
  endelse
;
  return, d
;
end
