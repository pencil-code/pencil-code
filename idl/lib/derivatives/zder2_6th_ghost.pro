;;
;;  $Id$
;;
;;  Second derivative d^2 / dz^2
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function zder2,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
;AB: the following should not be correct
; if (coord_system ne 'cartesian') then $
;     message, "zder2_6th_ghost: not yet implemented for coord_system='" + coord_system + "'"
;
;  Calculate fmx, fmy, and fmz, based on the input array size.
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'zder2_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
  d = make_array(size=s)
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  l1 = nghostx & l2 = fmx-nghostx-1
  m1 = nghosty & m2 = fmy-nghosty-1
  n1 = nghostz & n2 = fmz-nghostz-1
;
;  Check for degenerate case (no z-derivative)
;
  if (ldegenerated[2] or (fmz eq 1)) then return, d
;
  if (lequidist[2]) then begin
    fdz = dz_1[n1]^2/180.
  endif else begin
    if (fmz ne mz) then $
        message, "zder2_6th_ghost: not implemented for subvolumes on a non-equidistant grid in z."
    fdz = 1./180.
  endelse
;
  d[l1:l2,m1:m2,n1:n2,*] = $
       (-490.*fdz)*f[l1:l2,m1:m2,n1:n2,*] $
      + (270.*fdz)*(f[l1:l2,m1:m2,n1-1:n2-1,*]+f[l1:l2,m1:m2,n1+1:n2+1,*]) $
      -  (27.*fdz)*(f[l1:l2,m1:m2,n1-2:n2-2,*]+f[l1:l2,m1:m2,n1+2:n2+2,*]) $
      +   (2.*fdz)*(f[l1:l2,m1:m2,n1-3:n2-3,*]+f[l1:l2,m1:m2,n1+3:n2+3,*])
;
  if (not lequidist[2]) then begin
    ; Nonuniform mesh correction:
    ; d2f/dz2 = zeta'^2*f" + zeta"*f', see also the manual.
    ; will also work on subvolumes like zder2(ss[10:16,20:26,*])
    df_dz = zder(f)
    for n = n1, n2 do begin
      d[l1:l2,m1:m2,n,*] *= dz_1[n]^2
      df_dz[l1:l2,m1:m2,n,*] *= dz_tilde[n]
    endfor
    d += df_dz
  endif
;
  if (coord_system eq 'spherical') then begin
    if ((fmx ne mx) or (fmy ne my)) then $
        message, "zder_6th_ghost: not implemented for x- or y-subvolumes in spherical coordinates."
    sin_y = sin(y)
    sin1th = 1./sin_y
    i_sin = where(abs(sin_y) lt 1e-5) ; sinth_min=1e-5
    if (i_sin[0] ne -1) then sin1th[i_sin] = 0.
    for l = l1, l2 do d[l,*,*,*] /= x[l]^2
    for m = m1, m2 do d[*,m,*,*] *= sin1th[m]^2
  endif
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
