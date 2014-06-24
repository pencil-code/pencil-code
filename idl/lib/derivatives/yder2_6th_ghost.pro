;;
;;  $Id$
;;
;;  Second derivative d^2 / dy^2
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function yder2,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
;     message, "yder2_6th_ghost: not yet implemented for coord_system='" + coord_system + "'"
;
;  Calculate fmx, fmy, and fmz, based on the input array size.
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'yder2_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
  d = make_array(size=s)
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  l1 = nghostx & l2 = fmx-nghostx-1
  m1 = nghosty & m2 = fmy-nghosty-1
  n1 = nghostz & n2 = fmz-nghostz-1
;
;  Check for degenerate case (no y-derivative)
;
  if (ldegenerated[1] or (fmy eq 1)) then return, d
;
  if (lequidist[1]) then begin
    fdy = dy_1[m1]^2/180.
  endif else begin
    if (fmy ne my) then $
        message, "yder2_6th_ghost: not implemented for subvolumes on a non-equidistant grid in y."
    fdy = 1./180.
  endelse
;
  d[l1:l2,m1:m2,n1:n2,*] = $
       (-490.*fdy)*f[l1:l2,m1:m2,n1:n2,*] $
      + (270.*fdy)*(f[l1:l2,m1-1:m2-1,n1:n2,*]+f[l1:l2,m1+1:m2+1,n1:n2,*]) $
      -  (27.*fdy)*(f[l1:l2,m1-2:m2-2,n1:n2,*]+f[l1:l2,m1+2:m2+2,n1:n2,*]) $
      +   (2.*fdy)*(f[l1:l2,m1-3:m2-3,n1:n2,*]+f[l1:l2,m1+3:m2+3,n1:n2,*])
;
  if (not lequidist[1]) then begin
    ; Nonuniform mesh correction:
    ; d2f/dy2  = f"*psi'^2 + psi"f', see also the manual.
    ; will also work on subvolumes like yder2(ss[10:16,*,20:26])
    df_dy = yder(f)
    for m = m1, m2 do begin
      d[l1:l2,m,n1:n2,*] *= dy_1[m]^2
      df_dy[l1:l2,m,n1:n2,*] *= dy_tilde[m]
    endfor
    d += df_dy
  endif
;
;  Set ghost zones.
;
  if (any (coord_system eq ['cylindric','spherical'])) then $
      for l = l1, l2 do d[l,*,*,*] /= x[l]^2
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
