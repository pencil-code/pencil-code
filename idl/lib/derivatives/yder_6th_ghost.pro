;;
;;  $Id$
;;
;;  First derivative df / dy
;;  - 6th-order
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function yder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
  common cdat_coords, coord_system
  common pc_precision, zero, one
;
;  Default values.
;
  default, one, 1.d0
  default, ghost, 0
;
;  Calculate fmx, fmy, and fmz, based on the input array size.
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'yder_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
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
    fdy = dy_1[m1]/60.
  endif else begin
    if (fmy ne my) then $
        message, "yder_6th_ghost: not implemented for subvolumes on a non-equidistant grid in y."
    fdy = one/60.
  endelse
;
  d[l1:l2,m1:m2,n1:n2,*] = $
         (45.*fdy)*(f[l1:l2,m1+1:m2+1,n1:n2,*]-f[l1:l2,m1-1:m2-1,n1:n2,*]) $
      -   (9.*fdy)*(f[l1:l2,m1+2:m2+2,n1:n2,*]-f[l1:l2,m1-2:m2-2,n1:n2,*]) $
      +      (fdy)*(f[l1:l2,m1+3:m2+3,n1:n2,*]-f[l1:l2,m1-3:m2-3,n1:n2,*])
;
  if (not lequidist[1]) then for m = m1, m2 do d[*,m,*,*] *= dy_1[m]
;
  if (any (coord_system eq ['cylindric','spherical'])) then $
      for l = l1, l2 do d[l,*,*,*] /= x[l]
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
