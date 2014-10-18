;;
;;  $Id$
;;
;;  Second derivative d^2 / dy dx =^= yder (xder (f))
;;  - 6th-order
;;  - with ghost cells
;;
function yderxder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
;AB: the following should not be correct
; if (coord_system ne 'cartesian') then $
;     message, "yderxder_6th_ghost: not yet implemented for coord_system='" + coord_system + "'"
;
;  Calculate fmx, fmy, and fmz, based on the input array size.
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'yderxder_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
  d = make_array(size=s)
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  l1 = nghostx & l2 = fmx-nghostx-1
  m1 = nghosty & m2 = fmy-nghosty-1
  n1 = nghostz & n2 = fmz-nghostz-1
;
;  Check for degenerate case (no xy-derivative)
;
  if (ldegenerated[0] or ldegenerated[1] or (fmx eq 1) or (fmy eq 1)) then return, d
;
;  Calculate d^2 / dy dx (f)
;
  fac = one/60.^2
  if (lequidist[0]) then begin
    if (fmx ne mx) then $
        message, "yderxder_6th_ghost: not implemented for x-subvolumes on a non-equidistant grid in x."
    fac *= dx_1[l1]
  endif
  if (lequidist[1]) then begin
    if (fmy ne my) then $
        message, "yderxder_6th_ghost: not implemented for y-subvolumes on a non-equidistant grid in y."
    fac *= dy_1[m1]
  endif
;
  d[l1:l2,m1:m2,n1:n2,*] = $
       (45.*fac)*( ( 45.*(f[l1+1:l2+1,m1+1:m2+1,n1:n2,*]-f[l1-1:l2-1,m1+1:m2+1,n1:n2,*])   $
                    - 9.*(f[l1+2:l2+2,m1+1:m2+1,n1:n2,*]-f[l1-2:l2-2,m1+1:m2+1,n1:n2,*])   $
                    +    (f[l1+3:l2+3,m1+1:m2+1,n1:n2,*]-f[l1-3:l2-3,m1+1:m2+1,n1:n2,*]))  $
                  -( 45.*(f[l1+1:l2+1,m1-1:m2-1,n1:n2,*]-f[l1-1:l2-1,m1-1:m2-1,n1:n2,*])   $
                    - 9.*(f[l1+2:l2+2,m1-1:m2-1,n1:n2,*]-f[l1-2:l2-2,m1-1:m2-1,n1:n2,*])   $
                    +    (f[l1+3:l2+3,m1-1:m2-1,n1:n2,*]-f[l1-3:l2-3,m1-1:m2-1,n1:n2,*]))) $
      - (9.*fac)*( ( 45.*(f[l1+1:l2+1,m1+2:m2+2,n1:n2,*]-f[l1-1:l2-1,m1+2:m2+2,n1:n2,*])   $
                    - 9.*(f[l1+2:l2+2,m1+2:m2+2,n1:n2,*]-f[l1-2:l2-2,m1+2:m2+2,n1:n2,*])   $
                    +    (f[l1+3:l2+3,m1+2:m2+2,n1:n2,*]-f[l1-3:l2-3,m1+2:m2+2,n1:n2,*]))  $
                  -( 45.*(f[l1+1:l2+1,m1-2:m2-2,n1:n2,*]-f[l1-1:l2-1,m1-2:m2-2,n1:n2,*])   $
                    - 9.*(f[l1+2:l2+2,m1-2:m2-2,n1:n2,*]-f[l1-2:l2-2,m1-2:m2-2,n1:n2,*])   $
                    +    (f[l1+3:l2+3,m1-2:m2-2,n1:n2,*]-f[l1-3:l2-3,m1-2:m2-2,n1:n2,*]))) $
      +    (fac)*( ( 45.*(f[l1+1:l2+1,m1+3:m2+3,n1:n2,*]-f[l1-1:l2-1,m1+3:m2+3,n1:n2,*])   $
                    - 9.*(f[l1+2:l2+2,m1+3:m2+3,n1:n2,*]-f[l1-2:l2-2,m1+3:m2+3,n1:n2,*])   $
                    +    (f[l1+3:l2+3,m1+3:m2+3,n1:n2,*]-f[l1-3:l2-3,m1+3:m2+3,n1:n2,*]))  $
                  -( 45.*(f[l1+1:l2+1,m1-3:m2-3,n1:n2,*]-f[l1-1:l2-1,m1-3:m2-3,n1:n2,*])   $
                    - 9.*(f[l1+2:l2+2,m1-3:m2-3,n1:n2,*]-f[l1-2:l2-2,m1-3:m2-3,n1:n2,*])   $
                    +    (f[l1+3:l2+3,m1-3:m2-3,n1:n2,*]-f[l1-3:l2-3,m1-3:m2-3,n1:n2,*])))
;
  if (not lequidist[0]) then for l = l1, l2 do d[l,*,*,*] *= dx_1[l]
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
