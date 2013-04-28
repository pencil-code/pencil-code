;;
;;  $Id$
;;
;;  Second derivative d^2 / dx dz
;;  - 6th-order
;;  - with ghost cells
;;
function xderzder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
  if (coord_system ne 'cartesian') then $
      message, "xderzder_6th_ghost: not yet implemented for coord_system='" + coord_system + "'"
;
;  Calculate fmx, fmy, and fmz, based on the input array size.
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'xderzder_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
  d = make_array(size=s)
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  l1 = nghostx & l2 = fmx-nghostx-1
  m1 = nghosty & m2 = fmy-nghosty-1
  n1 = nghostz & n2 = fmz-nghostz-1
;
;  Check for degenerate case (no xz-derivative)
;
  if (ldegenerated[0] or ldegenerated[2] or (fmx eq 1) or (fmz eq 1)) then return, d
;
;  Calculate d2f/dxdz.
;
  fac = 1./60.^2
  if (lequidist[0]) then begin
    if (fmx ne mx) then $
        message, "xderzder_6th_ghost: not implemented for x-subvolumes on a non-equidistant grid in x."
    fac *= dx_1[l1]
  endif
  if (lequidist[2]) then begin
    if (fmz ne mz) then $
        message, "xderzder_6th_ghost: not implemented for z-subvolumes on a non-equidistant grid in z."
    fac *= dz_1[n1]
  endif
;
  d[l1:l2,m1:m2,n1:n2,*] = $
       (45.*fac)*( ( 45.*(f[l1+1:l2+1,m1:m2,n1+1:n2+1,*]-f[l1-1:l2-1,m1:m2,n1+1:n2+1,*])   $
                    - 9.*(f[l1+2:l2+2,m1:m2,n1+1:n2+1,*]-f[l1-2:l2-2,m1:m2,n1+1:n2+1,*])   $
                    +    (f[l1+3:l2+3,m1:m2,n1+1:n2+1,*]-f[l1-3:l2-3,m1:m2,n1+1:n2+1,*]))  $
                  -( 45.*(f[l1+1:l2+1,m1:m2,n1-1:n2-1,*]-f[l1-1:l2-1,m1:m2,n1-1:n2-1,*])   $
                    - 9.*(f[l1+2:l2+2,m1:m2,n1-1:n2-1,*]-f[l1-2:l2-2,m1:m2,n1-1:n2-1,*])   $
                    +    (f[l1+3:l2+3,m1:m2,n1-1:n2-1,*]-f[l1-3:l2-3,m1:m2,n1-1:n2-1,*]))) $
      - (9.*fac)*( ( 45.*(f[l1+1:l2+1,m1:m2,n1+2:n2+2,*]-f[l1-1:l2-1,m1:m2,n1+2:n2+2,*])   $
                    - 9.*(f[l1+2:l2+2,m1:m2,n1+2:n2+2,*]-f[l1-2:l2-2,m1:m2,n1+2:n2+2,*])   $
                    +    (f[l1+3:l2+3,m1:m2,n1+2:n2+2,*]-f[l1-3:l2-3,m1:m2,n1+2:n2+2,*]))  $
                  -( 45.*(f[l1+1:l2+1,m1:m2,n1-2:n2-2,*]-f[l1-1:l2-1,m1:m2,n1-2:n2-2,*])   $
                    - 9.*(f[l1+2:l2+2,m1:m2,n1-2:n2-2,*]-f[l1-2:l2-2,m1:m2,n1-2:n2-2,*])   $
                    +    (f[l1+3:l2+3,m1:m2,n1-2:n2-2,*]-f[l1-3:l2-3,m1:m2,n1-2:n2-2,*]))) $
      +    (fac)*( ( 45.*(f[l1+1:l2+1,m1:m2,n1+3:n2+3,*]-f[l1-1:l2-1,m1:m2,n1+3:n2+3,*])   $
                    - 9.*(f[l1+2:l2+2,m1:m2,n1+3:n2+3,*]-f[l1-2:l2-2,m1:m2,n1+3:n2+3,*])   $
                    +    (f[l1+3:l2+3,m1:m2,n1+3:n2+3,*]-f[l1-3:l2-3,m1:m2,n1+3:n2+3,*]))  $
                  -( 45.*(f[l1+1:l2+1,m1:m2,n1-3:n2-3,*]-f[l1-1:l2-1,m1:m2,n1-3:n2-3,*])   $
                    - 9.*(f[l1+2:l2+2,m1:m2,n1-3:n2-3,*]-f[l1-2:l2-2,m1:m2,n1-3:n2-3,*])   $
                    +    (f[l1+3:l2+3,m1:m2,n1-3:n2-3,*]-f[l1-3:l2-3,m1:m2,n1-3:n2-3,*])))
;
  if (not lequidist[0]) then for l = l1, l2 do d[l,*,*,*] *= dx_1[l]
  if (not lequidist[2]) then for n = n1, n2 do d[*,*,n,*] *= dz_1[n]
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
