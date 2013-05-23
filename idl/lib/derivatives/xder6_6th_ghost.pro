;;
;;  $Id$
;;
;;  Sixth derivative d^6 / dx^6
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;
function xder6,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t,ignoredx=ignoredx
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
  default, ignoredx, 0
;
  if (coord_system ne 'cartesian') then $
      message, "xder6_6th_ghost: not yet implemented for coord_system='" + coord_system + "'"
;
;  Calculate fmx, fmy, and fmz, based on the input array size.
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'xder6_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
  d = make_array(size=s)
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  l1 = nghostx & l2 = fmx-nghostx-1
  m1 = nghosty & m2 = fmy-nghosty-1
  n1 = nghostz & n2 = fmz-nghostz-1
;
;  Check for degenerate case (no x-derivative)
;
  if (ldegenerated[0] or (fmx eq 1)) then return, d
;
  if (ignoredx) then begin
    fdx=1.
  endif else begin
    if (lequidist[0]) then begin
      fdx = dx_1[l1]^6
    endif else message, "xder6_6th_ghost: non-equidistant grid in x only with ignoredx=1."
  endelse
;
  d[l1:l2,m1:m2,n1:n2,*] = $
        (-20.*fdx)*f[l1:l2,m1:m2,n1:n2,*] $
      +  (15.*fdx)*(f[l1-1:l2-1,m1:m2,n1:n2,*]+f[l1+1:l2+1,m1:m2,n1:n2,*]) $
      -   (6.*fdx)*(f[l1-2:l2-2,m1:m2,n1:n2,*]+f[l1+2:l2+2,m1:m2,n1:n2,*]) $
      +      (fdx)*(f[l1-3:l2-3,m1:m2,n1:n2,*]+f[l1+3:l2+3,m1:m2,n1:n2,*])
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
