;;
;;  $Id$
;;
;;  First derivative d / dz
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function zder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
;  Calculate fmx, fmy, and fmz, based on the input array size.
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'zder_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
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
    fdz = dz_1[n1]/60.
  endif else begin
    if (fmz ne mz) then $
        message, "zder_6th_ghost: not implemented for subvolumes on a non-equidistant grid in z."
    fdz = 1./60.
  endelse
;
  d[l1:l2,m1:m2,n1:n2,*] = $
         (45.*fdz)*(f[l1:l2,m1:m2,n1+1:n2+1,*]-f[l1:l2,m1:m2,n1-1:n2-1,*]) $
      -   (9.*fdz)*(f[l1:l2,m1:m2,n1+2:n2+2,*]-f[l1:l2,m1:m2,n1-2:n2-2,*]) $
      +      (fdz)*(f[l1:l2,m1:m2,n1+3:n2+3,*]-f[l1:l2,m1:m2,n1-3:n2-3,*])
;
  if (not lequidist[2]) then for n = n1, n2 do d[*,*,n,*] *= dz_1[n]
;
  if (coord_system eq 'spherical') then begin
    if ((fmx ne mx) or (fmy ne my)) then $
        message, "zder_6th_ghost: not implemented for x- or y-subvolumes in spherical coordinates."
    sin_y = sin(y)
    sin1th = 1./sin_y
    i_sin = where(abs(sin_y) lt 1e-5) ; sinth_min=1e-5
    if (i_sin[0] ne -1) then sin1th[i_sin] = 0.
    for l = l1, l2 do d[l,*,*,*] /= x[l]
    for m = m1, m2 do d[*,m,*,*] *= sin1th[m]
  endif
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
