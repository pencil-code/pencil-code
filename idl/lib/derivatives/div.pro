;;
;;  $Id$
;;
;;  Calculate divergence of vector field.
;;
function div,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, nx, ny, nz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
  s = size(f)
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  l1 = nghostx & l2 = fmx-nghostx-1
  m1 = nghosty & m2 = fmy-nghosty-1
  n1 = nghostz & n2 = fmz-nghostz-1
;
  w = xder(f[*,*,*,0]) + yder(f[*,*,*,1]) + zder(f[*,*,*,2])
;
  if (coord_system eq 'cylindric') then begin
    for l = l1, l2 do w[l,*,*] += f[l,*,*,0]/x[l]
  endif else if (coord_system eq 'spherical') then begin
    sin_y = sin(y[m1:m2])
    cotth = cos(y[m1:m2])/sin_y
    i_sin = where(abs(sin_y) lt 1e-5) ; sinth_min=1e-5
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
    corr = spread (cotth,0,nx) * spread(1.0/x[l1:l2],1,ny)
    for n = n1, n2 do w[l1:l2,m1:m2,n] += (2*f[l1:l2,m1:m2,n,0]+f[l1:l2,m1:m2,n,1])*corr
  endif
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
