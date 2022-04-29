;;
;;  $Id$
;;
;;  Calculate the curl of a 3-D vector field
;;
function curlx,f
  curlx_pro, f, crlx
  return, crlx
end
pro curlx_pro,f,curlx
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
  curlx = yder(f[*,*,*,2]) - zder(f[*,*,*,1])
;
  if (coord_system eq 'spherical') then begin
    sin_y = sin(y[m1:m2])
    cotth = cos(y[m1:m2])/sin_y
    i_sin = where(abs(sin_y) lt 1e-5) ; sinth_min=1e-5
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
    corr = spread(cotth,0,nx) * spread(1.0/x[l1:l2],1,ny)
    for n = n1, n2 do curlx[l1:l2,m1:m2,n] += f[l1:l2,m1:m2,n,2]*corr
  endif
;
end
;;
function curly,f
  curly_pro, f, crly
  return, crly
end
pro curly_pro,f,curly
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
  curly = zder(f[*,*,*,0]) - xder(f[*,*,*,2])
;
  if (coord_system eq 'spherical') then begin
    corr = spread(1.0/x[l1:l2],1,ny)
    for n = n1, n2 do curly[l1:l2,m1:m2,n] -= f[l1:l2,m1:m2,n,2]*corr
  end
;
end
;;
function curlz,f
  curlz_pro, f, crlz
  return, crlz
end
pro curlz_pro,f,curlz
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
  curlz = xder(f[*,*,*,1]) - yder(f[*,*,*,0])
;
  if (any (coord_system eq ['cylindric','spherical'])) then begin
    corr = spread(1.0/x[l1:l2],1,ny)
    for n = n1, n2 do curlz[l1:l2,m1:m2,n] += f[l1:l2,m1:m2,n,1]*corr
  end
;
end
;;
function curl,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
;  Default values.
;
  default, ghost, 0
;
  if ((size(f))[0] lt 4) then message, "curl is only implemented for 4D arrays."
  w = make_array(size=size(f))
;
  curlx_pro, f, crlx
  w[*,*,*,0]=crlx
  curly_pro, f, crlx
  w[*,*,*,1]=crlx
  curlz_pro, f, crlx
  w[*,*,*,2]=crlx & undefine, crlx
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
