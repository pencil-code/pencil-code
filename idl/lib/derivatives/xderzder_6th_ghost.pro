;;
;;  $Id$
;;
;;  Second derivative d^2 / dx dz =^= xder (zder (f))
;;  - 6th-order
;;  - with ghost cells
;;
function xderzder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
  d = zderxder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
  if (coord_system eq 'spherical') then begin
    d -= spread(1.0/x,[1,2],[my,mz]) * zder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
  endif
;
  return, d
;
end
