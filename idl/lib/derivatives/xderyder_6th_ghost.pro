;;
;;  $Id$
;;
;;  Second derivative d2f/dydx
;;  - 6th-order
;;  - with ghost cells
;;
function xderyder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
  tmp=yderxder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
  if (coord_system eq 'spherical') then begin
    r1=spread(1.0/x,[1,2],[my,mz])
    scr=yder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
    tmp=tmp-scr*r1
  endif
;
  return, tmp
;
end
