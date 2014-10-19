;;
;;  $Id$
;;
;;  Third derivative d^3 / dx dy dz =^= xder (yder (zder (f)))
;;  - 6th-order (9-point stencil) for inner grid cells
;;  - 4th-order (7-point stencil) for non-periodic boundary grid cells
;;
function yderzderxder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
  d = zderyderxder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
  if (coord_system eq 'spherical') then begin
    message, "yderzderxder_6th_ghost: not implemented for spherical grids."
  endif
;
  return, d
;
end
