;;
;;  $Id$
;;
;;  Second derivative d2f/dzdy
;;  - 6th-order
;;  - with ghost cells
;;
function yderzder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
  tmp=zderyder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
  if (coord_system eq 'spherical') then begin
    sin_y = sin(y)
    cotth = cos(y)/sin_y
    i_sin = where(abs(sin_y) lt 1e-5) ; sinth_min=1e-5
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
    r1_cotth=spread((1.0/x)#cotth,2,mz)
    scr=zder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
    tmp=tmp-scr*r1_cotth
  endif
;
  return, tmp
;
end
