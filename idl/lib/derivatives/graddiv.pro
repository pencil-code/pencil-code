;;
;;  $Id$
;;
;;  Calculate gradient of the divergence of a vector.
;;
function graddiv,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
  s = size(f)
  if (s[0] ne 4) then message, "grad is only implemented for 4D arrays."
  w = make_array(size=s)
;
  w[*,*,*,0]=   xder2(f[*,*,*,0])+xderyder(f[*,*,*,1])+xderzder(f[*,*,*,2])
  w[*,*,*,1]=yderxder(f[*,*,*,0])+   yder2(f[*,*,*,1])+yderzder(f[*,*,*,2])
  w[*,*,*,2]=zderxder(f[*,*,*,0])+zderyder(f[*,*,*,1])+   zder2(f[*,*,*,2])
;
  if (coord_system eq 'spherical') then begin
    sin_y = sin(y[m1:m2])
    cotth = cos(y[m1:m2])/sin_y
    i_sin = where(abs(sin_y) lt 1e-5) ; sinth_min=1e-5
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
    r1= spread(1.0/x[l1:l2],1,ny)
    r2= spread(1.0/x[l1:l2]^2,1,ny)
    r1_cotth= spread (cotth,0,nx) * spread(1.0/x[l1:l2],1,ny)
    r2_sinth2= 1./(spread (sin_y,0,nx) * spread(x[l1:l2],1,ny))^2
    if (i_sin[0] ne -1) then r2_sinth2[i_sin]=0.
    fxx=xder(f[*,*,*,0])
    fyx=xder(f[*,*,*,1])
    fxy=yder(f[*,*,*,0])
    fxz=zder(f[*,*,*,0])
    fyy=yder(f[*,*,*,1])
    fyz=zder(f[*,*,*,1])
    for n = n1, n2 do begin
      w[l1:l2,m1:m2,n,0] += 2.*fxx[l1:l2,m1:m2,n]*r1 $
                              +fyx[l1:l2,m1:m2,n]*r1_cotth $
                           -2.*f[l1:l2,m1:m2,n,0]*r2 $
                              -f[l1:l2,m1:m2,n,1]*r1*r1_cotth
      w[l1:l2,m1:m2,n,1] += 2.*fxy[l1:l2,m1:m2,n]*r1 $
                              +fyy[l1:l2,m1:m2,n]*r1_cotth $
                              -f[l1:l2,m1:m2,n,1]*r2_sinth2
      w[l1:l2,m1:m2,n,2] += 2.*fxz[l1:l2,m1:m2,n]*r1 $
                              +fyz[l1:l2,m1:m2,n]*r1_cotth
    endfor
  endif
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
