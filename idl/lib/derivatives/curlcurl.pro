;;
;;  $Id$
;;
;;  Calculate two consecutive curls of a 3-D vector field.
;;  This routine is faster and more memory efficient than calling
;;  "graddiv" and "del2" and subtracting them.
;;
;;  06-Oct-2014/Bourdin.KIS: coded
;;  10-jun-2016/MR: corrected to be in agreement with PC
;;
function curlcurl,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2, HIDDEN
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
  if (s[0] lt 4) then message, "curlcurl is only implemented for 4D arrays."
  w = make_array(size=s)
;
  ; The following block is optimized, please leave it as it is,
  ; unless you know exactly what you are doing. [Bourdin.KIS]
;
  tmp_add = yder2(f[*,*,*,0]) + zder2(f[*,*,*,0])
  w[*,*,*,0] = yderxder(f[*,*,*,1]) + zderxder(f[*,*,*,2]) - tmp_add
;
  tmp_add = xder2(f[*,*,*,1])
  tmp_add += zder2(f[*,*,*,1])
  w[*,*,*,1] = xderyder(f[*,*,*,0]) + zderyder(f[*,*,*,2]) - tmp_add
;
  tmp_add = xder2(f[*,*,*,2])
  tmp_add += yder2(f[*,*,*,2])
  w[*,*,*,2] = xderzder(f[*,*,*,0]) + yderzder(f[*,*,*,1]) - tmp_add
;
  if (coord_system eq 'cylindrical') then $
    print, 'curlcurl: Warning -- not fully implemented for cylindrical coordinates!'
  if (coord_system eq 'spherical') then begin
;
    sin_y = sin(y[m1:m2])
    cotth = cos(y[m1:m2])/sin_y
    i_sin = where(abs(sin_y) lt 1e-5)       ; sinth_min=1e-5
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
;
    r1 = spread(1.0/x[l1:l2],1,ny)          ; 1/r
;
    r1_cotth = spread (cotth,0,nx) * r1     ; cot(th)/r
    r2_sinth2 = 1./(spread (sin_y,0,nx) * spread(x[l1:l2],1,ny))^2 ; 1/(r sin(th))^2
    if (i_sin[0] ne -1) then r2_sinth2[i_sin]=0.
;
;  needs further optimization
;
    fxy=yder(f[*,*,*,0])
    fxz=zder(f[*,*,*,0])
    fyx=xder(f[*,*,*,1])
    fyy=yder(f[*,*,*,1])
    fzx=xder(f[*,*,*,2])
    fzy=yder(f[*,*,*,2])
    fzz=zder(f[*,*,*,2])
;
    for n = n1, n2 do begin

      w[l1:l2,m1:m2,n,0] +=  r1*      (fyy[l1:l2,m1:m2,n] + fzz[l1:l2,m1:m2,n]) $
                           + r1_cotth*(fyx[l1:l2,m1:m2,n] - fxy[l1:l2,m1:m2,n] + r1*f[l1:l2,m1:m2,n,1])

      w[l1:l2,m1:m2,n,1] +=  r1*     (fxy[l1:l2,m1:m2,n] - 2.*fyx[l1:l2,m1:m2,n]) $
                           + r1_cotth*fzz[l1:l2,m1:m2,n]

      w[l1:l2,m1:m2,n,2] +=  r1*     (fxz[l1:l2,m1:m2,n] - 2.*fzx[l1:l2,m1:m2,n]) $
                           - r1_cotth*fzy[l1:l2,m1:m2,n] + r2_sinth2*f[l1:l2,m1:m2,n,2] 
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
