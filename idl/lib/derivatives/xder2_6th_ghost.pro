;;
;;  $Id$
;;
;;  Second derivative d^2 / dx^2
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function xder2,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
  common cdat_coords, coord_system
  common pc_precision, zero, one
;
;  Default values.
;
  default, one, 1.d0
  default, ghost, 0
;
;  Calculate fmx, fmy, and fmz, based on the input array size.
;
  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'xder2_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
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
  if (lequidist[1]) then begin
    fdx = dx_1[l1]^2/180.
  endif else begin
    if (fmx ne mx) then $
        message, "xder2_6th_ghost: not implemented for subvolumes on a non-equidistant grid in x."
    fdx = one/180.
  endelse
;
  d[l1:l2,m1:m2,n1:n2,*] = $
       (-490.*fdx)*f[l1:l2,m1:m2,n1:n2,*] $
      + (270.*fdx)*(f[l1-1:l2-1,m1:m2,n1:n2,*]+f[l1+1:l2+1,m1:m2,n1:n2,*]) $
      -  (27.*fdx)*(f[l1-2:l2-2,m1:m2,n1:n2,*]+f[l1+2:l2+2,m1:m2,n1:n2,*]) $
      +   (2.*fdx)*(f[l1-3:l2-3,m1:m2,n1:n2,*]+f[l1+3:l2+3,m1:m2,n1:n2,*])
;
  if (not lequidist[0]) then begin
    ; Nonuniform mesh correction:
    ; d2f/dx2  = f"*xi'^2 + xi"f', see also the manual.
    ; will also work on subvolumes like xder2(ss[*,10:16,20:26])
    df_dx = xder(f)
    for l = l1, l2 do begin
      d[l,m1:m2,n1:n2,*] *= dx_1[l]^2
      df_dx[l,m1:m2,n1:n2,*] *= dx_tilde[l]
    endfor
    d += df_dx
  endif
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
