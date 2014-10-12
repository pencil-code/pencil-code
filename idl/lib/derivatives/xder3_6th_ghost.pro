;;
;;  $Id$
;;
;;  Third derivative d^3 / dx^3
;;  - 6th-order (9-point stencil) for inner grid cells
;;  - 4th-order (7-point stencil) for non-periodic boundary grid cells
;
; 12-Oct-2014/Bourdin.KIS: coded
;
function xder3,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
      message, 'xder3_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
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
    fdx = dx_1[l1]^3/240.
  endif else begin
    if (fmx ne mx) then $
        message, "xder3_6th_ghost: not implemented for subvolumes on a non-equidistant grid in x."
    fdx = 1./240.
  endelse
;
  if (lperi[1]) then begin
    d[l1:l2,m1:m2,n1:n2,*] = $
          (488.*fdx)*(f[l1-1:l2-1,m1:m2,n1:n2,*]-f[l1+1:l2+1,m1:m2,n1:n2,*]) $
        - (338.*fdx)*(f[l1-2:l2-2,m1:m2,n1:n2,*]-f[l1+2:l2+2,m1:m2,n1:n2,*]) $
        +  (72.*fdx)*(f[l1-3:l2-3,m1:m2,n1:n2,*]-f[l1+3:l2+3,m1:m2,n1:n2,*]) $
        -   (7.*fdx)*(shift (f[l1:l2,m1:m2,n1:n2,*], -4, 0, 0, 0)-shift (f[l1:l2,m1:m2,n1:n2,*], +4, 0, 0, 0))
  end else begin
    ; inner grid points (6th order)
    d[l1+1:l2-1,m1:m2,n1:n2,*] = $
          (488.*fdx)*(f[l1-0:l2-2,m1:m2,n1:n2,*]-f[l1+2:l2+0,m1:m2,n1:n2,*]) $
        - (338.*fdx)*(f[l1-1:l2-3,m1:m2,n1:n2,*]-f[l1+3:l2+1,m1:m2,n1:n2,*]) $
        +  (72.*fdx)*(f[l1-2:l2-4,m1:m2,n1:n2,*]-f[l1+4:l2+2,m1:m2,n1:n2,*]) $
        -   (7.*fdx)*(f[l1-3:l2-5,m1:m2,n1:n2,*]-f[l1+5:l2+3,m1:m2,n1:n2,*])
    ; boundary grid points (4th order)
;    −49/8 	29 	−461/8 	62 	−307/8 	13 	−15/8
    d[l1,m1:m2,n1:n2,*] = $
         (14880.*fdx)*f[l1,m1:m2,n1:n2,*] $
        -(13830.*fdx)*f[l1-1,m1:m2,n1:n2,*]-(9210.*fdx)*f[l1+1,m1:m2,n1:n2,*] $
        + (6960.*fdx)*f[l1-2,m1:m2,n1:n2,*]+(3120.*fdx)*f[l1+2,m1:m2,n1:n2,*] $
        - (1470.*fdx)*f[l1-3,m1:m2,n1:n2,*]- (450.*fdx)*f[l1+3,m1:m2,n1:n2,*]
    d[l2,m1:m2,n1:n2,*] = $
         (14880.*fdx)*f[l2,m1:m2,n1:n2,*] $
        -(13830.*fdx)*f[l1-1,m1:m2,n1:n2,*]-(9210.*fdx)*f[l1+1,m1:m2,n1:n2,*] $
        + (6960.*fdx)*f[l1-2,m1:m2,n1:n2,*]+(3120.*fdx)*f[l1+2,m1:m2,n1:n2,*] $
        - (1470.*fdx)*f[l1-3,m1:m2,n1:n2,*]- (450.*fdx)*f[l1+3,m1:m2,n1:n2,*]
  end
;
  if (not lequidist[0]) then begin
    ; Nonuniform mesh correction:
    ; d3f(xi)/dx3 = f'''*xi'^3 + 3*f"*xi'*xi" + f'*xi'''
    ; will also work on subvolumes like xder3(ss[*,10:16,20:26])
    df_dx = xder(f)
    d2f_dx2 = xder2(f)
    d[l1:l2,m1:m2,n1:n2,*] *= dx_1[l1:l2]^3
    message, "xder3: 'd2xi_dx2' is not yet implemented, needs more than 'dx_tilde'..."
    df_dx[l1:l2,m1:m2,n1:n2,*] *= d2xi_dx2[l1:l2]
    d2f_dx2[l1:l2,m1:m2,n1:n2,*] *= dx_tilde[l1:l2]
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
