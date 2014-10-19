;;
;;  $Id$
;;
;;  Third derivative d^3 / dz dy dx =^= zder (yder (xder (f)))
;;  - 6th-order (9-point stencil) for inner grid cells
;;  - 4th-order (7-point stencil) for non-periodic boundary grid cells
;
; 18-Oct-2014/Bourdin.KIS: coded
;
function zderyderxder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
      message, 'zderyderxder_6th_ghost: not implemented for '+strtrim(s[0],2)+'-D arrays'
  d = make_array(size=s)
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  l1 = nghostx & l2 = fmx-nghostx-1
  m1 = nghosty & m2 = fmy-nghosty-1
  n1 = nghostz & n2 = fmz-nghostz-1
;
;  Check for degenerate case (no xyz-derivative)
;
  if (any (ldegenerated) or (fmx eq 1) or (fmy eq 1) or (fmz eq 1)) then return, d
;
;  Calculate d^3 / dz dy dx (f)
;
  fac = one/60.^2
  if (lequidist[0]) then begin
    fac *= dx_1[l1]
  end else begin
    if (fmx ne mx) then $
        message, "zderyderxder_6th_ghost: not implemented for x-subvolumes on a non-equidistant grid in x."
  end
  if (lequidist[1]) then begin
    fac *= dy_1[m1]
  end else begin
    if (fmy ne my) then $
        message, "zderyderxder_6th_ghost: not implemented for y-subvolumes on a non-equidistant grid in y."
  end
  if (lequidist[2]) then begin
    fac *= dz_1[n1]
  end else begin
    if (fmz ne mz) then $
        message, "zderyderxder_6th_ghost: not implemented for z-subvolumes on a non-equidistant grid in z."
  end
;
; Differentiation scheme:
; d[l,m,n] = fac*( 45*(yderxder (f[l,m,n+1]) - yderxder (f[l,m,n-1]))
;                 - 9*(yderxder (f[l,m,n+2]) - yderxder (f[l,m,n-2]))
;                 +   (yderxder (f[l,m,n+3]) - yderxder (f[l,m,n-3])) )
; with yderxder (f) as:
; d[l,m,n] = fac*( 45*(xder (f[l,m+1,n]) - xder (f[l,m-1,n]))
;                 - 9*(xder (f[l,m+2,n]) - xder (f[l,m-2,n]))
;                 +   (xder (f[l,m+3,n]) - xder (f[l,m-3,n])) )
;
  d[l1:l2,m1:m2,n1:n2,*] = $
       (91125*fac)*(f[l1+1:l2+1,m1+1:m2+1,n1+1:n2+1,*]-f[l1-1:l2-1,m1+1:m2+1,n1+1:n2+1,*]) $
      -(18225*fac)*(f[l1+2:l2+2,m1+1:m2+1,n1+1:n2+1,*]-f[l1-2:l2-2,m1+1:m2+1,n1+1:n2+1,*]) $
      + (2025*fac)*(f[l1+3:l2+3,m1+1:m2+1,n1+1:n2+1,*]-f[l1-3:l2-3,m1+1:m2+1,n1+1:n2+1,*]) $
      -(91125*fac)*(f[l1+1:l2+1,m1-1:m2-1,n1+1:n2+1,*]-f[l1-1:l2-1,m1-1:m2-1,n1+1:n2+1,*]) $
      +(18225*fac)*(f[l1+2:l2+2,m1-1:m2-1,n1+1:n2+1,*]-f[l1-2:l2-2,m1-1:m2-1,n1+1:n2+1,*]) $
      - (2025*fac)*(f[l1+3:l2+3,m1-1:m2-1,n1+1:n2+1,*]-f[l1-3:l2-3,m1-1:m2-1,n1+1:n2+1,*]) $
      -(18225*fac)*(f[l1+1:l2+1,m1+2:m2+2,n1+1:n2+1,*]-f[l1-1:l2-1,m1+2:m2+2,n1+1:n2+1,*]) $
      + (3645*fac)*(f[l1+2:l2+2,m1+2:m2+2,n1+1:n2+1,*]-f[l1-2:l2-2,m1+2:m2+2,n1+1:n2+1,*]) $
      -  (405*fac)*(f[l1+3:l2+3,m1+2:m2+2,n1+1:n2+1,*]-f[l1-3:l2-3,m1+2:m2+2,n1+1:n2+1,*]) $
      +(18225*fac)*(f[l1+1:l2+1,m1-2:m2-2,n1+1:n2+1,*]-f[l1-1:l2-1,m1-2:m2-2,n1+1:n2+1,*]) $
      - (3645*fac)*(f[l1+2:l2+2,m1-2:m2-2,n1+1:n2+1,*]-f[l1-2:l2-2,m1-2:m2-2,n1+1:n2+1,*]) $
      +  (405*fac)*(f[l1+3:l2+3,m1-2:m2-2,n1+1:n2+1,*]-f[l1-3:l2-3,m1-2:m2-2,n1+1:n2+1,*]) $
      + (2025*fac)*(f[l1+1:l2+1,m1+3:m2+3,n1+1:n2+1,*]-f[l1-1:l2-1,m1+3:m2+3,n1+1:n2+1,*]) $
      -  (405*fac)*(f[l1+2:l2+2,m1+3:m2+3,n1+1:n2+1,*]-f[l1-2:l2-2,m1+3:m2+3,n1+1:n2+1,*]) $
      +   (45*fac)*(f[l1+3:l2+3,m1+3:m2+3,n1+1:n2+1,*]-f[l1-3:l2-3,m1+3:m2+3,n1+1:n2+1,*]) $
      - (2025*fac)*(f[l1+1:l2+1,m1-3:m2-3,n1+1:n2+1,*]-f[l1-1:l2-1,m1-3:m2-3,n1+1:n2+1,*]) $
      +  (405*fac)*(f[l1+2:l2+2,m1-3:m2-3,n1+1:n2+1,*]-f[l1-2:l2-2,m1-3:m2-3,n1+1:n2+1,*]) $
      -   (45*fac)*(f[l1+3:l2+3,m1-3:m2-3,n1+1:n2+1,*]-f[l1-3:l2-3,m1-3:m2-3,n1+1:n2+1,*]) $
      -(91125*fac)*(f[l1+1:l2+1,m1+1:m2+1,n1-1:n2-1,*]-f[l1-1:l2-1,m1+1:m2+1,n1-1:n2-1,*]) $
      +(18225*fac)*(f[l1+2:l2+2,m1+1:m2+1,n1-1:n2-1,*]-f[l1-2:l2-2,m1+1:m2+1,n1-1:n2-1,*]) $
      - (2025*fac)*(f[l1+3:l2+3,m1+1:m2+1,n1-1:n2-1,*]-f[l1-3:l2-3,m1+1:m2+1,n1-1:n2-1,*]) $
      +(91125*fac)*(f[l1+1:l2+1,m1-1:m2-1,n1-1:n2-1,*]-f[l1-1:l2-1,m1-1:m2-1,n1-1:n2-1,*]) $
      -(18225*fac)*(f[l1+2:l2+2,m1-1:m2-1,n1-1:n2-1,*]-f[l1-2:l2-2,m1-1:m2-1,n1-1:n2-1,*]) $
      + (2025*fac)*(f[l1+3:l2+3,m1-1:m2-1,n1-1:n2-1,*]-f[l1-3:l2-3,m1-1:m2-1,n1-1:n2-1,*]) $
      +(18225*fac)*(f[l1+1:l2+1,m1+2:m2+2,n1-1:n2-1,*]-f[l1-1:l2-1,m1+2:m2+2,n1-1:n2-1,*]) $
      - (3645*fac)*(f[l1+2:l2+2,m1+2:m2+2,n1-1:n2-1,*]-f[l1-2:l2-2,m1+2:m2+2,n1-1:n2-1,*]) $
      +  (405*fac)*(f[l1+3:l2+3,m1+2:m2+2,n1-1:n2-1,*]-f[l1-3:l2-3,m1+2:m2+2,n1-1:n2-1,*]) $
      -(18225*fac)*(f[l1+1:l2+1,m1-2:m2-2,n1-1:n2-1,*]-f[l1-1:l2-1,m1-2:m2-2,n1-1:n2-1,*]) $
      + (3645*fac)*(f[l1+2:l2+2,m1-2:m2-2,n1-1:n2-1,*]-f[l1-2:l2-2,m1-2:m2-2,n1-1:n2-1,*]) $
      -  (405*fac)*(f[l1+3:l2+3,m1-2:m2-2,n1-1:n2-1,*]-f[l1-3:l2-3,m1-2:m2-2,n1-1:n2-1,*]) $
      - (2025*fac)*(f[l1+1:l2+1,m1+3:m2+3,n1-1:n2-1,*]-f[l1-1:l2-1,m1+3:m2+3,n1-1:n2-1,*]) $
      +  (405*fac)*(f[l1+2:l2+2,m1+3:m2+3,n1-1:n2-1,*]-f[l1-2:l2-2,m1+3:m2+3,n1-1:n2-1,*]) $
      -   (45*fac)*(f[l1+3:l2+3,m1+3:m2+3,n1-1:n2-1,*]-f[l1-3:l2-3,m1+3:m2+3,n1-1:n2-1,*]) $
      + (2025*fac)*(f[l1+1:l2+1,m1-3:m2-3,n1-1:n2-1,*]-f[l1-1:l2-1,m1-3:m2-3,n1-1:n2-1,*]) $
      -  (405*fac)*(f[l1+2:l2+2,m1-3:m2-3,n1-1:n2-1,*]-f[l1-2:l2-2,m1-3:m2-3,n1-1:n2-1,*]) $
      +   (45*fac)*(f[l1+3:l2+3,m1-3:m2-3,n1-1:n2-1,*]-f[l1-3:l2-3,m1-3:m2-3,n1-1:n2-1,*]) $
      -(18225*fac)*(f[l1+1:l2+1,m1+1:m2+1,n1+1:n2+2,*]-f[l1-1:l2-1,m1+1:m2+1,n1+1:n2+2,*]) $
      + (3645*fac)*(f[l1+2:l2+2,m1+1:m2+1,n1+1:n2+2,*]-f[l1-2:l2-2,m1+1:m2+1,n1+1:n2+2,*]) $
      -  (405*fac)*(f[l1+3:l2+3,m1+1:m2+1,n1+1:n2+2,*]-f[l1-3:l2-3,m1+1:m2+1,n1+1:n2+2,*]) $
      +(18225*fac)*(f[l1+1:l2+1,m1-1:m2-1,n1+1:n2+2,*]-f[l1-1:l2-1,m1-1:m2-1,n1+1:n2+2,*]) $
      - (3645*fac)*(f[l1+2:l2+2,m1-1:m2-1,n1+1:n2+2,*]-f[l1-2:l2-2,m1-1:m2-1,n1+1:n2+2,*]) $
      +  (405*fac)*(f[l1+3:l2+3,m1-1:m2-1,n1+1:n2+2,*]-f[l1-3:l2-3,m1-1:m2-1,n1+1:n2+2,*]) $
      + (3645*fac)*(f[l1+1:l2+1,m1+2:m2+2,n1+1:n2+2,*]-f[l1-1:l2-1,m1+2:m2+2,n1+1:n2+2,*]) $
      -  (729*fac)*(f[l1+2:l2+2,m1+2:m2+2,n1+1:n2+2,*]-f[l1-2:l2-2,m1+2:m2+2,n1+1:n2+2,*]) $
      +   (81*fac)*(f[l1+3:l2+3,m1+2:m2+2,n1+1:n2+2,*]-f[l1-3:l2-3,m1+2:m2+2,n1+1:n2+2,*]) $
      - (3645*fac)*(f[l1+1:l2+1,m1-2:m2-2,n1+1:n2+2,*]-f[l1-1:l2-1,m1-2:m2-2,n1+1:n2+2,*]) $
      +  (729*fac)*(f[l1+2:l2+2,m1-2:m2-2,n1+1:n2+2,*]-f[l1-2:l2-2,m1-2:m2-2,n1+1:n2+2,*]) $
      -   (81*fac)*(f[l1+3:l2+3,m1-2:m2-2,n1+1:n2+2,*]-f[l1-3:l2-3,m1-2:m2-2,n1+1:n2+2,*]) $
      -  (405*fac)*(f[l1+1:l2+1,m1+3:m2+3,n1+1:n2+2,*]-f[l1-1:l2-1,m1+3:m2+3,n1+1:n2+2,*]) $
      +   (81*fac)*(f[l1+2:l2+2,m1+3:m2+3,n1+1:n2+2,*]-f[l1-2:l2-2,m1+3:m2+3,n1+1:n2+2,*]) $
      -    (9*fac)*(f[l1+3:l2+3,m1+3:m2+3,n1+1:n2+2,*]-f[l1-3:l2-3,m1+3:m2+3,n1+1:n2+2,*]) $
      +  (405*fac)*(f[l1+1:l2+1,m1-3:m2-3,n1+1:n2+2,*]-f[l1-1:l2-1,m1-3:m2-3,n1+1:n2+2,*]) $
      -   (81*fac)*(f[l1+2:l2+2,m1-3:m2-3,n1+1:n2+2,*]-f[l1-2:l2-2,m1-3:m2-3,n1+1:n2+2,*]) $
      +    (9*fac)*(f[l1+3:l2+3,m1-3:m2-3,n1+1:n2+2,*]-f[l1-3:l2-3,m1-3:m2-3,n1+1:n2+2,*]) $
      +(18225*fac)*(f[l1+1:l2+1,m1+1:m2+1,n1-1:n2-2,*]-f[l1-1:l2-1,m1+1:m2+1,n1-1:n2-2,*]) $
      - (3645*fac)*(f[l1+2:l2+2,m1+1:m2+1,n1-1:n2-2,*]-f[l1-2:l2-2,m1+1:m2+1,n1-1:n2-2,*]) $
      +  (405*fac)*(f[l1+3:l2+3,m1+1:m2+1,n1-1:n2-2,*]-f[l1-3:l2-3,m1+1:m2+1,n1-1:n2-2,*]) $
      -(18225*fac)*(f[l1+1:l2+1,m1-1:m2-1,n1-1:n2-2,*]-f[l1-1:l2-1,m1-1:m2-1,n1-1:n2-2,*]) $
      + (3645*fac)*(f[l1+2:l2+2,m1-1:m2-1,n1-1:n2-2,*]-f[l1-2:l2-2,m1-1:m2-1,n1-1:n2-2,*]) $
      -  (405*fac)*(f[l1+3:l2+3,m1-1:m2-1,n1-1:n2-2,*]-f[l1-3:l2-3,m1-1:m2-1,n1-1:n2-2,*]) $
      - (3645*fac)*(f[l1+1:l2+1,m1+2:m2+2,n1-1:n2-2,*]-f[l1-1:l2-1,m1+2:m2+2,n1-1:n2-2,*]) $
      +  (729*fac)*(f[l1+2:l2+2,m1+2:m2+2,n1-1:n2-2,*]-f[l1-2:l2-2,m1+2:m2+2,n1-1:n2-2,*]) $
      -   (81*fac)*(f[l1+3:l2+3,m1+2:m2+2,n1-1:n2-2,*]-f[l1-3:l2-3,m1+2:m2+2,n1-1:n2-2,*]) $
      + (3645*fac)*(f[l1+1:l2+1,m1-2:m2-2,n1-1:n2-2,*]-f[l1-1:l2-1,m1-2:m2-2,n1-1:n2-2,*]) $
      -  (729*fac)*(f[l1+2:l2+2,m1-2:m2-2,n1-1:n2-2,*]-f[l1-2:l2-2,m1-2:m2-2,n1-1:n2-2,*]) $
      +   (81*fac)*(f[l1+3:l2+3,m1-2:m2-2,n1-1:n2-2,*]-f[l1-3:l2-3,m1-2:m2-2,n1-1:n2-2,*]) $
      +  (405*fac)*(f[l1+1:l2+1,m1+3:m2+3,n1-1:n2-2,*]-f[l1-1:l2-1,m1+3:m2+3,n1-1:n2-2,*]) $
      -   (81*fac)*(f[l1+2:l2+2,m1+3:m2+3,n1-1:n2-2,*]-f[l1-2:l2-2,m1+3:m2+3,n1-1:n2-2,*]) $
      +    (9*fac)*(f[l1+3:l2+3,m1+3:m2+3,n1-1:n2-2,*]-f[l1-3:l2-3,m1+3:m2+3,n1-1:n2-2,*]) $
      -  (405*fac)*(f[l1+1:l2+1,m1-3:m2-3,n1-1:n2-2,*]-f[l1-1:l2-1,m1-3:m2-3,n1-1:n2-2,*]) $
      +   (81*fac)*(f[l1+2:l2+2,m1-3:m2-3,n1-1:n2-2,*]-f[l1-2:l2-2,m1-3:m2-3,n1-1:n2-2,*]) $
      -    (9*fac)*(f[l1+3:l2+3,m1-3:m2-3,n1-1:n2-2,*]-f[l1-3:l2-3,m1-3:m2-3,n1-1:n2-2,*]) $
      + (2025*fac)*(f[l1+1:l2+1,m1+1:m2+1,n1+1:n2+3,*]-f[l1-1:l2-1,m1+1:m2+1,n1+1:n2+3,*]) $
      -  (405*fac)*(f[l1+2:l2+2,m1+1:m2+1,n1+1:n2+3,*]-f[l1-2:l2-2,m1+1:m2+1,n1+1:n2+3,*]) $
      +   (45*fac)*(f[l1+3:l2+3,m1+1:m2+1,n1+1:n2+3,*]-f[l1-3:l2-3,m1+1:m2+1,n1+1:n2+3,*]) $
      - (2025*fac)*(f[l1+1:l2+1,m1-1:m2-1,n1+1:n2+3,*]-f[l1-1:l2-1,m1-1:m2-1,n1+1:n2+3,*]) $
      +  (405*fac)*(f[l1+2:l2+2,m1-1:m2-1,n1+1:n2+3,*]-f[l1-2:l2-2,m1-1:m2-1,n1+1:n2+3,*]) $
      -   (45*fac)*(f[l1+3:l2+3,m1-1:m2-1,n1+1:n2+3,*]-f[l1-3:l2-3,m1-1:m2-1,n1+1:n2+3,*]) $
      -  (405*fac)*(f[l1+1:l2+1,m1+2:m2+2,n1+1:n2+3,*]-f[l1-1:l2-1,m1+2:m2+2,n1+1:n2+3,*]) $
      +   (81*fac)*(f[l1+2:l2+2,m1+2:m2+2,n1+1:n2+3,*]-f[l1-2:l2-2,m1+2:m2+2,n1+1:n2+3,*]) $
      -    (9*fac)*(f[l1+3:l2+3,m1+2:m2+2,n1+1:n2+3,*]-f[l1-3:l2-3,m1+2:m2+2,n1+1:n2+3,*]) $
      +  (405*fac)*(f[l1+1:l2+1,m1-2:m2-2,n1+1:n2+3,*]-f[l1-1:l2-1,m1-2:m2-2,n1+1:n2+3,*]) $
      -   (81*fac)*(f[l1+2:l2+2,m1-2:m2-2,n1+1:n2+3,*]-f[l1-2:l2-2,m1-2:m2-2,n1+1:n2+3,*]) $
      +    (9*fac)*(f[l1+3:l2+3,m1-2:m2-2,n1+1:n2+3,*]-f[l1-3:l2-3,m1-2:m2-2,n1+1:n2+3,*]) $
      +   (45*fac)*(f[l1+1:l2+1,m1+3:m2+3,n1+1:n2+3,*]-f[l1-1:l2-1,m1+3:m2+3,n1+1:n2+3,*]) $
      -    (9*fac)*(f[l1+2:l2+2,m1+3:m2+3,n1+1:n2+3,*]-f[l1-2:l2-2,m1+3:m2+3,n1+1:n2+3,*]) $
      +      (fac)*(f[l1+3:l2+3,m1+3:m2+3,n1+1:n2+3,*]-f[l1-3:l2-3,m1+3:m2+3,n1+1:n2+3,*]) $
      -   (45*fac)*(f[l1+1:l2+1,m1-3:m2-3,n1+1:n2+3,*]-f[l1-1:l2-1,m1-3:m2-3,n1+1:n2+3,*]) $
      +    (9*fac)*(f[l1+2:l2+2,m1-3:m2-3,n1+1:n2+3,*]-f[l1-2:l2-2,m1-3:m2-3,n1+1:n2+3,*]) $
      -      (fac)*(f[l1+3:l2+3,m1-3:m2-3,n1+1:n2+3,*]-f[l1-3:l2-3,m1-3:m2-3,n1+1:n2+3,*]) $
      - (2025*fac)*(f[l1+1:l2+1,m1+1:m2+1,n1-1:n2-3,*]-f[l1-1:l2-1,m1+1:m2+1,n1-1:n2-3,*]) $
      +  (405*fac)*(f[l1+2:l2+2,m1+1:m2+1,n1-1:n2-3,*]-f[l1-2:l2-2,m1+1:m2+1,n1-1:n2-3,*]) $
      -   (45*fac)*(f[l1+3:l2+3,m1+1:m2+1,n1-1:n2-3,*]-f[l1-3:l2-3,m1+1:m2+1,n1-1:n2-3,*]) $
      + (2025*fac)*(f[l1+1:l2+1,m1-1:m2-1,n1-1:n2-3,*]-f[l1-1:l2-1,m1-1:m2-1,n1-1:n2-3,*]) $
      -  (405*fac)*(f[l1+2:l2+2,m1-1:m2-1,n1-1:n2-3,*]-f[l1-2:l2-2,m1-1:m2-1,n1-1:n2-3,*]) $
      +   (45*fac)*(f[l1+3:l2+3,m1-1:m2-1,n1-1:n2-3,*]-f[l1-3:l2-3,m1-1:m2-1,n1-1:n2-3,*]) $
      +  (405*fac)*(f[l1+1:l2+1,m1+2:m2+2,n1-1:n2-3,*]-f[l1-1:l2-1,m1+2:m2+2,n1-1:n2-3,*]) $
      -   (81*fac)*(f[l1+2:l2+2,m1+2:m2+2,n1-1:n2-3,*]-f[l1-2:l2-2,m1+2:m2+2,n1-1:n2-3,*]) $
      +    (9*fac)*(f[l1+3:l2+3,m1+2:m2+2,n1-1:n2-3,*]-f[l1-3:l2-3,m1+2:m2+2,n1-1:n2-3,*]) $
      -  (405*fac)*(f[l1+1:l2+1,m1-2:m2-2,n1-1:n2-3,*]-f[l1-1:l2-1,m1-2:m2-2,n1-1:n2-3,*]) $
      +   (81*fac)*(f[l1+2:l2+2,m1-2:m2-2,n1-1:n2-3,*]-f[l1-2:l2-2,m1-2:m2-2,n1-1:n2-3,*]) $
      -    (9*fac)*(f[l1+3:l2+3,m1-2:m2-2,n1-1:n2-3,*]-f[l1-3:l2-3,m1-2:m2-2,n1-1:n2-3,*]) $
      -   (45*fac)*(f[l1+1:l2+1,m1+3:m2+3,n1-1:n2-3,*]-f[l1-1:l2-1,m1+3:m2+3,n1-1:n2-3,*]) $
      +    (9*fac)*(f[l1+2:l2+2,m1+3:m2+3,n1-1:n2-3,*]-f[l1-2:l2-2,m1+3:m2+3,n1-1:n2-3,*]) $
      -      (fac)*(f[l1+3:l2+3,m1+3:m2+3,n1-1:n2-3,*]-f[l1-3:l2-3,m1+3:m2+3,n1-1:n2-3,*]) $
      +   (45*fac)*(f[l1+1:l2+1,m1-3:m2-3,n1-1:n2-3,*]-f[l1-1:l2-1,m1-3:m2-3,n1-1:n2-3,*]) $
      -    (9*fac)*(f[l1+2:l2+2,m1-3:m2-3,n1-1:n2-3,*]-f[l1-2:l2-2,m1-3:m2-3,n1-1:n2-3,*]) $
      +      (fac)*(f[l1+3:l2+3,m1-3:m2-3,n1-1:n2-3,*]-f[l1-3:l2-3,m1-3:m2-3,n1-1:n2-3,*])
;
  if (not lequidist[0]) then for l = l1, l2 do d[l,*,*,*] *= dx_1[l]
  if (not lequidist[1]) then for m = m1, m2 do d[*,m,*,*] *= dy_1[m]
  if (not lequidist[2]) then for n = n1, n2 do d[*,*,n,*] *= dz_1[n]
;
  if (any (coord_system eq ['cylindric','spherical'])) then $
      message, "zderyderxder_6th_ghost: not implemented for cylindrical or spherical grids."
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
