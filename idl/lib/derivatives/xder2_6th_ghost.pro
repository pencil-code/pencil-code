;;
;;  $Id$
;;
;;  Second derivative d^2/dx^2
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function xder2,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
;
;  Default values.
;
  default, ghost, 0
;
;  Assume nghost=3 for now.
;
  default, nghost, 3
;
;  Calculate mx, my, and mz, based on the input array size.
;
  s=size(f) & d=make_array(size=s)
  mx=s[1] & my=s[2] & mz=s[3]
;
;  Check for degenerate case (no x-derivative)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (mx eq 1) then return, fltarr(mx,my,mz)
;
  l1=nghost & l2=mx-nghost-1
  m1=nghost & m2=my-nghost-1
  n1=nghost & n2=mz-nghost-1
;
  nx = mx - 2*nghost
  ny = my - 2*nghost
  nz = mz - 2*nghost
;
  if (lequidist[0]) then begin
    dx2=1./(180.*(x[4]-x[3])^2)
  endif else begin
    dx2=dx_1[l1:l2]^2/180.
;
;  Nonuniform mesh correction.
;  d2f/dx2  = f"*xi'^2 + xi"f', see also the manual.
;
    d1=xder(f)
  endelse
;
  if (s[0] eq 2) then begin
    if (not ldegenerated[0]) then begin
      if (not lequidist[0]) then begin
        dx2 =    spread(dx2,     1,ny)
        dd  = d1*spread(dx_tilde,1,my)
      endif
      d[l1:l2,m1:m2]=dx2* $
          (-490.*f[l1:l2,m1:m2] $
           +270.*(f[l1-1:l2-1,m1:m2]+f[l1+1:l2+1,m1:m2]) $
            -27.*(f[l1-2:l2-2,m1:m2]+f[l1+2:l2+2,m1:m2]) $
             +2.*(f[l1-3:l2-3,m1:m2]+f[l1+3:l2+3,m1:m2]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 3) then begin
    if (not ldegenerated[0]) then begin
      if (not lequidist[0]) then begin
        dx2 =    spread(dx2,     [1,2],[ny,nz])
        dd  = d1*spread(dx_tilde,[1,2],[my,mz])
      endif
      d[l1:l2,m1:m2,n1:n2]=dx2* $
          (-490.*f[l1:l2,m1:m2,n1:n2] $
           +270.*(f[l1-1:l2-1,m1:m2,n1:n2]+f[l1+1:l2+1,m1:m2,n1:n2]) $
            -27.*(f[l1-2:l2-2,m1:m2,n1:n2]+f[l1+2:l2+2,m1:m2,n1:n2]) $
             +2.*(f[l1-3:l2-3,m1:m2,n1:n2]+f[l1+3:l2+3,m1:m2,n1:n2]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[0]) then begin
      if (not lequidist[0]) then begin
        dx2 =    spread(dx2,     [1,2,3],[ny,nz,s[4]])
        dd  = d1*spread(dx_tilde,[1,2,3],[my,mz,s[4]])
      endif
      d[l1:l2,m1:m2,n1:n2,*]=dx2* $
          (-490.*f[l1:l2,m1:m2,n1:n2,*] $
           +270.*(f[l1-1:l2-1,m1:m2,n1:n2,*]+f[l1+1:l2+1,m1:m2,n1:n2,*]) $
            -27.*(f[l1-2:l2-2,m1:m2,n1:n2,*]+f[l1+2:l2+2,m1:m2,n1:n2,*]) $
             +2.*(f[l1-3:l2-3,m1:m2,n1:n2,*]+f[l1+3:l2+3,m1:m2,n1:n2,*]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: xder2_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
  endelse
;
;  Apply correction only for nonuniform mesh.
;
  if (not lequidist[0]) then d=d+dd
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
