;;
;;  $Id$
;;
;;  First derivative d/dx
;;  - 6th-order
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function xder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
    dx2=1./(60.*(x[4]-x[3]))
  endif else begin
    dx2=dx_1[l1:l2]/60.
  endelse
;
  if (s[0] eq 2) then begin
;
    if (not ldegenerated[0]) then begin
      if (not lequidist[0]) then dx2=spread(dx2,1,ny)
      d[l1:l2,m1:m2]=dx2* $
          ( +45.*(f[l1+1:l2+1,m1:m2]-f[l1-1:l2-1,m1:m2]) $
             -9.*(f[l1+2:l2+2,m1:m2]-f[l1-2:l2-2,m1:m2]) $
                +(f[l1+3:l2+3,m1:m2]-f[l1-3:l2-3,m1:m2]) )
    endif else d[l1:l2,*]=0.
;
  endif else if (s[0] eq 3) then begin
    if (not ldegenerated[0]) then begin
      if (not lequidist[0]) then dx2=spread(dx2,[1,2],[ny,nz])
      d[l1:l2,m1:m2,n1:n2]=dx2* $
          ( +45.*(f[l1+1:l2+1,m1:m2,n1:n2]-f[l1-1:l2-1,m1:m2,n1:n2]) $
             -9.*(f[l1+2:l2+2,m1:m2,n1:n2]-f[l1-2:l2-2,m1:m2,n1:n2]) $
                +(f[l1+3:l2+3,m1:m2,n1:n2]-f[l1-3:l2-3,m1:m2,n1:n2]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[0]) then begin
      if (not lequidist[0]) then dx2=spread(dx2,[1,2,3],[ny,nz,s[4]])
      d[l1:l2,m1:m2,n1:n2,*]=dx2* $
          ( +45.*(f[l1+1:l2+1,m1:m2,n1:n2,*]-f[l1-1:l2-1,m1:m2,n1:n2,*]) $
             -9.*(f[l1+2:l2+2,m1:m2,n1:n2,*]-f[l1-2:l2-2,m1:m2,n1:n2,*]) $
                +(f[l1+3:l2+3,m1:m2,n1:n2,*]-f[l1-3:l2-3,m1:m2,n1:n2,*]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: xder_6th_ghost not implemented for ', $
        strtrim(s[0],2), '-D arrays'
  endelse
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
