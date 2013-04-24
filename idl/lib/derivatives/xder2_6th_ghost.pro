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
    dx2=replicate(1./(180.*(x[4]-x[3])^2),nx)
  endif else begin
    dx2=dx_1[l1:l2]^2/180.
  endelse
;
  if (s[0] eq 2) then begin
    if (not ldegenerated[0]) then begin
      for m=m1,m2 do begin
        d[l1:l2,m]=dx2* $
            (-490.*f[l1:l2,m] $
             +270.*(f[l1-1:l2-1,m]+f[l1+1:l2+1,m]) $
              -27.*(f[l1-2:l2-2,m]+f[l1+2:l2+2,m]) $
               +2.*(f[l1-3:l2-3,m]+f[l1+3:l2+3,m]) )
      endfor   
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 3) then begin
    if (not ldegenerated[0]) then begin
      for m=m1,m2 do begin
        for n=n1,n2 do begin
          d[l1:l2,m,n]=dx2* $
              (-490.*f[l1:l2,m,n] $
               +270.*(f[l1-1:l2-1,m,n]+f[l1+1:l2+1,m,n]) $
                -27.*(f[l1-2:l2-2,m,n]+f[l1+2:l2+2,m,n]) $
                 +2.*(f[l1-3:l2-3,m,n]+f[l1+3:l2+3,m,n]) )
        endfor
      endfor  
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[0]) then begin
      for m=m1,m2 do begin
        for n=n1,n2 do begin
          for p=0,s[4]-1 do begin
            d[l1:l2,m,n,p]=dx2* $
                (-490.*f[l1:l2,m,n,p] $
                 +270.*(f[l1-1:l2-1,m,n,p]+f[l1+1:l2+1,m,n,p]) $
                  -27.*(f[l1-2:l2-2,m,n,p]+f[l1+2:l2+2,m,n,p]) $
                   +2.*(f[l1-3:l2-3,m,n,p]+f[l1+3:l2+3,m,n,p]) )
          endfor
        endfor
      endfor  
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: xder2_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
  endelse
;
;  Apply correction for nonuniform mesh.
;  d2f/dx2  = f"*xi'^2 + xi"f', see also the manual.
;
  if (not lequidist[0]) then $
     d=nonuniform_mesh_correction_x(d,f)
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
