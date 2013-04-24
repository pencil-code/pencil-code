;;
;;  $Id: xder_6th_ghost.pro 18721 2012-05-08 23:38:47Z Bourdin.KIS $
;;
;;  Correction for nonequidistant grid for the x direction.
;;  d2f/dx2  = f"*xi'^2 + xi"f', see also the manual.
;;  - Adapted from xder_6th_ghost
;;
function nonuniform_mesh_correction_x,d,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
  s=size(f) 
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
    dx2=replicate(1./(60.*(x[4]-x[3])),nx)
  endif else begin
    dx2=dx_1[l1:l2]/60.
  endelse
;
  if (s[0] eq 2) then begin
;
    if (not ldegenerated[0]) then begin
      for m=m1,m2 do begin
        d[l1:l2,m]=d[l1:l2,m] + dx2*dx_tilde* $
            ( +45.*(f[l1+1:l2+1,m]-f[l1-1:l2-1,m]) $
               -9.*(f[l1+2:l2+2,m]-f[l1-2:l2-2,m]) $
                  +(f[l1+3:l2+3,m]-f[l1-3:l2-3,m]) )
      endfor
    endif else d[l1:l2,*]=0.
;
  endif else if (s[0] eq 3) then begin
    if (not ldegenerated[0]) then begin
      for m=m1,m2 do begin
        for n=n1,n2 do begin
          d[l1:l2,m,n]=d[l1:l2,m,n] + dx2*dx_tilde* $
              ( +45.*(f[l1+1:l2+1,m,n]-f[l1-1:l2-1,m,n]) $
                 -9.*(f[l1+2:l2+2,m,n]-f[l1-2:l2-2,m,n]) $
                    +(f[l1+3:l2+3,m,n]-f[l1-3:l2-3,m,n]) )
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
            d[l1:l2,m,n,p]=d[l1:l2,m,n,p] + dx2*dx_tilde* $
                ( +45.*(f[l1+1:l2+1,m,n,p]-f[l1-1:l2-1,m,n,p]) $
                   -9.*(f[l1+2:l2+2,m,n,p]-f[l1-2:l2-2,m,n,p]) $
                      +(f[l1+3:l2+3,m,n,p]-f[l1-3:l2-3,m,n,p]) )
          endfor
        endfor
      endfor
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: nonuniform_mesh_correction_x not implemented for ', $
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
