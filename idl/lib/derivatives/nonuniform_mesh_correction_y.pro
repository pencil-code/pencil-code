;;
;;  $Id: yder_6th_ghost.pro 18721 2012-05-08 23:38:47Z Bourdin.KIS $
;;
;;  Correction for nonequidistant grid for the y direction.
;;  d2f/dy2  = f"*psi'^2 + psi"f', see also the manual.
;;  - Adapted from yder_6th_ghost
;;
function nonuniform_mesh_correction_y,d,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
  common cdat_coords,coord_system
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
;  Check for degenerate case (no y-derivative)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (my eq 1) then return, fltarr(mx,my,mz)
;
  l1=nghost & l2=mx-nghost-1
  m1=nghost & m2=my-nghost-1
  n1=nghost & n2=mz-nghost-1
;
  nx = mx - 2*nghost
  ny = my - 2*nghost
  nz = mz - 2*nghost
;
  if (lequidist[1]) then begin
    dy2=replicate(1./(60.*(y[4]-y[3])),ny)
  endif else begin
    dy2=dy_1[m1:m2]/60.
  endelse
;
  if (s[0] eq 3) then begin
    if (not ldegenerated[1]) then begin
      ; will also work on slices like yder(ss[10,*,n1:n2])
      for l=l1,l2 do begin
        for n=n1,n2 do begin
          d[l,m1:m2,n]=d[l,m1:m2,n] + dy2*dy_tilde* $
              ( +45.*(f[l,m1+1:m2+1,n]-f[l,m1-1:m2-1,n]) $
                 -9.*(f[l,m1+2:m2+2,n]-f[l,m1-2:m2-2,n]) $
                    +(f[l,m1+3:m2+3,n]-f[l,m1-3:m2-3,n]) )
        endfor
      endfor
      if (coord_system ne 'cartesian') then begin
        for m=m1,m2 do begin
          for n=n1,n2 do d[l1:l2,m,n]=d[l1:l2,m,n]/x[l1:l2]
        endfor  
      endif
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[1]) then begin
      ; will also work on slices like yder(uu[10,*,*,*,])
      for l=l1,l2 do begin
        for n=n1,n2 do begin
          for p=0,s[4]-1 do begin   
            d[l,m1:m2,n,p]=d[l,m1:m2,n,p] + dy2*dy_tilde* $
                ( +45.*(f[l,m1+1:m2+1,n,p]-f[l,m1-1:m2-1,n,p]) $
                   -9.*(f[l,m1+2:m2+2,n,p]-f[l,m1-2:m2-2,n,p]) $
                      +(f[l,m1+3:m2+3,n,p]-f[l,m1-3:m2-3,n,p]) )
          endfor
        endfor
      endfor
      if (coord_system ne 'cartesian') then $
        for m=m1,m2 do begin
          for n=n1,n2 do begin
            for p=0,s[4]-1 do d[l1:l2,m,n,p]=d[l1:l2,m,n,p]/x[l1:l2]
          endfor
        endfor  
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: yder_6th_ghost not implemented for ', $
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
