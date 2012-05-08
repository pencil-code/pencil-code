;;
;;  $Id$
;;
;;  First derivative d/dy
;;  - 6th-order
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function yder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
  s=size(f) & d=make_array(size=s)
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
  xx=spread(x,[1,2],[my,mz])
;
  if (lequidist[1]) then begin
    dy2=1./(60.*(y[4]-y[3]))
  endif else begin
    dy2=dy_1[m1:m2]/60.
  endelse
;
  if (s[0] eq 3) then begin
    if (not ldegenerated[1]) then begin
      if (not lequidist[1]) then dy2=spread(dy2,[0,2],[nx,nz])
      ; will also work on slices like yder(ss[10,*,n1:n2])
      d[l1:l2,m1:m2,n1:n2]=dy2* $
          ( +45.*(f[l1:l2,m1+1:m2+1,n1:n2]-f[l1:l2,m1-1:m2-1,n1:n2]) $
             -9.*(f[l1:l2,m1+2:m2+2,n1:n2]-f[l1:l2,m1-2:m2-2,n1:n2]) $
                +(f[l1:l2,m1+3:m2+3,n1:n2]-f[l1:l2,m1-3:m2-3,n1:n2]) )
      if (coord_system ne 'cartesian') then d=d/xx
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[1]) then begin
      if (not lequidist[1]) then dy2=spread(dy2,[0,2,3],[nx,nz,s[4]])
      ; will also work on slices like yder(uu[10,*,*,*,])
      d[l1:l2,m1:m2,n1:n2,*]=dy2* $
          ( +45.*(f[l1:l2,m1+1:m2+1,n1:n2,*]-f[l1:l2,m1-1:m2-1,n1:n2,*]) $
             -9.*(f[l1:l2,m1+2:m2+2,n1:n2,*]-f[l1:l2,m1-2:m2-2,n1:n2,*]) $
                +(f[l1:l2,m1+3:m2+3,n1:n2,*]-f[l1:l2,m1-3:m2-3,n1:n2,*]) )
      if (coord_system ne 'cartesian') then $
         for i=0,s[4]-1 do d[*,*,*,i]=d[*,*,*,i]/xx
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
