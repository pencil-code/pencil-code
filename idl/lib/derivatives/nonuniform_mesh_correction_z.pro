;;
;;  $Id: zder_6th_ghost.pro 18721 2012-05-08 23:38:47Z Bourdin.KIS $
;;
;;  Correction for nonequidistant grid for the z direction. 
;;  d2f/dz2 = zeta'^2*f" + zeta"*f', see also the manual.
;;  - Adapted from zder_6th_ghost
;;
function nonuniform_mesh_correction_z,d,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
;  Check for degenerate case (no z-derivative)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (mz eq 1) then return, fltarr(mx,my,mz)
;
  l1=nghost & l2=mx-nghost-1
  m1=nghost & m2=my-nghost-1
  n1=nghost & n2=mz-nghost-1
;
  nx = mx - 2*nghost
  ny = my - 2*nghost
  nz = mz - 2*nghost
;
  sin1th=1./sin(y)
  i_sin=where(abs(sin(y)) lt 1e-5) ;sinth_min=1e-5
  if (i_sin[0] ne -1) then sin1th[i_sin]=0.
;
  if (lequidist[2]) then begin
    dz2=replicate(1./(60.*(z[4]-z[3])),nz) 
  endif else begin
    dz2=dz_1[n1:n2]/60.
  endelse
;
  if (s[0] eq 2) then begin
    if (n2 gt n1) then begin
      for l=l1,l2 do begin 
        d[l,n1:n2]=d[l,n1:n2] + dz2*dz_tilde* $
            ( +45.*(f[l,n1+1:n2+1]-f[l,n1-1:n2-1]) $
               -9.*(f[l,n1+2:n2+2]-f[l,n1-2:n2-2]) $
                  +(f[l,n1+3:n2+3]-f[l,n1-3:n2-3]) )
      endfor
      if (coord_system eq 'spherical') then begin
        for l=l1,l2 do begin   
          d[l,n1:n2]=d[l,n1:n2]/x[l]*sin1th[nghost]
        endfor  
      endif
    endif else d[l1:l2,n1:n2]=0.
;
  endif else if (s[0] eq 3) then begin
    if (not ldegenerated[2]) then begin
      ; will also work on slices like zder(ss[10,20,*])
      for l=l1,l2 do begin 
        for m=m1,m2 do begin 
          d[l,m,n1:n2]=d[l,m,n1:n2] + dz2*dz_tilde* $
              ( +45.*(f[l,m,n1+1:n2+1]-f[l,m,n1-1:n2-1]) $
                 -9.*(f[l,m,n1+2:n2+2]-f[l,m,n1-2:n2-2]) $
                    +(f[l,m,n1+3:n2+3]-f[l,m,n1-3:n2-3]) )
        endfor
      endfor 
      if (coord_system eq 'spherical') then begin 
        for l=l1,l2 do begin 
          for m=m1,m2 do begin 
            d[l,m,n1:n2]=d[l,m,n1:n2]/x[l]*sin1th[m]
          endfor
        endfor 
      endif
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[2]) then begin
      ; will also work on slices like zder(uu[10,20,*,*])
      for l=l1,l2 do begin 
        for m=m1,m2 do begin 
          for p=0,s[4]-1 do begin 
             d[l,m,n1:n2,p]=d[l,m,n1:n2,p] + dz2*dz_tilde* $
                 ( +45.*(f[l,m,n1+1:n2+1,p]-f[l,m,n1-1:n2-1,p]) $
                    -9.*(f[l,m,n1+2:n2+2,p]-f[l,m,n1-2:n2-2,p]) $
                       +(f[l,m,n1+3:n2+3,p]-f[l,m,n1-3:n2-3,p]) )
          endfor
        endfor 
      endfor
      if (coord_system eq 'spherical') then begin
        for l=l1,l2 do begin 
          for m=m1,m2 do begin 
            d[l,m,n1:n2,*]=d[l,m,n1:n2,*]/x[l]*sin1th[m]
          endfor
        endfor 
      endif   
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: nonuniform_mesh_correction_z not implemented for ', $
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
