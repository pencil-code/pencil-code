;;
;;  $Id$
;;
;;  First derivative d/dz
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function zder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
  common cdat_coords,coord_system
;
;  Default values.
;
  default, ghost, 0
;
;  Calculate nx, ny, and nz, based on the input array size.
;
  s=size(f) & d=make_array(size=s)
  nx=s[1] & ny=s[2] & nz=s[3]
;
  xx=spread(x,[1,2],[ny,nz])
  yy=spread(y,[0,2],[nx,nz])
  sin1th=1./sin(yy)
  i_sin=where(abs(sin(yy)) lt 1e-5) ;sinth_min=1e-5
  if (i_sin[0] ne -1) then sin1th[i_sin]=0.
;

; 26-jun-2007/dintrans: 2-D case only means (x,z) for the moment
  if (s[0] eq 2) then nz=s[2]
;
;  Check for degenerate case (no z-extension).
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (nz eq 1) then return,fltarr(nx,ny,nz)
;
;  Determine location of ghost zones, assume nghost=3 for now.
;
  l1=3 & l2=nx-4 & m1=3 & m2=ny-4 & n1=3 & n2=nz-4
;
  if (lequidist[2]) then begin
    dz2=1./(60.*(z[4]-z[3])) 
  endif else begin
    dz2=dz_1[n1:n2]/60.
  endelse
;
  if (s[0] eq 2) then begin
    if (n2 gt n1) then begin
      if (lequidist[2] eq 0) then dz2=spread(dz2,0,s[2])
      d[l1:l2,n1:n2]=dz2* $
          ( +45.*(f[l1:l2,n1+1:n2+1]-f[l1:l2,n1-1:n2-1]) $
             -9.*(f[l1:l2,n1+2:n2+2]-f[l1:l2,n1-2:n2-2]) $
                +(f[l1:l2,n1+3:n2+3]-f[l1:l2,n1-3:n2-3]) )
      if (coord_system eq 'spherical') then d=d/xx*sin1th
    endif else d[l1:l2,n1:n2]=0.
;
  endif else if (s[0] eq 3) then begin
    if (n2 gt n1) then begin
      if (lequidist[2] eq 0) then dz2=spread(dz2,[0,0],[s[2],s[1]])
      ; will also work on slices like zder(ss[10,20,*])
      d[l1:l2,m1:m2,n1:n2]=dz2* $
          ( +45.*(f[l1:l2,m1:m2,n1+1:n2+1]-f[l1:l2,m1:m2,n1-1:n2-1]) $
             -9.*(f[l1:l2,m1:m2,n1+2:n2+2]-f[l1:l2,m1:m2,n1-2:n2-2]) $
                +(f[l1:l2,m1:m2,n1+3:n2+3]-f[l1:l2,m1:m2,n1-3:n2-3]) )
      if (coord_system eq 'spherical') then d=d/xx*sin1th
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (n2 gt n1) then begin
      if (lequidist[2] eq 0) then dz2=spread(dz2,[0,0,3],[s[2],s[1],s[4]])
      ; will also work on slices like zder(uu[10,20,*,*])
      d[l1:l2,m1:m2,n1:n2,*]=dz2* $
          ( +45.*(f[l1:l2,m1:m2,n1+1:n2+1,*]-f[l1:l2,m1:m2,n1-1:n2-1,*]) $
             -9.*(f[l1:l2,m1:m2,n1+2:n2+2,*]-f[l1:l2,m1:m2,n1-2:n2-2,*]) $
                +(f[l1:l2,m1:m2,n1+3:n2+3,*]-f[l1:l2,m1:m2,n1-3:n2-3,*]) )
      if (coord_system eq 'spherical') then $
        d[l1:l2,m1:m2,*,0:s[4]-1]=d[l1:l2,m1:m2,*,0:s[4]-1]/xx*sin1th
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: zder_6th_ghost not implemented for ', $
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
