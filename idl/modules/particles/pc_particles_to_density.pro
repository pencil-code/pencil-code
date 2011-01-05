;;
;;  $Id$
;;
;;  Convert positions of particles to a grid density field.
;;
;;  Author: Anders Johansen
;;
function pc_particles_to_density, xxp, x, y, z, $
    cic=cic, tsc=tsc, fine=fine, ghost=ghost, normalize=normalize, $
    density=density, datadir=datadir, quiet=quiet
COMMON pc_precision, zero, one
;;
;;  Set default values.
;;
default, cic, 0
default, tsc, 0
default, fine, 1
default, ghost, 0
default, normalize, 0
default, density, 1
default, quiet, 0
;;
;;  Read parameters.
;;
pc_read_param, obj=param, datadir=datadir, /quiet
;;
;;  Set real precision.
;;
pc_set_precision, datadir=datadir, /quiet
;;
;;  Read rhopswarm from data directory.
;;
pc_read_const, obj=cst, datadir=datadir, /quiet
;
npar=0L
npar=n_elements(xxp[*,0])
nx=n_elements(x)
ny=n_elements(y)
nz=n_elements(z)
;
if (nx gt 1) then dx=x[1]-x[0] else dx=1.0*one
if (ny gt 1) then dy=y[1]-y[0] else dy=1.0*one
if (nz gt 1) then dz=z[1]-z[0] else dz=1.0*one
dx_1=1.0d/dx   & dy_1=1.0d/dy   & dz_1=1.0d/dz
dx_2=1.0d/dx^2 & dy_2=1.0d/dy^2 & dz_2=1.0d/dz^2
;;
;;  Set interpolation scheme
;;
interpolation_scheme='ngp'
if (cic) then interpolation_scheme='cic'
if (tsc) then interpolation_scheme='tsc'
;;
;;  The CIC and TSC schemes work with ghost cells, so if x, y, z are given
;;  without ghost zones, add the ghost zones automatically.
;;
if (not ghost) then begin
  mx=nx+6 & l1=3 & l2=l1+nx-1
  x2=fltarr(mx)*one
  x2[l1:l2]=x
  for l=l1-1,   0,-1 do x2[l]=x2[l+1]-dx
  for l=l2+1,mx-1,+1 do x2[l]=x2[l-1]+dx
  x=x2
;
  my=ny+6 & m1=3 & m2=m1+ny-1
  y2=fltarr(my)*one
  y2[m1:m2]=y
  for m=m1-1,   0,-1 do y2[m]=y2[m+1]-dy
  for m=m2+1,my-1,+1 do y2[m]=y2[m-1]+dy
  y=y2
;
  mz=nz+6 & n1=3 & n2=n1+nz-1
  z2=fltarr(mz)*one
  z2[n1:n2]=z
  for n=n1-1,   0,-1 do z2[n]=z2[n+1]-dz
  for n=n2+1,mz-1,+1 do z2[n]=z2[n-1]+dz
  z=z2
endif else begin
  mx=nx & nx=mx-6 & l1=3 & l2=l1+nx-1
  my=ny & ny=my-6 & m1=3 & m2=m1+ny-1
  mz=nz & nz=mz-6 & n1=3 & n2=n1+nz-1
endelse
;;
;;  Possible to map the particles on a finer grid.
;;
if (fine gt 1) then begin
;
  x0=param.xyz0[0] & y0=param.xyz0[1] & z0=param.xyz0[2]
  x1=param.xyz1[0] & y1=param.xyz1[1] & z1=param.xyz1[2]
;
  if (nx gt 1) then begin
    nx=fine*nx
    mx=nx+6
    dx=dx/fine
    x=fltarr(nx)
    for i=0,mx-1 do x[i]=(i-3)*dx+dx/2
  endif
;
  if (ny gt 1) then begin
    ny=fine*ny
    my=ny+6
    dy=dy/fine
    y=fltarr(ny)
    y[0]=y0+dy/2
    for i=0,my-1 do y[i]=(i-3)*dy+dy/2
  endif
;
  if (nz gt 1) then begin
    nz=fine*nz
    mz=nz+6
    dz=dz/fine
    z=fltarr(nz)
    z[0]=z0+dz/2
    for i=0,mz-1 do z[i]=(i-3)*dz+dz/2
  endif
;
  dx_1=1.0d/dx   & dy_1=1.0d/dy   & dz_1=1.0d/dz
  dx_2=1.0d/dx^2 & dy_2=1.0d/dy^2 & dz_2=1.0d/dz^2
;
endif
;;
;;  Define density array.
;;
np=fltarr(mx,my,mz)*one
;;
;;  Three different ways to assign particle density to the grid are
;;  implemented:   (see the book by Hockney & Eastwood)
;;    0. NGP (Nearest Grid Point)
;;    1. CIC (Cloud In Cell)
;;    2. TSC (Triangular Shaped Cloud)
;;
case interpolation_scheme of
;;
;;  Assign particle density to the grid using the zeroth order NGP method.
;;
  'ngp': begin
;
    if (not quiet) then print, 'Assigning density using NGP method.'
;
    for k=0L,npar-1 do begin
;  
      ix = round((xxp[k,0]-x[0])*dx_1)
      iy = round((xxp[k,1]-y[0])*dy_1)
      iz = round((xxp[k,2]-z[0])*dz_1)
      if (ix eq l2+1) then ix=ix-1
      if (iy eq m2+1) then iy=iy-1
      if (iz eq n2+1) then iz=iz-1
      if (ix eq l1-1) then ix=ix+1
      if (iy eq m1-1) then iy=iy+1
      if (iz eq n1-1) then iz=iz+1
;;
;;  Particles are assigned to the nearest grid point.
;;
      np[ix,iy,iz]=np[ix,iy,iz]+1.0
;
    endfor ; loop over particles
;
  end ; 'ngp'
;;
;;  Assign particle density to the grid using the first order CIC method.
;;
  'cic': begin
;
    if (not quiet) then print, 'Assigning density using CIC method.'
;
    for k=0L,npar-1 do begin
;;  Find nearest grid point     
      ix0 = round((xxp[k,0]-x[0])*dx_1)
      iy0 = round((xxp[k,1]-y[0])*dy_1)
      iz0 = round((xxp[k,2]-z[0])*dz_1)
      if (ix0 eq l2+1) then ix0=ix0-1
      if (iy0 eq m2+1) then iy0=iy0-1
      if (iz0 eq n2+1) then iz0=iz0-1
      if (ix0 eq l1-1) then ix0=ix0+1
      if (iy0 eq m1-1) then iy0=iy0+1
      if (iz0 eq n1-1) then iz0=iz0+1
;;  Find lower grid point in surrounding grid points.        
      if ( (x[ix0] gt xxp[k,0]) and (nx ne 1) ) then ix0=ix0-1
      if ( (y[iy0] gt xxp[k,1]) and (ny ne 1) ) then iy0=iy0-1
      if ( (z[iz0] gt xxp[k,2]) and (nz ne 1) ) then iz0=iz0-1
;;  Don't assign density to degenerate directions. 
      ix1=ix0 & if (nx ne 1) then ix1=ix0+1
      iy1=iy0 & if (ny ne 1) then iy1=iy0+1
      iz1=iz0 & if (nz ne 1) then iz1=iz0+1
;;  Calculate weight of each particle on the grid.
      for ixx=ix0,ix1 do begin & for iyy=iy0,iy1 do begin & for izz=iz0,iz1 do begin
        weight=1.0
        if (nx ne 1) then weight=weight*( 1.0d - abs(xxp[k,0]-x[ixx])*dx_1 )
        if (ny ne 1) then weight=weight*( 1.0d - abs(xxp[k,1]-y[iyy])*dy_1 )
        if (nz ne 1) then weight=weight*( 1.0d - abs(xxp[k,2]-z[izz])*dz_1 )
        np[ixx,iyy,izz]=np[ixx,iyy,izz]+weight
      endfor & endfor & endfor
; 
    endfor ; end loop over particles
;
  end ; 'cic'
;;
;;  Assign particle density to the grid using the second order TSC method.
;;
  'tsc': begin
;
    if (not quiet) then print, 'Assigning density using TSC method.'
;
    for k=0L,npar-1 do begin
;;  Find nearest grid point     
      ix0=l1 & iy0=m1 & iz0=n1
      if (nx ne 1) then begin
        ix0 = round((xxp[k,0]-x[0])*dx_1)
        if (ix0 eq l2+1) then ix0=ix0-1
        if (ix0 eq l1-1) then ix0=ix0+1
;;  Each particle affects its nearest grid point and the two neighbours of that
;;  grid point in all directions.
        ixx0=ix0-1 & ixx1=ix0+1
      endif else begin
        ixx0=ix0 & ixx1=ix0
      endelse
      if (ny ne 1) then begin
        iy0 = round((xxp[k,1]-y[0])*dy_1)
        if (iy0 eq m2+1) then iy0=iy0-1
        if (iy0 eq m1-1) then iy0=iy0+1
        iyy0=iy0-1 & iyy1=iy0+1
      endif else begin
        iyy0=iy0 & iyy1=iy0
      endelse
      if (nz ne 1) then begin
        iz0 = round((xxp[k,2]-z[0])*dz_1)
        if (iz0 eq n2+1) then iz0=iz0-1
        if (iz0 eq n1-1) then iz0=iz0+1
        izz0=iz0-1 & izz1=iz0+1
      endif else begin
        izz0=iz0 & izz1=iz0
      endelse
;;  Calculate weight of each particle on the grid.
      for ixx=ixx0,ixx1 do begin & for iyy=iyy0,iyy1 do begin & for izz=izz0,   izz1 do begin
        if ( ((ixx-ix0) eq -1) or ((ixx-ix0) eq +1) ) then begin
          weight_x = 1.125d - 1.5d*abs(xxp[k,0]-x[ixx])  *dx_1 + $
                              0.5d*abs(xxp[k,0]-x[ixx])^2*dx_2
        endif else begin
          if (nx ne 1) then $
          weight_x = 0.75d  -         (xxp[k,0]-x[ixx])^2*dx_2
        endelse
        if ( ((iyy-iy0) eq -1) or ((iyy-iy0) eq +1) ) then begin
          weight_y = 1.125d - 1.5d*abs(xxp[k,1]-y[iyy])  *dy_1 + $
                              0.5d*abs(xxp[k,1]-y[iyy])^2*dy_2
        endif else begin
          if (ny ne 1) then $
          weight_y = 0.75d  -         (xxp[k,1]-y[iyy])^2*dy_2
        endelse
        if ( ((izz-iz0) eq -1) or ((izz-iz0) eq +1) ) then begin
          weight_z = 1.125d - 1.5d*abs(xxp[k,2]-z[izz])  *dz_1 + $
                              0.5d*abs(xxp[k,2]-z[izz])^2*dz_2
        endif else begin
          if (nz ne 1) then $
          weight_z = 0.75d  -         (xxp[k,2]-z[izz])^2*dz_2
        endelse
;
        weight=1.0
        if (nx ne 1) then weight=weight*weight_x
        if (ny ne 1) then weight=weight*weight_y
        if (nz ne 1) then weight=weight*weight_z
;
        np[ixx,iyy,izz]=np[ixx,iyy,izz]+weight
      endfor & endfor & endfor
; 
    endfor ; end loop over particles
;
  end ; 'tsc'
;
endcase
;;
;;  Fold density from ghost cells into main array.
;;
if (cic or tsc) then begin
;
  if (param.lshear) then begin
    print, 'pc_particles_to_density: folding not implemented for shear (yet)'
    stop
  endif
;
  if (nz ne 1) then begin
    np[l1-1:l2+1,m1-1:m2+1,n1]= $
        np[l1-1:l2+1,m1-1:m2+1,n1] + np[l1-1:l2+1,m1-1:m2+1,n2+1]
    np[l1-1:l2+1,m1-1:m2+1,n2]= $
        np[l1-1:l2+1,m1-1:m2+1,n2] + np[l1-1:l2+1,m1-1:m2+1,n1-1]
  endif
;
  if (ny ne 1) then begin
    np[l1-1:l2+1,m1,n1:n2]= $
        np[l1-1:l2+1,m1,n1:n2,*] + np[l1-1:l2+1,m2+1,n1:n2]
    np[l1-1:l2+1,m2,n1:n2]= $
        np[l1-1:l2+1,m2,n1:n2,*] + np[l1-1:l2+1,m1-1,n1:n2]
  endif
;
  if (nx ne 1) then begin
    np[l1,m1:m2,n1:n2]=np[l1,m1:m2,n1:n2] + np[l2+1,m1:m2,n1:n2]
    np[l2,m1:m2,n1:n2]=np[l2,m1:m2,n1:n2] + np[l1-1,m1:m2,n1:n2]
  endif
endif
;;
;;  Trim the arrays of ghost zones.
;;
x=x[l1:l2]
y=y[m1:m2]
z=z[n1:n2]
np=np[l1:l2,m1:m2,n1:n2]
;;
;;  Convert from count to density.
;;
if (density) then np=np*cst.rhop_swarm
;;
;;  Normalize the density to have an average of one.
;;
if (normalize) then np=np/mean(np)
;;
;;  Purge missing directions from np before returning it.
;;
return, reform(np)
;
end
