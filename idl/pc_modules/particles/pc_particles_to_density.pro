;
;  $Id: pc_particles_to_density.pro,v 1.23 2007-08-03 13:33:49 ajohan Exp $
;
;  Convert positions of particles to a grid density field.
;
;  Author: Anders Johansen
;
function pc_particles_to_density, xxp, x, y, z, $
    cic=cic, tsc=tsc, fine=fine, normalize=normalize, $
    datadir=datadir, vvp=vvp, vprms=vprms, lscalar_rms=lscalar_rms, $
    quiet=quiet

COMMON pc_precision, zero, one
;
;  Set default values.
;
default, cic, 0
default, tsc, 0
default, lscalar_rms, 0
default, quiet, 0
default, normalize, 0
;
;  Read run parameters and define auxiliary parameters.
;
pc_set_precision, datadir=datadir
pc_read_param, obj=par, datadir=datadir, /quiet
pc_read_dim  , obj=dim, datadir=datadir, /quiet
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
;
x0=par.xyz0[0] & y0=par.xyz0[1] & z0=par.xyz0[2]
x1=par.xyz1[0] & y1=par.xyz1[1] & z1=par.xyz1[2]
;
;  Set interpolation scheme
;
interpolation_scheme='ngp'
if (cic) then interpolation_scheme='cic'
if (tsc) then interpolation_scheme='tsc'
;
;  Possible to map the particles on a finer grid.
;
default, fine, 1

if (fine gt 1) then begin
;
  if (nx gt 1) then begin
    nx=fine*nx
    dx=dx/fine
    x=fltarr(nx)*one
    x[0]=x0+dx/2
    for i=1,nx-1 do begin
      x[i]=x[0]+i*dx
    endfor
  endif
;
  if (ny gt 1) then begin
    ny=fine*ny
    dy=dy/fine
    y=fltarr(ny)*one
    y[0]=y0+dy/2
    for i=1,ny-1 do begin
      y[i]=y[0]+i*dy
    endfor
  endif
;
  if (nz gt 1) then begin
    nz=fine*nz
    dz=dz/fine
    z=fltarr(nz)*one
    z[0]=z0+dz/2
    for i=1,nz-1 do begin
      z[i]=z[0]+i*dz
    endfor
  endif
;
  dx_1=1.0d/dx   & dy_1=1.0d/dy   & dz_1=1.0d/dz
  dx_2=1.0d/dx^2 & dy_2=1.0d/dy^2 & dz_2=1.0d/dz^2
;
endif
;  The CIC and TSC schemes work with ghost cells, so if x, y, z are given
;  without ghost zones, add the ghost zones automatically.
if (cic or tsc) then begin
  nghost=1

  mx=nx+2*nghost
  l1=nghost & l2=l1+nx-1
  my=ny+2*nghost
  m1=nghost & m2=m1+ny-1
  mz=nz+2*nghost
  n1=nghost & n2=n1+nz-1

  if (nx ne 1) then begin
    x2=fltarr(mx)*one
    x2[l1:l2]=x
    for l=l1-1,   0,-1 do x2[l]=x2[l+1]-dx
    for l=l2+1,mx-1,+1 do x2[l]=x2[l-1]+dx
    x=x2
    nx=mx
  endif
  if (ny ne 1) then begin
    y2=fltarr(my)*one
    y2[m1:m2]=y
    for m=m1-1,   0,-1 do y2[m]=y2[m+1]-dy
    for m=m2+1,my-1,+1 do y2[m]=y2[m-1]+dy
    y=y2
    ny=my
  endif
  if (nz ne 1) then begin
    z2=fltarr(mz)*one
    z2[n1:n2]=z
    for n=n1-1,   0,-1 do z2[n]=z2[n+1]-dz
    for n=n2+1,mz-1,+1 do z2[n]=z2[n-1]+dz
    z=z2
    nz=mz
  endif
endif else begin
  nghost=0
  mx=nx
  l1=0 & l2=l1+nx-1
  my=ny
  m1=0 & m2=m1+ny-1
  mz=nz+2*nghost
  n1=0 & n2=n1+nz-1
endelse
;
;  Define density and velocity dispersion arrays.
;
np=fltarr(mx,my,mz)*one

if (arg_present(vprms)) then begin
  if (n_elements(vvp) eq 0) then begin
    print, 'pc_particles_to_density: ERROR - '+ $
        'vvp must be supplied to calculate vprms!'
    return, 1
  endif
  vvpm=fltarr(nx,ny,nz,3)*one
  if (lscalar_rms) then begin
    vprms=fltarr(nx,ny,nz)*one
  endif else begin
    vvp2m=fltarr(nx,ny,nz,3)*one
    vprms=fltarr(nx,ny,nz,3)*one
  endelse
endif
;
;  Three different ways to assign particle density to the grid are implemented:
;  (see the book by Hockney & Eastwood)
;    0. NGP (Nearest Grid Point)
;    1. CIC (Cloud In Cell)
;    2. TSC (Triangular Shaped Cloud)
;
case interpolation_scheme of
;
;  Assign particle density to the grid using the zeroth order NGP method.
;
  'ngp': begin
;
    if (not quiet) then print, 'Assigning density using NGP method.'
;  Nearest grid point is stored in array.
    ineargrid=intarr(npar,3)
;
    for k=0L,npar-1 do begin
;  
      ix = round((xxp[k,0]-x[0])*dx_1)
      iy = round((xxp[k,1]-y[0])*dy_1)
      iz = round((xxp[k,2]-z[0])*dz_1)
      if (ix gt l2) then ix=l2
      if (ix lt l1) then ix=l1
      if (iy gt m2) then iy=m2
      if (iy lt m1) then iy=m1
      if (iz gt n2) then iz=n2
      if (iz lt n1) then iz=n1
      ineargrid[k,*]=[ix,iy,iz]
;
;  Particles are assigned to the nearest grid point.
;
      np[ix,iy,iz]=np[ix,iy,iz]+1.0*one
;  Mean velocity (needed for velocity dispersion).
      if (arg_present(vprms)) then begin
        vvpm[ix,iy,iz,0]  = vvpm[ix,iy,iz,0] + vvp[k,0]
        vvpm[ix,iy,iz,1]  = vvpm[ix,iy,iz,1] + vvp[k,1]
        vvpm[ix,iy,iz,2]  = vvpm[ix,iy,iz,2] + vvp[k,2]
;        if (not lscalar_rms) then begin
          vvp2m[ix,iy,iz,0] = vvp2m[ix,iy,iz,0] + vvp[k,0]^2
          vvp2m[ix,iy,iz,1] = vvp2m[ix,iy,iz,1] + vvp[k,1]^2
          vvp2m[ix,iy,iz,2] = vvp2m[ix,iy,iz,2] + vvp[k,2]^2
;        endif
      endif
;
    endfor ; loop over particles
;
;  Calculate velocity dispersions.
;
    if (arg_present(vprms)) then begin
      ii=array_indices(np,where(np gt 1.0))
      ii3=intarr(n_elements(ii[0,*]))
;  Normalize sums to get average.
      for j=0,2 do begin
        ii3[*]=j
        vvpm[ii[0,*],ii[1,*],ii[2,*],ii3] = $
            vvpm[ii[0,*],ii[1,*],ii[2,*],ii3]/np[ii[0,*],ii[1,*],ii[2,*]]
;        if (not lscalar_rms) then $
            vvp2m[ii[0,*],ii[1,*],ii[2,*],ii3] = $
            vvp2m[ii[0,*],ii[1,*],ii[2,*],ii3]/np[ii[0,*],ii[1,*],ii[2,*]]
      endfor
;  Scalar velocity dispersion.
;      if (lscalar_rms) then begin
;        for k=0L,npar-1 do begin
;          ix=ineargrid[k,0] & iy=ineargrid[k,1] & iz=ineargrid[k,2]
;          vprms[ix,iy,iz]=vprms[ix,iy,iz]+ $
;              total((vvp[k,*]-vvpm[ix,iy,iz,*])^2,2)
;        endfor
;        ii=where(np gt 1.0)
;        vprms[ii] = sqrt(vprms[ii]/np[ii])
;      endif else begin
;  Vector velocity dispersion.
        vprms=vvp2m-vvpm^2
;        i0=where(vprms lt 0.0)
;        if (i0[0] ne -1) then vprms[i0]=0.0
        vprms=sqrt(vprms)
;      endelse
    endif
;
  end ; 'ngp'
;
;  Assign particle density to the grid using the first order CIC method.
;
  'cic': begin
;
    if (not quiet) then print, 'Assigning density using CIC method.'
;
    for k=0L,npar-1 do begin
;  Find nearest grid point     
      ix0 = round((xxp[k,0]-x[0])*dx_1)
      iy0 = round((xxp[k,1]-y[0])*dy_1)
      iz0 = round((xxp[k,2]-z[0])*dz_1)
      if (ix0 gt l2) then ix0=l2
      if (ix0 lt l1) then ix0=l1
      if (iy0 gt m2) then iy0=m2
      if (iy0 lt m1) then iy0=m1
      if (iz0 gt n2) then iz0=n2
      if (iz0 lt n1) then iz0=n1
;  Find lower grid point in surrounding grid points.        
      if ( (x[ix0] gt xxp[k,0]) and (nx ne 1) ) then ix0=ix0-1
      if ( (y[iy0] gt xxp[k,1]) and (ny ne 1) ) then iy0=iy0-1
      if ( (z[iz0] gt xxp[k,2]) and (nz ne 1) ) then iz0=iz0-1
;  Don't assign density to degenerate directions. 
      ix1=ix0 & if (nx ne 1) then ix1=ix0+1
      iy1=iy0 & if (ny ne 1) then iy1=iy0+1
      iz1=iz0 & if (nz ne 1) then iz1=iz0+1
;  Calculate weight of each particle on the grid.
      for ixx=ix0,ix1 do begin & for iyy=iy0,iy1 do begin & for izz=iz0,iz1 do begin
        weight=1.0d
        if (nx ne 1) then weight=weight*( 1.0d - abs(xxp[k,0]-x[ixx])*dx_1 )
        if (ny ne 1) then weight=weight*( 1.0d - abs(xxp[k,1]-y[iyy])*dy_1 )
        if (nz ne 1) then weight=weight*( 1.0d - abs(xxp[k,2]-z[izz])*dz_1 )
        np[ixx,iyy,izz]=np[ixx,iyy,izz]+weight
      endfor & endfor & endfor
; 
    endfor ; end loop over particles
;
  end ; 'cic'
;
;  Assign particle density to the grid using the second order TSC method.
;
  'tsc': begin
;
    if (not quiet) then print, 'Assigning density using TSC method.'
;
    for k=0L,npar-1 do begin
;  Find nearest grid point     
      ix0=l1 & iy0=m1 & iz0=n1
      if (nx ne 1) then begin
        ix0 = round((xxp[k,0]-x[0])*dx_1)
        if (ix0 gt l2) then ix0=l2
        if (ix0 lt l1) then ix0=l1
;  Each particle affects its nearest grid point and the two neighbours of that
;  grid point in all directions.
        ixx0=ix0-1 & ixx1=ix0+1
      endif else begin
        ixx0=ix0 & ixx1=ix0
      endelse
      if (ny ne 1) then begin
        iy0 = round((xxp[k,1]-y[0])*dy_1)
        if (iy0 gt m2) then iy0=m2
        if (iy0 lt m1) then iy0=m1
        iyy0=iy0-1 & iyy1=iy0+1
      endif else begin
        iyy0=iy0 & iyy1=iy0
      endelse
      if (nz ne 1) then begin
        iz0 = round((xxp[k,2]-z[0])*dz_1)
        if (iz0 gt n2) then iz0=n2
        if (iz0 lt n1) then iz0=n1
        izz0=iz0-1 & izz1=iz0+1
      endif else begin
        izz0=iz0 & izz1=iz0
      endelse
;  Calculate weight of each particle on the grid.
      for ixx=ixx0,ixx1 do begin & for iyy=iyy0,iyy1 do begin & for izz=izz0,izz1 do begin
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
        weight=1.0d
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

endcase
;
;  Fold density from ghost cells into main array.
;
if (cic or tsc) then begin
;
  if (nz ne 1) then begin
    np[l1-1:l2+1,m1-1:m2+1,n1]= $
        np[l1-1:l2+1,m1-1:m2+1,n1] + np[l1-1:l2+1,m1-1:m2+1,n2+1]
    np[l1-1:l2+1,m1-1:m2+1,n2]= $
        np[l1-1:l2+1,m1-1:m2+1,n2] + np[l1-1:l2+1,m1-1:m2+1,n1-1]
  endif
;
  if (ny ne 1) then begin
    np[l1-1:l2+1,m1,n1:n2]=np[l1-1:l2+1,m1,n1:n2] + np[l1-1:l2+1,m2+1,n1:n2]
    np[l1-1:l2+1,m2,n1:n2]=np[l1-1:l2+1,m2,n1:n2] + np[l1-1:l2+1,m1-1,n1:n2]
  endif
;
  if (nx ne 1) then begin
    np[l1,m1:m2,n1:n2]=np[l1,m1:m2,n1:n2] + np[l2+1,m1:m2,n1:n2]
    np[l2,m1:m2,n1:n2]=np[l2,m1:m2,n1:n2] + np[l1-1,m1:m2,n1:n2]
  endif
;
;  Trim the array of ghost zones.
;
  np=np[l1:l2,m1:m2,n1:n2]
endif
;
;  Normalize the density to have an average of one.
;
if (normalize) then np=np/mean(np)
;
;  Purge missing directions from np before returning it.
;
return, reform(np)
;
end
