;
;  $Id$
;
;  Convert positions and velocities of particles to a grid velocity field.
;
;  Author: Anders Johansen
;
function pc_particles_to_velocity, xxp, vvp, x, y, z, $
    cic=cic, tsc=tsc, datadir=datadir, quiet=quiet

COMMON pc_precision, zero, one
;
;  Set default values.
;
default, cic, 0
default, tsc, 0
default, quiet, 0
;
;  Read run parameters and define auxiliary parameters.
;
pc_set_precision, datadir=datadir
pc_read_param, obj=par, datadir=datadir, /quiet
pc_read_dim  , obj=dim, datadir=datadir, /quiet
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
;
x0=par.xyz0[0] & y0=par.xyz0[1] & z0=par.xyz0[2]
x1=par.xyz1[0] & y1=par.xyz1[1] & z1=par.xyz1[2]
l1=dim.l1 & l2=dim.l2 & mx=dim.mx
m1=dim.m1 & m2=dim.m2 & my=dim.my
n1=dim.n1 & n2=dim.n2 & mz=dim.mz
;
;  Set interpolation scheme
;
interpolation_scheme='ngp'
if (cic) then interpolation_scheme='cic'
if (tsc) then interpolation_scheme='tsc'
;  The CIC and TSC schemes work with ghost cells, so if x, y, z are given
;  without ghost zones, add the ghost zones automatically.
if (cic or tsc) then begin
  if (nx ne mx) then begin
    x2=fltarr(mx)*one
    x2[l1:l2]=x
    for l=l1-1,   0,-1 do x2[l]=x2[l+1]-dx
    for l=l2+1,mx-1,+1 do x2[l]=x2[l-1]+dx
    x=x2
    nx=mx
  endif
  if (ny ne my) then begin
    y2=fltarr(my)*one
    y2[m1:m2]=y
    for m=m1-1,   0,-1 do y2[m]=y2[m+1]-dy
    for m=m2+1,my-1,+1 do y2[m]=y2[m-1]+dy
    y=y2
    ny=my
  endif
  if (nz ne mz) then begin
    z2=fltarr(mz)*one
    z2[n1:n2]=z
    for n=n1-1,   0,-1 do z2[n]=z2[n+1]-dz
    for n=n2+1,mz-1,+1 do z2[n]=z2[n-1]+dz
    z=z2
    nz=mz
  endif
endif
;
;  Possible to map the particles on a finer grid.
;
default, fine, 1

if (fine gt 1) then begin
;
  if (nx gt 1) then begin
    nx=fine*nx
    dx=dx/fine
    x=fltarr(nx)
    x[0]=x0+dx/2
    for i=1,nx-1 do begin
      x[i]=x[0]+i*dx
    endfor
  endif
;
  if (ny gt 1) then begin
    ny=fine*ny
    dy=dy/fine
    y=fltarr(ny)
    y[0]=y0+dy/2
    for i=1,ny-1 do begin
      y[i]=y[0]+i*dy
    endfor
  endif
;
  if (nz gt 1) then begin
    nz=fine*nz
    dz=dz/fine
    z=fltarr(nz)
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
;
;  Define velocity and density array.
;
ww=fltarr(nx,ny,nz,3)*one
rhop=fltarr(nx,ny,nz)*one
;
;  Three different ways to assign particle velocity to the grid are implemented:
;  (see the book by Hockney & Eastwood)
;    0. NGP (Nearest Grid Point)
;    1. CIC (Cloud In Cell)
;    2. TSC (Triangular Shaped Cloud)
;
case interpolation_scheme of
;
;  Assign particle velocity to the grid using the zeroth order NGP method.
;
  'ngp': begin
;
    if (not quiet) then print, 'Assigning velocity using NGP method.'
;
    for k=0L,npar-1 do begin
;  
      ix = round((xxp[k,0]-x[0])*dx_1)
      iy = round((xxp[k,1]-y[0])*dy_1)
      iz = round((xxp[k,2]-z[0])*dz_1)
      if (ix eq nx) then ix=ix-1
      if (iy eq ny) then iy=iy-1
      if (iz eq nz) then iz=iz-1
      if (ix eq -1) then ix=0
      if (iy eq -1) then iy=0
      if (iz eq -1) then iz=0
;
;  Particles are assigned to the nearest grid point.
;
      ww[ix,iy,iz,*]=ww[ix,iy,iz,*]+cst.rhop_tilde*vvp[k,*]
      rhop[ix,iy,iz]=rhop[ix,iy,iz]+cst.rhop_tilde
;
    endfor ; loop over particles
;
  end ; 'ngp'
;
;  Assign particle velocity to the grid using the first order CIC method.
;
  'cic': begin
;
    if (not quiet) then print, 'Assigning velocity using CIC method.'
;
    for k=0L,npar-1 do begin
;  Find nearest grid point     
      ix0 = round((xxp[k,0]-x[0])*dx_1)
      iy0 = round((xxp[k,1]-y[0])*dy_1)
      iz0 = round((xxp[k,2]-z[0])*dz_1)
      if (ix0 eq nx) then ix0=ix0-1
      if (iy0 eq ny) then iy0=iy0-1
      if (iz0 eq nz) then iz0=iz0-1
      if (ix0 eq -1) then ix0=0
      if (iy0 eq -1) then iy0=0
      if (iz0 eq -1) then iz0=0
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
        weight=cst.rhop_tilde
        if (nx ne 1) then weight=weight*( 1.0d - abs(xxp[k,0]-x[ixx])*dx_1 )
        if (ny ne 1) then weight=weight*( 1.0d - abs(xxp[k,1]-y[iyy])*dy_1 )
        if (nz ne 1) then weight=weight*( 1.0d - abs(xxp[k,2]-z[izz])*dz_1 )
        ww[ixx,iyy,izz,*]=ww[ixx,iyy,izz,*]+weight*vvp[k,*]
        rhop[ixx,iyy,izz]=rhop[ixx,iyy,izz]+weight
      endfor & endfor & endfor
; 
    endfor ; end loop over particles
;
  end ; 'cic'
;
;  Assign particle velocity to the grid using the second order TSC method.
;
  'tsc': begin
;
    if (not quiet) then print, 'Assigning velocity using TSC method.'
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
        weight=cst.rhop_tilde
        if (nx ne 1) then weight=weight*weight_x
        if (ny ne 1) then weight=weight*weight_y
        if (nz ne 1) then weight=weight*weight_z
;
        ww[ixx,iyy,izz,*]=ww[ixx,iyy,izz,*]+weight*vvp[k,*]
        rhop[ixx,iyy,izz]=rhop[ixx,iyy,izz]+weight
      endfor & endfor & endfor
; 
    endfor ; end loop over particles
;
  end ; 'tsc'
;
endcase
;
;  Fold velocity from ghost cells into main array.
;
if (cic or tsc) then begin
;
  if (nz ne 1) then begin
    ww[l1-1:l2+1,m1-1:m2+1,n1,*]= $
        ww[l1-1:l2+1,m1-1:m2+1,n1,*] + ww[l1-1:l2+1,m1-1:m2+1,n2+1,*]
    ww[l1-1:l2+1,m1-1:m2+1,n2,*]= $
        ww[l1-1:l2+1,m1-1:m2+1,n2,*] + ww[l1-1:l2+1,m1-1:m2+1,n1-1,*]
    rhop[l1-1:l2+1,m1-1:m2+1,n1]= $
        rhop[l1-1:l2+1,m1-1:m2+1,n1] + rhop[l1-1:l2+1,m1-1:m2+1,n2+1]
    rhop[l1-1:l2+1,m1-1:m2+1,n2]= $
        rhop[l1-1:l2+1,m1-1:m2+1,n2] + rhop[l1-1:l2+1,m1-1:m2+1,n1-1]
  endif
;
  if (ny ne 1) then begin
    ww[l1-1:l2+1,m1,n1:n2,*]= $
        ww[l1-1:l2+1,m1,n1:n2,*] + ww[l1-1:l2+1,m2+1,n1:n2,*]
    ww[l1-1:l2+1,m2,n1:n2,*]= $
        ww[l1-1:l2+1,m2,n1:n2,*] + ww[l1-1:l2+1,m1-1,n1:n2,*]
    rhop[l1-1:l2+1,m1,n1:n2]= $
        rhop[l1-1:l2+1,m1,n1:n2,*] + rhop[l1-1:l2+1,m2+1,n1:n2]
    rhop[l1-1:l2+1,m2,n1:n2]= $
        rhop[l1-1:l2+1,m2,n1:n2,*] + rhop[l1-1:l2+1,m1-1,n1:n2]
  endif
;
  if (nx ne 1) then begin
    ww[l1,m1:m2,n1:n2,*]=ww[l1,m1:m2,n1:n2,*] + ww[l2+1,m1:m2,n1:n2,*]
    ww[l2,m1:m2,n1:n2,*]=ww[l2,m1:m2,n1:n2,*] + ww[l1-1,m1:m2,n1:n2,*]
    rhop[l1,m1:m2,n1:n2]=rhop[l1,m1:m2,n1:n2] + rhop[l2+1,m1:m2,n1:n2]
    rhop[l2,m1:m2,n1:n2]=rhop[l2,m1:m2,n1:n2] + rhop[l1-1,m1:m2,n1:n2]
  endif
endif
;
;  Normalize momentum sum with total density.
;
for iz=0,nz-1 do begin & for iy=0,ny-1 do begin & for ix=0,nx-1 do begin
  if (rhop[ix,iy,iz] ne 0.0) then ww[ix,iy,iz,*]=ww[ix,iy,iz,*]/rhop[ix,iy,iz]
endfor & endfor & endfor
;
;  Trim the array of ghost zones.
;
if (cic or tsc) then ww=ww[l1:l2,m1:m2,n1:n2,*]
;
;  Purge missing directions from ww before returning it.
;
return, reform(ww)
;
end
