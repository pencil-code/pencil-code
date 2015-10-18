;
;  $Id$
;
;  Calculate the gas velocity at the position of particles by interpolation from
;  nearest grid points.
;
;  Author: Anders Johansen
;
function pc_gas_velocity_at_particle, xxp, uu, x, y, z, $
    cic=cic, tsc=tsc, datadir=datadir, quiet=quiet, par=par, dim=dim, $
    ixmin=ixmin, ixmax=ixmax, iymin=iymin, iymax=iymax, $
    izmin=izmin, izmax=izmax

common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
;  Set default values.
;
default, cic, 0
default, tsc, 0
default, quiet, 0
default, ixmin, -1
default, ixmax, -1
default, iymin, -1
default, iymax, -1
default, izmin, -1
default, izmax, -1
;
;  Read run parameters and define auxiliary parameters.
;
if (not arg_present(dim)) then pc_read_dim, obj=dim, datadir=datadir, /quiet
if (not arg_present(par)) then pc_read_param, obj=par, dim=dim, datadir=datadir, /quiet
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
;  The CIC and TSC schemes work only with ghost cells in x, y, z and uu.
if (cic or tsc) then begin
  sizeuu=size(uu)
  if (nx ne mx or ny ne my or nz ne mz) then begin
    print, 'ERROR: x, y, z arrays must contain ghost zones for cic and tsc interpolation'
    stop
  endif
  if (sizeuu[1] ne mx or sizeuu[2] ne my or sizeuu[3] ne mz) then begin
    print, 'ERROR: uu must contain ghost zones for cic and tsc interpolation'
    stop
  endif
endif
;
;  Define gas velocity array.
;
uup=fltarr(npar,3)*one
;
;  Three different ways to interpolate gas velocity to the positions of
;  particles.
;  (see the book by Hockney & Eastwood)
;    0. NGP (Nearest Grid Point)
;    1. CIC (Cloud In Cell)
;    2. TSC (Triangular Shaped Cloud)
;
case interpolation_scheme of
;
;  Assign gas velocity to particles using the zeroth order NGP method.
;
  'ngp': begin
;
    if (not quiet) then print, 'Assigning gas velocity using NGP method.'
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
;  Gas velocity assigned from the nearest grid point.
;
      uup[k,*]=uu[ix,iy,iz,*]
;
    endfor ; loop over particles
;
  end ; 'ngp'
;
;  Assign gas velocity to particles using the first order CIC method.
;
  'cic': begin
;
    if (not quiet) then print, 'Assigning gas velocity using CIC method.'
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
;  Calculate weight of each grid point.
      for ixx=ix0,ix1 do begin & for iyy=iy0,iy1 do begin & for izz=iz0,iz1 do begin
        weight=1.0d
        if (nx ne 1) then weight=weight*( 1.0d - abs(xxp[k,0]-x[ixx])*dx_1 )
        if (ny ne 1) then weight=weight*( 1.0d - abs(xxp[k,1]-y[iyy])*dy_1 )
        if (nz ne 1) then weight=weight*( 1.0d - abs(xxp[k,2]-z[izz])*dz_1 )
        uup[k,*]=uup[k,*]+weight*uu[ixx,iyy,izz,*]
      endfor & endfor & endfor
; 
    endfor ; end loop over particles
;
  end ; 'cic'
;
;  Assign gas velocity to particles using the second order TSC method.
;
  'tsc': begin
;
    if (not quiet) then print, 'Assigning gas velocity using TSC method.'
;
    for k=0L,npar-1 do begin
;  Find nearest grid point     
      ix0=l1 & iy0=m1 & iz0=n1
      if (nx ne 1) then ix0 = round((xxp[k,0]-x[0])*dx_1)
      if (ny ne 1) then iy0 = round((xxp[k,1]-y[0])*dy_1)
      if (nz ne 1) then iz0 = round((xxp[k,2]-z[0])*dz_1)
      if (ix0 eq nx) then ix0=ix0-1
      if (iy0 eq ny) then iy0=iy0-1
      if (iz0 eq nz) then iz0=iz0-1
      if (ix0 eq -1) then ix0=0
      if (iy0 eq -1) then iy0=0
      if (iz0 eq -1) then iz0=0
;  Because IDL is slow, allow the user to specify only a range of grid points
;  to consider.
      lcontinue=1
      if (ixmin gt -1 and ixmax gt -1 and $
          iymin gt -1 and iymax gt -1 and $
          izmin gt -1 and izmax gt -1) then begin
        if (ix0 ge ixmin and ix0 le ixmax and $
            iy0 ge iymin and iy0 le iymax and $
            iz0 ge izmin and iz0 le izmax) then begin
          lcontinue=1
        endif else begin
          lcontinue=0
        endelse
      endif
;
      if (lcontinue) then begin
;  Each particle is affected by the nearest grid point and the two neighbours
;  of that grid point in all directions.
        ixx0=ix0-1 & ixx1=ix0+1
        iyy0=iy0-1 & iyy1=iy0+1
        izz0=iz0-1 & izz1=iz0+1
        if (nx eq 1) then begin
          ixx0=ix0 & ixx1=ix0
        endif
        if (ny eq 1) then begin
          iyy0=iy0 & iyy1=iy0
        endif
        if (nz eq 1) then begin
          izz0=iz0 & izz1=iz0
        endif
;  Calculate weight of each grid point at the particle.
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
          weight=1.0
          if (nx ne 1) then weight=weight*weight_x
          if (ny ne 1) then weight=weight*weight_y
          if (nz ne 1) then weight=weight*weight_z
;
          uup[k,*]=uup[k,*]+weight*uu[ixx,iyy,izz,*]
        endfor & endfor & endfor
      endif
; 
    endfor ; end loop over particles
;
  end ; 'tsc'

endcase
;
;  Return gas velocity.
;
return, uup
;
end
