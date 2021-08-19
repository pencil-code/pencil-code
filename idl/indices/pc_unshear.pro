;;
;;  $Id$
;;
;;  Transform from sheared to unsheared frame by shifting data along the
;;  y-direction to match last purely periodic state.
;;
;;  This routine is distinct from idl/modules/shear/pc_unshift.
;;  This routine also allows for unshifting in Fourier space.
;;
;;  Author : Anders Johansen
;;  Date   : 06/06/2008
;;
;;  Input:
;;          a : trimmed 2-D or 3-D array of scalar, vector, or matrix data
;;     deltay : azimuthal distance between periodic points at inner and outer
;;              radial edge of the box. May be calculated from
;;                deltay=-Sshear*Lx*t, where Lx=xyz0(0)-xyz1(0)
;;          x : trimmed array of x coordinates
;;      param : structure containing initial conditions, used if deltay is not
;;              provided
;;          t : time of snapshot, used if deltay is not provided
;;    datadir : directory for reading param, used if param is not provided
;;
function pc_unshear, a, deltay=deltay, xax=xax, x0=x0, Lx=Lx, Ly=Ly, $
    param=param, t=t, interpolation_type=interpolation_type, datadir=datadir, $
    nowrap=nowrap
;
; Default values.
;
default, x0, 0.0d
default, interpolation_type, 'fourier'
default, datadir, 'data'
;
; If deltay is not provided, the program can read the relevant information
; from param.nml.
;
if (not is_defined(deltay) or not keyword_set(Lx) or not keyword_set(Ly) ) then begin
  if (not keyword_set(param)) then pc_read_param, obj=param, datadir=datadir, /quiet
  if (not is_defined(t) and not is_defined(deltay)) then begin
    print, 'pc_unshear: must provide t when deltay is not provided'
    return, a
  endif
  if not keyword_set(Lx) then Lx=param.Lxyz[0]
  if not keyword_set(Ly) then Ly=param.Lxyz[1]
  if not is_defined(deltay) then deltay=-param.Sshear*Lx*t
endif
;
; xax must always be provided.
;
if (not keyword_set(xax)) then begin
  print, 'pc_unshear: must provide 1-D array of x-coordinates'
  return, a
endif
;
; Check if a has correct size.
;
s=size(a)
;
if (s[0] lt 2) then begin
  print, 'error: pc_unshear only works for 2-D and 3-D arrays'
  return, a
endif
;
nx=s[1]
ny=s[2]
if (s[0] gt 2) then nz=s[3] else nz=1
if (s[0] gt 3) then nv=s[4] else nv=1
if (s[0] gt 4) then nm=s[5] else nm=1
;
if (n_elements(xax) ne nx) then begin
  print, 'error: pc_unshear: xax must have same number of elements as a[*,0,0]'
  return, a
endif
;
;  The array must be trimmed of ghost cells for the interpolation to be
;  meaningful.
;
;if ( (a[l1,m1,n1] eq a[l1,ny-3,3]) $
;     or (total(a[0:2,0:2,*]) eq 0.0 and max(a) ne 0.0) ) then begin
;  print, 'error: pc_unshear: it seems that a still contains ghost cells.'
;  print, '                   You must trim a before unshearing!'
;  return, a
;endif
;
;  Define wavenumbers based on the Nyquist frequency.
;
k0=2*!dpi/Ly
ky_fft=k0*double([indgen(ny/2+1),-reverse(indgen(ny/2-1)+1)])
;
;  Different interpolation schemes.
;
if (interpolation_type eq 'linear') then begin
;
;  Linear interpolation.
;
  dy=Ly/ny
  for ix=0,nx-1 do begin
    pencil_y=reform(a[ix,*,*,*,*])
;
;  Compute the amount of shift at each x.
;  By default, wrap back each time S*t is an integer,
;  unless /nowrap is given as keyword.
;
    if keyword_set(nowrap) then begin
      deltay_x=deltay*(xax[ix]-x0)/Lx
    endif else begin
      deltay_x=(deltay mod Ly)*(xax[ix]-x0)/Lx
    endelse
;
    deltay_x_int=fix(deltay_x/dy)
    deltay_x_fra=deltay_x/dy-deltay_x_int
;
;  Integer displacement is done as a simple circular shift.
;
    pencil_y=shift(pencil_y,deltay_x_int)
;
;  Fractional displacement by linear interpolation.
;
    pencil_y_ghost=dblarr(ny+2)
    m1=1 & m2=m1+ny-1
    pencil_y_ghost[m1:m2]=pencil_y
    pencil_y_ghost[0]=pencil_y_ghost[m2]
    pencil_y_ghost[m2+1]=pencil_y_ghost[m1]
    for iy=m1,m2 do begin
      if (deltay_x_fra gt 0.0) then begin
        pencil_y[iy-1,*,*,*,*]=pencil_y_ghost[iy,*,*,*,*]+(pencil_y_ghost[iy+1,*,*,*,*]-pencil_y_ghost[iy,*,*,*,*])*deltay_x_fra
      endif else begin
        pencil_y[iy-1,*,*,*,*]=pencil_y_ghost[iy,*,*,*,*]+(pencil_y_ghost[iy,*,*,*,*]-pencil_y_ghost[iy-1,*,*,*,*])*deltay_x_fra
      endelse
    endfor
    a[ix,*,*,*,*]=pencil_y
  endfor
endif else if (interpolation_type eq 'sixth') then begin
;
;  Sixth order polynomial interpolation.
;
  dy=Ly/ny
  for ix=0,nx-1 do begin
;
;  Compute the amount of shift at each x.
;  By default, wrap back each time S*t is an integer,
;  unless /nowrap is given as keyword.
;
    if keyword_set(nowrap) then begin
      deltay_x=deltay*(xax[ix]-x0)/Lx
    endif else begin
      deltay_x=(deltay mod Ly)*(xax[ix]-x0)/Lx
    endelse
;
    deltay_x_int=fix(deltay_x/dy)
    deltay_x_fra=deltay_x/dy-deltay_x_int

    pencil_y=reform(a[ix,*])
;
;  Integer displacement is done as a simple circular shift.
;
    pencil_y=shift(pencil_y,deltay_x_int)
;
;  Fractional displacement by interpolation.
;    
    pencil_y_ghost=dblarr(ny+6)
    m1=3 & m2=m1+ny-1
    pencil_y_ghost[m1:m2]=pencil_y
    pencil_y_ghost[0:2]=pencil_y_ghost[m2-2:m2]
    pencil_y_ghost[m2+1:m2+3]=pencil_y_ghost[m1:m1+2]
    frak=abs(deltay_x_fra)
    c1 = -          (frak+1.)*frak*(frak-1.)*(frak-2.)*(frak-3.)/120.
    c2 = +(frak+2.)          *frak*(frak-1.)*(frak-2.)*(frak-3.)/24.
    c3 = -(frak+2.)*(frak+1.)     *(frak-1.)*(frak-2.)*(frak-3.)/12.
    c4 = +(frak+2.)*(frak+1.)*frak          *(frak-2.)*(frak-3.)/12.
    c5 = -(frak+2.)*(frak+1.)*frak*(frak-1.)          *(frak-3.)/24.
    c6 = +(frak+2.)*(frak+1.)*frak*(frak-1.)*(frak-2.)          /120.
    if (deltay_x_fra lt 0.0) then begin
      for iy=m1,m2 do begin
        pencil_y[iy-3]=c1*pencil_y_ghost[iy-2] + $
                       c2*pencil_y_ghost[iy-1] + $
                       c3*pencil_y_ghost[iy]   + $
                       c4*pencil_y_ghost[iy+1] + $
                       c5*pencil_y_ghost[iy+2] + $
                       c6*pencil_y_ghost[iy+3]
      endfor
    endif else begin
      for iy=m1,m2 do begin
        pencil_y[iy-3]=c1*pencil_y_ghost[iy+2] + $
                       c2*pencil_y_ghost[iy+1] + $
                       c3*pencil_y_ghost[iy]   + $
                       c4*pencil_y_ghost[iy-1] + $
                       c5*pencil_y_ghost[iy-2] + $
                       c6*pencil_y_ghost[iy-3]
      endfor
    endelse
    a[ix,*,*,*,*]=pencil_y
  endfor
endif else if (interpolation_type eq 'fourier') then begin
;
;  Fourier transform along y and interpolate in Fourier space.
;
  for ix=0,nx-1 do begin
;
;  Define complex shift array.
;  By default, wrap back each time S*t is an integer,
;  unless /nowrap is given as keyword.
;
    if keyword_set(nowrap) then begin
      deltay_x=deltay*(xax[ix]-x0)/Lx
    endif else begin
      deltay_x=(deltay mod Ly)*(xax[ix]-x0)/Lx
    endelse
;
;  Fourier transform along y.
;
    plane_yz=reform(a[ix,*,*,*,*])
    plane_yz_ky=fft(plane_yz,dim=1)
;
;  Shift plane by the amount deltay_x in Fourier space.
;
    plane_yz_ky=plane_yz_ky*spread(exp(complex(0.0d,-ky_fft*deltay_x)),[1,2,3],[nz,nv,nm])
;
;  Transform back to real space.
;
    plane_yz=fft(plane_yz_ky,dim=1,/inverse)
    a[ix,*,*,*,*]=plane_yz
;
  endfor
;
endif else begin
  print, 'pc_unshear: no such interpolation type ', interpolation_type
  return, a
endelse
;
return, a
;
end
