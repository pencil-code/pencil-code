;;
;;  Transform from sheared to unsheared frame by shifting data along the
;;  y-direction to match last purely periodic state.
;;
;;  Input:
;;         a : 2-D or 3-D array of scalar, vector, or matrix data.
;;    deltay : azimuthal distance between periodic points at inner and outer
;;             edge of the box. May be calculated from
;;             deltay=-qshear*Omega*Lx*t
;;         x : array of x coordinates
;;
;;
;;  Author : Anders Johansen
;;  Date   : 06/06/2008
;;
function pc_unshear, a, deltay=deltay, x=x, Lx=Lx, Ly=Ly
;
; Default values.
;
default, deltay, 0.0d
default, Lx, 2*!dpi
default, Ly, 2*!dpi
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
if (n_elements(x) ne nx) then begin
  print, 'error: pc_unshear: x must have same number of elements as a[*,0,0]'
  return, a
endif
;
;  Define wavenumbers based on the Nyquist frequency.
;
k0=2*!dpi/Ly
ky_fft=k0*double([indgen(ny/2+1),-reverse(indgen(ny/2-1)+1)])
;
;  Fourier transform along y and interpolate in Fourier space.
;
for ix=0,nx-1 do begin
;
;  Define complex shift array.
;
  deltay_x=-(deltay mod (2*!dpi))*x[ix]/Lx
;
;  Fourier transform along y.
;
  plane_yz=reform(a[ix,*,*,*,*])
  plane_yz_ky=fft(plane_yz,dim=1)
;
;  Shift plane by the amount deltay_x in Fourier space.
;
  plane_yz_ky=plane_yz_ky*spread(exp(complex(0.0d,ky_fft*deltay_x)),[1,2,3],[nz,nv,nm])
;
;  Transform back to real space.
;
  plane_yz=fft(plane_yz_ky,dim=1,/inverse)
  a[ix,*,*,*]=plane_yz
;
endfor
;
return, a
;
end
