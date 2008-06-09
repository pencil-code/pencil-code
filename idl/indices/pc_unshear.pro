;;
;;  $Id: pc_unshear.pro,v 1.2 2008-06-09 14:50:09 ajohan Exp $
;;
;;  Transform from sheared to unsheared frame by shifting data along the
;;  y-direction to match last purely periodic state.
;;
;;  Author : Anders Johansen
;;  Date   : 06/06/2008
;;
;;  Input:
;;          a : trimmed 2-D or 3-D array of scalar, vector, or matrix data
;;     deltay : azimuthal distance between periodic points at inner and outer
;;              radial edge of the box. May be calculated from
;;                deltay=-Sshear*Lx*t
;;          x : trimmed array of x coordinates
;;      param : structure containing initial conditions, used if deltay is not
;;              provided
;;          t : time of snapshot, used if deltay is not provided
;;    datadir : directory for reading param, used if param is not provided
;;
function pc_unshear, a, deltay=deltay, x=x, Lx=Lx, Ly=Ly, $
    param=param, t=t, datadir=datadir
;
; If deltay is not provided, the program can read the relevant information
; from param.
;
if (not keyword_set(deltay)) then begin
  if (not keyword_set(param)) then begin
    if (keyword_set(datadir)) then begin
      pc_read_param, obj=param, datadir=datadir, /quiet
    endif else begin
      print, 'pc_unshear: must provide either param or datadir or {deltay,x,Lx,Ly}'
      return, a
    endelse
  endif
  if (not keyword_set(t)) then begin
    print, 'pc_unshear: must provide t when deltay is not provided'
    return, a
  endif
  Lx=param.Lxyz[0]
  Ly=param.Lxyz[1]
  deltay=-param.Sshear*Lx*t
endif
;
; x must always be provided.
;
if (not keyword_set(x)) then begin
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
if (n_elements(x) ne nx) then begin
  print, 'error: pc_unshear: x must have same number of elements as a[*,0,0]'
  return, a
endif
;
;  The array must be trimmed of ghost cells for the interpolation to be
;  meaningful.
;
if ( (a[3,3,3] eq a[3,3,ny-3]) $
     or (total(a[0:2,0:2,*]) eq 0.0 and max(a) ne 0.0) ) then begin
  print, 'error: pc_unshear: it seems that a still contains ghost cells.'
  print, '                   You must trim a before unshearing!'
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
  a[ix,*,*,*,*]=plane_yz
;
endfor
;
return, a
;
end
