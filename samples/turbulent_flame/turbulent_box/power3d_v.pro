PRO power3d_v,aa,oplt=oplt,title=title,xtitle=xtitle,ytitle=ytitle,$
  size=size,spectrum=spectrum,wavenumbers=wavenumbers,$
  average=average,debug=debug,noplot=noplot,half=half
;
;  Calculate 3d power spectrum.  The resulting spectrum satisfies
;  total(a^2)=wavenumber(1)*total(spectrum).
;
;  Only handles the common case where all sides of the periodic
;  domain are equal (=size).
;
;  if /half is set, then total(spectrum)=<bb^2>/2, otherwise there is a factor of 2
;
if n_elements(title ) eq 0 then title =''
if n_elements(xtitle) eq 0 then xtitle='Wavenumber'
if n_elements(ytitle) eq 0 then ytitle='Power'

a1=aa(*,*,*,0)
a2=aa(*,*,*,1)
a3=aa(*,*,*,2)
;
s=size(a1) & nx=s[1] & ny=s[2] & nz=s[3]
if n_elements(size) eq 0 then size=2.*!pi
k0=2.*!pi/size

if not keyword_set(noplot) then print,'doing fft ...'
ta1=fft(a1,-1)                                    ; complex 3D FFT
ta1=shift(ta1,nx/2,ny/2,nz/2)                     ; shift to corner
ta2=fft(a2,-1)                                    ; complex 3D FFT
ta2=shift(ta2,nx/2,ny/2,nz/2)                     ; shift to corner
ta3=fft(a3,-1)                                    ; complex 3D FFT
ta3=shift(ta3,nx/2,ny/2,nz/2)                     ; shift to corner

x1=indgen(nx)-nx/2
y1=indgen(ny)-ny/2
z1=indgen(nz)-nz/2
x1=reform(x1,nx,1,1)
y1=reform(y1,1,ny,1)
z1=reform(z1,1,1,nz)
x1=rebin(x1,nx,ny,nz)
y1=rebin(y1,nx,ny,nz)
z1=rebin(z1,nx,ny,nz)
rad=sqrt(x1^2+y1^2+z1^2)                        ; wave vector radius

imax=min([nx,ny,nz]/2)                       ; cut corners
;if keyword_set(nocorners) then $
; imax=min([nx,ny,nz]/2) $                      ; cut corners
;else $
; imax=fix(max(rad))                            ; include corners
spectrum=fltarr(imax)

if not keyword_set(noplot) then print,'looping over shells'
for i=0,imax-1 do begin                         ; cut the corners?
  if keyword_set(debug) then print,i
  ii=where((rad ge i-.5) and (rad lt i+.5))     ; indices into shell
  tmp=float(ta1(ii)*conj(ta1(ii))+$
            ta2(ii)*conj(ta2(ii))+$
            ta3(ii)*conj(ta3(ii)))                ; power of those points
  if keyword_set(half) then tmp=tmp/2.
  if keyword_set(average) then begin
    spectrum(i)=aver(tmp)*((i+1)^3-i^3)*4.*!pi  ; average*volume
  end else $
    spectrum(i)=total(tmp)                      ; sum over shell
endfor

wavenumbers=indgen(imax)*k0                     ; dimensional wavenos
spectrum=spectrum/k0                            ; normalize to k0

if not keyword_set(noplot) then begin
  if n_elements(oplt) eq 0 then begin
    plot_oo,wavenumbers(1:*),spectrum(1:*),$      ; skip DC in log-log
      title=title,xtitle=xtitle,ytitle=ytitle
  endif else begin
    oplot,wavenumbers(1:*),spectrum(1:*),$
      line=oplt
  endelse
endif

return
end
