;
; $Id$
;
; Script to check the implementation of hyperviscosity.
;
; The velocity field is set as a wave with independent wavenumber and amplitude
; for each component. This scripts reads in the value of the hyperviscosity and
; compares with the analytical expectation.
;
pc_read_param, obj=par, /quiet
pc_read_dim, obj=dim, /quiet
nx=dim.nxgrid & ny=dim.nygrid & nz=dim.nzgrid
;
; Define wavenumbers as IDL sees them.
;
kkx=[indgen(nx/2+1),indgen(nx/2-1)-nx/2+1]
kky=[indgen(ny/2+1),indgen(ny/2-1)-ny/2+1]
kkz=[indgen(nz/2+1),indgen(nz/2-1)-nz/2+1]
;
; Read snapshot (at it=1).
;
pc_read_var, obj=ff, /trimall, /quiet
;
; Amplitudes from start.in.
;
ampl_ux=par.ampl_ux[0]
ampl_uy=par.ampl_uy[0]
ampl_uz=par.ampl_uz[0]
ampl_lnrho=par.ampllnrho[0]
;
; Wavenumbers from start.in.
;
kx_ux=par.kx_ux[0]
ky_ux=par.ky_ux[0]
kz_ux=par.kz_ux[0]
kx_uy=par.kx_uy[0]
ky_uy=par.ky_uy[0]
kz_uy=par.kz_uy[0]
kx_uz=par.kx_uz[0]
ky_uz=par.ky_uz[0]
kz_uz=par.kz_uz[0]
kx_lnrho=par.kx_lnrho[0]
ky_lnrho=par.ky_lnrho[0]
kz_lnrho=par.kz_lnrho[0]
;
; Integer version of wavenumbers.
;
ikx_ux=fix(kx_ux)
iky_ux=fix(ky_ux)
ikz_ux=fix(kz_ux)
ikx_uy=fix(kx_uy)
iky_uy=fix(ky_uy)
ikz_uy=fix(kz_uy)
ikx_uz=fix(kx_uz)
iky_uz=fix(ky_uz)
ikz_uz=fix(kz_uz)
ikx_lnrho=fix(kx_lnrho)
iky_lnrho=fix(ky_lnrho)
ikz_lnrho=fix(kz_lnrho)
;
; Tell the user what you found in start.in.
;
print, 'The following wavenumbers and amplitudes were read from start.in:'
print, ' kx_ux,  ky_ux,  kz_ux, ampl_ux = ', ampl_ux, kx_ux, ky_ux, kz_ux, $
    format='(A,4f10.3)'
print, 'ikx_ux, iky_ux, ikz_ux          =       ', ikx_ux, iky_ux, ikz_ux, $
    format='(A,3i10)'
print, ' kx_uy,  ky_uy,  kz_uy, ampl_uy = ', ampl_uy, kx_uy, ky_uy, kz_uy, $
    format='(A,4f10.3)'
print, 'ikx_uy, iky_uy, ikz_uy          =       ', ikx_uy, iky_uy, ikz_uy, $
    format='(A,3i10)'
print, ' kx_uz,  ky_uz,  kz_uz, ampl_uz = ', ampl_uz, kx_uz, ky_uz, kz_uz, $
    format='(A,4f10.3)'
print, 'ikx_uz, iky_uz, ikz_uz          =       ', ikx_uz, iky_uz, ikz_uz, $
    format='(A,3i10)'
print, ''

print, 'I am now going to compare the expected value of the hyperviscosity '
print, 'with what the code has calculated.'
print, ''
;
; Fourier transform of hyperviscosity (without nu).
;
hyvx_k=fft(ff.hyv[*,*,*,0])
hyvy_k=fft(ff.hyv[*,*,*,1])
hyvz_k=fft(ff.hyv[*,*,*,2])
;
; Amplitude of hyv_x at (kx_ux,ky_ux,kz_ux)
;
print, 'Analytical amplitude of hyv_x at (kx,ky,kz)=(', $
    kx_ux, ky_ux, kz_ux, ' ) = ', $
    abs( (-kx_ux^2-ky_ux^2-kz_ux^2)^3*ampl_ux + $
    1/3.0*(-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
    (-1.0)*kx_ux*kx_ux*ampl_ux ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_ux)
      iky=ny/2-j*(ny/2-iky_ux)
      ikz=nz/2-k*(nz/2-ikz_ux)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvx_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvx_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_x at (kx_uy,ky_uy,kz_uy)
;
print, 'Analytical amplitude of hyv_x at (kx,ky,kz)=(', $
    kx_uy, ky_uy, kz_uy, ' ) = ', $
    abs( 1/3.0*(-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
    (-1.0)*kx_uy*ky_uy*ampl_uy ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_uy)
      iky=ny/2-j*(ny/2-iky_uy)
      ikz=nz/2-k*(nz/2-ikz_uy)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvx_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvx_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_x at (kx_uz,ky_uz,kz_uz)
;
print, 'Analytical amplitude of hyv_x at (kx,ky,kz)=(', $
    kx_uz, ky_uz, kz_uz, ' ) = ', $
    abs( (1/3.)*(-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
    (-1.0)*kx_uz*kz_uz*ampl_uz ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_uz)
      iky=ny/2-j*(ny/2-iky_uz)
      ikz=nz/2-k*(nz/2-ikz_uz)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvx_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvx_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_x at (kx_ux+kx_lnrho,ky_ux+ky_lnrho,kz_ux+kz_lnrho)
;
print, 'Analytical amplitude of hyv_x at (kx,ky,kz)=(', $
    kx_ux+kx_lnrho, ky_ux+ky_lnrho, kz_ux+kz_lnrho, ' ) = ', $
    abs( (-kx_ux^2-ky_ux^2-kz_ux^2)^2*4/3.* $
                            kx_ux*ampl_ux*(-1.0)*kx_lnrho*ampl_lnrho/2 + $
         (-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
                            ky_ux*ampl_ux*(-1.0)*ky_lnrho*ampl_lnrho/2 + $
         (-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
                            kz_ux*ampl_ux*(-1.0)*kz_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_ux+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_ux+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_ux+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvx_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvx_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_x at (kx_uy+kx_lnrho,ky_uy+ky_lnrho,kz_uy+kz_lnrho)
;
print, 'Analytical amplitude of hyv_x at (kx,ky,kz)=(', $
    kx_uy+kx_lnrho, ky_uy+ky_lnrho, kz_uy+kz_lnrho, ' ) = ', $
    abs( (-kx_uy^2-ky_uy^2-kz_uy^2)^2*(-2/3.)* $
                            ky_uy*ampl_uy*(-1.0)*kx_lnrho*ampl_lnrho/2 + $
         (-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
                            kx_uy*ampl_uy*(-1.0)*ky_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_uy+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_uy+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_uy+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvx_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvx_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_x at (kx_uz+kx_lnrho,ky_uz+ky_lnrho,kz_uz+kz_lnrho)
;
print, 'Analytical amplitude of hyv_x at (kx,ky,kz)=(', $
    kx_uz+kx_lnrho, ky_uz+ky_lnrho, kz_uz+kz_lnrho, ' ) = ', $
    abs( (-kx_uz^2-ky_uz^2-kz_uz^2)^2*(-2/3.)* $
                            kz_uz*ampl_uz*(-1.0)*kx_lnrho*ampl_lnrho/2 + $
         (-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
                            kx_uz*ampl_uz*(-1.0)*kz_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_uz+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_uz+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_uz+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvx_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvx_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_y at (kx_ux,ky_ux,kz_ux)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_ux, ky_ux, kz_ux, ' ) = ', $
    abs( 1/3.0*(-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
    (-1.0)*ky_ux*kx_ux*ampl_ux ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_ux)
      iky=ny/2-j*(ny/2-iky_ux)
      ikz=nz/2-k*(nz/2-ikz_ux)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvy_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvy_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_y at (kx_uy,ky_uy,kz_uy)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_uy, ky_uy, kz_uy, ' ) = ', $
    abs( (-kx_uy^2-ky_uy^2-kz_uy^2)^3*ampl_uy + $
    1/3.0*(-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
    (-1.0)*ky_uy*ky_uy*ampl_uy ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_uy)
      iky=ny/2-j*(ny/2-iky_uy)
      ikz=nz/2-k*(nz/2-ikz_uy)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvy_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvy_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_y at (kx_uz,ky_uz,kz_uz)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_uz, ky_uz, kz_uz, ' ) = ', $
    abs( 1/3.0*(-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
    (-1.0)*ky_uz*kz_uz*ampl_uz ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_uz)
      iky=ny/2-j*(ny/2-iky_uz)
      ikz=nz/2-k*(nz/2-ikz_uz)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvy_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvy_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_y at (kx_ux+kx_lnrho,ky_ux+ky_lnrho,kz_ux+kz_lnrho)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_ux+kx_lnrho, ky_ux+ky_lnrho, kz_ux+kz_lnrho, ' ) = ', $
    abs( (-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
                            ky_ux*ampl_ux*(-1.0)*kx_lnrho*ampl_lnrho/2 + $
    (-2/3.)*(-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
                            kx_ux*ampl_ux*(-1.0)*ky_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_ux+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_ux+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_ux+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvy_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvy_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_y at (kx_uy+kx_lnrho,ky_uy+ky_lnrho,kz_uy+kz_lnrho)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_uy+kx_lnrho, ky_uy+ky_lnrho, kz_uy+kz_lnrho, ' ) = ', $
    abs( (-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
                            kx_uy*ampl_uy*(-1.0)*kx_lnrho*ampl_lnrho/2 + $
    (+4/3.)*(-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
                            ky_uy*ampl_uy*(-1.0)*ky_lnrho*ampl_lnrho/2 + $
    (-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
                            kz_uy*ampl_uy*(-1.0)*kz_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_uy+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_uy+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_uy+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvy_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvy_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_y at (kx_uz+kx_lnrho,ky_uz+ky_lnrho,kz_uz+kz_lnrho)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_uz+kx_lnrho, ky_uz+ky_lnrho, kz_uz+kz_lnrho, ' ) = ', $
    abs( (-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
                            ky_uz*ampl_uz*(-1.0)*kz_lnrho*ampl_lnrho/2 + $
    (-2/3.)*(-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
                            kz_uz*ampl_uz*(-1.0)*ky_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_uz+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_uz+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_uz+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvy_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvy_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_z at (kx_ux,ky_ux,kz_ux)
;
print, 'Analytical amplitude of hyv_z at (kx,ky,kz)=(', $
    kx_ux, ky_ux, kz_ux, ' ) = ', $
    abs( 1/3.0*(-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
    (-1.0)*kz_ux*kx_ux*ampl_ux ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_ux)
      iky=ny/2-j*(ny/2-iky_ux)
      ikz=nz/2-k*(nz/2-ikz_ux)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvz_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvz_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_z at (kx_uy,ky_uy,kz_uy)
;
print, 'Analytical amplitude of hyv_z at (kx,ky,kz)=(', $
    kx_uy, ky_uy, kz_uy, ' ) = ', $
    abs( 1/3.0*(-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
    (-1.0)*kz_uy*ky_uy*ampl_uy ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_uy)
      iky=ny/2-j*(ny/2-iky_uy)
      ikz=nz/2-k*(nz/2-ikz_uy)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvz_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvz_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_z at (kx_uy,ky_uy,kz_uy)
;
print, 'Analytical amplitude of hyv_z at (kx,ky,kz)=(', $
    kx_uz, ky_uz, kz_uz, ' ) = ', $
    abs( (-kx_uz^2-ky_uz^2-kz_uz^2)^3*ampl_uz + $
    1/3.0*(-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
    (-1.0)*kz_uz*kz_uz*ampl_uz ), format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-ikx_uz)
      iky=ny/2-j*(ny/2-iky_uz)
      ikz=nz/2-k*(nz/2-ikz_uz)
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvz_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvz_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_z at (kx_ux+kx_lnrho,ky_ux+ky_lnrho,kz_ux+kz_lnrho)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_ux+kx_lnrho, ky_ux+ky_lnrho, kz_ux+kz_lnrho, ' ) = ', $
    abs( (-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
                            kz_ux*ampl_ux*(-1.0)*kx_lnrho*ampl_lnrho/2 + $
    (-2/3.)*(-kx_ux^2-ky_ux^2-kz_ux^2)^2* $
                            kx_ux*ampl_ux*(-1.0)*kz_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_ux+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_ux+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_ux+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvz_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvz_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_z at (kx_uy+kx_lnrho,ky_uy+ky_lnrho,kz_uy+kz_lnrho)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_uy+kx_lnrho, ky_uy+ky_lnrho, kz_uy+kz_lnrho, ' ) = ', $
    abs( (-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
                            kz_uy*ampl_uy*(-1.0)*ky_lnrho*ampl_lnrho/2 + $
    (-2/3.)*(-kx_uy^2-ky_uy^2-kz_uy^2)^2* $
                            ky_uy*ampl_uy*(-1.0)*kz_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_uy+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_uy+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_uy+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvz_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvz_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''
;
; Amplitude of hyv_z at (kx_uy+kx_lnrho,ky_uy+ky_lnrho,kz_uy+kz_lnrho)
;
print, 'Analytical amplitude of hyv_y at (kx,ky,kz)=(', $
    kx_uy+kx_lnrho, ky_uy+ky_lnrho, kz_uy+kz_lnrho, ' ) = ', $
    abs( (-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
                            kx_uz*ampl_uz*(-1.0)*kx_lnrho*ampl_lnrho/2 + $
         (-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
                            ky_uz*ampl_uz*(-1.0)*ky_lnrho*ampl_lnrho/2 + $
    (+4/3.)*(-kx_uz^2-ky_uz^2-kz_uz^2)^2* $
                            kz_uz*ampl_uz*(-1.0)*kz_lnrho*ampl_lnrho/2 ), $
    format='(A,3f5.2,A,f10.3)'
tot=0.0
for i=+1,-1,-2 do begin
  for j=+1,-1,-2 do begin
    for k=+1,-1,-2 do begin
      ikx=nx/2-i*(nx/2-(ikx_uz+ikx_lnrho))
      iky=ny/2-j*(ny/2-(iky_uz+iky_lnrho))
      ikz=nz/2-k*(nz/2-(ikz_uz+ikz_lnrho))
      print, 'Numerical amplitude at (ikx,iky,ikz)=(', $
          kkx[ikx], kky[iky], kkz[ikz], ')=', $
          abs(hyvz_k[ikx,iky,ikz]), format='(A,3i5,A,f10.3)'
      tot=tot+abs(hyvz_k[ikx,iky,ikz])
    endfor
  endfor
endfor
print, 'Total amplitude = ', tot, format='(A65,f10.3)'
print, ''

end
