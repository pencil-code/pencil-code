;
;  Calculate power spectrum of variable gg along x, y and z directions.
;
pro power_snapshot, gg, g_x=g_x, g_y=g_y, g_z=g_z, plot=plot

nx=n_elements(gg[*,0,0])
ny=n_elements(reform(gg[0,*,0]))
nz=n_elements(reform(gg[0,0,*]))

g_x=fltarr(nx/2)
g_y=fltarr(ny/2)
g_z=fltarr(nz/2)

for m=0,ny-1 do begin
  for n=0,nz-1 do begin
    AA=fft(gg[*,m,n],-1)
    for j=1,nx/2-1 do begin
      g_x[j] = g_x[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
    endfor
  endfor
endfor

for l=0,nx-1 do begin
  for n=0,nz-1 do begin
    AA=fft(gg[l,*,n],-1)
    for j=1,ny/2-1 do begin
      g_y[j] = g_y[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
    endfor
  endfor
endfor

for l=0,nx-1 do begin
  for m=0,ny-1 do begin
    AA=fft(gg[l,m,*],-1)
    for j=1,nz/2-1 do begin
      g_z[j] = g_z[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
    endfor
  endfor
endfor

g_x=g_x/(ny*nz)
g_y=g_y/(nx*nz)
g_z=g_z/(nx*ny)

if (plot) then begin
  pmin=min([g_x[1:nx/2-1],g_y[1:ny/2-1],g_z[1:nz/2-1]])
  pmax=max([g_x[1:nx/2-1],g_y[1:ny/2-1],g_z[1:nz/2-1]])

  linestyles=[0,1,2]
  
  plot, g_x, xtitle='k', ytitle='|g(k)|', $
      xrange=[1.0,nx/2], $
      yrange=[10.0^floor(alog10(pmin)),10.0^ceil(alog10(pmax))], $
      /xlog, /ylog, $
      linestyle=linestyles[0]
  oplot, g_y, linestyle=linestyles[1]
  oplot, g_z, linestyle=linestyles[2]
  legend, ['k!Dx!N','k!Dy!N','k!Dz!N'], linestyle=linestyles, /bottom
endif


end
