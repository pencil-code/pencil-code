pro power_snapshot, gg=gg, gx_x=gx_x, gx_y=gx_y, gx_z=gx_z

nx=n_elements(gg[*,0,0])
ny=n_elements(reform(gg[0,*,0]))
nz=n_elements(reform(gg[0,0,*]))

gx_x=fltarr(nx/2)
gx_y=fltarr(ny/2)
gx_z=fltarr(nz/2)

for m=0,ny-1 do begin
  for n=0,nz-1 do begin
    AA=fft(gg[*,m,n],-1)
    for j=1,nx/2-1 do begin
      gx_x[j] = gx_x[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
    endfor
  endfor
endfor

for l=0,nx-1 do begin
  for n=0,nz-1 do begin
    AA=fft(gg[l,*,n],-1)
    for j=1,ny/2-1 do begin
      gx_y[j] = gx_y[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
    endfor
  endfor
endfor

for l=0,nx-1 do begin
  for m=0,ny-1 do begin
    AA=fft(gg[l,m,*],-1)
    for j=1,nz/2-1 do begin
      gx_z[j] = gx_z[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
    endfor
  endfor
endfor

gx_x=gx_x/(ny*nz)
gx_y=gx_y/(nx*nz)
gx_z=gx_z/(nx*ny)


end
