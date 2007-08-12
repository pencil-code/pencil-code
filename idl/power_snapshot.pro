;
;  Calculate power spectrum of variable gg along x, y and z directions.
;
pro power_snapshot, gg=gg, g_x=g_x, g_y=g_y, g_z=g_z, $
                    plot=plot,doubleprec=doubleprec,  $
                    ps=ps, filename=filename, nolegend=nolegend

default, plot, 1
default, nolegend, 0

zero=0.
if keyword_set(doubleprec) then begin
  zero=0.D0
  gg=double(gg)
endif
default, ps, 0

sizegg=size(gg)

nx=1 & ny=1 & nz=1 & g_x=1.0 & g_y=1.0 & g_z=1.0

nx=sizegg[1]
g_x=fltarr(nx/2)*zero
if (sizegg[0] ge 2) then begin
  ny=sizegg[2]
  g_y=fltarr(ny/2)*zero
  if (sizegg[0] ge 3) then begin
    nz=sizegg[3]
    g_z=fltarr(nz/2)*zero
  endif
endif

for m=0,ny-1 do begin
  for n=0,nz-1 do begin
    AA=fft(reform(gg[*,m,n]),-1)
    for j=0,nx/2-1 do begin
      if keyword_set(doubleprec) then begin
        g_x[j] = g_x[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
      endif else begin
        g_x[j] = g_x[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
      endelse
    endfor
  endfor
endfor

for l=0,nx-1 do begin
  for n=0,nz-1 do begin
    AA=fft(reform(gg[l,*,n]),-1)
    for j=0,ny/2-1 do begin
      if keyword_set(doubleprec) then begin
        g_y[j] = g_y[j] + 2*sqrt(double(AA[j])^2+double(imaginary(AA[j]))^2)
      endif else begin
        g_y[j] = g_y[j] + 2*sqrt(double(AA[j])^2+double(imaginary(AA[j]))^2)
      endelse
    endfor
  endfor
endfor

for l=0,nx-1 do begin
  for m=0,ny-1 do begin
    AA=fft(reform(gg[l,m,*]),-1)
    for j=0,nz/2-1 do begin
      if keyword_set(doubleprec) then begin
        g_z[j] = g_z[j] + 2*sqrt(double(AA[j])^2+imaginary(AA[j])^2)
      endif else begin
        g_z[j] = g_z[j] + 2*sqrt(float(AA[j])^2+imaginary(AA[j])^2)
      endelse
    endfor
  endfor
endfor


g_x=g_x/(ny*nz)
g_y=g_y/(nx*nz)
g_z=g_z/(nx*ny)

if (plot) then begin
  pmin=min(g_x[1:nx/2-1])
  pmax=max(g_x[1:nx/2-1])
  if (ny gt 1) then begin
    pmin=min([pmin,g_y[1:ny/2-1]])
    pmax=max([pmax,g_y[1:ny/2-1]])
  endif
  if (nz gt 1) then begin
    pmin=min([pmin,g_z[1:nz/2-1]])
    pmax=max([pmax,g_z[1:nz/2-1]])
  endif

  linestyles=[0,1,2]

  if (ps) then begin
    default, filename, 'power.eps'
    set_plot, 'ps'
    device, /encapsulated, color=color, xsize=8.7, ysize=8.0, $
        font_size=11, filename=filename
    !p.font=-1
    !p.charsize=1.0
    thick=3
    !p.charthick=thick & !p.thick=thick & !x.thick=thick & !y.thick=thick
  endif
  
  plot, g_x, xtitle='k/k0', ytitle='|g(k)|', $
      xrange=[1.0,nx/2], $
      yrange=[10.0^floor(alog10(pmin)),10.0^ceil(alog10(pmax))], $
      /xlog, /ylog, $
      linestyle=linestyles[0]
  if (ny gt 1) then oplot, g_y, linestyle=linestyles[1]
  if (nz gt 1) then oplot, g_z, linestyle=linestyles[2]
;  legend, ['k!Dx!N','k!Dy!N','k!Dz!N'], linestyle=linestyles, /bottom
  if (not nolegend) then legend, ['1','2','3'], linestyle=linestyles, /bottom

  if (ps) then begin
    device, /close
    set_plot, 'x'
  endif
endif


end
