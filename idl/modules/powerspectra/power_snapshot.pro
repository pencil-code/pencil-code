;;
;; $Id$
;;
;; Calculate energy spectrum of 3-D cube.
;;
;;   var           : 3-D variable to make power spectrum of
;;   eks           : shell-integrated spectrum
;;   vkx           : 1-D averaged spectrum <f(kx,y,z)>_(y,z)
;;   vky           : 1-D averaged spectrum <f(x,ky,z)>_(x,z)
;;   vkz           : 1-D averaged spectrum <f(x,y,kz)>_(x,y)
;;   ks            : scalar shell wavenumber
;;   kx, ky, kz    : scalar wavenumber in x, y and z
;;   k0x, k0y, k0z : unit wavenumbers
;;   Lx, Ly, Lz    : box dimensions, gives unit wavenumbers k0=2*!pi/L
;;   nks           : number of shells to sum over
;;   deltak        : radius difference of shells
;;
pro power_snapshot, var, eks=eks, vkx=vkx, vky=vky, vkz=vkz, $
    nshell=nshell, ks=ks, kx=kx, ky=ky, kz=kz, $
    k0x=k0x, k0y=k0y, k0z=k0z, Lx=Lx, Ly=Ly, Lz=Lz, deltak=deltak, nks=nks, $
    double=double, plot=plot, ps=ps, filename=filename, $
    nolegend=nolegend, quiet=quiet, debug=debug
;
;  Default values.
;
default, plot, 0
default, ps, 0
default, nolegend, 0
default, quiet, 0
default, debug, 0
;
zero=0.0
one =1.0
if (keyword_set(double)) then begin
  zero=0.0d
  one =1.0d
  var=double(var)
endif
;
sizevar=size(var)
;
if (sizevar[0] ge 1) then nx=sizevar[1] else nx=1
if (sizevar[0] ge 2) then ny=sizevar[2] else ny=1
if (sizevar[0] ge 3) then nz=sizevar[3] else nz=1
;
;  1-D spectrum in x-direction
;
if (arg_present(vkx)) then begin
  vkx=fltarr(nx/2)*zero
  for n=0,nz-1 do begin & for m=0,ny-1 do begin
    AA=fft(reform(var[*,m,n]),-1)
    for j=0,nx/2-1 do vkx[j] = vkx[j] + 2*abs(AA[j])
  endfor & endfor
  vkx=vkx/(ny*nz)
endif
;
;  1-D spectrum in y-direction
;
if (arg_present(vky)) then begin
  vky=fltarr(ny/2)*zero
  for n=0,nz-1 do begin & for l=0,nx-1 do begin
    AA=fft(reform(var[l,*,n]),-1)
    for j=0,ny/2-1 do vky[j] = vky[j] + 2*abs(AA[j])
  endfor & endfor
  vky=vky/(nx*nz)
endif
;
;  1-D spectrum in z-direction
;
if (arg_present(vkz)) then begin
  vkz=fltarr(nz/2)*zero
  for m=0,ny-1 do begin & for l=0,nx-1 do begin
    AA=fft(reform(var[l,m,*]),-1)
    for j=0,nz/2-1 do vkz[j] = vkz[j] + 2*abs(AA[j])
  endfor & endfor
  vkz=vkz/(nx*ny)
endif
;
;  3-D shell-integrated spectrum.
;
if (arg_present(eks)) then begin
;
;  Define directional wavenumbers kx, ky, kz.
;
  if (n_elements(k0x) eq 0) then k0x=one else k0x=k0x*one
  if (n_elements(k0y) eq 0) then k0y=one else k0y=k0y*one
  if (n_elements(k0z) eq 0) then k0z=one else k0z=k0z*one
  if (n_elements(Lx) ne 0)  then k0x=2*!pi/Lx*one
  if (n_elements(Ly) ne 0)  then k0y=2*!pi/Ly*one
  if (n_elements(Lz) ne 0)  then k0z=2*!pi/Lz*one
  kx=[indgen(nx/2+1),-reverse(indgen(nx/2-1)+1)]*k0x
  ky=[indgen(ny/2+1),-reverse(indgen(ny/2-1)+1)]*k0y
  kz=[indgen(nz/2+1),-reverse(indgen(nz/2-1)+1)]*k0z
;
;  Define scalar wavenumber over which to bin.
;
  default, deltak, 1.0*one
  if (n_elements(nks) eq 0) then begin
    ks=indgen(nx/2)*one
  endif else begin
    ks=indgen(nks)*deltak
  endelse
  if (not quiet) then begin
    print, 'k0x, k0y, k0z=', k0x, k0y, k0z
    print, 'Going to sum in shells of radius ks=', ks
  endif
;
;  Define array to hold shell summed power and number of elements in each shell.
;  
  eks=fltarr(n_elements(ks))
  nshell=lonarr(n_elements(ks))
;
;  Transform to wavenumber space.
;
  vkk=fft(var)
;
;  Loop over all vector k, calculate |k| and assign to proper shell.
;
  for ikz=0,nz-1 do begin & for iky=0,ny-1 do begin & for ikx=0,nx-1 do begin
    k=round(sqrt(kx[ikx]^2+ky[iky]^2+kz[ikz]^2)/deltak)
    if (debug) then print, ikx, kx[ikx], iky, ky[iky], ikz, kz[ikz], k
    if (k lt n_elements(ks)) then begin
      eks[k]=eks[k]+abs(vkk[ikx,iky,ikz])^2
      nshell[k]=nshell[k]+1
    endif else begin
      if (debug) then print, '** Mode not allocated to any shell'
    endelse
  endfor & endfor & endfor
;
;  Print total number of modes that were allocated to a shell.
;
  if (not quiet) then begin
    print, 'total(nshell)=', total(nshell), ' / total(modes)=', n_elements(vkk)
  endif
endif
;
;  Make simple plot if requested.
;
if (plot) then begin
  pmin=min(vkx[1:nx/2-1])
  pmax=max(vkx[1:nx/2-1])
  if (ny gt 1) then begin
    pmin=min([pmin,vky[1:ny/2-1]])
    pmax=max([pmax,vky[1:ny/2-1]])
  endif
  if (nz gt 1) then begin
    pmin=min([pmin,vkz[1:nz/2-1]])
    pmax=max([pmax,vkz[1:nz/2-1]])
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
  
  plot, vkx, xtitle='k/k0', ytitle='|g(k)|', $
      xrange=[1.0,nx/2], $
      yrange=[10.0^floor(alog10(pmin)),10.0^ceil(alog10(pmax))], $
      /xlog, /ylog, $
      linestyle=linestyles[0]
  if (ny gt 1) then oplot, vky, linestyle=linestyles[1]
  if (nz gt 1) then oplot, vkz, linestyle=linestyles[2]
;  legend, ['k!Dx!N','k!Dy!N','k!Dz!N'], linestyle=linestyles, /bottom
  if (not nolegend) then legend, ['1','2','3'], linestyle=linestyles, /bottom

  if (ps) then begin
    device, /close
    set_plot, 'x'
  endif
endif
;
end
