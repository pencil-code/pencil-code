;
;  $Id$
;
;  Output particles positions in ascii file.
;
;  Author: Anders Johansen
;
pro pc_particles_to_ascii, xxp, filename=filename, npar=npar, $
    xshift=xshift, yshift=yshift, zshift=zshift, deltay=deltay, $
    lwrite_tauf=lwrite_tauf, param=param, datadir=datadir
;
;  Default values.
;
default, filename, './particles.dat'
default, npar, n_elements(xxp[*,0])
default, xshift, 0.0
default, yshift, 0.0
default, zshift, 0.0
default, deltay, 0.0
default, lwrite_tauf, 0
default, datadir, './data'
;
;  Define box size
;
if (n_elements(param) eq 0) then pc_read_param, obj=param, datadir=datadir, /quiet
x0=param.xyz0[0] & x1=param.xyz1[0] & Lx=param.Lxyz[0]
y0=param.xyz0[1] & y1=param.xyz1[1] & Ly=param.Lxyz[1]
z0=param.xyz0[2] & z1=param.xyz1[2] & Lz=param.Lxyz[2]
;
;  Read friction time from parameter files.
;
if (lwrite_tauf) then begin
  pc_read_pdim, obj=pdim, datadir=datadir, /quiet
  if (max(param.tausp_species) eq 0.0) then begin
;
;  Single friction time.
;
    ipar_fence_species=0
  endif else begin
;
;  Multiple friction times.
;
    npar_species=n_elements(param.tausp_species)
    npar_per_species=pdim.npar/npar_species
    ipar_fence_species=lonarr(npar_species)
    ipar_fence_species[0]=npar_per_species-1
    for jspec=1,npar_species-1 do begin
      ipar_fence_species[jspec]=ipar_fence_species[jspec-1]+npar_per_species
    endfor
  endelse
endif
;
;  Open file for writing.
;
close, 1
openw, 1, filename
;
;  Loop over particles.
;
for ipar=0L,npar-1 do begin
  xxpar=xxp[ipar,*]
;
;  Shift in the x-direction.
;
  if (xshift ne 0.0) then begin
    xxpar[0]=xxpar[0]+xshift
    if (xxpar[0] lt x0) then begin
      xxpar[0]=xxpar[0]+Lx
      if (deltay ne 0.0) then begin
        xxpar[1]=xxpar[1]-deltay
        if (xxpar[1] lt y0) then xxpar[1]=xxpar[1]+Ly
        if (xxpar[1] gt y1) then xxpar[1]=xxpar[1]-Ly
      endif
    endif
    if (xxpar[0] gt x1) then begin
      xxpar[0]=xxpar[0]-Lx
      if (deltay ne 0.0) then begin
        xxpar[1]=xxpar[1]+deltay
        if (xxpar[1] lt y0) then xxpar[1]=xxpar[1]+Ly
        if (xxpar[1] gt y1) then xxpar[1]=xxpar[1]-Ly
      endif
    endif
  endif
;
;  Shift in the y-direction.
;
  if (yshift ne 0.0) then begin
    xxpar[1]=xxpar[1]+yshift
    if (xxpar[1] lt y0) then xxpar[1]=xxpar[1]+Ly
    if (xxpar[1] gt y1) then xxpar[1]=xxpar[1]-Ly
  endif
;
;  Shift in the z-direction.
;
  if (zshift ne 0.0) then begin
    xxpar[2]=xxpar[2]+zshift
    if (xxpar[2] lt z0) then xxpar[2]=xxpar[2]+Lz
    if (xxpar[2] gt z1) then xxpar[2]=xxpar[2]-Lz
  endif
;
;  Write particle friction time to file as well.
;
  if (lwrite_tauf) then begin
    if (max(param.tausp_species) eq 0.0) then begin
      tauf=param.tausp
    endif else begin
      ispec=0
      while (ipar gt ipar_fence_species[ispec]) do begin
        ispec=ispec+1
      endwhile
      tauf=param.tausp_species[ispec]
    endelse
    printf, 1, xxpar, tauf, format='(4f9.4)'
  endif else begin
    printf, 1, xxpar, format='(3f9.4)'
  endelse
endfor
;
;  Close file.
;
close, 1
;
end
