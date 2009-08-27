pro write_particle_series, n0=n0, n1=n1, npar=npar, imgdir=imgdir, $
    random=random, seed=seed, sym=sym

default, n0, 0
default, n1, 10
default, npar, 1000
default, imgdir, '.'
default, random, 0
default, seed, 1
default, sym, 3

for i=n0,n1 do begin
  zeros=''
  if (i lt 100) then begin
    zeros='0'
  endif
  if (i lt 10) then begin
    zeros='00'
  endif
  varfile='PVAR'+strtrim(i,2)
  pc_read_pvar, obj=fp, varfile=varfile, /quiet
 
  if (i eq n0) then begin 
    if (random) then begin
      ii=randomi_norep_aj(seed,npar,n_elements(fp.xx[*,0]))
    endif else begin
      ii=indgen(npar,/long)
    endelse
  endif

  filename='PVAR'+zeros+strtrim(i,2)+'.eps'
  pc_plot_par, fp.xx[ii,*], color=0, sym=sym, $
      filename=filename, imgdir=imgdir, /ps
endfor


end
