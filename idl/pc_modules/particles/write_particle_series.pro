pro write_particle_series, n0=n0, n1=n1, npar=npar, imgdir=imgdir

default, n0, 0
default, n1, 10
default, npar, 1000
default, imgdir, '.'


for i=n0,n1 do begin
  if (i lt 100) then begin
    zeros='0'
  endif
  if (i lt 10) then begin
    zeros='00'
  endif
  varfile='PVAR'+strtrim(i,2)
  pc_read_pvar, obj=fp, varfile=varfile, /quiet
  filename='PVAR'+zeros+strtrim(i,2)+'.eps'
  pc_plot_par, fp.xx[0:npar-1,*], color=0, $
      filename=filename, imgdir=imgdir, /ps
endfor


end
