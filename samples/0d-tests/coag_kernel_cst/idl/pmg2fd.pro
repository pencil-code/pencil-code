; $Id$
pc_read_param, obj=par
pc_read_var, obj=ff, /trimall

!p.charsize=1.5

ndustspec = n_elements(ff.nd[0,*])
md0  = par.md0
nd00 = 1.0
nd0  = 1.0
deltamd = par.deltamd

md    = fltarr(ndustspec)
nd    = fltarr(ndustspec)
fd    = fltarr(ndustspec)
fd_an = fltarr(ndustspec)
nd_an = fltarr(ndustspec)
md[0] = md0
for i=1,ndustspec-1 do begin
  md[i]=md[0]*deltamd^i
endfor

for idust=0,ndustspec-1 do begin
  sdust = strtrim(string(idust),2)
  string = 'nd['+sdust+'] = nd'+sdust+'(3,3,3)'
  res = execute(string)
endfor      
fd0 = nd00
;
; Calculate fd
;
for idust=1,ndustspec-2 do begin
  fd[idust] = nd[idust]/(0.5*(md[idust+1]-md[idust-1]))
endfor

eta = par.dkern_cst*nd00*ff.t
fr = 1/(1.+eta/2)
fd0_an = nd0/md(0)
for imd=0,ndustspec-1 do begin
  k = md[imd]/md[0]
  nd_an[imd] = nd00*fr^2*(1-fr)^(k-1)
  fd_an[imd] = nd_an[imd]/md[0]
endfor

plot, alog10(md/md[0]), md/md[0]*fd/fd0*md/md[0], $
      yrange=[0.0,0.6], xtitle='log(!8m!6!dd!n/!8m!6!d0!n)', $
      ytitle='!8m!6!dd!u2!n!8f!6!dd!n/(!8m!6!d0!u2!n!8f!6!d0!n)', psym=10
oplot, alog10(md/md[0]), md/md[0]*fd_an/fd0_an*md/md[0]


end
