; $Id$
pc_read_var, obj=ff, /trimall

!p.charsize=1.5
int=0.
for k=0,ndustspec-1 do begin
  int = int+ff.nd(k)*ff.md(k)
endfor

fd_an = fltarr(ndustspec)

nd00 = par.nd0
mdave0 = par.mdave0
eta = par.dkern_cst*int*t
fr = exp(-eta)
for i=0,ndustspec-1 do begin
  k = md(i)/md(0)
  fd_an(i) = nd00*fr*exp( -md(i)/mdave0*(1-sqrt(1-fr))^2 ) / $
       ( 2*sqrt(!pi)*mdave0^(-0.5)*md(i)^1.5*(1-fr)^0.75 )
endfor

plot, alog10(md/md(0)), md/md(0)*fd*md/md(0), $
      yrange=[0.0,0.6], xtitle='log(!8m!6!dd!n/!8m!6!d0!n)', $
      ytitle='!8m!6!dd!u2!n!8f!6!dd!n/(!8m!6!d0!u2!n!8f!6!d0!n)', psym=10
oplot, alog10(md/md(0)), md/md(0)*fd_an*md/md(0)


end
