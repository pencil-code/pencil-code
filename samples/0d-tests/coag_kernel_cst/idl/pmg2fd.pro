; $Id: pmg2fd.pro,v 1.1 2003-12-30 13:32:26 ajohan Exp $
!p.charsize=1.5

ndustspec = n_elements(ind)
mg0  = par.mg0
nd00 = par.nd00
deltamg = par.deltamg

mg    = fltarr(ndustspec)
nd    = fltarr(ndustspec)
fd    = fltarr(ndustspec)
fd_an = fltarr(ndustspec)
nd_an = fltarr(ndustspec)
mg(0) = mg0
for i=1,ndustspec-1 do begin
  mg(i)=mg(0)*deltamg^i
endfor

for idust=0,ndustspec-1 do begin
  sdust = strtrim(string(idust),2)
  string = 'nd('+sdust+') = nd'+sdust+'(3,3,3)'
  res = execute(string)
endfor      
fd0 = nd00
;
; Calculate fd
;
for idust=1,ndustspec-2 do begin
  fd(idust) = nd(idust)/(0.5*(mg(idust+1)-mg(idust-1)))
endfor

eta = par.dkern_cst*nd00*t
fr = 1/(1.+eta/2)
fd0_an = nd0/mg(0)
for img=0,ndustspec-1 do begin
  k = mg(img)/mg(0)
  nd_an(img) = nd00*fr^2*(1-fr)^(k-1)
  fd_an(img) = nd_an(img)/mg(0)
endfor

plot, alog10(mg/mg(0)), mg/mg(0)*fd/fd0*mg/mg(0), $
      yrange=[0.0,0.6], xtitle='log(!8m!6!dg!n/!8m!6!d0!n)', $
      ytitle='!8m!6!dg!u2!n!8f!6!dd!n/(!8m!6!d0!u2!n!8f!6!d0!n)', psym=10
oplot, alog10(mg/mg(0)), mg/mg(0)*fd_an/fd0_an*mg/mg(0)


end
