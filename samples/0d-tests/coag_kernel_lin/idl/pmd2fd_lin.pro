; $Id: pmd2fd_lin.pro,v 1.2 2004-02-02 14:22:08 ajohan Exp $
!p.charsize=1.5

ndustspec = n_elements(ind)
md0  = par.md0
nd00 = par.nd00
mdave0 = par.mdave0
deltamd = par.deltamd

md      = fltarr(ndustspec)
mdplus  = fltarr(ndustspec)
mdminus = fltarr(ndustspec)
nd      = fltarr(ndustspec)
rhod    = fltarr(ndustspec)
fd      = fltarr(ndustspec)
fd_an   = fltarr(ndustspec)

for idust=0,ndustspec-1 do begin
  sdust = strtrim(string(idust),2)
  string = 'nd('+sdust+') = nd'+sdust+'(3,3,3)'
  res = execute(string)
  string = 'rhod('+sdust+') = rhod'+sdust+'(3,3,3)'
  res = execute(string)
endfor
;
; Calculate grain masses
;
for i=0,ndustspec-1 do begin
  mdminus(i) = md0*deltamd^i
  mdplus(i) = md0*deltamd^(i+1)
  md(i)=0.5*(mdplus(i)+mdminus(i))
endfor
k = where((nd ne 0.) and (rhod ne 0.))
md(k) = rhod(k)/nd(k)


int=0.
for k=0,ndustspec-1 do begin
  int = int+nd(k)*md(k)
endfor

fd0 = 1.
;
; Calculate fd
;
for i=0,ndustspec-1 do begin
  fd(i) = nd(i)/(mdplus(i)-mdminus(i))
endfor

eta = par.dkern_cst*int*t
fr = exp(-eta)
for i=0,ndustspec-1 do begin
  k = md(i)/md(0)
  fd_an(i) = nd00*fr*exp( -md(i)/mdave0*(1-sqrt(1-fr))^2 ) / $
       ( 2*sqrt(!pi)*mdave0^(-0.5)*md(i)^1.5*(1-fr)^0.75 )
endfor

plot, alog10(md/md(0)), md/md(0)*fd/fd0*md/md(0), $
      yrange=[0.0,0.6], xtitle='log(!8m!6!dd!n/!8m!6!d0!n)', $
      ytitle='!8m!6!dd!u2!n!8f!6!dd!n/(!8m!6!d0!u2!n!8f!6!d0!n)', psym=10
oplot, alog10(md/md(0)), md/md(0)*fd_an*md/md(0)


end
