; $Id: pzaver_sph.pro,v 1.24 2016/05/03 20:55:30 rei Exp $
;
@parameters
;
; t1=time of saturation
default, t1, 0
;
; Read grid to compute x and y
pc_read_grid,o=grid,/trimxyz
r=grid.x & theta=grid.y
x=r#sin(theta) & y=r#cos(theta)
;
pc_read_param,o=par2, /param
;
;
pc_read_zaver,o=o, variables=['alp11xy','alp12xy','alp13xy','alp21xy','alp22xy','alp23xy','alp31xy','alp32xy','alp33xy', $
'eta111xy','eta112xy','eta121xy','eta122xy','eta131xy', 'eta132xy','eta211xy','eta212xy','eta221xy','eta222xy','eta231xy',$
'eta232xy','eta311xy','eta312xy','eta321xy','eta322xy','eta331xy','eta332xy', $
'ux2mxy','uy2mxy','uz2mxy','uxmxy','uymxy','uzmxy','TTmxy','rhomxy','oumxy'], njump=njump
;
nr=n_elements(o.alp11xy(*,0,0)) & ntheta=n_elements(o.alp11xy(0,*,0))
nq=n_elements(o.alp11xy(0,0,*))
ttt=o.t

nt0=min(where(ttt ge t1))
ntu=nq-nt0
print, ntu, nq,nt0

fac2=1./ntu
;
aij=fltarr(3,3,nr,ntheta,ntu)
;
aij(0,0,*,*,*)=o.alp11xy[*,*,nt0:nq-1] & aij(0,1,*,*,*)=o.alp12xy[*,*,nt0:nq-1] & aij(0,2,*,*,*)=o.alp13xy[*,*,nt0:nq-1]
aij(1,0,*,*,*)=o.alp21xy[*,*,nt0:nq-1] & aij(1,1,*,*,*)=o.alp22xy[*,*,nt0:nq-1] & aij(1,2,*,*,*)=o.alp23xy[*,*,nt0:nq-1]
aij(2,0,*,*,*)=o.alp31xy[*,*,nt0:nq-1] & aij(2,1,*,*,*)=o.alp32xy[*,*,nt0:nq-1] & aij(2,2,*,*,*)=o.alp33xy[*,*,nt0:nq-1]
;
eijk=fltarr(3,3,2,nr,ntheta,ntu)
;
eijk(0,0,0,*,*,*)=o.eta111xy[*,*,nt0:nq-1] & eijk(0,0,1,*,*,*)=o.eta112xy[*,*,nt0:nq-1]
eijk(0,1,0,*,*,*)=o.eta121xy[*,*,nt0:nq-1] & eijk(0,1,1,*,*,*)=o.eta122xy[*,*,nt0:nq-1]
eijk(0,2,0,*,*,*)=o.eta131xy[*,*,nt0:nq-1] & eijk(0,2,1,*,*,*)=o.eta132xy[*,*,nt0:nq-1]
eijk(1,0,0,*,*,*)=o.eta211xy[*,*,nt0:nq-1] & eijk(1,0,1,*,*,*)=o.eta212xy[*,*,nt0:nq-1]
eijk(1,1,0,*,*,*)=o.eta221xy[*,*,nt0:nq-1] & eijk(1,1,1,*,*,*)=o.eta222xy[*,*,nt0:nq-1]
eijk(1,2,0,*,*,*)=o.eta231xy[*,*,nt0:nq-1] & eijk(1,2,1,*,*,*)=o.eta232xy[*,*,nt0:nq-1]
eijk(2,0,0,*,*,*)=o.eta311xy[*,*,nt0:nq-1] & eijk(2,0,1,*,*,*)=o.eta312xy[*,*,nt0:nq-1]
eijk(2,1,0,*,*,*)=o.eta321xy[*,*,nt0:nq-1] & eijk(2,1,1,*,*,*)=o.eta322xy[*,*,nt0:nq-1]
eijk(2,2,0,*,*,*)=o.eta331xy[*,*,nt0:nq-1] & eijk(2,2,1,*,*,*)=o.eta332xy[*,*,nt0:nq-1]
;
aijt=fltarr(3,3,nr,ntheta) & eijkt=fltarr(3,3,2,nr,ntheta)
;
; take out the zeros and the flollowing point by the resetting

; 1. find smallest zero number
ntest0=min(where(abs(aij[0,0,nr/4,ntheta/4,*]) lt 1e-5))
if (ntest0 ge 10) then ntest0=ntest0 mod 10
; 2. get the intervals
int_test=fix(par2.daainit/par2.d2davg)
; 3. construct a good array
print, 'int_test= ', int_test
print, 'ntest0= ', ntest0
good2=intarr(ntu)

for i=0,ntu/int_test do begin
  if (ntest0 eq 0) then begin
  good2[(int_test-2)*i:(int_test-2)*i+7]=indgen(8)+2+i*int_test
  endif else begin
    if(ntest0 ge 2) then begin
       good2[(int_test-2)*i:(int_test-2)*i+ntest0-2]=indgen(ntest0-1)+1+i*int_test
    endif
    if(ntest0 le 8) then begin
      good2[(int_test-2)*i+ntest0-1:(i+1)*(int_test-2)-1]=indgen(int_test-ntest0-1)+ntest0+2+i*int_test
    endif
  endelse
endfor

bet=where(good2 gt 0)
good=good2[bet]
help, good

fac1=1./n_elements(good)

;
; time averages
;
for i=0,2 do begin
  for j=0,2 do begin
    aijt(i,j,*,*)=fac1*total(aij(i,j,*,*,good),5)
  endfor
endfor
;
for i=0,2 do begin
  for j=0,2 do begin
    for k=0,1 do begin
      eijkt(i,j,k,*,*)=fac1*total(eijk(i,j,k,*,*,good),6)
    endfor
  endfor
endfor
;
; Compute normalizations
cv=0.6 & gamma=5./3 & alp_MLT=5./3.
;
urmsxy=sqrt(o.ux2mxy-o.uxmxy^2+o.uy2mxy-o.uymxy^2+o.uz2mxy-o.uzmxy^2)
ppmxy=fac1*total(o.rhomxy[*,*,nt0:nq-1]*o.TTmxy[*,*,nt0:nq-1]*cv*(gamma-1.),3)
r2sinth=r^2#sin(theta) & ppr=fltarr(nr) & Hp=fltarr(nr,ntheta)

for j=0, ntheta-1 do begin
   Hp[*,j]=-1./(deriv(r,alog(ppmxy[*,j])))
endfor
urmst=fac2*total(urmsxy(*,*,nt0:nq-1),3)
kinhel=fac2*total(o.oumxy(*,*,nt0:nq-1),3)
;
alp0=urmst/3. & etat0=fltarr(nr,ntheta)
for i=0,nr-1 do begin
  etat0(i,*)=urmst(i,*)*alp_MLT*Hp(i,*)/3.
endfor
;
rr=fltarr(nr,ntheta)
for i=0,ntheta-1 do begin 
  rr(*,i)=r(*)
endfor
;
for i=0,2 do $
  for j=0,2 do $
    eijkt(i,j,1,*,*) *= rr      ; MR: because of convention (8) in Schrinner et al. 1997, different from PC 
;
eijkt=-eijkt ; Sign difference w.r.t. Schrinner et al.?
;
; Compute alpha, gamma, beta, delta, and kappa
alpij=fltarr(3,3,nr,ntheta)
alpij(0,0,*,*)=aijt(0,0,*,*)-eijkt(0,1,1,*,*)/rr
alpij(0,1,*,*)=0.5*(aijt(0,1,*,*)+eijkt(0,0,1,*,*)/rr+aijt(1,0,*,*)-eijkt(1,1,1,*,*)/rr)
alpij(0,2,*,*)=0.5*(aijt(0,2,*,*)+aijt(2,0,*,*)-eijkt(2,1,1,*,*)/rr)
alpij(1,0,*,*)=alpij(0,1,*,*)
alpij(1,1,*,*)=aijt(1,1,*,*)+eijkt(1,0,1,*,*)/rr
alpij(1,2,*,*)=0.5*(aijt(1,2,*,*)+aijt(2,1,*,*)+eijkt(2,0,1,*,*)/rr)     ; MR: added factor 1/r 
alpij(2,0,*,*)=alpij(0,2,*,*)
alpij(2,1,*,*)=alpij(1,2,*,*)
alpij(2,2,*,*)=aijt(2,2,*,*)
;
; gamma
;
gami=fltarr(3,nr,ntheta)
gami(0,*,*)=-0.5*(aijt(1,2,*,*)-aijt(2,1,*,*)-eijkt(2,0,1,*,*)/rr)
gami(1,*,*)=-0.5*(aijt(2,0,*,*)-aijt(0,2,*,*)-eijkt(2,1,1,*,*)/rr)
gami(2,*,*)=-0.5*(aijt(0,1,*,*)-aijt(1,0,*,*)+eijkt(0,0,1,*,*)/rr+eijkt(1,1,1,*,*)/rr)   ; MR: reverted sign at eijkt(0,0,1,*,*)/rr
;
; beta
;
betij=fltarr(3,3,nr,ntheta)
betij(0,0,*,*)=-0.5*eijkt(0,2,1,*,*)
betij(0,1,*,*)=0.25*(eijkt(0,2,0,*,*)-eijkt(1,2,1,*,*))
betij(0,2,*,*)=0.25*(eijkt(0,0,1,*,*)-eijkt(2,2,1,*,*)-eijkt(0,1,0,*,*))
betij(1,0,*,*)=betij(0,1,*,*)
betij(1,1,*,*)=0.5*eijkt(1,2,0,*,*)
betij(1,2,*,*)=0.25*(eijkt(1,0,1,*,*)+eijkt(2,2,0,*,*)-eijkt(1,1,0,*,*))
betij(2,0,*,*)=betij(0,2,*,*)
betij(2,1,*,*)=betij(1,2,*,*)
betij(2,2,*,*)=0.5*(eijkt(2,0,1,*,*)-eijkt(2,1,0,*,*))
;
; delta
;
deli=fltarr(3,nr,ntheta)
deli(0,*,*)=0.25*(eijkt(1,1,0,*,*)-eijkt(1,0,1,*,*)+eijkt(2,2,0,*,*))
deli(1,*,*)=0.25*(eijkt(0,0,1,*,*)-eijkt(0,1,0,*,*)+eijkt(2,2,1,*,*))
deli(2,*,*)=-0.25*(eijkt(0,2,0,*,*)+eijkt(1,2,1,*,*))
;
; kappa: a 3-rang tensor
;
kapijk=fltarr(3,3,3,nr,ntheta)
for i=0,2 do begin
   kapijk(i,0,0,*,*)=-eijkt(i,0,0,*,*)
   kapijk(i,0,1,*,*)=-0.5*(eijkt(i,0,1,*,*)+eijkt(i,1,0,*,*))
   kapijk(i,0,2,*,*)=-0.5*eijkt(i,2,0,*,*)
   kapijk(i,1,0,*,*)=kapijk(i,0,1,*,*)
   kapijk(i,1,1,*,*)=-eijkt(i,1,1,*,*)
   kapijk(i,1,2,*,*)=-0.5*eijkt(i,2,1,*,*)
   kapijk(i,2,0,*,*)=kapijk(i,0,2,*,*)
   kapijk(i,2,1,*,*)=kapijk(i,1,2,*,*)
   kapijk(i,2,2,*,*)=1e-9*etat0
endfor
;
; Plots
;
nomit=5 & nl=20 & lev=-1.+2.*findgen(nl)/float(nl-1)
;
!x.margin=[.5,.5] & !y.margin=[0.,2.] & !p.charsize=2.
;
x0=0.09 & x1=0.11 & y0=0.75 & y1=0.89 & dx0=.333 & dy0=0.333
;
alp_lab=[['!7a!8!drr!n!6','!7a!8!d!7h!8r!n!6','!7a!8!d!7u!8r!n!6'],['!7a!8!dr!7h!n!6','!7a!8!d!7hh!n!6','!7a!8!d!7uh!n!6'],['!7a!8!dr!7u!n!6','!7a!8!d!7hu!n!6','!7a!8!d!7uu!n!6']]
;
if !d.name eq 'PS' then begin
  device,file='alpij.eps',xsize=16,ysize=28,yoffset=3,/color,/encapsul,BITS=8
  !p.charthick=2 & !p.thick=2 & !x.thick=2 & !y.thick=2
end
;
if !d.name eq 'X' then window,0,xsize=400,ysize=700
!p.multi=[0,3,3]
;
for i=0,2 do begin
  for j=0,2 do begin
    fac=max(abs(alpij(i,j,*,nomit:ntheta-nomit)/alp0))
    contour,alpij(i,j,*,*)/alp0,x,y,xs=4,ys=4,/fi,/iso,lev=fac*lev,tit=alp_lab(i,j)
    colorbar,range=[min(fac*lev),max(fac*lev)],pos=[x0+i*dx0,y0-j*dy0,x1+i*dx0,y1-j*dy0],/vert,ytickformat='(F7.1)',yticks=2,ytickv=[min(fac*lev),0.,max(fac*lev)],yaxis=0,char=2.5
  endfor
endfor
if !d.name eq 'PS' then device,/close_file
;
; gamma
;
gam_lab=['!7c!8!dr!n!6','!7c!8!d!7h!n!6','!7c!8!d!7u!n!6']
;
if !d.name eq 'PS' then begin
  device,file='gami.eps',xsize=5.333,ysize=28,yoffset=3,/color,/encapsul,BITS=8
  !p.charthick=2 & !p.thick=2 & !x.thick=2 & !y.thick=2
end
;
if !d.name eq 'X' then window,1,xsize=133,ysize=700
!p.multi=[0,1,3]
;
for j=0,2 do begin
  i=0
  fac=max(abs(gami(j,*,nomit:ntheta-nomit)/alp0))
  contour,gami(j,*,*)/alp0,x,y,xs=4,ys=4,/fi,/iso,lev=fac*lev,tit=gam_lab(j)
  colorbar,range=[min(fac*lev),max(fac*lev)],pos=[x0+0.12,y0-j*dy0,x1+0.16,y1-j*dy0],/vert,ytickformat='(F7.1)',yticks=2,ytickv=[min(fac*lev),0.,max(fac*lev)],yaxis=0,char=2.5
endfor
if !d.name eq 'PS' then device,/close_file
;
;contour,gami(0,*,*)/alp0,x,y,/fi,/iso,lev=fac*lev,tit='!7c!8!dr!n!6'
;contour,gami(1,*,*)/alp0,x,y,/fi,/iso,lev=fac*lev,tit='!7c!8!d!7h!n!6'
;contour,gami(2,*,*)/alp0,x,y,/fi,/iso,lev=fac*lev,tit='!7c!8!d!7u!n!6'
;
bet_lab=[['!7b!8!drr!n!6','!7b!8!d!7h!8r!n!6','!7b!8!d!7u!8r!n!6'],['!7b!8!dr!7h!n!6','!7b!8!d!7hh!n!6','!7b!8!d!7uh!n!6'],['!7b!8!dr!7u!n!6','!7b!8!d!7hu!n!6','!7b!8!d!7uu!n!6']]
;
if !d.name eq 'PS' then begin
  device,file='betij.eps',xsize=16,ysize=28,yoffset=3,/color,/encapsul,BITS=8
  !p.charthick=2 & !p.thick=2 & !x.thick=2 & !y.thick=2
end
;
if !d.name eq 'X' then window,2,xsize=400,ysize=700
!p.multi=[0,3,3]
;
for i=0,2 do begin
  for j=0,2 do begin
    fac=max(abs(betij(i,j,*,nomit:ntheta-nomit)/etat0))
    contour,betij(i,j,*,*)/etat0,x,y,xs=4,ys=4,/fi,/iso,lev=fac*lev,tit=bet_lab(i,j)
    colorbar,range=[min(fac*lev),max(fac*lev)],pos=[x0+i*dx0,y0-j*dy0,x1+i*dx0,y1-j*dy0],/vert,ytickformat='(F7.1)',yticks=2,ytickv=[min(fac*lev),0.,max(fac*lev)],yaxis=0,char=2.5
  endfor
endfor
;
if !d.name eq 'PS' then device,/close_file
;
; delta
;
del_lab=['!7d!8!dr!n!6','!7d!8!d!7h!n!6','!7d!8!d!7u!n!6']
;
if !d.name eq 'PS' then begin
  device,file='deli.eps',xsize=5.333,ysize=28,yoffset=3,/color,/encapsul,BITS=8
  !p.charthick=2 & !p.thick=2 & !x.thick=2 & !y.thick=2
end
;
if !d.name eq 'X' then window,3,xsize=133,ysize=700
!p.multi=[0,1,3]
;
for j=0,2 do begin
  i=0
  fac=max(abs(deli(j,*,nomit:ntheta-nomit)/etat0))
  contour,deli(j,*,*)/etat0,x,y,xs=4,ys=4,/fi,/iso,lev=fac*lev,tit=del_lab(j)
  colorbar,range=[min(fac*lev),max(fac*lev)],pos=[x0+0.12,y0-j*dy0,x1+0.16,y1-j*dy0],/vert,ytickformat='(F7.1)',yticks=2,ytickv=[min(fac*lev),0.,max(fac*lev)],yaxis=0,char=2.5
endfor
if !d.name eq 'PS' then device,/close_file
;
;contour,deli(0,*,*)/etat0,x,y,/fi,/iso,lev=fac*lev,tit='!7d!8!dr!n!6'
;contour,deli(1,*,*)/etat0,x,y,/fi,/iso,lev=fac*lev,tit='!7d!8!d!7h!n!6'
;contour,deli(2,*,*)/etat0,x,y,/fi,/iso,lev=fac*lev,tit='!7d!8!d!7u!n!6'
;
; kappa
;
kap_lab=[[['!7j!8!drr!8r!n!6','!7j!8!d!7h!8r!8r!n!6','!7j!8!d!7u!8r!8r!n!6'],['!7j!8!dr!7h!8r!n!6','!7j!8!d!7hh!8r!n!6','!7j!8!d!7uh!8r!n!6'],['!7j!8!dr!7u!8r!n!6','!7j!8!d!7hu!8r!n!6','!7j!8!d!7uu!8r!n!6']],[['!7j!8!drr!7h!n!6','!7j!8!d!7h!8r!7h!n!6','!7j!8!d!7u!8r!7h!n!6'],['!7j!8!dr!7h!7h!n!6','!7j!8!d!7hh!7h!n!6','!7j!8!d!7uh!7h!n!6'],['!7j!8!dr!7u!7h!n!6','!7j!8!d!7hu!7h!n!6','!7j!8!d!7uu!7h!n!6']],[['!7j!8!drr!7u!n!6','!7j!8!d!7h!8r!7u!n!6','!7j!8!d!7u!8r!7u!n!6'],['!7j!8!dr!7h!7u!n!6','!7j!8!d!7hh!7u!n!6','!7j!8!d!7uh!7u!n!6'],['!7j!8!dr!7u!7u!n!6','!7j!8!d!7hu!7u!n!6','!7j!8!d!7uu!7u!n!6']]]

x0=0.030 & x1=0.0366 & y0=0.75 & y1=0.89 & dx0=.111 & dy0=0.333
;
if !d.name eq 'PS' then begin
  device,file='kapijk.eps',xsize=48,ysize=28,yoffset=3,/color,/encapsul,BITS=8
  !p.charthick=2 & !p.thick=2 & !x.thick=2 & !y.thick=2
end
;
if !d.name eq 'X' then window,4,xsize=1200,ysize=700
!p.multi=[0,9,3]
;
for j=0,2 do begin
   for i=0,2 do begin
     for k=0,2 do begin
        fac=max(abs(kapijk(i,j,k,*,nomit:ntheta-nomit)/etat0))
        contour,kapijk(i,j,k,*,*)/etat0,x,y,xs=4,ys=4,/fi,/iso,lev=fac*lev,tit=kap_lab(i,j,k)
        colorbar,range=[min(fac*lev),max(fac*lev)],pos=[x0+i*3*dx0+k*dx0,y0-j*dy0,x1+i*3*dx0+k*dx0,y1-j*dy0],/vert,ytickformat='(F7.1)',yticks=2,ytickv=[min(fac*lev),0.,max(fac*lev)],yaxis=0,char=2.5
     endfor
  endfor
endfor
;
if !d.name eq 'PS' then device,/close_file
;
!p.multi=0
;
save,file='coefs.sav',x,y,urmst,Hp,alp0,etat0,alpij,gami,betij,deli,kapijk,kinhel
;
end
