;$Id$
;
;  Calculates pdf of all 4 components of alpha_ij, and saves the
;  pdf results for the full array of all 4 components and the fit.
;
;  normalized with respect to etat0, taken from parameters.pro,
;  and this in turn from running kinshear_Re14
;
pc_read_ts,o=o
;
;  set minimum time after which averaging begins
;
t2=1e30
spawn,'touch parameters.pro'
@parameters
@data/index.pro
@data/testfield_info.dat
default,run,''
default,npdf,31
;
;  introduce abbreviations
;
tt=o.t
urms=o.urms
nt=n_elements(tt)
alp=fltarr(nt,2,2)&alpm=fltarr(2,2)&alprms=fltarr(2,2)&alprms_err=fltarr(2,2)
eta=fltarr(nt,2,2)&etam=fltarr(2,2)&etarms=fltarr(2,2)&etarms_err=fltarr(2,2)
xpdf_alp=fltarr(npdf,2,2) & ypdf_alp=fltarr(npdf,2,2)
xpdf_eta=fltarr(npdf,2,2) & ypdf_eta=fltarr(npdf,2,2)
amax=1.6
;
;  alpha tensor
;
alp(*,0,0)=o.alp11
alp(*,1,0)=o.alp21
;
;  eta tensor as in EMF_i = ... -eta_ij J_j
;  so these new etas are also referred to as eta^*
;
eta(*,0,1)=-o.eta11
eta(*,1,1)=-o.eta21
;
;  read extra fields (if itestfield eq 'B11-B22'
;  as opposed to just itestfield eq 'B11-B21')
;
if itestfield eq 'B11-B22' then begin
  alp(*,0,1)=o.alp12
  alp(*,1,1)=o.alp22
  eta(*,0,0)=+o.eta12
  eta(*,1,0)=+o.eta22
endif
;
;  range of time where to do the analysis
;
tmin=min(tt)
tmax=max(tt)
default,t1,(tmin+tmax)/2.
good=where(tt gt t1 and tt lt t2)
kf=5.
;
;  if itestfield eq 'B11-B21' then use special index limits
;
!x.title='!6'
if itestfield eq 'B11-B21' or itestfield eq 'B11-B21+B=0' then begin
  !p.multi=[0,2,2]
  !p.charsize=1.6
  index_max=0
  index_min=1
endif else begin
  !p.multi=[0,2,4]
  !p.charsize=2.6
  index_max=1
  index_min=0
endelse
;
;  give modified alpmax values in parameters.pro file
;
ialpcount=0
for i=0,1 do begin
for j=0,index_max do begin
  !p.title='!7a!6!d'+str(i)+str(j)+'!n'
  alpij=alp(*,i,j)
  pdf,alpij(good)/(etat0*kf),xx,yy,n=npdf
  plot,xx,yy,ps=10
  xpdf_alp(*,i,j)=xx
  ypdf_alp(*,i,j)=yy
  if ialpcount eq 0 then alpgood=alpij(good) else alpgood=[alpgood,alpij(good)]
  ialpcount=ialpcount+1
endfor
endfor
print,'alpmax=',amax
;
for i=0,1 do begin
for j=index_min,1 do begin
  !p.title='!7g!6!d'+str(i)+str(j)+'!n'
  pdf,eta(*,i,j)/etat0,xx,yy,n=npdf
  plot,xx,yy,ps=10
  xpdf_eta(*,i,j)=xx
  ypdf_eta(*,i,j)=yy
endfor
endfor
print,'etamax=',amax
;
wait,.1
!p.multi=0
print,max(abs(alpgood/(etat0*kf)))
amax=0.4
pdf,alpgood/(etat0*kf),xx,yy,n=21,fmin=-amax,fmax=amax
plot_io,xx,yy,ps=10,yr=[.001,3.]
p=moments0(xx,yy,fit=fit,xfit=xfit)
oplot,xfit,fit
;
save,file='alphaeta_pdf.sav',xx,yy,xfit,fit,xpdf_alp,ypdf_alp,xpdf_eta,ypdf_eta
END
