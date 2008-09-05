;$Id$
;
;  This routine plots alpha and eta results from time series.
;  The rms value of the fluctuations (departure from mean)
;  are being written to the file alphaeta_rms.pro, and then
;  added to cvs (this would not work if one is off-line).
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
;
;  introduce abbreviations
;
tt=o.t
urms=o.urms
nt=n_elements(tt)
alp=fltarr(nt,2,2)&alpm=fltarr(2,2)&alprms=fltarr(2,2)&alprms_err=fltarr(2,2)
eta=fltarr(nt,2,2)&etam=fltarr(2,2)&etarms=fltarr(2,2)&etarms_err=fltarr(2,2)
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
  ;!p.multi=[0,2,4]
  !p.multi=[0,2,2]
  ;!p.charsize=2.6
  !p.charsize=1.6
  index_max=1
  index_min=0
endelse
;
;  give modified alpmax values in parameters.pro file
;
for i=0,1 do begin
for j=0,index_max do begin
  !p.title='!7a!6!d'+str(i)+str(j)+'!n'
  pc_fluct_stat,tt,/plo,alp(*,i,j),fm,frms,frms_err,good=good
  alpm(i,j)=fm & alprms(i,j)=frms & alprms_err(i,j)=frms_err
endfor
endfor
;
for i=0,1 do begin
for j=index_min,1 do begin
  !p.title='!7g!6!d'+str(i)+str(j)+'!n'
  pc_fluct_stat,tt,/plo,eta(*,i,j),fm,frms,frms_err,good=good
  etam(i,j)=fm & etarms(i,j)=frms & etarms_err(i,j)=frms_err
endfor
endfor
;
print,alprms
print
print,etarms
print
;
if itestfield eq 'B11-B22' then begin
  alprms_all=sqrt(.25*total(alprms^2))
  etarms_all=sqrt(.25*total(etarms^2))
  etatrms_all=sqrt(.5*etarms(0,0)^2+etarms(1,1)^2)
  eta21rms_all=etarms(1,0)
  eta12rms_all=etarms(0,1)
  alprms_err_all=sqrt(.25*total(alprms_err^2))
  etarms_err_all=sqrt(.25*total(etarms_err^2))
  etatrms_err_all=sqrt(.5*etarms_err(0,0)^2+etarms_err(1,1)^2)
  eta21rms_err_all=etarms_err(1,0)
  eta12rms_err_all=etarms_err(0,1)
endif else begin
  alprms_all=sqrt(.5*total(alprms(*,0)^2))
  etarms_all=sqrt(.5*total(etarms(*,1)^2))
  etatrms_all=etarms(1,1)
  eta21rms_all=etarms(1,0)
  eta12rms_all=etarms(1,0)
  alprms_err_all=sqrt(.5*total(alprms_err(*,0)^2))
  etarms_err_all=sqrt(.5*total(etarms_err(*,1)^2))
  etatrms_err_all=etarms_err(1,1)
  eta21rms_err_all=etarms_err(1)
  eta12rms_err_all=0.
endelse
;
fo='(a,e8.2)'
openw,1,'alphaeta_rms.pro'
printf,1,'alprms_all=',alprms_all,fo=fo
printf,1,'etarms_all=',etarms_all,fo=fo
printf,1,'etatrms_all=',etatrms_all,fo=fo
printf,1,'eta21rms_all=',eta21rms_all,fo=fo
printf,1,'eta12rms_all=',eta12rms_all,fo=fo
printf,1,'alprms_err_all=',alprms_err_all,fo=fo
printf,1,'etarms_err_all=',etarms_err_all,fo=fo
printf,1,'etatrms_err_all=',etatrms_err_all,fo=fo
printf,1,'eta21rms_err_all=',eta21rms_err_all,fo=fo
printf,1,'eta12rms_err_all=',eta12rms_err_all,fo=fo
close,1
;
spawn,'cat alphaeta_rms.pro'
spawn,'cvs add alphaeta_rms.pro'
END
