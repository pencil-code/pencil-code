;
; $Id$
;
;  preliminary version for reading point distribution function (pdf) files
;  (at the time of check-in on horseshoe I don't remember details
;  of when this routine is really used, but since horseshoe will
;  shut down I better check this in as it is...)
;
pc_read_dim,nprocx=nprocx,nprocy=nprocy,nprocz=nprocz,/quiet
nproc=nprocx*nprocy*nprocz
print,'nproc=',nproc
;
t=0. & n_pdf=0L & pdf_dx=0. & pdf_min=0. & pdf_mean=0. & pdf_rms=0.
;
;  open all files
;  note: this method only works on up to 128 processors
;
for iproc=0,nproc-1 do begin
  lun=iproc+1
  filename='data/proc'+str(iproc)+'/pdf_rhocc.dat'
  filename='data/proc'+str(iproc)+'/pdf_rhocc.dat'
  filename='data/proc'+str(iproc)+'/pdf_lncc.dat'
  filename='data/proc'+str(iproc)+'/pdf_lngcc.dat'
  filename='data/proc'+str(iproc)+'/pdf_gcc.dat'
  filename='data/proc'+str(iproc)+'/pdf_cc.dat'
  close,lun
  openr,lun,filename
endfor
;
nt=1000
for it=0,nt-1 do begin
for iproc=0,nproc-1 do begin
  lun=iproc+1
  readf,lun,t,n_pdf,pdf_dx,pdf_min,pdf_mean,pdf_rms
  print,t,n_pdf,pdf_dx,pdf_min,pdf_mean,pdf_rms
  ;
  pdf_xx=(findgen(n_pdf)+.5)*pdf_dx+pdf_min
  pdf_yy1=fltarr(n_pdf)
  readf,lun,pdf_yy1
  ;
  ;  accumulate
  ;
  if iproc eq 0 then pdf_yy=pdf_yy1 else pdf_yy=pdf_yy+pdf_yy1
  ;
  ;plot,pdf_xx,pdf_yy,yr=[0,5000],ps=10
  ;wait,.1
endfor
;
;  normalize
;
;if t gt 850. then begin
  norm=total(pdf_yy)*pdf_dx
  plot,pdf_xx,pdf_yy/norm,yr=[0,1]*1.4,ps=10,xr=[-1,1]*3
  wait,.05
;  stop
;endif
endfor
;
for iproc=0,nproc-1 do begin
  close,lun
endfor
END
