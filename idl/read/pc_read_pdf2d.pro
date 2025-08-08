;
; $Id$
;
;  Reading two-dimensional probability distribution function (pdf) files 
;
t=0. & n_pdf=0L & pdf_dx=0. & pdf_min=0. & pdf_mean=0. & pdf_rms=0. & lun=123
;
;  open all files
;  note: this method only works on up to 128 processors
;
filename='data/pdf2d_FI_mixfrac.dat'
close,lun
openr,lun,filename
;
nt=1
for it=0,nt-1 do begin
   readf,lun,t, n_pdf, n_pdfy, pdf_dx, pdf_dy, pdf_max, pdf_min, pdfy_max, pdfy_min, pdf_mean, pdf_rms
   print,t, n_pdf, n_pdfy, pdf_dx, pdf_dy, pdf_max, pdf_min, pdfy_max, pdfy_min, pdf_mean, pdf_rms
   pdf_yy=fltarr(n_pdf, n_pdfy)
   readf,lun,pdf_yy
endfor
close,lun
;
; Plot results
;
tvscl,pdf_yy
;
END
