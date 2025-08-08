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
; Count number of times
;
nt=0
while not EOF(lun) do begin
   readf,lun,t, n_pdf, n_pdfy, pdf_dx, pdf_dy, pdf_max, pdf_min, pdfy_max, pdfy_min, pdf_mean, pdf_rms
   pdf_yy=fltarr(n_pdf, n_pdfy)
   readf,lun,pdf_yy
   nt=nt+1
endwhile
close,lun
;
; Store data in arrays of the correct size
;
openr,lun,filename
arr_t=fltarr(nt)
arr_pdf_dx=fltarr(nt)
arr_pdf_dy=fltarr(nt)
arr_pdf_min=fltarr(nt)
arr_pdf_max=fltarr(nt)
arr_pdfy_min=fltarr(nt)
arr_pdfy_max=fltarr(nt)
arr_pdf_yy=fltarr(n_pdf, n_pdfy,nt)
pdf_yy=fltarr(n_pdf, n_pdfy)
;
print,"Reading ",nt," times..."
nt=0
while not EOF(lun) do begin
   readf,lun,t, n_pdf, n_pdfy, pdf_dx, pdf_dy, pdf_max, pdf_min, pdfy_max, pdfy_min, pdf_mean, pdf_rms
   readf,lun,pdf_yy
;
   arr_t[nt]          = t
   arr_pdf_dx[nt]     = pdf_dx
   arr_pdf_dy[nt]     = pdf_dy
   arr_pdf_min[nt]    = pdf_min
   arr_pdf_max[nt]    = pdf_max
   arr_pdfy_min[nt]   = pdfy_min
   arr_pdfy_max[nt]   = pdfy_max
   arr_pdf_yy[*,*,nt] = pdf_yy
;   
   print,"t=",t
   nt=nt+1
endwhile
close,lun

x=fltarr(n_pdf)
y=fltarr(n_pdfy)
print,"arr_t=",arr_t
;
; Plot results
;
for it=0,nt-1 do begin
   ;
   ; Modify pdf array for better visualization
   ;
   yy=arr_pdf_yy(*,*,it)
   good=where(yy eq 0)
   yy(good)=0.001
   xx=alog(yy)
   ;
   ; Make x and y arrays to set ranges for the two axis
   ;
   x(0)=arr_pdf_min[it]
   for i=1,n_pdf-1 do begin
      x(i)=x(i-1)+arr_pdf_dx[it]
   end
   y(0)=arr_pdfy_min[it]
   for i=1,n_pdfy-1 do begin
      y(i)=y(i-1)+arr_pdf_dy[it]
   end
   ;
   ; Visualize
   ;
   contour,xx-min(xx),x,y,/fill,nlev=40,xtit="mixture fraction",ytit="flame index",charsize=2.0
   ;
   wait,1.5
endfor


;
END
