pro pc_read_pdf2d,pdf,x,y,noplot=noplot,datadir=datadir,arr_t=arr_t,arr_pdf_dx=arr_pdf_dx,arr_pdf_dy=arr_pdf_dy,arr_pdf_min=arr_pdf_min,arr_pdf_max=arr_pdf_max,arr_pdfy_min=arr_pdfy_min,arr_pdfy_max=arr_pdfy_max
;
; $Id$
  lnoplot=0 & if (keyword_set(noplot)) then lnoplot=1
  print,"lnoplot=",lnoplot
;
;  Reading two-dimensional probability distribution function (pdf) files 
;
t=0. & n_pdf=0L & pdf_dx=0. & pdf_min=0. & pdf_mean=0. & pdf_rms=0. & lun=123
;
;  open all files
;  note: this method only works on up to 128 processors
;
if (keyword_set(datadir)) then begin
   datadir=datadir
endif else begin
   datadir="data"
endelse
filename=datadir+'/pdf2d_FI_mixfrac.dat'
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
pdf=fltarr(n_pdf, n_pdfy,nt)
pdf_yy=fltarr(n_pdf, n_pdfy)
x=fltarr(n_pdf,nt)
y=fltarr(n_pdfy,nt)
;
print,"Reading ",nt," times..."
it=0
while not EOF(lun) do begin
   readf,lun,t, n_pdf, n_pdfy, pdf_dx, pdf_dy, pdf_max, pdf_min, pdfy_max, pdfy_min, pdf_mean, pdf_rms
   readf,lun,pdf_yy
;
   arr_t[it]          = t
   arr_pdf_dx[it]     = pdf_dx
   arr_pdf_dy[it]     = pdf_dy
   arr_pdf_min[it]    = pdf_min
   arr_pdf_max[it]    = pdf_max
   arr_pdfy_min[it]   = pdfy_min
   arr_pdfy_max[it]   = pdfy_max
   pdf[*,*,it] = pdf_yy
   ;
   ; Make x and y arrays to set ranges for the two axis
   ;
   x(0,it)=arr_pdf_min[it]
   for i=1,n_pdf-1 do begin
      x(i,it)=x(i-1,it)+arr_pdf_dx[it]
   end
   y(0,it)=arr_pdfy_min[it]
   for i=1,n_pdfy-1 do begin
      y(i,it)=y(i-1,it)+arr_pdf_dy[it]
   end
;   
   it=it+1
endwhile
close,lun
;
; Plot results
;
if (not lnoplot) then begin
for it=0,nt-1 do begin
   ;
   ; Modify pdf array for better visualization
   ;
   yy=pdf(*,*,it)
   good=where(yy eq 0)
   yy(good)=0.001
   xx=alog(yy)
   ;
   ; Visualize
   ;
      contour,xx-min(xx),x(*,it),y(*,it),/fill,nlev=40,xtit="mixture fraction",ytit="flame index",charsize=2.0
      wait,1.5
endfor
endif


;
END
