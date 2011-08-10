FUNCTION inputs, file, DOUBLE=double, ONE=one,_EXTRA=extra
;; Read scalar field from binary data file. If file does not exist,
;; return NaNs.
;; Keywords:
;;   DOUBLE  -- flag for double precision data
;;   ONE     -- if 1.0, use single precision, if 1.0D0, use double precision
;;              Use with start.pro's variable ONE like this:
;;                var=inputs('var.dat',ONE=ONE)
;; All other keywords (e.g. /SWAP_ENDIAN) are passed on to the OPENR
;; statement.
  common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
  ;
  default, ONE, 1.0
  if (keyword_set(double)) then ONE=1.D0
  ;
  field=fltarr(nx,ny,nz)*ONE
  if (file_test(file)) then begin
    openr,lun,file,/f77,/get_lun,_EXTRA=extra
    readu,lun,field
    close,lun
    free_lun,lun
  endif else begin
    message,/informational,"No such file: '" + file + "'"
    field = field*!VALUES.F_NAN
  endelse
  ;
  return, field

END
