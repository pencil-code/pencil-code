FUNCTION inputs,file, _EXTRA=extra
;; Read scalar field from binary data file. If file does not exist,
;; return NaNs.
;; Any keywords (e.g. /SWAP_ENDIAN) are passed on to the OPENR
;; statement.
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
;
field=fltarr(nx,ny,nz)
if ((findfile(file))[0] ne '') then begin
  close,1
  openr,1,file,/f77,_EXTRA=extra
  readu,1,field
  close,1
endif else begin
  message,/informational,"No such file: '" + file + "'"
  field = field*!VALUES.F_NAN
endelse
;
return, field
;
END
