; $Id: pc_read_pdim.pro,v 1.1 2005-01-09 16:56:36 ajohan Exp $
;
;   Read stuff from pdim.dat
;
pro pc_read_pdim, npar=npar, mpvar=mpvar, object=object, datadir=datadir, $
    PRINT=PRINT, QUIET=QUIET
COMPILE_OPT IDL2,HIDDEN

; Default data directory
default,datadir,'data'

;
; Initialize / set default returns for ALL variables
;
npar=0L
mpvar=0L

; Get a unit number
GET_LUN, file

filename=datadir+'/pdim.dat'

; Check for existence and read the data
dummy=findfile(filename, COUNT=found)
if (found gt 0) then begin
  IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' + filename + '...'

  openr,file,filename
  readf,file,npar,mpvar
  FREE_LUN,file
endif else begin
  FREE_LUN,file
  message, 'ERROR: cannot find file ' + filename
endelse

; Build structure of all the variables
object = CREATE_STRUCT(name=object, ['npar','mpvar'], npar, mpvar)

; If requested print a summary
if (keyword_set(PRINT)) then begin
  print, '     (npar,mpvar) = (',npar,',',mpvar,')'
endif

end
