;
; $Id$
;
;   Read particle dimension data.
;
pro pc_read_spdim,nspar=nspar,mspvar=mspvar,object=object,datadir=datadir, $
    PRINT=PRINT, QUIET=QUIET
COMPILE_OPT IDL2,HIDDEN
;
; Default data directory
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
;
; Initialize / set default returns for ALL variables
;
nspar=0L
mspvar=0L

;
; Get a unit number
;
GET_LUN, file

filename=datadir+'/spdim.dat'
;
; Check for existence and read the data
;
if (file_test(filename)) then begin
  IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' + filename + '...'
  openr,file,filename
  readf,file,nspar,mspvar
  FREE_LUN,file
endif else begin
  FREE_LUN,file
  message, 'ERROR: cannot find file ' + filename
endelse
;
; Build structure of all the variables
;
object = CREATE_STRUCT(name=objectname, ['nspar','mspvar'], nspar, mspvar)
;
; If requested print a summary
;
if (keyword_set(PRINT)) then begin
  print, '     (nspar,mspvar) = (',nspar,',',mspvar,')'
endif

end
