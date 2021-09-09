;
; $Id$
;
;   Read particle dimension data.
;
pro pc_read_qdim,nqpar=nqpar,mqvar=mqvar,object=object,datadir=datadir, $
    PRINT=PRINT, QUIET=QUIET
COMPILE_OPT IDL2,HIDDEN
;
; Default data directory
;
datadir = pc_get_datadir(datadir)
;
; Initialize / set default returns for ALL variables
;
nqpar=0L
mqvar=0L

;
; Get a unit number
;
GET_LUN, file

filename=datadir+'/qdim.dat'
;
; Check for existence and read the data
;
if (file_test(filename)) then begin
  if is_valid(object,'QDIM',filename) then return
  IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' + filename + '...'
  openr,file,filename
  readf,file,nqpar,mqvar
  FREE_LUN,file
endif else begin
  FREE_LUN,file
  message, 'ERROR: cannot find file ' + filename
endelse
;
; Build structure of all the variables
;
object = CREATE_STRUCT(name='PC_QDIM:'+strtrim(filename,2), ['nqpar','mqvar'], nqpar, mqvar)
;
; If requested print a summary
;
if (keyword_set(PRINT)) then begin
  print, '     (nqpar,mqvar) = (',nqpar,',',mqvar,')'
endif

end
