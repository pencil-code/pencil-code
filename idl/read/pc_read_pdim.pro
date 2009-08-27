;
; $Id$
;
;   Read particle dimension data.
;
pro pc_read_pdim, npar=npar, mpvar=mpvar, object=object, datadir=datadir, $
    PRINT=PRINT, QUIET=QUIET
COMPILE_OPT IDL2,HIDDEN
;
; Default data directory
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
;
; Initialize / set default returns for ALL variables
;
npar=0L
mpvar=0L
npar_stalk=0L
;
; Get a unit number
;
GET_LUN, file

filename=datadir+'/pdim.dat'
;
; Check for existence and read the data
;
dummy=findfile(filename, COUNT=found)
if (found gt 0) then begin
  IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' + filename + '...'
  openr,file,filename
  readf,file,npar,mpvar,npar_stalk
  FREE_LUN,file
endif else begin
  FREE_LUN,file
  message, 'ERROR: cannot find file ' + filename
endelse
;
; Build structure of all the variables
;
object = CREATE_STRUCT(name=objectname, ['npar','mpvar','npar_stalk'], npar, mpvar, npar_stalk)
;
; If requested print a summary
;
if (keyword_set(PRINT)) then begin
  print, '     (npar,mpvar,npar_stalk) = (',npar,',',mpvar,',',npar_stalk,')'
endif

end
