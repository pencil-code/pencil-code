;;
;; $Id$
;;
;;  Read particle dimension data.
;;
pro pc_read_pdim, npar=npar, mpvar=mpvar, object=object, datadir=datadir, $
    print=print, quiet=quiet
compile_opt IDL2,HIDDEN
;
; Default data directory.
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
;
; Initialize all variables.
;
npar=0L
mpvar=0L
npar_stalk=0L
;
; Check for existence and read the data.
;
filename=datadir+'/pdim.dat'
if (file_test(filename)) then begin
  if (not keyword_set(quiet)) then print, 'Reading ' + filename + '...'
  get_lun, file
  openr, file, filename
  readf, file, npar, mpvar, npar_stalk
  free_lun, file
endif else begin
  message, 'ERROR: cannot find file ' + filename
endelse
;
; Build structure of all the variables.
;
object = create_struct(name=objectname, ['npar','mpvar','npar_stalk'], npar, mpvar, npar_stalk)
;
; Print a summary if requested.
;
if (keyword_set(print)) then begin
  print, '     (npar,mpvar,npar_stalk) = (',npar,',',mpvar,',',npar_stalk,')'
endif
;
end
