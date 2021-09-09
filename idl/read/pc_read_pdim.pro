;;
;; $Id$
;;
;;  Read particle dimension data.
;;
pro pc_read_pdim, npar=npar, mpvar=mpvar, mpaux=mpaux, object=object, datadir=datadir, $
    print=print, quiet=quiet
compile_opt IDL2,HIDDEN
;
; Default data directory.
;
datadir = pc_get_datadir(datadir)
;
; Initialize all variables.
;
npar=0L
mpvar=0L
npar_stalk=0L
mpaux=0L
;
; Check for existence and read the data.
;
filename=datadir+'/pdim.dat'
if (file_test(filename)) then begin
  if is_valid(object,'PDIM',filename) then return
  if (not keyword_set(quiet)) then print, 'Reading ' + filename + '...'
  get_lun, file
  openr, file, filename
  readf, file, npar, mpvar, npar_stalk, mpaux
  close, file
  free_lun, file
endif else $
  message, 'ERROR: cannot find file ' + filename
;
; Build structure of all the variables.
;
object = create_struct(name='PC_PDIM:'+strtrim(filename,2), ['npar','mpvar','npar_stalk','mpaux'], npar, mpvar, npar_stalk, mpaux)
;
; Print a summary if requested.
;
if (keyword_set(print)) then $
  print, '     (npar,mpvar,npar_stalk,mpaux) = (',npar,',',mpvar,',',npar_stalk,',',mpaux,')'
;
end
