; $Id: pc_read_param.pro,v 1.11 2006-06-04 18:47:52 ajohan Exp $
;
;   Read param.nml
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;  $Date: 2006-06-04 18:47:52 $
;  $Revision: 1.11 $
;
;  27-nov-02/tony: coded mostly from Wolgang's start.pro
;
;  REQUIRES: external 'nl2idl' perl script (WD)
;  
;  
function param
COMPILE_OPT IDL2,HIDDEN
; Dummy to keep IDL from complaining. The real param() routine will be
; compiled below
end
;
pro pc_read_param, object=object, dim=dim, datadir=datadir, $
                   param2=param2, PRINT=PRINT,QUIET=QUIET,HELP=HELP
COMPILE_OPT IDL2,HIDDEN
  COMMON pc_precision, zero, one
; If no meaningful parameters are given show some help!
  IF ( keyword_set(HELP) ) THEN BEGIN
    print, "Usage: "
    print, ""
    print, "pc_read_param, object=object,"
    print, "               datadir=datadir, proc=proc,"
    print, "               /PRINT, /QUIET, /HELP,"
    print, "               /param2"
    print, ""
    print, "Returns the parameters of a Pencil-Code run."
    print, "Returns zeros and empty in all variables on failure."
    print, ""
    print, "   datadir: specify the root data directory. Default is './data'        [string]"
    print, ""
    print, "   object: optional structure in which to return all the above as tags  [struct]"
    print, ""
    print, "   /param2: for reading param2.nml"
    print, "   /PRINT: instruction to print all variables to standard output"
    print, "   /QUIET: instruction not to print any 'helpful' information"
    print, "   /HELP: display this usage information, and exit"
    print
    return
  ENDIF

; Default data directory

default, datadir, 'data'
if (n_elements(dim) eq 0) then pc_read_dim, datadir=datadir, object=dim, $
    quiet=quiet
pc_set_precision,dim=dim,quiet=quiet   ; check precision is set
precision=dim.precision


; Build the full path and filename
if keyword_set(param2) then begin
  filename=datadir+'/param2.nml'
endif else begin
  filename=datadir+'/param.nml'
endelse

;If double precision, force input from params.nml to be doubles
if ((precision eq 'S') or (precision eq 's')) then begin ; single precision
  nl2idl_d_opt = ''
endif else if ((precision eq 'D') or (precision eq 'd')) then begin ; double precision
  nl2idl_d_opt = '-d'
endif  ;; Erroneous case already caught by pc_set_precision

; Check for existance and read the data
dummy = findfile(filename, COUNT=found)
if (found gt 0) then begin
    IF ( not keyword_set(QUIET) ) THEN print, 'Reading ' + filename + '...'
    spawn, 'for d in . $TMPDIR $TMP /tmp /var/tmp; do if [ -d $d -a -w $d ]; then echo $d; fi; done', /SH, result
    if (strlen(result[0])) le 0 then begin
      message, "Can't find writeable directory for temporary files"
    endif else begin
      tmpdir = result[0]
    endelse
    tmpfile = tmpdir+'/param.pro'
    ;; Write content of param.nml to temporary file:
    spawn, '$PENCIL_HOME/bin/nl2idl '+nl2idl_d_opt+' -m '+filename+'> ' $
         + tmpfile , result
    spawn, "sed -i -e 's/,$/, \$/g' param.pro", result
    ;; Compile that file. Should be easy, but is incredibly awkward, as
    ;; there is no way in IDL to compile a given file at run-time
    ;; outside the command line:
    ;; Save old path and pwd
    _path = !path
    cd, tmpdir, CURRENT=_pwd
    !path = '.:'
    resolve_routine, 'param', /IS_FUNCTION
    ;; Restore old path and pwd
    !path = _path & cd, _pwd
    ;; Delete temporary file
    ; file_delete, tmpfile      ; not in IDL <= 5.3
    spawn, 'rm -f '+tmpfile, /SH
    object = param()
                                ;
                                ;
endif else begin
    message, 'Warning: cannot find file '+ filename
endelse

; If requested print a summary
if keyword_set(PRINT) then begin
    print, 'For GLOBAL calculation domain:'
    print, '    NO SUMMARY INFORMATION CONFIGURED - edit pc_read_param.pro'
endif

end
