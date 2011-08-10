;  $Id$
;
;  Read startup parameters
;
;forward_function safe_get_tag

pfile = datatopdir+'/'+'param2.nml'
if (file_test(pfile)) then begin

  if (quiet le 2) then print, 'Reading param2.nml..'
  spawn, 'for d in . $TMPDIR $TMP /tmp /var/tmp; do if [ -d $d -a -w $d ]; then echo $d; fi; done', /SH, result
  if (strlen(result[0])) le 0 then begin
    message, "Can't find writeable directory for temporary files"
  endif else begin
    tmpdir = result[0]
  endelse
  tmpfile = tmpdir+'/param2.pro'
  ;; Write content of param2.nml to temporary file:
  spawn, '$PENCIL_HOME/bin/nl2idl '+nl2idl_d_opt+' -f param2 -m '+datatopdir+'/param2.nml > ' $
         + tmpfile , result
  ;; Compile that file. Should be easy, but is incredibly awkward, as
  ;; there is no way in IDL to compile a given file at run-time
  ;; outside the command line:
  ;; Save old path and pwd
  _path = !path
  cd, tmpdir, CURRENT=_pwd
  !path = '.:'
  resolve_routine, 'param2', /IS_FUNCTION
  ;; Restore old path and pwd
  !path = _path & cd, _pwd
  ;; Delete temporary file
  ; file_delete, tmpfile      ; not in IDL <= 5.3
  spawn, 'rm -f '+tmpfile, /SH
  par2 = param2()

  ;; Abbreviate some frequently used parameters
  default, STRUCT=par2, ['nu','cs0'], 0.

  if (lhydro) then begin
    nu=par2.nu
  endif

  if (ldensity) then begin
    cs0=par2.cs0
  endif
  if (lentropy) then begin
    hcond0=par2.hcond0 & hcond1=par2.hcond1 & hcond2=par2.hcond2
    luminosity=par2.luminosity & wheat=par2.wheat
    cool=par2.cool & wcool=par2.wcool
    Fbot=par2.Fbot
  endif
  if (lmagnetic) then begin
    eta=par2.eta
    b_ext=par2.b_ext
  endif
  if (leos_fixed_ionization) then begin
    yH0=par2.yH0
  endif
endif else begin
  if (quiet le 4) then print, 'Note: the file ', pfile,' does not yet exist.'
  par2={lwrite_aux:0L}
endelse
