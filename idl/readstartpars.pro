;  $Id: readstartpars.pro,v 1.9 2003-08-15 18:38:22 brandenb Exp $
;
;  Read startup parameters
;
pfile = datatopdir+'/'+'param2.nml'
dummy = findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  if (quiet le 2) then print, 'Reading param2.nml..'
  spawn, '$PENCIL_HOME/bin/nl2idl -1 -m '+datatopdir+'/param2.nml', result
  res = flatten_strings(result)
  ;; For people with an unclean shell: remove everything up to the
  ;; opening brace:
  brace = strpos(res,'{')
  if (brace lt 0) then message, 'TROUBLE: no brace found in <'+res+'>'
  if ((brace ne 0) and (quiet le 4)) then begin
    print, "Your shell produces output when it shouldn't; you'd better"
    print, "fix your prompt."
    print, "Trying to clean up the mess.."
    res = strmid(res,brace)
  endif
  ;; Execute the resulting line
  if (execute('par2 = '+res) ne 1) then $
      message, 'There was a problem with param.nml', /INFO
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
endif else begin
  if (quiet le 4) then print, 'Note: the file ', pfile,' does not yet exist.'
  par2={lwrite_aux:0L}
endelse
