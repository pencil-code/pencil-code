;  $Id: readstartpars.pro,v 1.11 2003-12-30 13:17:54 dobler Exp $
;
;  Read startup parameters
;
pfile = datatopdir+'/'+'param2.nml'
dummy = findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  if (quiet le 2) then print, 'Reading param2.nml..'
  spawn, '$PENCIL_HOME/bin/nl2idl -1 -m -M '+strtrim(maxtags,2) $
         + ' '+datatopdir+'/param2.nml', result
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
  ;; Execute the resulting line(s); format is `{ block1 } { block2 } ..'.
  ;; Need to add each block in its own execute() call in order to
  ;; remain below the limit of structure tags that can be set within
  ;; one execute() statement (typically 398 or 570).
  par2 = { dummy: 0 }            ; initialize structure for appending
  brace = strpos(res,'}')
  iblock = 0
  repeat begin
    iblock = iblock+1
    print, "  param2.nml: block "+strtrim(iblock,2)
    block = strmid(res,0,brace+1)
    if (execute('par2 = create_struct(par2,'+block+')') ne 1) then $
        message, 'There was a problem with param2.nml', /INFO
    res = strmid(res,brace+1)
    brace = strpos(res,'}')
  endrep until (brace lt 0)
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
  if (lionization_fixed) then begin
    yH0=par2.yH0
  endif
endif else begin
  if (quiet le 4) then print, 'Note: the file ', pfile,' does not yet exist.'
  par2={lwrite_aux:0L}
endelse
