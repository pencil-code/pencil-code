; $Id: r.pro,v 1.45 2003-04-09 14:50:56 theine Exp $

;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced on one processor
;;; You should have run `start.pro' once before.
;;; $Id: r.pro,v 1.45 2003-04-09 14:50:56 theine Exp $

function param2
; Dummy to keep IDL from complaining. The real param() routine will be
; compiled below
end

;
;  read data
;
if ((n_elements(started) le 0) or (n_elements(read_all) gt 0)) then begin
  message, "You need to run start.pro first: use `.rnew start'"
endif
undefine, read_all
;
default, datadir, 'data'
default, varfile, 'var.dat'
;
if (lhydro)     then uu    = fltarr(mx,my,mz,3)*one
if (ldensity)   then lnrho = fltarr(mx,my,mz  )*one
if (lentropy)   then ss    = fltarr(mx,my,mz  )*one
if (lmagnetic)  then aa    = fltarr(mx,my,mz,3)*one
if (lradiation) then ff    = fltarr(mx,my,mz,3)*one
if (lradiation) then ee    = fltarr(mx,my,mz  )*one
if (lpscalar )  then lncc  = fltarr(mx,my,mz  )*one
;
;  Read startup parameters
;
pfile = datatopdir+'/'+'param2.nml'
dummy = findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  print, 'Reading param2.nml..'
  spawn, '$PENCIL_HOME/bin/nl2idl -1 -m '+datatopdir+'/param2.nml', result
  res = flatten_strings(result)
  ;; For people with an unclean shell: remove everything up to the
  ;; opening brace:
  brace = strpos(res,'{')
  if (brace lt 0) then message, 'TROUBLE: no brace found in <'+res+'>'
  if (brace ne 0) then begin
    print, "Your shell produces output when it shouldn't; you'd better"
    print, "fix your prompt."
    print, "Trying to clean up the mess.."
    res = strmid(res,brace)
  endif
  ;; Execute the resulting line
  if (execute('par2 = '+res) ne 1) then $
      message, 'There was a problem with param.nml', /INFO
  if (lhydro) then begin
    cs0=par2.cs0 & nu=par2.nu
  endif
  if (lentropy) then begin
    hcond0=par2.hcond0 & hcond1=par2.hcond1 & hcond2=par2.hcond2
    luminosity=par2.luminosity & wheat=par2.wheat
    cool=par2.cool & wcool=par2.wcool
    Fbot=par2.Fbot
  endif
endif else begin
  print, 'Warning: cannot find file ', pfile
endelse

;
;  Read data
;
close,1
openr,1, datadir+'/'+varfile, /F77
  ;
  if iuu ne 0 and ilnrho ne 0 and ient ne 0 and iaa ne 0 then begin
    print,'MHD with entropy'
    readu,1,uu,lnrho,ss,aa
  end else if iuu ne 0 and ilnrho ne 0 and ient eq 0 and iaa ne 0 then begin
    print,'hydro without entropy, but with magnetic field'
    readu,1,uu,lnrho,aa
  end else if iuu ne 0 and ilnrho ne 0 and ient ne 0 and ie ne 0 then begin
    print,'hydro with entropy, density and radiation'
    readu,1,uu,lnrho,ss,ee,ff
  end else if iuu ne 0 and ilnrho ne 0 and ient ne 0 and iaa eq 0 then begin
    print,'hydro with entropy, but no magnetic field'
    readu,1,uu,lnrho,ss
  end else if iuu ne 0 and ilnrho ne 0 and ilncc ne 0 and iaa eq 0 then begin
    print,'hydro with entropy, but no magnetic field'
    readu,1,uu,lnrho,lncc
  end else if iuu ne 0 and ilnrho ne 0 and ient eq 0 and iaa eq 0 then begin
    print,'hydro with no entropy and no magnetic field'
    readu,1,uu,lnrho
  end else if iuu ne 0 and ilnrho eq 0 and ient eq 0 and iaa eq 0 then begin
    print,'just velocity (Burgers)'
    readu,1,uu
  end else if iuu eq 0 and ilnrho eq 0 and ient eq 0 and iaa ne 0 then begin
    print,'just magnetic field (kinematic)'
    readu,1,aa
  end else if iuu eq 0 and ilnrho eq 0 and ient eq 0 and iaa eq 0 and ilncc ne 0 then begin
    print,'just passive scalar (no field nor hydro)'
    readu,1,lncc
  end else if iuu eq 0 and ilnrho ne 0 and ient eq 0 and iaa eq 0 then begin
    print,'just density (probably just good for tests)'
    readu,1,lnrho
  end else if iuu eq 0 and ilnrho eq 0 and ient eq 0 and iaa eq 0 and ie ne 0 then begin
    print,'just radiation'
    readu,1,ee,ff
  end else begin
    print,'not prepared...'
  end
  ;
if (lshear) then begin
  readu,1, t, x, y, z, dx, dy, dz, deltay
end else begin
  readu,1, t, x, y, z, dx, dy, dz
end
;
xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
rr = sqrt(xx^2+yy^2+zz^2)
;
;  Summarise data
;
xyz = ['x', 'y', 'z']
fmt = '(A,4G15.6)'
print, '  var            minval         maxval          mean           rms'
;AB:  does not work ok for kinematic case. I suggest to use only f.
;WD: Ought to work now (the namelist lphysics need to know about the
;    relevant logicals). What is f ?
;
if (lhydro) then $
    for j=0,2 do $
      print, FORMAT=fmt, 'uu_'+xyz[j]+'   =', $
      minmax(uu(*,*,*,j)), mean(uu(*,*,*,j),/DOUBLE), rms(uu(*,*,*,j),/DOUBLE)
if (ldensity) then $
    print, FORMAT=fmt, 'lnrho  =', $
    minmax(lnrho), mean(lnrho,/DOUBLE), rms(lnrho,/DOUBLE)
if (lentropy) then $
    print, FORMAT=fmt, 'ss     =', $
      minmax(ss), mean(ss,/DOUBLE), rms(ss,/DOUBLE)
if (lpscalar) then $
    print, FORMAT=fmt, 'lncc   =', $
      minmax(lncc), mean(lncc,/DOUBLE), rms(lncc,/DOUBLE)
if (lradiation) then $
    for j=0,2 do $
      print, FORMAT=fmt, 'ff_'+xyz[j]+'   =', $
      minmax(ff(*,*,*,j)), mean(ff(*,*,*,j),/DOUBLE), rms(ff(*,*,*,j),/DOUBLE)
if (lradiation) then $
    print, FORMAT=fmt, 'ee     =', $
      minmax(ee), mean(ee,/DOUBLE), rms(ee,/DOUBLE)
if (lmagnetic) then begin
    for j=0,2 do $
      print, FORMAT=fmt, 'aa_'+xyz[j]+'   =', $
      minmax(aa(*,*,*,j)), mean(aa(*,*,*,j),/DOUBLE), rms(aa(*,*,*,j),/DOUBLE)
    if (cpar gt 0) then begin
      eta=par2.eta
    end
end
;
print,'t = ',t
;
if (par.lradiation ne 0) then begin
  if (par.output_Qrad) then begin
    Qrad=fltarr(mx,my,mz)*one
    Srad=fltarr(mx,my,mz)*one
    kappa=fltarr(mx,my,mz)*one
    TT=fltarr(mx,my,mz)*one
    readu,1,Qrad,Srad,kappa,TT
    print, FORMAT=fmt, 'Qrad   =', $
      minmax(Qrad), mean(Qrad,/DOUBLE), rms(Qrad,/DOUBLE)
    print, FORMAT=fmt, 'Srad   =', $
      minmax(Srad), mean(Srad,/DOUBLE), rms(Srad,/DOUBLE)
    print, FORMAT=fmt, 'kappa   =', $
      minmax(kappa), mean(kappa,/DOUBLE), rms(kappa,/DOUBLE)
    print, FORMAT=fmt, 'TT   =', $
      minmax(TT), mean(TT,/DOUBLE), rms(TT,/DOUBLE)
  end
end
;
if (par.lionization ne 0) then begin
  if (par.output_yH) then begin
    yyH=fltarr(mx,my,mz)*one
    readu,1,yyH
    print, FORMAT=fmt, 'yyH   =', $
      minmax(yyH), mean(yyH,/DOUBLE), rms(yyH,/DOUBLE)
  end
end
;
close,1
;
END

; End of file
