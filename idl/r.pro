; $Id: r.pro,v 1.49 2003-06-19 21:42:09 mee Exp $

;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced on one processor
;;; You should have run `start.pro' once before.
;;; $Id: r.pro,v 1.49 2003-06-19 21:42:09 mee Exp $

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
endif else begin
  print, 'Warning: cannot find file ', pfile
endelse

;
;  Read data
;
if (par2.lwrite_aux ne 0) then totalvars=nvar+naux else totalvars=nvar
varcontent=REPLICATE({varcontent, variable:'UNKNOWN', idlvar:'dummy', idlinit:'fltarr(mx,my,mz)*one', skip:0},totalvars+1)

;
; Declare ALL variables that may occur
;
varcontent[iuu].variable = 'Velocity (uu)'
varcontent[iuu].idlvar  = 'uu'
varcontent[iuu].idlinit = 'fltarr(mx,my,mz,3)*one'
varcontent[iuu].skip  = 2

varcontent[ilnrho].variable = 'Log density (lnrho)'
varcontent[ilnrho].idlvar  = 'lnrho'
varcontent[ilnrho].idlinit = 'fltarr(mx,my,mz)*one'

varcontent[ient].variable = 'Entropy (ss)'
varcontent[ient].idlvar   = 'ss'
varcontent[ient].idlinit  = 'fltarr(mx,my,mz)*one'

varcontent[iaa].variable = 'Magnetic vector potential (aa)'
varcontent[iaa].idlvar   = 'aa'
varcontent[iaa].idlinit  = 'fltarr(mx,my,mz,3)*one'
varcontent[iaa].skip  = 2

varcontent[ifx].variable = 'Radiation vector ?something? (ff)'
varcontent[ifx].idlvar   = 'ff'
varcontent[ifx].idlinit  = 'fltarr(mx,my,mz,3)*one'
varcontent[ifx].skip  = 2

varcontent[ie].variable = 'Radiation scalar ?something? (ee)'
varcontent[ie].idlvar   = 'ee'
varcontent[ie].idlinit  = 'fltarr(mx,my,mz)*one'

varcontent[ilncc].variable = 'Log passive scalar (lncc)'
varcontent[ilncc].idlvar   = 'lncc'
varcontent[ilncc].idlinit  = 'fltarr(mx,my,mz)*one'

if (par2.lwrite_aux ne 0) then begin
    varcontent[iQrad].variable = 'Radiation (Qrad)'
    varcontent[iQrad].idlvar   = 'Qrad'
    varcontent[iQrad].idlinit  = 'fltarr(mx,my,mz)*one'
    
    varcontent[iSrad].variable = 'Radiation (Srad)'
    varcontent[iSrad].idlvar   = 'Srad'
    varcontent[iSrad].idlinit  = 'fltarr(mx,my,mz)*one'
    
    varcontent[ikappa].variable = 'Radiation (kappa)'
    varcontent[ikappa].idlvar   = 'kappa'
    varcontent[ikappa].idlinit  = 'fltarr(mx,my,mz)*one'
    
; May need special condition as can be maux or mvar variable?
    varcontent[iTT].variable = 'Temperature (TT)'
    varcontent[iTT].idlvar   = 'TT'
    varcontent[iTT].idlinit  = 'fltarr(mx,my,mz)*one'

    varcontent[iyH].variable = 'Hydrogen ionization fraction (yyH)'
    varcontent[iyH].idlvar   = 'yyH'
    varcontent[iyH].idlinit  = 'fltarr(mx,my,mz)*one'

    varcontent[ishock].variable = 'Shock characteristic (nu_shock)'
    varcontent[ishock].idlvar   = 'nu_shock'
    varcontent[ishock].idlinit  = 'fltarr(mx,my,mz)*one'
end

varcontent[0].variable = 'UNKNOWN'
varcontent[0].idlvar   = 'dummy'
varcontent[0].idlinit  = '0.'
varcontent[0].skip  = 0

; Prepare for read
res=''
content=''
for i=1,totalvars do begin
  res     = res + ',' + varcontent[i].idlvar
  content = content + ', ' + varcontent[i].variable
  ; Initialise variable
  if (varcontent[i].variable eq 'UNKNOWN') then $
           message, 'Unknown variable at position ' + str(i)  $
                                    + ' needs declaring in r.pro', /INFO   
  if (execute(varcontent[i].idlvar+'='+varcontent[i].idlinit,0) ne 1) then $
           message, 'Error initialising ' + varcontent[i].variable $
                                    +' - '+ varcontent[i].idlvar, /INFO
;If it's a vector quantity skip the required number of elements
  i=i+varcontent[i].skip
end

dummy=0.

content = strmid(content,2)
print,'File contains: '+content

close,1
openr,1, datadir+'/'+varfile, /F77
  if (execute('readu,1'+res) ne 1) then $
           message, 'Error reading: ' + 'readu,1'+res
 
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
fmt = '(A9,A,4G15.6)'
print, '  var             minval         maxval          mean           rms'
;AB: Does not work ok for kinematic case. I suggest to use only f.
;WD: Ought to work now (the namelist lphysics need to know about the
;    relevant logicals). What is f ?
;
for i=1,totalvars do begin
  if (varcontent[i].skip eq 2) then begin
      for j=0,2 do begin
          cmd = "print, FORMAT=fmt,strmid('"+varcontent[i].idlvar+"_'+xyz["+str(j)+"]+'        ',0,8),'=', " $
            + "minmax("+varcontent[i].idlvar+"(*,*,*,"+str(j)+")), " $
            + "mean("+varcontent[i].idlvar+"(*,*,*,"+str(j)+"),/DOUBLE), " $
            + "rms("+varcontent[i].idlvar+"(*,*,*,"+str(j)+"),/DOUBLE)"
          if (execute(cmd,1) ne 1) then $
                          message, 'Error printing stats for ' + varcontent[i].variable         
      end
  end else begin
      cmd = "print, FORMAT=fmt,strmid('"+varcontent[i].idlvar+"        ',0,8),'=', " $
        + "minmax("+varcontent[i].idlvar+"(*,*,*)), " $
        + "mean("+varcontent[i].idlvar+"(*,*,*),/DOUBLE), " $
        + "rms("+varcontent[i].idlvar+"(*,*,*),/DOUBLE)
      if (execute(cmd,1) ne 1) then $
        message, 'Error printing stats for ' + varcontent[i].variable         
  end
  i=i+varcontent[i].skip
end


if (lmagnetic) then begin
    if (cpar gt 0) then begin
      eta=par2.eta
    end
end
;
print,'t = ',t
;
close,1
;
END

; End of file
