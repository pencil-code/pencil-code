; $Id: r.pro,v 1.21 2002-06-17 21:03:43 dobler Exp $

;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced on one processor
;;; You should have run `start.pro' once before.

;
;  read data
;
if ((n_elements(started) le 0) or (n_elements(read_all) gt 0)) then begin
  message, "You need to run start.pro first: use `.rnew start'"
endif
undefine, read_all
;
default, datadir, 'tmp'
default, file, 'var.dat'
;
if (iuu ne 0)    then uu    = fltarr(mx,my,mz,3)*one
if (ilnrho ne 0) then lnrho = fltarr(mx,my,mz  )*one
if (ient ne 0)   then ss = fltarr(mx,my,mz  )*one
if (iaa ne 0)    then aa = fltarr(mx,my,mz,3)*one
;
;  Read startup parameters
;
pfile=datatopdir+'/'+'param2.nml'
dummy=findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  print, 'Generating and reading param2.nml..'
  spawn, '../../../bin/nl2idl -f param2 tmp/param2.nml > tmp/param2.pro'
  resolve_routine, 'param2', /IS_FUNCTION, /COMPILE_FULL_FILE
  par2=param2()
  if (lhydro) then begin
    cs0=par2.cs0 & nu=par2.nu
;  cs0=1. & nu=0.
  endif
  if (lentropy) then begin
    hcond0=par2.hcond0 & hcond1=par2.hcond1
    hcond2=par2.hcond2 & whcond=par2.whcond
    cheat=par2.cheat & wheat=par2.wheat
    cool=par2.cool & wcool=par2.wcool
    Fheat=par2.Fheat
  endif
endif else begin
  print, 'Warning: cannot find file ', pfile
endelse

;
;  Read data
;
close,1
openr,1, datadir+'/'+file, /F77
  ;
  if iuu ne 0 and ilnrho ne 0 and ient ne 0 and iaa ne 0 then begin
    print,'MHD with entropy'
    readu,1,uu,lnrho,ss,aa
  end else if iuu ne 0 and ilnrho ne 0 and ient eq 0 and iaa ne 0 then begin
    print,'hydro without entropy, but with magnetic field'
    readu,1,uu,lnrho,aa
  end else if iuu ne 0 and ilnrho ne 0 and ient ne 0 and iaa eq 0 then begin
    print,'hydro with entropy, but no magnetic field'
    readu,1,uu,lnrho,ss
  end else if iuu ne 0 and ilnrho ne 0 and ient eq 0 and iaa eq 0 then begin
    print,'hydro with no entropy and no magnetic field'
    readu,1,uu,lnrho
  end else if iuu ne 0 and ilnrho eq 0 and ient eq 0 and iaa eq 0 then begin
    print,'just velocity (Burgers)'
    readu,1,uu
  end else if iuu eq 0 and ilnrho eq 0 and ient eq 0 and iaa ne 0 then begin
    print,'just magnetic ffield (kinematic)'
    readu,1,aa
  end else begin
    print,'not prepared...'
  end
  ;
readu,1, t, x, y, z
close,1
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
if (lmagnetic) then $
    for j=0,2 do $
      print, FORMAT=fmt, 'aa_'+xyz[j]+'   =', $
      minmax(aa(*,*,*,j)), mean(aa(*,*,*,j),/DOUBLE), rms(aa(*,*,*,j),/DOUBLE)
;
print,'t = ',t
;
END

; End of file
