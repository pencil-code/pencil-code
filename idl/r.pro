; $Id: r.pro,v 1.18 2002-06-13 15:55:49 brandenb Exp $

;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced on one processor
;;; Assumes you have run `start.pro' once before.

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
;if (cpar gt 0) then begin
;  print, 'Reading param2.nml..'
;  spawn, '../../../bin/nl2idl tmp/param2.nml > tmp/param2.pro'
;  @tmp/param2.pro
; cs0=par.cs0 & nu=par.nu
;  cs0=1. & nu=0.
; hcond0=par.hcond0 & hcond1=par.hcond1
; hcond2=par.hcond2 & whcond=par.whcond
; cheat=par.cheat & wheat=par.wheat
; cool=par.cool & wcool=par.wcool
; Fheat=par.Fheat
;endif else begin
;  print, 'Warning: cannot find file ', pfile
;endelse

;
;  Read data
;
;AB: the following is not quite save, nvar=7 could mean other things...
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
    print,'just magnetic ffield (kinematic)
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
print, ' var        minval         maxval            mean           rms'
;AB:  does not work ok for kinematic case. I suggest to use only f.
;for j=0,2 do $
;    print, FORMAT=fmt, 'uu_'+xyz[j]+' =', $
;    minmax(uu(*,*,*,j)), mean(uu(*,*,*,j),/DOUBLE), rms(uu(*,*,*,j),/DOUBLE)
;print, FORMAT=fmt, 'lnrho  =', $
;    minmax(lnrho), mean(lnrho,/DOUBLE), rms(lnrho,/DOUBLE)
;if (lentropy) then $
;    print, FORMAT=fmt, 'ss  =', $
;      minmax(ss), mean(ss,/DOUBLE), rms(ss,/DOUBLE)
;if (lmagnetic) then $
;    for j=0,2 do $
;      print, FORMAT=fmt, 'aa_'+xyz[j]+' =', $
;      minmax(aa(*,*,*,j)), mean(aa(*,*,*,j),/DOUBLE), rms(aa(*,*,*,j),/DOUBLE)
;
print,'t = ',t
;
END

; End of file
