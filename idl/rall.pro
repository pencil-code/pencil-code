;;;;;;;;;;;;;;;;;;;;
;;;   rall.pro   ;;;
;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   09-Sep-2001
;;;  $Id: rall.pro,v 1.18 2002-08-11 04:00:11 brandenb Exp $
;;;
;;;  Description:
;;;   Read data from all processors and combine them into one array
;;;   for each variable.
;;;   Overwrites nx, ny, etc, from start.pro, thus you need to run
;;;   start.pro after this if you want to continue working on data
;;;   from individual processors (e.g. if you want to run r.pro).

function param2
; Dummy to keep IDL from complaining. The real param() routine will be
; compiled below
end

;
;  need to run start first: check whether this has been done
;
if (n_elements(started) le 0) then begin
  message, "You need to run start.pro first: use `.rnew start'"
endif

;
;  Read startup parameters
;
pfile=datatopdir+'/'+'param2.nml'
dummy=findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  print, 'Generating and reading param2.nml..'
  spawn, '$PENCIL_HOME/bin/nl2idl -f param2 -m tmp/param2.nml > tmp/param2.pro'
  resolve_routine, 'param2', /IS_FUNCTION
  par2=param2()
  if (lhydro) then begin
    cs0=par2.cs0 & nu=par2.nu
;  cs0=1. & nu=0.
  endif
  if (lentropy) then begin
    hcond0=par2.hcond0 & hcond1=par2.hcond1 & hcond2=par2.hcond2
    Luminosity=par2.Luminosity & wheat=par2.wheat
    cool=par2.cool & wcool=par2.wcool
    Fbot=par2.Fbot
  endif
endif else begin
  print, 'Warning: cannot find file ', pfile
endelse

;
;  read global sizes
;
nprocx=0L & nprocy=0L & nprocz=0L
close,1
openr,1,datatopdir+'/'+'dim.dat'
readf,1,mx,my,mz,nvar
readf,1,prec
readf,1,nghostx,nghosty,nghostz
readf,1,nprocx,nprocy,nprocz
close,1
;
ncpus = nprocx*nprocy*nprocz
;
;  read local sizes
;
datadir=datatopdir+'/proc0'
mxloc=0L & myloc=0L & mzloc=0L
;
close,1
openr,1,datadir+'/'+'dim.dat'
readf,1,mxloc,myloc,mzloc
close,1
;
nxloc=mxloc-2*nghostx
nyloc=myloc-2*nghosty
nzloc=mzloc-2*nghostz
;
;  read data
;
default, dataopdir, 'tmp'
default, file, 'var.dat'
;
x = fltarr(mx) & y = fltarr(my) & z = fltarr(mz)
xloc = fltarr(mxloc) & yloc = fltarr(myloc) & zloc = fltarr(mzloc)
if (lhydro) then begin
  uu    = fltarr(mx,my,mz,3)*one
  uu_loc = fltarr(mxloc,myloc,mzloc,3)*one
endif
if (ldensity) then begin
  lnrho = fltarr(mx,my,mz)*one
  lnrho_loc = fltarr(mxloc,myloc,mzloc)*one
endif
if (lentropy ) then begin
  ss = fltarr(mx,my,mz)*one
  ss_loc = fltarr(mxloc,myloc,mzloc)*one
endif
if (lmagnetic) then begin
  aa = fltarr(mx,my,mz,3)*one
  aa_loc = fltarr(mxloc,myloc,mzloc,3)*one
endif
;
for i=0,ncpus-1 do begin        ; read data from individual files
  tag='proc'+str(i)
  datadir=datatopdir+'/proc'+strtrim(i,2)
  ; read processor position
  dummy=''
  ipx=0L &ipy=0L &ipz=0L
  close,1
  openr,1,datadir+'/dim.dat'
  readf,1, dummy
  readf,1, dummy
  readf,1, dummy
  readf,1, ipx,ipy,ipz
  ; read data
  close,1
  openr,1, datadir+'/'+file, /F77
    ;
    if iuu ne 0 and ilnrho ne 0 and ient ne 0 and iaa ne 0 then begin
      id='MHD with entropy'
      readu,1,uu_loc,lnrho_loc,ss_loc,aa_loc
    end else if iuu ne 0 and ilnrho ne 0 and ient eq 0 and iaa ne 0 then begin
      id='hydro without entropy, but with magnetic field'
      readu,1,uu_loc,lnrho_loc,aa_loc
    end else if iuu ne 0 and ilnrho ne 0 and ient ne 0 and iaa eq 0 then begin
      id='hydro with entropy, but no magnetic field'
      readu,1,uu_loc,lnrho_loc,ss_loc
    end else if iuu ne 0 and ilnrho ne 0 and ient eq 0 and iaa eq 0 then begin
      id='hydro with no entropy and no magnetic field'
      readu,1,uu_loc,lnrho_loc
    end else if iuu ne 0 and ilnrho eq 0 and ient eq 0 and iaa eq 0 then begin
      id='just velocity (Burgers)'
      readu,1,uu_loc
    end else if iuu eq 0 and ilnrho eq 0 and ient eq 0 and iaa ne 0 then begin
      id='just magnetic ffield (kinematic)'
      readu,1,aa_loc
    end else begin
      id='not prepared...'
    end
    if (i eq 0) then begin
      print, id
      print, FORMAT='(A,$)', "Reading: "
    endif
    print, FORMAT='(A," ",$)', tag
    ;
  readu,1, t, xloc, yloc, zloc
  close,1
  ;
  ;  Don't overwrite ghost zones of processor to the left (and
  ;  accordingly in y and z direction makes a difference on the
  ;  diagonals)
  ;
  if (ipx eq 0) then begin
    i0x=ipx*nxloc & i1x=i0x+mxloc-1
    i0xloc=0 & i1xloc=mxloc-1
  endif else begin
    i0x=ipx*nxloc+nghostx & i1x=i0x+mxloc-1-nghostx
    i0xloc=nghostx & i1xloc=mxloc-1
  endelse
  ;
  if (ipy eq 0) then begin
    i0y=ipy*nyloc & i1y=i0y+myloc-1
    i0yloc=0 & i1yloc=myloc-1
  endif else begin
    i0y=ipy*nyloc+nghosty & i1y=i0y+myloc-1-nghosty
    i0yloc=nghosty & i1yloc=myloc-1
  endelse
  ;
  if (ipz eq 0) then begin
    i0z=ipz*nzloc & i1z=i0z+mzloc-1
    i0zloc=0 & i1zloc=mzloc-1
  endif else begin
    i0z=ipz*nzloc+nghostz & i1z=i0z+mzloc-1-nghostz
    i0zloc=nghostz & i1zloc=mzloc-1
  endelse
  ;
  x[i0x:i1x] = xloc[i0xloc:i1xloc]
  y[i0y:i1y] = yloc[i0yloc:i1yloc]
  z[i0z:i1z] = zloc[i0zloc:i1zloc]
  uu [i0x:i1x,i0y:i1y,i0z:i1z,*] =  $
      uu_loc [i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*]
  lnrho[i0x:i1x,i0y:i1y,i0z:i1z]   = $
      lnrho_loc[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc]
  if (ient ne 0) then ss[i0x:i1x,i0y:i1y,i0z:i1z]   = $
      ss_loc[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc]
  if (iaa ne 0) then aa [i0x:i1x,i0y:i1y,i0z:i1z,*] =  $
      aa_loc [i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*]
endfor
print
;
xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
;
;  Summarise data
;
xyz = ['x', 'y', 'z']
fmt = '(A,4G15.6)'
print, ' var        minval         maxval            mean           rms'
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

print,'t = ',t

; reset datadir to more reasonable default
datadir=datatopdir+'/proc0'
;
;  reset boundary values for (full) physical domain (not sub-domain)
;
l1=3 & l2=mx-4
m1=3 & m2=my-4
n1=3 & n2=mz-4
;
;  fix z3=ztop which was local top in start.pro
;
if (lgravz) then begin
  ztop=z[n2] & z3=ztop
endif
;
read_all = 1                    ; marker for r.pro

end

; End of file
