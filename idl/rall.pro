;;;;;;;;;;;;;;;;;;;;
;;;   rall.pro   ;;;
;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   09-Sep-2001
;;;
;;;  Description:
;;;   Read data from all processors and combine them into one array
;;;   for each variable.
;;;   Overwrites nx, ny, etc, from start.pro, thus you need to run
;;;   start.pro after this if you want to continue working on data
;;;   from individual processors (e.g. if you want to run r.pro).

;
;  need to run start first: check whether this has been done
;
if (n_elements(started) le 0) then begin
  message, "You need to run start.pro first: use `.rnew start'"
endif
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
uu    = fltarr(mx,my,mz,3)*one
lnrho = fltarr(mx,my,mz)*one
if (lentropy) then ss = fltarr(mx,my,mz)*one
if (lmagnetic) then aa = fltarr(mx,my,mz,3)*one
;
xloc = fltarr(mxloc) & yloc = fltarr(myloc) & zloc = fltarr(mzloc)
uu_loc = fltarr(mxloc,myloc,mzloc,3)*one
lnrho_loc = fltarr(mxloc,myloc,mzloc)*one
if (lentropy) then ss_loc = fltarr(mxloc,myloc,mzloc)*one
if (lmagnetic) then aa_loc = fltarr(mxloc,myloc,mzloc,3)*one
;
for i=0,ncpus-1 do begin        ; read data from individual files
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
  if (lentropy and lmagnetic) then readu,1, uu_loc, lnrho_loc, ss_loc, aa_loc
  if ((not lentropy) and lmagnetic) then readu,1, uu_loc, lnrho_loc, aa_loc
  if (lentropy and not lmagnetic)   then readu,1, uu_loc, lnrho_loc, ss_loc
  if ((not lentropy) and (not lmagnetic)) then readu,1, uu_loc, lnrho_loc
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
  if (lentropy) then ss[i0x:i1x,i0y:i1y,i0z:i1z]   = $
      ss_loc[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc]
  if (lmagnetic) then aa [i0x:i1x,i0y:i1y,i0z:i1z,*] =  $
      aa_loc [i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*]
endfor
;
xx = spread(x, [1,2], [my,mz])
yy = spread(x, [0,2], [mx,mz])
zz = spread(x, [0,1], [mx,my])

print, ' var        minval       maxval'
for j=0,2 do print, 'uu_', strtrim(j,2), ' =', minmax(uu(*,*,*,j))
print, 'lnrho  =', minmax(lnrho)
if (lentropy) then print, 'ss  =', minmax(ss)
if (lmagnetic) then $
    for j=0,2 do print, 'aa_', strtrim(j,2), ' =', minmax(aa(*,*,*,j))
;
print,'t = ',t

; reset datadir to more reasonable default
datadir=datatopdir+'/proc0'

read_all = 1                    ; marker for r.pro

end

; End of file
