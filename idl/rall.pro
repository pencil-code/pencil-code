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
;;;   Overwrites nx, ny, etc, thus you need to run start.pro after
;;;   this if you want to continue working on data from individual
;;;   processors (e.g. if you want to run r.pro).

;
;  read global sizes
;
nprocx=0L & nprocy=0L & nprocz=0L
close,1
openr,1,datatopdir+'/'+'dim.dat'
readf,1,nx,ny,nz,nw
readf,1,prec
readf,1,nghx,nghy,nghz
readf,1,nprocx,nprocy,nprocz
close,1
;
nxtot = nx+2*nghx
nytot = ny+2*nghy
nztot = nz+2*nghz
ncpus = nprocx*nprocy*nprocz
;
;  read local sizes
;
datadir=datatopdir+'/proc0'
nxloc=0L & nyloc=0L & nzloc=0L
;
close,1
openr,1,datadir+'/'+'dim.dat'
readf,1,nxloc,nyloc,nzloc
close,1
nxloctot=nxloc+2*nghx
nyloctot=nyloc+2*nghy
nzloctot=nzloc+2*nghz
;
;  read data
;
default, dataopdir, 'tmp'
default, file, 'var.dat'
;
x = fltarr(nxtot) & y = fltarr(nytot) & z = fltarr(nztot)
uu = fltarr(nxtot,nytot,nztot,3)*one
lam = fltarr(nxtot,nytot,nztot)*one
aa = fltarr(nxtot,nytot,nztot,3)*one
;
xloc = fltarr(nxloctot) & yloc = fltarr(nyloctot) & zloc = fltarr(nzloctot)
uu_loc = fltarr(nxloctot,nyloctot,nzloctot,3)*one
lam_loc = fltarr(nxloctot,nyloctot,nzloctot)*one
aa_loc = fltarr(nxloctot,nyloctot,nzloctot,3)*one
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
  readu,1, uu_loc, lam_loc, aa_loc
  readu,1, t, xloc, yloc, zloc
  close,1
  ;
  ;  Don't overwrite ghost zones of processor to the left (and
  ;  accordingly in y and z direction makes a difference on the
  ;  diagonals)
  ;
  if (ipx eq 0) then begin
    i0x=ipx*nxloc & i1x=i0x+nxloctot-1
    i0xloc=0 & i1xloc=nxloctot-1
  endif else begin
    i0x=ipx*nxloc+nghx & i1x=i0x+nxloctot-1-nghx
    i0xloc=nghx & i1xloc=nxloctot-1
  endelse
  ;
  if (ipy eq 0) then begin
    i0y=ipy*nyloc & i1y=i0y+nyloctot-1
    i0yloc=0 & i1yloc=nyloctot-1
  endif else begin
    i0y=ipy*nyloc+nghy & i1y=i0y+nyloctot-1-nghy
    i0yloc=nghy & i1yloc=nyloctot-1
  endelse
  ;
  if (ipz eq 0) then begin
    i0z=ipz*nzloc & i1z=i0z+nzloctot-1
    i0zloc=0 & i1zloc=nzloctot-1
  endif else begin
    i0z=ipz*nzloc+nghz & i1z=i0z+nzloctot-1-nghz
    i0zloc=nghz & i1zloc=nzloctot-1
  endelse
  ;
  x[i0x:i1x] = xloc[i0xloc:i1xloc]
  y[i0y:i1y] = yloc[i0yloc:i1yloc]
  z[i0z:i1z] = zloc[i0zloc:i1zloc]
  uu [i0x:i1x,i0y:i1y,i0z:i1z,*] =  $
      uu_loc [i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*]
  lam[i0x:i1x,i0y:i1y,i0z:i1z]   = $
      lam_loc[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc]
  aa [i0x:i1x,i0y:i1y,i0z:i1z,*] =  $
      aa_loc [i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*]
endfor
;
xx = spread(x, [1,2], [nytot,nztot])
yy = spread(x, [0,2], [nxtot,nztot])
zz = spread(x, [0,1], [nxtot,nytot])

print, ' var        minval       maxval'
for j=0,2 do print, 'uu_', strtrim(j,2), ' =', minmax(uu(*,*,*,j))
print, 'lam  =', minmax(lam)
for j=0,2 do print, 'aa_', strtrim(j,2), ' =', minmax(aa(*,*,*,j))
;
print,'t = ',t

; reset datadir to more reasonable default
datadir=datatopdir+'/proc0'
end

; End of file
