;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_read_all_dummy.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   09-Sep-2001
;;;  $Id$
;;;
;;;  Description:
;;;   Dummy routine for reading data from all processors.
;;;   Useful when using derivative routines on data from save files.
;;;   Overwrites nx, ny, etc, from start.pro, thus you need to run
;;;   start.pro after this if you want to continue working on data
;;;   from individual processors (e.g. if you want to run r.pro).

function param2
COMPILE_OPT IDL2,HIDDEN
; Dummy to keep IDL from complaining. The real param2() routine will be
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
@readstartpars

;
;  read global sizes
;
nprocx=0L & nprocy=0L & nprocz=0L
close,1
openr,1,datatopdir+'/'+dimfile
readf,1,mx,my,mz,nvar,naux
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
openr,1,datadir+'/'+dimfile
readf,1,mxloc,myloc,mzloc
close,1
;
nxloc=mxloc-2*nghostx
nyloc=myloc-2*nghosty
nzloc=mzloc-2*nghostz
;
;  read data
;
default, dataopdir, 'data'
default, varfile, 'var.dat'
default, ldustvelocity, 0
default, ldustdensity, 0
;
dx=2.*!pi/nx
dy=2.*!pi/ny
dz=2.*!pi/nz
;
x = dx*findgen(mx) & y = dy*findgen(my) & z = dz*findgen(mz)
xloc = fltarr(mxloc) & yloc = fltarr(myloc) & zloc = fltarr(mzloc)

;
;  reset boundary values and nx,ny,nz for (full) physical domain (not
;  sub-domain)
;
l1=3 & l2=mx-4
m1=3 & m2=my-4
n1=3 & n2=mz-4
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
;
;  fix z3=ztop which was local top in start.pro
;
if (lgravz) then begin
  ztop=z[n2] & z3=ztop
endif
;
;  free memory
;
undefine, uu_loc
undefine, lnrho_loc
undefine, ss_loc
undefine, aa_loc
undefine, lncc_loc
;
read_all = 1                    ; marker for r.pro

end

; End of file
