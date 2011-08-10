;;;;;;;;;;;;;;;;;;;
;;;  start.pro  ;;;
;;;;;;;;;;;;;;;;;;;
;;;
;;; Initialise coordinate arrays, detect precision and dimensions.
;;; Typically run only once before running `r.pro' and other
;;; plotting/analysing scripts.
;;; $Id$

common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
;
@xder_6th_ghost
@yder_6th_ghost
@zder_6th_ghost
@xder2_6th_ghost
@yder2_6th_ghost
@zder2_6th_ghost
;
default, proc, 0
if (n_elements(datatopdir) eq 0) then datatopdir=pc_get_datadir()
default, datafile, 'var.dat'
datadir=datatopdir+'/proc'+str(proc)
;
;  Read the dimensions and precision (single or double) from dim.dat
;
mx=0L & my=0L & mz=0L & nvar=0L
prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
openr,lun,datadir+'/'+'dim.dat',/get_lun
readf,lun,mx,my,mz,nvar
readf,lun,prec
readf,lun,nghostx,nghosty,nghostz
close,lun
free_lun,lun
;
mw=mx*my*mz  ;(this must be calculated; its not in dim.dat)
prec = (strtrim(prec,2))        ; drop leading zeros
prec = strmid(prec,0,1)
if ((prec eq 'S') or (prec eq 's')) then begin
  one = 1.e0
endif else if ((prec eq 'D') or (prec eq 'd')) then begin
  one = 1.D0
endif else begin
  print, "prec = `", prec, "' makes no sense to me"
  STOP
endelse
zero = 0*one
;
;  the following files contain the positions of variables in f
;
@data/index
print,'nname=',nname
;
;  Read grid
;
t=zero
x=fltarr(mx)*one & y=fltarr(my)*one & z=fltarr(mz)*one
dx=zero & dy=zero & dz=zero & dxyz=zero
gridfile=datadir+'/'+'grid.dat'
if (file_test(gridfile)) then begin
  print, 'Reading grid.dat..'
  openr,lun, gridfile, /F77, /get_lun
  readu,lun, t,x,y,z
  readu,lun, dx,dy,dz
  close,lun
  free_lun,lun
endif else begin
  print, 'Warning: cannot find file ', gridfile
endelse
;
;print,'calculating xx,yy,zz (comment this out if there is not enough memory)'
;xx = spread(x, [1,2], [my,mz])
;yy = spread(y, [0,2], [mx,mz])
;zz = spread(z, [0,1], [mx,my])
;
;  set boundary values for physical (sub)domain
;
l1=3 & l2=mx-4
m1=3 & m2=my-4
n1=3 & n2=mz-4
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
;
print, '..done'
;
started=1
END
