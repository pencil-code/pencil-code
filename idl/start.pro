;;;;;;;;;;;;;;;;;;;
;;;  start.pro  ;;;
;;;;;;;;;;;;;;;;;;;
;;;
;;; Initialise coordinate arrays, detect precision and dimensions.
;;; Typically run only once before running `r.pro' and other
;;; plotting/analysing scripts.

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
default, datatopdir, 'tmp'
default, file, 'var.dat'
datadir=datatopdir+'/proc'+str(proc)
;
;  Read the dimensions and precision (single or double) from dim.dat
;
mx=0L & my=0L & mz=0L & nvar=0L
prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
close,1
openr,1,datadir+'/'+'dim.dat'
readf,1,mx,my,mz,nvar
readf,1,prec
readf,1,nghostx,nghosty,nghostz
close,1
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
;  Read startup parameters
;
pfile=datatopdir+'/'+'param.nml'
dummy=findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  print, 'Reading param.nml..'
  spawn, '../../bin/nl2idl tmp/param.nml > tmp/param.pro'
  @tmp/param.pro
  x0=par.xyz0[0] & y0=par.xyz0[1] & z0=par.xyz0[2]
  Lx=par.Lxyz[0] & Ly=par.Lxyz[1] & Lz=par.Lxyz[2]
  cs0=par.cs0 & rho0=par.rho0
  gamma=par.gamma & gamma1=gamma-1.
  lgravz=par.lgravz & lgravr = par.lgravr
  lforcing=par.lforcing
  lentropy=par.lentropy
  lmagnetic=par.lmagnetic
  if (lgravz) then begin
    z1=par.z1 & z2=par.z2
    ztop=par.ztop & z3=ztop
    gravz=par.gravz
  endif
  if (lentropy) then begin
    hcond0=par.hcond0 & hcond1=par.hcond1
    hcond2=par.hcond2 & whcond=par.whcond
    mpoly0=par.mpoly0 & mpoly1=par.mpoly1
    mpoly2=par.mpoly2 & isothtop=par.isothtop
  endif
endif else begin
  print, 'Warning: cannot find file ', pfile
endelse
;
;  Read grid
;
t=zero
x=fltarr(mx)*one & y=fltarr(my)*one & z=fltarr(mz)*one
dx=zero &  dy=zero &  dz=zero & dxyz=zero
gfile=datadir+'/'+'grid.dat'
dummy=findfile(gfile, COUNT=cgrid)
if (cgrid gt 0) then begin
  print, 'Reading grid.dat..'
  openr,1, gfile, /F77
  readu,1, t,x,y,z
  readu,1, dx,dy,dz
  close,1
endif else begin
  print, 'Warning: cannot find file ', gfile
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
print, '..done'
;
started=1
END
