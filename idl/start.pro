;;;;;;;;;;;;;;;;;;;
;;;  start.pro  ;;;
;;;;;;;;;;;;;;;;;;;
;;;
;;; Initialise coordinate arrays, detect precision and dimensions.
;;; Typically run only once before running `r.pro' and other
;;; plotting/analysing scripts.
;;; $Id: start.pro,v 1.37 2002-08-11 04:00:11 brandenb Exp $

function param
; Dummy to keep IDL from complaining. The real param() routine will be
; compiled below
  message, $
      "This dummy  function should never be called" $
      + "-- make sure you have `tmp' in your !path ."
end

common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
;
@xder_6th_ghost
@yder_6th_ghost
@zder_6th_ghost
@xder2_6th_ghost
@yder2_6th_ghost
@zder2_6th_ghost
;
;  The following avoids a mysterious bug when using esrg_legend later
;  (box size was wrong, because lenstr(['a','b']) would be wrong,
;  because xyout would write all letters onto one posision) ..
;
@lenstr
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
;  the following files contain the positions of variables in f
;
@tmp/index
print,'nname=',nname
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
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
;
;  Read startup parameters
;
pfile=datatopdir+'/'+'param.nml'
dummy=findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  print, 'Reading param.nml..'
  spawn, '$PENCIL_HOME/bin/nl2idl -f param -m tmp/param.nml > tmp/param.pro'
  resolve_routine, 'param', /IS_FUNCTION
  par=param()
  x0=par.xyz0[0] & y0=par.xyz0[1] & z0=par.xyz0[2]
  Lx=par.Lxyz[0] & Ly=par.Lxyz[1] & Lz=par.Lxyz[2]
  ;
  lhydro    = par.lhydro
  lgravz    = par.lgravz
  lgravr    = par.lgravr
  ldensity  = par.ldensity
  lforcing  = par.lforcing
  lentropy  = par.lentropy
  lmagnetic = par.lmagnetic
  lradiation= par.lradiation
  lpscalar  = par.lpscalar
  ;
  if (ldensity) then begin
    cs0=par.cs0 & rho0=par.rho0
    gamma=par.gamma & gamma1=gamma-1.
    cs20 = cs0^2 & lnrho0 = alog(rho0)
  endif
  ;
  if (lgravz) then begin
    z1=par.z1 & z2=par.z2
    zref=par.zref
    gravz=par.gravz
    ztop=z[n2] & z3=ztop
  endif
  ;
  if (lentropy) then begin
    mpoly0=par.mpoly0 & mpoly1=par.mpoly1
    mpoly2=par.mpoly2 & isothtop=par.isothtop
  endif
endif else begin
  print, 'Warning: cannot find file ', pfile
endelse

;
print, '..done'
;
started=1
END
