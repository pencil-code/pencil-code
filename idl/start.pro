;;;;;;;;;;;;;;;;;;;
;;;  start.pro  ;;;
;;;;;;;;;;;;;;;;;;;
;;;
;;; Initialise coordinate arrays, detect precision and dimensions.
;;; Typically run only once before running `r.pro' and other
;;; plotting/analysing scripts.

common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
;
@xder_6th_ghost
@yder_6th_ghost
@zder_6th_ghost
@xder2_6th_ghost
@yder2_6th_ghost
@zder2_6th_ghost
;
default, datatopdir, 'tmp'
default, datadir, datatopdir+'/proc0'
default, file, 'var.dat'
;
;  Read the dimensions and precision (single or double) from dim.dat
;
nx=0L & ny=0L & nz=0L & nvar=0L
prec=''
nghx=0L & nghy=0L & nghz=0L
;
close,1
openr,1,datadir+'/'+'dim.dat'
readf,1,nx,ny,nz,nvar
readf,1,prec
readf,1,nghx,nghy,nghz
close,1
;
nw=nx*ny*nz  ;(this must be calculated; its not in dim.dat)
;nxtot = nx+2*nghx
;nytot = ny+2*nghy
;nztot = nz+2*nghz
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
x0=(y0=(z0=zero))
Lx=(Ly=(Lz=zero))
cs0=(gamma=(gamma1=zero))
gravz=(rho0=(grads0=zero))
z1=(z2=(z3=zero))
hcond0=(hcond1=(hcond2=(whcond=zero)))
mpoly0=(mpoly1=(mpoly2=zero))
isothtop=0L
bx_ext=(by_ext=(bz_ext=zero))
lgravz=(lgravr=0L)
;
pfile=datatopdir+'/'+'param.dat'
dummy=findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  print, 'Reading param.dat..'
  openr,1, pfile, /F77
  readu,1, x0,y0,z0,Lx,Ly,Lz
  readu,1, cs0,gamma,gamma1
  readu,1, gravz,rho0,grads0
  readu,1, z1,z2,ztop
  readu,1, hcond0,hcond1,hcond2,whcond
  readu,1, mpoly0,mpoly1,mpoly2,isothtop
  readu,1, lgravz,lgravr
  close,1
endif else begin
  print, 'Warning: cannot find file ', pfile
endelse
;
;  Read grid
;
t=zero
x=fltarr(nx)*one & y=fltarr(ny)*one & z=fltarr(nz)*one
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
;  set boundary values for physical (sub)domain
;
l1=3 & l2=nx-4
m1=3 & m2=ny-4
n1=3 & n2=nz-4
;
print, '..done'
;
started=1
END
