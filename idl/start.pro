;;;;;;;;;;;;;;;;;;;
;;;  start.pro  ;;;
;;;;;;;;;;;;;;;;;;;
;;;
;;; Initialise coordinate arrays, detect precision and dimensions.
;;; Typically run only once before running `r.pro' and other
;;; plotting/analysing scripts.

common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
;
@xder_6th_bc
@yder_6th_bc
@zder_6th_bc
@xder2_6th_bc
@yder2_6th_bc
@zder2_6th_bc
;
default, datatopdir, 'tmp'
default, datadir, datatopdir+'/proc0'
default, file, 'var.dat'
;
;  Read the dimensions and precision (single or double) from dim.dat
;
nx=0L & ny=0L & nz=0L & nw=0L
prec=''
nghx=0L & nghy=0L & nghz=0L
;
close,1
openr,1,datadir+'/'+'dim.dat'
readf,1,nx,ny,nz,nw
readf,1,prec
readf,1,nghx,nghy,nghz
close,1
;
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
Lx=zero & Ly=zero & Lz=zero
cs0=zero & gamma=zero & gamma1=zero
bx_ext=zero & by_ext=zero & bz_ext=zero
pfile=datatopdir+'/'+'param.dat'
dummy=findfile(pfile, COUNT=cpar)
if (cpar gt 0) then begin
  print, 'Reading grid.dat..'
  openr,1, datatopdir+'/'+'param.dat', /F77
  readu,1, Lx,Ly,Lz
  readu,1, cs0,gamma,gamma1
  close,1
endif else begin
  print, 'Warning: cannot find file ', pfile
endelse
;
;  Read grid
;
xx=fltarr(nx,ny,nz)*one & yy=xx & zz=xx
t=zero
x=fltarr(nx)*one & y=fltarr(ny)*one & z=fltarr(nz)*one
Lx=zero &  Ly=zero &  Lz=zero
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
print, '..done'
;
END

; End of file
