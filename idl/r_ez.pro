;  $Id: r_ez.pro,v 1.5 2002-08-11 04:00:11 brandenb Exp $
;
;  this routine is an easy-to-use replacement of r.pro,
;  which can sometimes fail (if the perl script fails,
;  for example).
;
common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
;
;  compile the right derivative routines
;
@xder_6th_ghost
@yder_6th_ghost
@zder_6th_ghost
@xder2_6th_ghost
@yder2_6th_ghost
@zder2_6th_ghost
;
mx=1 & my=1 & mz=1 & nvar=1
openr,1,'tmp/dim.dat'
readf,1,mx,my,mz,mvar
close,1
print,mx,my,mz,mvar
;
dir='tmp/proc0/'
file='var.dat'
f=fltarr(mx,my,mz,mvar)
x=fltarr(mx)
y=fltarr(my)
z=fltarr(mz)
openr,1,dir+file,/f77
readu,1,f
readu,1,t,x,y,z,dx,dy,dz
close,1
print,'t=',t
 
@tmp/index

;  calculate bb
if iax ne 0 then bb=curl(f(*,*,*,iax-1:iaz-1))
END
