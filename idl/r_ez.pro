;  $Id: r_ez.pro,v 1.11 2003-03-18 22:29:19 brandenb Exp $
;
;  this routine is an easy-to-use replacement of r.pro,
;  which can sometimes fail (if nml2idl fails, for example).
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
mx=1 & my=1 & mz=1 & mvar=1
openr,1,datatopdir+'/dim.dat'
readf,1,mx,my,mz,mvar
close,1
print,mx,my,mz,mvar
;
default, datadir, datatopdir+'/proc0/'
default, varfile,'var.dat'
f=fltarr(mx,my,mz,mvar)
x=fltarr(mx)
y=fltarr(my)
z=fltarr(mz)
openr,1,datadir+'/'+varfile,/f77
readu,1,f
if (lshear) then readu,1,t,x,y,z,dx,dy,dz,deltay else readu,1,t,x,y,z,dx,dy,dz
close,1
;
for ivar=0,mvar-1 do print,ivar,min(f(*,*,*,ivar)),max(f(*,*,*,ivar))
print,'t=',t
 
@data/index

;  calculate bb
if iax ne 0 then bb=curl(f(*,*,*,iax-1:iaz-1))
END
