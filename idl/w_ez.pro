;  $Id: w_ez.pro,v 1.3 2002-11-27 11:25:01 brandenb Exp $
;
;  this routine is an easy-to-use write routine
;  This assumes that the data have been read in with r_ez.pro
;
common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
;
openw, 1, datadir+'/'+varfile,/F77
writeu, 1, f
if (lshear) then writeu, 1, t, x,y,z, dx,dy,dz, deltay else writeu, 1, t, x,y,z, dx,dy,dz
close, 1
print, 't=',t
;
END
