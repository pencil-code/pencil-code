;  $Id: w_ez.pro,v 1.1 2002-08-21 03:33:22 brandenb Exp $
;
;  this routine is an easy-to-use write routine
;  This assumes that the data have been read in with r_ez.pro
;
common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
;
openw,1,dir+file,/f77
writeu,1,f
writeu,1,t,x,y,z,dx,dy,dz,deltay
close,1
print,'t=',t
;
END
